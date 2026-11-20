"""Stream a field animation to FFmpeg without temporary PNG files."""

import argparse
import shutil
import subprocess
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
from tqdm import tqdm
from results import load_run, load_fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fps", type=int, default=60)
    parser.add_argument(
        "--stride", type=int, default=1, help="Render every Nth saved frame"
    )
    parser.add_argument(
        "--max-frames", type=int, help="Limit frames for a quick preview"
    )
    parser.add_argument("--codec", choices=["libx264", "h264_nvenc"], default="libx264")
    args = parser.parse_args()
    if (
        args.fps < 1
        or args.stride < 1
        or (args.max_frames is not None and args.max_frames < 1)
    ):
        parser.error("fps, stride, and max-frames must be positive")
    path, meta = load_run(args.run)
    fields = load_fields(path, meta)
    output = args.output or path / "fields.mp4"
    if output.exists():
        parser.error(f"Output already exists: {output}")
    x = (np.arange(fields.shape[1]) + 1) * meta["dx"]
    y = (np.arange(fields.shape[2]) + 1) * meta["dy"]
    xx, yy = np.meshgrid(x, y, indexing="ij")
    fig = plt.figure(figsize=(10, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="Ez (V/m)", zlim=(-1, 1))
    ax.view_init(elev=25, azim=-140)
    indices = list(range(0, len(fields), args.stride))
    if args.max_frames:
        indices = indices[: args.max_frames]
    candidates = dict.fromkeys(
        filter(None, [shutil.which("ffmpeg"), "/usr/bin/ffmpeg"])
    )
    for candidate in candidates:
        if not Path(candidate).is_file():
            continue
        encoders = subprocess.run(
            [candidate, "-hide_banner", "-encoders"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        if any(
            len(line.split()) > 1 and line.split()[1] == args.codec
            for line in encoders.splitlines()
        ):
            matplotlib.rcParams["animation.ffmpeg_path"] = candidate
            break
    else:
        parser.error(
            f"No FFmpeg installation provides {args.codec}; install a suitable build or choose --codec h264_nvenc"
        )
    writer = FFMpegWriter(
        fps=args.fps, codec=args.codec, extra_args=["-pix_fmt", "yuv420p"]
    )
    surface = None
    try:
        with writer.saving(fig, str(output), dpi=100):
            for index in tqdm(indices, desc="Rendering"):
                if surface is not None:
                    surface.remove()
                surface = ax.plot_surface(
                    xx, yy, fields[index], cmap="viridis", vmin=-1, vmax=1
                )
                step = index * meta["sample_every"]
                ax.set_title(
                    f'Ez under the patch | step {step} | t = {step * meta["dt"] * 1e12:.3f} ps'
                )
                writer.grab_frame()
    finally:
        plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
