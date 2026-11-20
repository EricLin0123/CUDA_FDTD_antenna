"""Build and run the CUDA antenna and strip reference in an isolated directory."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="New run directory (must not exist)"
    )
    args = parser.parse_args()
    subprocess.run(["make", "main"], cwd=ROOT, check=True)
    output = (
        args.output or ROOT / "runs" / datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    ).resolve()
    output.mkdir(parents=True, exist_ok=False)
    metadata = dict(
        format_version=1,
        status="running",
        grid=[60, 100, 16],
        steps=8000,
        sample_every=5,
        slice_z=2,
        dt=0.441e-12,
        dx=0.389e-3,
        dy=0.400e-3,
        dz=0.265e-3,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
    manifest = output / "metadata.json"
    manifest.write_text(json.dumps(metadata, indent=2) + "\n")
    try:
        for name, executable in [("antenna", "fdtd"), ("strip", "fdtd-strip")]:
            folder = output / name
            folder.mkdir()
            with (folder / "simulation.log").open("w") as log:
                subprocess.run(
                    [str(ROOT / "build" / executable)],
                    cwd=folder,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
    except BaseException:
        metadata["status"] = "failed"
        raise
    else:
        metadata["status"] = "complete"
    finally:
        manifest.write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"Run directory: {output}", flush=True)
    print(f'Plots: {sys.executable} {ROOT / "plotEz.py"} {output}')
    print(f'Video: {sys.executable} {ROOT / "plotCPU.py"} {output}')


if __name__ == "__main__":
    main()
