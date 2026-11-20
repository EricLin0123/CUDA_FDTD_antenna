"""Save waveform and reflection plots for an antenna/strip simulation pair."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from results import load_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    path, meta = load_run(args.run)
    antenna = np.loadtxt(path / 'antenna' / 'Incident.txt') / 3
    strip = np.loadtxt(path / 'strip' / 'Incident_strip.txt') / 3
    gaussian = np.loadtxt(path / 'antenna' / 'Gauss.txt')
    if any(a.shape != (meta['steps'],) or not np.isfinite(a).all()
           for a in (antenna, strip, gaussian)):
        raise ValueError('Waveforms must be finite and match the run step count')
    reflected = antenna - strip
    time = np.arange(meta['steps']) * meta['dt'] * 1e12
    fig, ax = plt.subplots()
    for data, label in [(antenna, 'Antenna / 3'), (strip, 'Strip / 3'),
                        (reflected, 'Net reflected'), (gaussian, 'Gaussian source')]:
        ax.plot(time, data, label=label)
    ax.set(xlabel='Time (ps)', ylabel='Amplitude', title='Feed waveforms')
    ax.grid(True)
    ax.legend()
    fig.savefig(path / 'waveforms.png', dpi=150)
    plt.close(fig)

    # Preserve the original source-normalized calculation and /3 scaling.
    # Mask frequencies where the source has negligible energy.
    source = np.abs(np.fft.rfft(gaussian))
    response = np.abs(np.fft.rfft(reflected))
    valid = source > source.max() * 1e-8
    db = np.full_like(source, np.nan)
    db[valid] = 20 * np.log10(np.maximum(response[valid] / source[valid],
                                         np.finfo(float).tiny))
    frequency = np.fft.rfftfreq(meta['steps'], meta['dt']) / 1e9
    fig, ax = plt.subplots()
    ax.plot(frequency, db)
    ax.set(xlabel='Frequency (GHz)', ylabel='Source-normalized reflection (dB)',
           title='Patch reflection (original normalization)', xlim=(0, 20), ylim=(-40, 5))
    ax.grid(True)
    fig.savefig(path / 'reflection.png', dpi=150)
    plt.close(fig)
    print(f'Saved {path / "waveforms.png"} and {path / "reflection.png"}')


if __name__ == '__main__':
    main()
