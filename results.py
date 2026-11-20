"""Shared validation for simulation results."""
import json
from pathlib import Path
import numpy as np


def load_run(path):
    path = Path(path)
    metadata = json.loads((path / 'metadata.json').read_text())
    if metadata.get('format_version') != 1 or metadata.get('status') != 'complete':
        raise ValueError('Expected a complete run with format_version 1')
    return path, metadata


def load_fields(path, metadata):
    fields = np.load(path / 'antenna' / 'fields.npy', mmap_mode='r', allow_pickle=False)
    expected = ((metadata['steps'] + metadata['sample_every'] - 1)
                // metadata['sample_every'], *metadata['grid'][:2])
    if fields.shape != expected or fields.dtype.kind != 'f':
        raise ValueError(f'Invalid field array: {fields.shape}, expected {expected}')
    return fields
