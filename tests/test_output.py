"""Check the C++/NumPy boundary and reject incomplete simulation results."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from results import load_run, load_fields

ROOT = Path(__file__).resolve().parents[1]


class OutputTests(unittest.TestCase):
    def test_cpp_writer_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            source = path / 'writer.cpp'
            source.write_text('''#include "output.h"
int main() {
    FieldOutput output("fields.npy", 2, 2, 3);
    double a[] = {0., 1., -2., 1e-12, 1e12, 3.141592653589793};
    output.write(a, 6);
    output.write(a, 6);
    output.close();
}
''')
            subprocess.run(['g++', '-std=c++14', '-I', str(ROOT), str(source),
                            '-o', str(path / 'writer')], check=True)
            subprocess.run([str(path / 'writer')], cwd=path, check=True)
            actual = np.load(path / 'fields.npy', mmap_mode='r')
            self.assertEqual(actual.shape, (2, 2, 3))
            np.testing.assert_array_equal(actual[0].ravel(),
                                          [0., 1., -2., 1e-12, 1e12, np.pi])
            np.testing.assert_array_equal(actual[0], actual[1])

    def test_failed_run_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / 'metadata.json').write_text(json.dumps(
                {'format_version': 1, 'status': 'failed'}))
            with self.assertRaises(ValueError):
                load_run(directory)

    def test_wrong_field_shape_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / 'antenna').mkdir()
            np.save(path / 'antenna/fields.npy', np.zeros((1, 2, 3)))
            with self.assertRaises(ValueError):
                load_fields(path, {'steps': 10, 'sample_every': 5, 'grid': [2, 3, 4]})


if __name__ == '__main__':
    unittest.main()
