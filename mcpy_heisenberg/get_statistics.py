import argparse

import h5py

from mcpy_heisenberg import Result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'))

    args = parser.parse_args()

    with h5py.File(args.file, 'r') as f:
        results = Result(f)

        print(results.T, results.H, *results.get_statistics())
