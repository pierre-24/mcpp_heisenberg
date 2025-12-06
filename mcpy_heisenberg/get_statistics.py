import argparse

import h5py

from results import Result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'), nargs='*')

    args = parser.parse_args()

    for fx in args.file:
        with h5py.File(fx, 'r') as f:
            results = Result(f)

            print(results.T, results.H, *results.get_statistics())
