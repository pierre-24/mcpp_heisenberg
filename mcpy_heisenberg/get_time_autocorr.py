import argparse
import h5py

import numpy as np
from matplotlib import pyplot as plt

from results import Result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'))
    parser.add_argument('-n', type=int, default=1, help='On what to autocorr (0=energy, 1=spins)')
    parser.add_argument('-N', type=int, default=200, help='Number of time series')

    args = parser.parse_args()

    with h5py.File(args.file, 'r') as f:
        results = Result(f)

        autocov = Result.sample_autocovariance(results.aggregated_data()[:, args.n], args.N)

        # rough estimate of t
        t = np.where(autocov < np.exp(-1))[0]
        print("estimated t =", t[0] if len(t) > 0 else '-1')

        plt.plot(np.arange(0, args.N), autocov)
        plt.show()
