import argparse
import h5py

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from results import Result

def cfit(t, tau: float):
    return np.exp(-t / tau)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'))
    parser.add_argument('-n', type=int, default=1, help='On what to autocorr (0=energy, 1=spins)')
    parser.add_argument('-N', type=int, default=100, help='Number of time series')
    parser.add_argument('-o', '--output')

    args = parser.parse_args()

    with h5py.File(args.file, 'r') as f:
        results = Result(f)

        autocov = Result.sample_autocovariance(results.aggregated_data()[:, args.n], args.N)
        T = np.arange(0, args.N)

        r = curve_fit(cfit, T, autocov)
        print('τ =', r[0][0])

        figure = plt.figure(figsize=(8, 5))
        ax = figure.subplots()

        ax.plot(T, cfit(T, *r[0]), '--', color='grey')
        ax.plot(T, autocov, 'o-')

        ax.set(xlabel='$τ$', ylabel='$\\langle s_i(t)s_j(t+τ)\\rangle$')
        plt.tight_layout()

        if args.output:
            figure.savefig(args.output)
        else:
            plt.show()
