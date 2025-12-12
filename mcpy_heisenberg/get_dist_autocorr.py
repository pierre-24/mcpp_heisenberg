import argparse
import h5py

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from results import Result

def cfit(x, c1, c2, zeta):
    return c1 + c2 * np.exp(-x / zeta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'))
    parser.add_argument('-R', type=float, default=10.0, help='Maximum distance')
    parser.add_argument('-o', '--output')

    args = parser.parse_args()

    with h5py.File(args.file, 'r') as f:
        results = Result(f)

        autocov = results.get_dist_autocorr(args.R)

        r = curve_fit(cfit, autocov[:, 0], autocov[:, 1])
        print(r[0])

        X = np.linspace(0, args.R, 100)

        figure = plt.figure(figsize=(8, 5))
        ax = figure.subplots()

        ax.plot(X, cfit(X, *r[0]), '--', color='grey')
        ax.plot(autocov[:, 0], autocov[:, 1], 'o-')

        ax.set(xlabel='$r$', ylabel='$\\langle s_i(0)s_j(r)\\rangle$')
        plt.tight_layout()

        if args.output:
            figure.savefig(args.output)
        else:
            plt.show()
