import argparse
import h5py

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

from results import Result

def cfit(x, c1, c2, zeta):
    return np.exp(-x / zeta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='H5 file', type=argparse.FileType('rb'))
    parser.add_argument('-R', type=float, default=10.0, help='Maximum distance')

    args = parser.parse_args()

    with h5py.File(args.file, 'r') as f:
        results = Result(f)

        autocov = results.get_dist_autocorr(args.R)

        r = curve_fit(cfit, autocov[:, 0], autocov[:, 1])
        print(r[0])

        X = np.linspace(0, args.R, 100)
        plt.plot(X, cfit(X, *r[0]), '--', color='grey', label='$\\zeta$={}'.format(r[0][2]))
        plt.plot(autocov[:, 0], autocov[:, 1], 'o-')
        plt.legend()
        plt.show()
