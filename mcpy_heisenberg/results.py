import numpy as np
import h5py

from numpy.typing import NDArray


class Result:
    def __init__(self, f: h5py.File):
        self.file = f

        self.kB, self.T, self.muB, self.H = self.file['results/parameters'][:]
        self.number_of_sites = self.file['geometry/positions'].shape[1]
        self.N = self.file['results/configs'].shape[0]

    def aggregated_data(self) -> NDArray[float]:
        return self.file['results/aggregated_data'][:]

    def get_statistics(self, skip: int = 0) -> tuple[float, float, float, float, float]:
        """Return <E>/N, C_V = (<E²>-<E>²)/(N * kB * T²), <M>/N, <|M|>/N, χ = (<M²>-<M>²)/(N*kB*T)"""

        aggs = self.aggregated_data()

        E = np.mean(aggs[skip:, 0])
        E2 = np.mean(aggs[skip:, 0] ** 2)
        Mabs = np.mean(np.abs(aggs[skip:, 1]))
        M = np.mean(aggs[skip:, 1])
        M2 = np.mean(aggs[skip:, 1] ** 2)

        return (
            E / self.number_of_sites,
            (E2 - E**2) / (self.number_of_sites * self.T ** 2 * self.kB),
            M / self.number_of_sites,
            Mabs / self.number_of_sites,
            (M2 - M**2) / (self.number_of_sites * self.T * self.kB)
        )

    @staticmethod
    def sample_autocovariance(x, tmax):
        """
        Compute the autocorrelation of the time series x for t = 0,1,...,tmax-1.

        From https://hef.ru.nl/~tbudd/mct/lectures/cluster_algorithms.html
        """

        x_shifted = x - np.mean(x)
        autocorr = np.array([np.dot(x_shifted[:len(x)-t],x_shifted[t:]) / (len(x) - t) for t in range(tmax)])
        return autocorr / autocorr[0]

    def get_dist_autocorr(self, r_max: float = 10.0, round_: int = 3):
        positions = self.file['geometry/positions'][:]
        lattice = self.file['geometry/lattice_vectors'][:]
        configs = self.file['results/configs']

        N = {}
        statistics = {}

        for i in range(self.number_of_sites):
            conf_i = configs[:, i]

            for j in range(i + 1, self.number_of_sites):
                # use PBC
                r = positions[:, j] - positions[:, i]
                for k in range(3):
                    if r[k] > .5:
                        r[k] -= 1
                    elif r[k] < -.5:
                        r[k] += 1

                # get norm
                norm = round(np.linalg.norm(sum(r[k] * lattice[:, k] for k in range(3))), round_)

                if norm <= r_max:
                    conf_j = configs[:, j]
                    if norm not in N:
                        N[norm] = 0
                        statistics[norm] = 0

                    N[norm] += 1
                    statistics[norm] += np.mean(conf_i * conf_j)

        return np.array([(0, 1)] + [(k, statistics[k] / N[k]) for k in sorted(N.keys())])
