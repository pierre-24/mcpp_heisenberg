import numpy as np
import h5py

from numpy.typing import NDArray


class Result:
    def __init__(self, f: h5py.File):
        self.file = f

        self.T, self.H = self.file['results/T&H'][:]
        self.number_of_sites = self.file['geometry/positions'].shape[1]
        self.N = self.file['results/configs'].shape[0]

    def aggregated_data(self) -> NDArray[float]:
        return self.file['results/aggregated_data'][:]

    def get_statistics(self) -> tuple[float, float, float, float]:
        """Return <E>/N, C_V = (<E²>-<E>²)/NT², <|M|>/N, χ = (<M²>-<M>²)/NT"""

        aggs = self.aggregated_data()

        E = np.mean(aggs[:, 0])
        E2 = np.mean(aggs[:, 0] ** 2)
        Mabs = np.mean(np.abs(aggs[:, 1]))
        M = np.mean(np.abs(aggs[:, 1]))
        M2 = np.mean(aggs[:, 1] ** 2)

        return (
            E / self.number_of_sites,
            (E2 - E**2) / self.number_of_sites / (self.T ** 2),
            Mabs / self.number_of_sites,
            (M2 - M**2) / self.number_of_sites / self.T
        )

    @staticmethod
    def sample_autocovariance(x, tmax):
        """
        Compute the autocorrelation of the time series x for t = 0,1,...,tmax-1.

        From https://hef.ru.nl/~tbudd/mct/lectures/cluster_algorithms.html
        """

        x_shifted = x - np.mean(x)
        autocorr = np.array([np.dot(x_shifted[:len(x)-t],x_shifted[t:]) / len(x) for t in range(tmax)])
        return autocorr / autocorr[0]
