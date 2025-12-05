"""
Just a set of scripts to post-analyze the results
"""

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
        """Return <E>/N, (<E²>-<E>²)/NT², <|M|>/N, (<M²>-<|M|>²)/N"""

        aggs = self.aggregated_data()

        return (
            np.mean(aggs[:, 0]) / self.number_of_sites,
            (np.mean(aggs[:, 0] ** 2) - np.mean(aggs[:, 0])**2) / self.number_of_sites / (self.T ** 2),
            np.mean(np.abs(aggs[:, 1])) / self.number_of_sites,
            (np.mean(aggs[:, 1] ** 2) - np.mean(np.abs(aggs[:, 1]))**2) / self.number_of_sites
        )
