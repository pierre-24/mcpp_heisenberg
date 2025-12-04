#include <mcpp_heisenberg/mc.hpp>

#include <vector>
#include <queue>
#include <set>

namespace mch {

void IsingMonteCarloRunner::sweep(double T, double H) {
  assert(T > 0);

  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.number_of_magnetic_sites(); ++ispin) {
    double dE = _hamiltonian.delta_energy(_spins, ispin);

    if ((dE < 0) || (dis(_rng) < exp(-dE / T))) {
      _spins.at(ispin) *= -1;
      _energy += dE;
    }
  }
}

uint64_t IsingMonteCarloRunner::cluster_update(double T, double H) {
  assert(T > 0);

  // pick a spin
  uint64_t i = std::uniform_int_distribution<>(0, _hamiltonian.number_of_magnetic_sites() - 1)(_rng);

  // select cluster
  std::uniform_real_distribution<> dis(.0, 1.);

  std::set<uint64_t> cluster;
  cluster.insert(i);

  std::queue<uint64_t> newly_added;
  newly_added.push(i);

  while (!newly_added.empty()) {
    uint64_t j = newly_added.front();
    newly_added.pop();

    for (auto& neighbor : _hamiltonian.neighbors(j)) {
      uint64_t k = neighbor.first;
      if (!cluster.contains(k)
          && (_spins.at(j) == _spins.at(k))
          && (dis(_rng) < (1 - exp(-2 * neighbor.second / T)))) {
        cluster.insert(k);
        newly_added.push(k);
      }
    }
  }

  // swap
  for (auto& ispin : cluster) {
    _spins.at(ispin) *= -1;
  }

  // re-compute energy
  _energy = _hamiltonian.energy(_spins);

  return cluster.size();
}

}  // namespace mch
