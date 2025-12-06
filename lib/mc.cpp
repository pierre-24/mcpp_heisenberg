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

    if ((dE < 0) || (dis(_rng) < exp(-dE / (_kB * T)))) {
      _spins.at(ispin) *= -1;
      _energy += dE;
    }
  }
}

uint64_t IsingMonteCarloRunner::cluster_update(double T, double H) {
  assert(T > 0);

  // pick a spin
  uint64_t i = std::uniform_int_distribution<uint64_t>(0, _hamiltonian.number_of_magnetic_sites() - 1)(_rng);
  double invspin = _spins.at(i);

  // select cluster
  std::uniform_real_distribution<> dis(.0, 1.);

  uint64_t cluster_size = 1;
  _spins.at(i) *= -1;

  std::queue<uint64_t> newly_added;
  newly_added.push(i);

  while (!newly_added.empty()) {
    uint64_t j = newly_added.front();
    newly_added.pop();

    for (auto& neighbor : _hamiltonian.neighbors(j)) {
      uint64_t k = neighbor.first;
      if ((_spins.at(k) == invspin)
          && (dis(_rng) < (1 - exp(-2 * neighbor.second / (_kB * T))))) {
        newly_added.push(k);
        _spins.at(k) *= -1;
        cluster_size++;
      }
    }
  }

  // re-compute energy
  _energy = _hamiltonian.energy(_spins);

  return cluster_size;
}

}  // namespace mch
