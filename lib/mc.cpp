#include <mcpp_heisenberg/mc.hpp>

#include <vector>
#include <queue>
#include <set>

namespace mch {

void IsingMonteCarloRunner::sweep(double T, double H) {
  assert(T > 0);
  assert(_initialized_energy);

  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.number_of_magnetic_sites(); ++ispin) {
    auto pair = _hamiltonian.P_i(_spins, ispin, _kB * T, -_spins.at(ispin), _muB * H);

    if (pair.first <= 0 || dis(_rng) < pair.second) {
      _spins.at(ispin) *= -1;
      _energy += pair.first;
    }
  }
}

}  // namespace mch
