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
    double dE = _hamiltonian.delta_energy(_spins, ispin, H / _muB);

    if ((dE < 0) || (dis(_rng) < exp(-dE / (_kB * T)))) {
      _spins.at(ispin) *= -1;
      _energy += dE;
    }
  }
}

}  // namespace mch
