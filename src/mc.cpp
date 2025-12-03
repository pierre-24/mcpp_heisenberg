#include <mcpp_heisenberg/mc.hpp>

namespace mch {

void MonteCarloRunner::sweep(double T, double H) {
  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.N(); ++ispin) {
    _spins.at(ispin) *= -1;
    double nE = _hamiltonian.energy(_spins);

    if ((nE < _energy) || (dis(_rng) < exp(-(nE - _energy) / T))) {
      _energy = nE;
    } else {
      _spins.at(ispin) *= -1;
    }
  }
}

}  // namespace mch
