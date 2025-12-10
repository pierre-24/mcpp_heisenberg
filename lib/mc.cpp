#include <mcpp_heisenberg/mc.hpp>

#include <vector>
#include <queue>
#include <set>

namespace mch {

void IsingMonteCarloRunner::sweep(double T, double H) {
  assert(T > 0);
  assert(this->_initialized_energy);

  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.number_of_magnetic_sites(); ++ispin) {
    auto pair = _hamiltonian.P_i(_spins, ispin, _kB * T, -_spins.at(ispin), _muB * H);

    if (pair.first <= 0 || dis(_rng) < pair.second) {
      _spins.at(ispin) *= -1;
      this->_energy += pair.first;
    }
  }
}

void QuantumIsingMonteCarloRunner::sweep(double T, double H) {
  assert(T > 0);
  assert(this->_initialized_energy);

  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.number_of_magnetic_sites(); ++ispin) {
    // pick a (new) value
    std::uniform_int_distribution<> dist_i(0, _possible_values.at(ispin).size() - 1);
    double new_value = _spins.at(ispin);
    while (new_value == _spins.at(ispin)) {
      new_value = _possible_values.at(ispin).at(dist_i(_rng));
    }

    auto pair = _hamiltonian.P_i(_spins, ispin, _kB * T, new_value, _muB * H);

    if (pair.first <= 0 || dis(_rng) < pair.second) {
      _spins.at(ispin) = new_value;
      this->_energy += pair.first;
    }
  }
}

}  // namespace mch
