#include <mcpp_heisenberg/hamiltonian.hpp>

#include <cassert>

namespace mch {

double Hamiltonian::energy(const arma::uvec &spins) const {
  assert(spins.n_rows == _N);

  double energy = .0;

  for (auto& pair : _pairs) {
    energy += pair.second * spins.at(pair.first.first) * spins.at(pair.first.second);
  }

  return -1. * energy;  // assume unique pairs
}

}  // namespace mch
