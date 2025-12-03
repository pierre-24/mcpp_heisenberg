#ifndef INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
#define INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_

#include <cmath>
#include <vector>
#include <cstdint>
#include <utility>

#include <mcpp_heisenberg/arma.hpp>

namespace mch {

/**
 * Check if a value is within `[v-Δ, v+Δ]`
 */
class Approx {
 protected:
  double _val;
  double _delta{0};

 public:
  explicit Approx(double val, double delta = .0): _val{val}, _delta{delta} {}

  bool operator==(double other) {
    return fabs(other - _val) < _delta;
  }
};

using jpair_t = std::pair<std::pair<uint64_t, uint64_t>, double>;

/**
 * Heisenberg Hamiltonian
 */
class Hamiltonian {
 protected:
  /// Number of magnetic site/spins
  uint64_t _N{0};
  /// (unique) interaction pairs
  std::vector<jpair_t> _pairs;
 public:
  Hamiltonian() = default;

  /// Create an Heisenberg hamiltonian, with `N` sites and `pairs` pair interactions
  Hamiltonian(uint64_t N, const std::vector<jpair_t>& pairs): _N{N}, _pairs(pairs) {}

  /// Get number of sites
  uint64_t N() const { return _N; }

  /// Compute the energy
  double energy(const arma::uvec& spins) const;
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
