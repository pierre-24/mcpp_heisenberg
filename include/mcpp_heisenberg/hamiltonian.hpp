#ifndef INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
#define INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_

#include <cmath>
#include <vector>
#include <cstdint>
#include <utility>
#include <string>

#include <mcpp_heisenberg/arma.hpp>
#include <mcpp_heisenberg/geometry.hpp>

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

  bool operator==(double other) const {
    return fabs(other - _val) < _delta;
  }
};

using jpair_t = std::pair<std::pair<uint64_t, uint64_t>, double>;

struct jpairdef_t {
  /// Magnetic sites
  std::string site_1, site_2;
  /// Distance criterion
  double distance;
  /// Coupling value to apply if matched
  double J;

  [[nodiscard]] bool match(const std::string& s1, const std::string& s2, double d, double delta = 1e-3) const {
    return ((s1 == site_1 && s2 == site_2) || (s1 == site_2 && s2 == site_1)) && (Approx(distance, delta) == d);
  }
};

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
  [[nodiscard]] double energy(const arma::vec& spins) const;

  static Hamiltonian from_geometry(
      const Geometry& geometry, const std::vector<std::string>& magnetic_sites, std::vector<jpairdef_t> pair_defs);
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
