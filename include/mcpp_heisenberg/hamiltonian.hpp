#ifndef INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
#define INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_

#include <cmath>
#include <vector>
#include <cstdint>
#include <utility>
#include <string>
#include <map>
#include <set>

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

using neighborlist_t = std::set<std::pair<uint64_t, double>>;

/**
 * Ising Hamiltonian
 */
class IsingHamiltonian {
 protected:
  /// Number of magnetic site/spins
  uint64_t _N{0};
  /// (unique) interaction pairs
  std::vector<jpair_t> _pairs;
  /// Neighbour list
  std::map<uint64_t, neighborlist_t> _neighbor_list;

 public:
  IsingHamiltonian() = default;

  /// Create an Heisenberg hamiltonian, with `N` sites and `pairs` pair interactions
  IsingHamiltonian(uint64_t N, const std::vector<jpair_t>& pairs): _N{N}, _pairs(pairs) {
    // build the neighbors list
    for (auto& pair : _pairs) {
      if (!_neighbor_list.contains(pair.first.first)) {
        _neighbor_list[pair.first.first] = neighborlist_t();
      }

      _neighbor_list[pair.first.first].insert({pair.first.second, pair.second});

      if (!_neighbor_list.contains(pair.first.second)) {
        _neighbor_list[pair.first.second] = neighborlist_t();
      }

      _neighbor_list[pair.first.second].insert({pair.first.first, pair.second});
    }
  }

  /// Get number of sites
  uint64_t N() const { return _N; }

  /// Compute the energy
  [[nodiscard]] double energy(const arma::vec& spins) const;

  /// Compute the change in energy due to flip of spin `i`
  [[nodiscard]] double delta_energy(const arma::vec& spins, uint64_t i) const;

  /// Get neighbors of i, a list of `(k, J_ik)`, where `k` is the index of the neighbor and `J_ik` is
  /// The coupling strength.
  neighborlist_t neighbors(uint64_t i) {
    assert(i < _N);

    return _neighbor_list[i];
  }

  /// Create IsingHamiltonian from a geometry
  static IsingHamiltonian from_geometry(const Geometry& geometry, std::vector<jpairdef_t> pair_defs);
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
