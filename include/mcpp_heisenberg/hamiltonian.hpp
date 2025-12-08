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
  uint64_t _n_magnetic_sites{0};
  /// (unique) interaction pairs
  std::vector<jpair_t> _pairs;
  /// Neighbour list
  std::map<uint64_t, neighborlist_t> _neighbor_list;
  /// Magnetic anisotropy (constant term)
  double _magnetic_anisotropy{.0};

 public:
  IsingHamiltonian() = default;

  /// Create an Heisenberg hamiltonian, with `N` sites and `pairs` pair interactions
  IsingHamiltonian(uint64_t N, const std::vector<jpair_t>& pairs): _n_magnetic_sites{N}, _pairs(pairs) {
    // build the neighbor list
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
  [[nodiscard]] uint64_t number_of_magnetic_sites() const { return _n_magnetic_sites; }

  /// Get number of pairs
  [[nodiscard]] uint64_t number_of_pairs() const { return _pairs.size(); }

  /// Set magnetic anisotropy
  void set_magnetic_anisotropy(double m) { _magnetic_anisotropy = m; }

  /// Compute the energy
  [[nodiscard]] double energy(const arma::vec& spins, double muBH = .0) const;

  /// Compute the change in energy due to flip of spin `i`
  [[nodiscard]] double delta_energy(const arma::vec& spins, uint64_t i, double muBH = .0) const;

  /// Compute (ΔE_i, P_i) associated to flipping spin i
  [[nodiscard]] std::pair<double, double> P_i(const arma::vec& spins, uint64_t i, double kBT, double muBH = .0) const;

  /// Get neighbors of i, a list of `(k, J_ik)`, where `k` is the index of the neighbor and `J_ik` is
  /// the coupling strength between these two.
  [[nodiscard]] neighborlist_t neighbors(uint64_t i) const {
    assert(i < _n_magnetic_sites);

    return _neighbor_list.at(i);
  }

  /// Create IsingHamiltonian from a geometry
  static IsingHamiltonian from_geometry(const Geometry& geometry, std::vector<jpairdef_t> pair_defs);

  /// Save in H5 group
  void to_h5_group(HighFive::Group& group) const;
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_HAMILTONIAN_HPP_
