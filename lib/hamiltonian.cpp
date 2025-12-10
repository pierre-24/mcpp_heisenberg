#include <mcpp_heisenberg/hamiltonian.hpp>

#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>

namespace mch {

double IsingHamiltonian::energy(const arma::vec& spins, double muBH) const {
  assert(spins.n_rows == _n_magnetic_sites);

  double energy = .0;

  for (auto& pair : _pairs) {
    energy += pair.second * spins.at(pair.first.first) * spins.at(pair.first.second);
  }

  return -1. * energy - muBH * arma::sum(spins) - arma::dot(_magnetic_anisotropies, arma::pow(spins, 2));
}

double IsingHamiltonian::delta_energy(const arma::vec& spins, uint64_t i, double target_value, double muBH) const {
  assert(i < _n_magnetic_sites);

  double dE = .0;

  if (_neighbor_list.contains(i)) {
    for (auto& neighbor : _neighbor_list.at(i)) {
      dE += neighbor.second * spins.at(neighbor.first);
    }
  }

  return (spins.at(i) - target_value) * (dE + muBH) +
         (pow(spins.at(i), 2) - pow(target_value, 2)) * _magnetic_anisotropies.at(i);
}

std::pair<double, double> IsingHamiltonian::P_i(
    const arma::vec& spins, uint64_t i, double kBT, double target_value, double muBH) const {
  double dE = delta_energy(spins, i, target_value, muBH);
  return {dE, dE <= 0 ? 1 : exp(-dE / kBT)};
}

IsingHamiltonian IsingHamiltonian::from_geometry(
    const Geometry& geometry, std::vector<jpairdef_t> pair_defs,
    const std::map<std::string, double>& magnetic_anisotropies) {
  elapsed::Chrono chrono;
  LOGI << "*> Create pair list from geometry (" << geometry.number_of_atoms() << " magnetic sites)";

  // prepare geometry
  auto& positions = geometry.positions();
  auto anisotropies = arma::vec(geometry.number_of_atoms(), arma::fill::value(.0));

  auto nx = 0;
  auto ion_types = std::vector<std::string>();
  for (auto& it : geometry.ions()) {
    for (uint64_t iatm = 0; iatm < it.second; ++iatm) {
      // list ion types
      ion_types.push_back(it.first);

      // make anisotropies
      if (magnetic_anisotropies.contains(it.first)) {
        anisotropies.subvec(nx, nx + it.second - 1).fill(magnetic_anisotropies.at(it.first));
      }
    }

    nx += it.second;
  }

  // make pairs
  std::vector<jpair_t> pairs;

  for (uint64_t iatm = 0; iatm < geometry.number_of_atoms(); ++iatm) {
    for (uint64_t jatm = iatm + 1; jatm < geometry.number_of_atoms(); ++jatm) {
      arma::vec diff = positions.row(jatm).t() - positions.row(iatm).t();

      // use PBC
      for (int kc = 0; kc < 3; ++kc) {
        if (diff.at(kc) < -.5) {
          diff.at(kc) += 1;
        } else if (diff.at(kc) > .5) {
          diff.at(kc) -= 1;
        }
      }

      // get norm using lattice vectors
      arma::vec tdiff(3);
      for (int kc = 0; kc < 3; ++kc) {
        tdiff += diff.at(kc) * geometry.lattice().row(kc).t();
      }

      double norm = arma::norm(tdiff);

      LOGV << iatm << "," << jatm << ":: " << norm;
      for (auto& criterion : pair_defs) {
        if (criterion.match(ion_types[iatm], ion_types[jatm], norm)) {  // TODO(pierre): delta
          pairs.push_back(std::make_pair(std::make_pair(iatm, jatm), criterion.J));
        }
      }
    }
  }

  LOGI << "<* Done with list, got " << pairs.size() << " pair interactions (took " << chrono.format() << ")";

  return { geometry.number_of_atoms(), pairs, anisotropies };
}

void IsingHamiltonian::to_h5_group(HighFive::Group& group) const {
  // save pairs
  std::vector<std::array<uint64_t, 2>> pairs(_pairs.size());
  std::transform(
      _pairs.cbegin(), _pairs.cend(), pairs.begin(),
      [](auto& t) { return std::array<uint64_t, 2>({t.first.first, t.first.second}); });

  auto dset_pairs = group.createDataSet("pairs", pairs);
  dset_pairs.write(pairs);

  // save J
  std::vector<double> Js(_pairs.size());
  std::transform(_pairs.cbegin(), _pairs.cend(), Js.begin(), [](auto& t) { return t.second; });

  auto dset_Js = group.createDataSet("J", Js);
  dset_Js.write(Js);

  // Save magnetic anisotropies
  group.createDataSet<double>(
      "magnetic_anisotropies", HighFive::DataSpace({_n_magnetic_sites}))
      .write_raw(_magnetic_anisotropies.memptr());
}

}  // namespace mch
