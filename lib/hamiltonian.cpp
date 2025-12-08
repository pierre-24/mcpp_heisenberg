#include <mcpp_heisenberg/hamiltonian.hpp>

#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

namespace mch {

double IsingHamiltonian::energy(const arma::vec& spins, double muBH) const {
  assert(spins.n_rows == _n_magnetic_sites);

  double energy = .0;

  for (auto& pair : _pairs) {
    energy += pair.second * spins.at(pair.first.first) * spins.at(pair.first.second);
  }

  return -1. * energy - muBH * arma::sum(spins) - _magnetic_anisotropy;  // assume unique pairs
}

double IsingHamiltonian::delta_energy(const arma::vec& spins, uint64_t i, double muBH) const {
  assert(i < _n_magnetic_sites);

  double dE = .0;

  for (auto& neighbor : _neighbor_list.at(i)) {
    dE += neighbor.second * spins.at(neighbor.first);
  }

  return 2 * spins.at(i) * (dE + muBH);
}

std::pair<double, double> IsingHamiltonian::P_i(const arma::vec& spins, uint64_t i, double kBT, double muBH) const {
  double dE = delta_energy(spins, i, muBH);
  return {dE, dE <= 0 ? 1 : exp(-dE / kBT)};
}

IsingHamiltonian IsingHamiltonian::from_geometry(const Geometry& geometry, std::vector<jpairdef_t> pair_defs) {
  elapsed::Chrono chrono;
  LOGI << "*> Create pair list from geometry (" << geometry.number_of_atoms() << " magnetic sites)";

  // prepare geometry
  auto& positions = geometry.positions();

  // list ion types
  auto ion_types = std::vector<std::string>();
  for (auto& it : geometry.ions()) {
    for (uint64_t iatm = 0; iatm < it.second; ++iatm) {
      ion_types.push_back(it.first);
    }
  }

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

  return {geometry.number_of_atoms(), pairs};
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
}

}  // namespace mch
