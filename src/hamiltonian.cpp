#include <mcpp_heisenberg/hamiltonian.hpp>

#include <cassert>
#include <string>
#include <vector>
#include <utility>

namespace mch {

double Hamiltonian::energy(const arma::uvec &spins) const {
  assert(spins.n_rows == _N);

  double energy = .0;

  for (auto& pair : _pairs) {
    energy += pair.second * spins.at(pair.first.first) * spins.at(pair.first.second);
  }

  return -1. * energy;  // assume unique pairs
}

Hamiltonian Hamiltonian::from_geometry(
const Geometry& geometry, const std::vector<std::string>& magnetic_sites, std::vector<jpairdef_t> pair_defs) {
  // prepare geometry
  auto& positions = geometry.positions();
  auto ion_types = std::vector<std::string>();
  uint64_t N = 0;

  for (auto& it : geometry.ions()) {
    // list ion types
    for (uint64_t iatm = 0; iatm < it.second; ++iatm) {
      ion_types.push_back(it.first);
    }

    // count number of magnetic site
    for (auto& s : magnetic_sites) {
      if (s == it.first) {
        N += it.second;
      }
    }
  }

  LOGD << "Got " << N << " magnetic sites";
  LOGD << "create pair list";

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

  LOGD << "Done with list, got " << pairs.size() << " pairs";

  return Hamiltonian(N, pairs);
}

}  // namespace mch
