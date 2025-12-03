#include <mcpp_heisenberg/mc.hpp>

#include <vector>

namespace mch {

void MonteCarloRunner::sweep(double T, double H) {
  assert(T > 0);

  LOGV << "Sweep at T=" << T << ", H=" << H;

  std::uniform_real_distribution<> dis(.0, 1.);

  for (uint64_t ispin = 0; ispin < _hamiltonian.N(); ++ispin) {
    double dE = _hamiltonian.delta_energy(_spins, ispin);

    if ((dE < 0) || (dis(_rng) < exp(-dE / T))) {
      _spins.at(ispin) *= -1;
      _energy += dE;
    }
  }

  _save_to_frame();
}

void MonteCarloRunner::save(HighFive::Group& group) {
  LOGD << "Save MC";

  std::vector<uint64_t> info = {_hamiltonian.N(), _frames.size()};
  group.createDataSet<uint64_t>("info", HighFive::DataSpace::From(info)).write(info);

  arma::vec energies(_frames.size());
  arma::mat configs(_hamiltonian.N(), _frames.size());  // HDF5 is row major

  uint64_t i = 0;
  for (auto& frame : _frames) {
    energies.at(i) = frame.first;
    configs.col(i) = frame.second;
    i++;
  }

  auto dset_configs = group.createDataSet<double>("configs", HighFive::DataSpace({ _frames.size(), _hamiltonian.N()}));
  dset_configs.write_raw(configs.memptr());

  auto dset_energies = group.createDataSet<double>("energies", HighFive::DataSpace({energies.n_rows}));
  dset_energies.write_raw(energies.memptr());
}

}  // namespace mch
