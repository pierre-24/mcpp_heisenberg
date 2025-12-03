#ifndef INCLUDE_MCPP_HEISENBERG_MC_HPP_
#define INCLUDE_MCPP_HEISENBERG_MC_HPP_

#include <random>
#include <vector>
#include <utility>

#include <highfive/highfive.hpp>

#include <mcpp_heisenberg/hamiltonian.hpp>

namespace mch {

class MonteCarloRunner {
 protected:
  /// Hamiltonian
  Hamiltonian _hamiltonian;
  /// Spin configuration
  arma::vec _spins;
  /// energy associated to `_spin`
  double _energy;
  /// Generator
  std::mt19937 _rng;

  /// Frames
  std::vector<std::pair<double, arma::vec>> _frames;

  /// Save current config in frame
  void _save_to_frame() { _frames.push_back(std::make_pair(_energy, _spins)); }

 public:
  explicit MonteCarloRunner(const Hamiltonian& hamiltonian): _hamiltonian{hamiltonian} {
    _spins.resize(hamiltonian.N());

    // initialize random generator
    std::random_device rd;
    _rng = std::mt19937(rd());

    // choose a random configuration
    auto spins = arma::vec(hamiltonian.N(), arma::fill::randu);
    spins.for_each([](double& e) { e = e >= .5 ? 1.0 : -1.0; });
    set_spins(spins);
  }

  /// Get current spin configuration
  [[nodiscard]] const arma::vec& spins() const { return _spins; }

  /// Set spin configuration
  void set_spins(const arma::vec& config) {
    assert(config.n_rows == _hamiltonian.N());
    _spins = config;
    _energy = _hamiltonian.energy(_spins);

    _save_to_frame();
  }

  /// Get current energy
  [[nodiscard]] double energy() const { return _energy; }

  /// Get number of magnetic sites
  [[nodiscard]] uint64_t N() { return _hamiltonian.N(); }

  /// Sweep over all spins (at a given temperature `T` and a given magnetic field `H`) and switch them if any
  void sweep(double T, double H = .0);

  /// do a cluster update
  void cluster_update(double T, double H = .0);

  /// Save results in `group`
  void save(HighFive::Group& group);
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_MC_HPP_
