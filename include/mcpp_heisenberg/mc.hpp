#ifndef INCLUDE_MCPP_HEISENBERG_MC_HPP_
#define INCLUDE_MCPP_HEISENBERG_MC_HPP_

#include <random>
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

 public:
  explicit MonteCarloRunner(const Hamiltonian& hamiltonian): _hamiltonian{hamiltonian} {
    _spins.resize(hamiltonian.N());

    // initialize random generator
    std::random_device rd;
    _rng = std::mt19937(rd());

    // choose a random configuration
    auto config = arma::randi<arma::vec>(hamiltonian.N(), arma::distr_param(0, 1));
    config.for_each([](double& val) { val = (val == 0) ? -1. : val; });
    set_spins(config);
  }

  /// Get current spin configuration
  [[nodiscard]] const arma::vec& spins() const { return _spins; }

  /// Set spin configuration
  void set_spins(const arma::vec& config) {
    assert(config.n_rows == _hamiltonian.N());
    _spins = config;
    _energy = _hamiltonian.energy(_spins);
  }

  /// Get current energy
  [[nodiscard]] double energy() const { return _energy; }

  /// Sweep over all spins (at a given temperature `T` and a given magnetic field `H`) and switch them if any
  void sweep(double T, double H = .0);
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_MC_HPP_
