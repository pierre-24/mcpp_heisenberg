#ifndef INCLUDE_MCPP_HEISENBERG_MC_HPP_
#define INCLUDE_MCPP_HEISENBERG_MC_HPP_

#include <random>
#include <vector>
#include <utility>

#include <highfive/highfive.hpp>

#include <mcpp_heisenberg/hamiltonian.hpp>

namespace mch {

class IsingMonteCarloRunner {
 protected:
  /// Value of Boltzmann constant
  double _kB {1.0};
  /// Value of muB
  double _muB{1.0};
  /// Ising Hamiltonian
  IsingHamiltonian _hamiltonian;
  /// Spin configuration
  arma::vec _spins;
  /// energy associated to `_spin`
  double _energy{.0};
  /// Generator
  std::mt19937 _rng;
  /// First
  bool _initialized_energy{false};

 public:
  IsingMonteCarloRunner() = default;

  explicit IsingMonteCarloRunner(
      const IsingHamiltonian& hamiltonian, const arma::vec& initial, double kB = 1.0, double muB = 1.0)
      : _kB{kB}, _muB{muB}, _hamiltonian{hamiltonian} {
    // initialize random generator
    std::random_device rd;
    _rng = std::mt19937(rd());

    LOGD << "initial config is " << initial << ", sum = " << arma::sum(initial);

    // choose a starting config
    set_spins(initial);
  }

  /// Get current spin configuration
  [[nodiscard]] const arma::vec& spins() const { return _spins; }

  /// Set spin configuration
  void set_spins(const arma::vec& config) {
    assert(config.n_rows == _hamiltonian.number_of_magnetic_sites());
    _spins = config;
  }

  /// Get current energy
  [[nodiscard]] double energy() const { return _energy; }

  /// Sweep over all spins (at a given temperature `T` and a given magnetic field `H`) and switch them if any
  void sweep(double T, double H = .0);

  /// reset energy, must be done if `H` or `T` change
  void reset_energy(double T, double H = .0) {
    LOGD << "reset energy for T=" << T << ", H=" << H;

    _energy = _hamiltonian.energy(_spins, _muB * H);
    _initialized_energy = true;
  }
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_MC_HPP_
