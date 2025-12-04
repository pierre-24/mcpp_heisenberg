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
  /// IsingHamiltonian
  IsingHamiltonian _hamiltonian;
  /// Spin configuration
  arma::vec _spins;
  /// energy associated to `_spin`
  double _energy;
  /// Generator
  std::mt19937 _rng;

  /// Frames
  std::vector<std::pair<double, arma::vec>> _frames;

 public:
  explicit IsingMonteCarloRunner(const IsingHamiltonian& hamiltonian): _hamiltonian{hamiltonian} {
    _spins.resize(hamiltonian.number_of_magnetic_sites());

    // initialize random generator
    std::random_device rd;
    _rng = std::mt19937(rd());

    // choose a random configuration
    auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::randu);
    spins.for_each([](double& e) { e = e >= .5 ? 1.0 : -1.0; });
    set_spins(spins);
  }

  /// Get current spin configuration
  [[nodiscard]] const arma::vec& spins() const { return _spins; }

  /// Set spin configuration
  void set_spins(const arma::vec& config) {
    assert(config.n_rows == _hamiltonian.number_of_magnetic_sites());
    _spins = config;
    _energy = _hamiltonian.energy(_spins);
  }

  /// Get current energy
  [[nodiscard]] double energy() const { return _energy; }

  /// Sweep over all spins (at a given temperature `T` and a given magnetic field `H`) and switch them if any
  void sweep(double T, double H = .0);

  /// do a cluster update
  uint64_t cluster_update(double T, double H = .0);
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_MC_HPP_
