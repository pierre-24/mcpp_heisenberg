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
  /// Value of Boltzmann constant
  double _kB {1.0};
  /// Value of muB
  double _muB{1.0};
  /// current energy
  double _energy{.0};
  /// First
  bool _initialized_energy{false};
  /// Generator
  std::mt19937 _rng;

 public:
  explicit MonteCarloRunner(double kB = 1.0, double muB = 1.0): _kB{kB}, _muB{muB} {
    // initialize random generator
    std::random_device rd;
    _rng = std::mt19937(rd());
  }

  virtual ~MonteCarloRunner() = default;

  /// Sweep over all spins (at a given temperature `T` and a given magnetic field `H`) and switch them if any
  virtual void sweep(double T, double H = .0)  = 0;

  /// reset energy, must be done if `H` or `T` change
  virtual void reset_energy(double T, double H = .0) = 0;

  /// Get current energy
  [[nodiscard]] double energy() const { return _energy; }

  /// Get current spin configuration
  [[nodiscard]] virtual const arma::vec& spins() const = 0;
};

class IsingMonteCarloRunner: public MonteCarloRunner{
 protected:
  /// Ising Hamiltonian
  IsingHamiltonian _hamiltonian;
  /// Current spin configuration
  arma::vec _spins;

 public:
  IsingMonteCarloRunner() = default;

  explicit IsingMonteCarloRunner(
      const IsingHamiltonian& hamiltonian, const arma::vec& initial, double kB = 1.0, double muB = 1.0)
      : MonteCarloRunner(kB, muB), _hamiltonian{hamiltonian} {
    assert(initial.n_rows == hamiltonian.number_of_magnetic_sites());

    LOGD << "initial config is " << initial << ", sum = " << arma::sum(initial);

    // choose a starting config
    set_spins(initial);
  }

  /// Get current spin configuration
  [[nodiscard]] const arma::vec& spins() const override { return _spins; }

  /// Set spin configuration
  void set_spins(const arma::vec& config) {
    assert(config.n_rows == _hamiltonian.number_of_magnetic_sites());
    _spins = config;
  }

  void sweep(double T, double H = .0) override;

  void reset_energy(double T, double H = .0) override {
    LOGD << "reset energy for T=" << T << ", H=" << H;

    this->_energy = _hamiltonian.energy(_spins, _muB * H);
    this->_initialized_energy = true;
  }
};

/// Account for {S_max, S_max - 1, ..., -S_max + 1, -S_max}
class QuantumIsingMonteCarloRunner : public IsingMonteCarloRunner{
 protected:
  std::vector<std::vector<double>> _possible_values;
 public:
  QuantumIsingMonteCarloRunner(
      const IsingHamiltonian& hamiltonian,
      const arma::vec& initial, const arma::vec& max_spins, double kB = 1.0, double muB = 1.0)
  : IsingMonteCarloRunner(hamiltonian, initial, kB, muB) {
    assert(max_spins.n_rows == initial.n_rows);

    for (auto& ms : max_spins) {
      uint64_t n = static_cast<uint64_t>(round(ms * 2));
      std::vector<double> possibilities(n + 1);
      for (uint64_t iv = 0; iv <= n; ++iv) {
        possibilities.at(iv) = ms - static_cast<double>(n - iv);
      }

      _possible_values.push_back(possibilities);
    }

    LOGV << "Spin values = " << _possible_values;
  }

  void sweep(double T, double H = .0) override;
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_MC_HPP_
