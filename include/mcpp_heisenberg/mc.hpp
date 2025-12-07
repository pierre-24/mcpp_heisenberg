#ifndef INCLUDE_MCPP_HEISENBERG_MC_HPP_
#define INCLUDE_MCPP_HEISENBERG_MC_HPP_

#include <random>
#include <vector>
#include <utility>

#include <highfive/highfive.hpp>

#include <mcpp_heisenberg/hamiltonian.hpp>

namespace mch {

class InitialConfig {
 protected:
  /// Number of sites
  uint64_t _number_of_sites{0};
 public:
  InitialConfig() = delete;

  explicit InitialConfig(uint64_t n): _number_of_sites{n} {}

  virtual arma::vec make() const = 0;
};

/// All spin ups
class FerriInitialConfig : public InitialConfig{
 protected:
  arma::vec _spin_values;
 public:
  FerriInitialConfig() = delete;

  explicit FerriInitialConfig(const arma::vec& spin_values)
      : InitialConfig(spin_values.n_rows), _spin_values{spin_values} {}

  [[nodiscard]] arma::vec make() const override {
    return _spin_values;
  }
};

/// Random configuration
class RandomInitialConfig: public FerriInitialConfig {
 public:
  RandomInitialConfig() = delete;

  explicit RandomInitialConfig(const arma::vec& spin_values) : FerriInitialConfig(spin_values) {}

  [[nodiscard]] arma::vec make() const override {
    auto config = FerriInitialConfig::make();

    // initialize random generator
    std::random_device rd;
    auto rng = std::mt19937(rd());
    std::bernoulli_distribution dis(0.5);

    config.for_each([&dis, &rng](auto& val) {
      if (dis(rng)) {
        val *= -1;
      }
    });

    return config;
  }
};

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

  /// Frames
  std::vector<std::pair<double, arma::vec>> _frames;

 public:
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
