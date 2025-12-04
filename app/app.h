#ifndef APP_APP_H_
#define APP_APP_H_

#include <vector>
#include <string>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include <toml++/toml.hpp>

namespace mch::app {

struct Parameters {
  /// Size of the supercell
  std::array<uint64_t, 3> supercell_size = {10, 10, 1};

  /// Magnetic sites
  std::vector<std::string> magnetic_sites = {"H"};

  /// Pair definitions
  std::vector<mch::jpairdef_t> pair_defs = {{"H", "H", 2.0, 1.0}};

  /// Temperature
  double T = 3.0;

  /// Magnetic field
  double H = 0.0;

  /// Number of steps
  uint64_t N = 20000;

  /// Update using TOML
  void update(toml::table& input);

  /// Print
  void print(std::ostream& stream) const;
};

struct Simulation {
  mch::Geometry geometry;
  mch::IsingHamiltonian hamiltonian;
  mch::IsingMonteCarloRunner runner;
};

/// Set log level
void set_log_level(plog::Severity default_ = plog::Severity::info);

/// Read TOML
void read_toml_input(const std::string& input_file, Parameters& parameters);

/// Prepare simulation
Simulation prepare_simulation(const Parameters& parameters, const std::string& geometry_file);

}  // namespace mch::app

#endif  // APP_APP_H_
