#ifndef APP_APP_H_
#define APP_APP_H_

#include <vector>
#include <string>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include <toml++/toml.hpp>

namespace mch::app {

/// Type of update at each MC step
enum MCStepType {
  /// Sweep through all spins and update one by one
  Sweep,

  /// Cluster update
  Cluster
};

struct Parameters {
  /// Size of the supercell
  std::array<uint64_t, 3> supercell_size = {1, 1, 1};

  /// Magnetic sites
  std::vector<std::string> magnetic_sites;

  /// Pair definitions
  std::vector<mch::jpairdef_t> pair_defs;

  /// Temperature
  double T = 0.1;

  /// Magnetic field
  double H = 0.0;

  /// Number of steps
  uint64_t N = 1000;

  /// Type of update
  MCStepType step_type = Sweep;

  /// Save interval
  uint64_t save_interval = 1000;

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

/// Save simulation to H5
void save_simulation(HighFive::File& file, const Simulation& simulation);

}  // namespace mch::app

#endif  // APP_APP_H_
