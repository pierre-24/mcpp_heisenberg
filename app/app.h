#ifndef APP_APP_H_
#define APP_APP_H_

#include <vector>
#include <string>
#include <utility>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include <toml++/toml.hpp>

namespace mch::app {

/// Starting configuration
enum StartConfig {
  /// Ferrimagnetic (all up) config
  Ferri,
  /// All down config
  FerriDown,
  /// Random config
  Random
};

struct Parameters {
  /// Size of the supercell
  std::array<uint64_t, 3> supercell_size = {1, 1, 1};

  /// Magnetic sites
  std::vector<std::string> magnetic_sites;

  /// Pair definitions
  std::vector<mch::jpairdef_t> pair_defs;

  /// Starting confi
  StartConfig start_config = Ferri;

  /// Boltzmann constant
  double kB = 1.0;

  /// Bohr magneton
  double muB = 1.0;

  /// Temperature
  double T = 0.1;

  /// Magnetic field
  double H = 0.0;

  /// Number of steps
  uint64_t N = 1000;

  /// Save interval
  uint64_t save_interval = 1000;

  /// Deflate level
  uint32_t deflate_level = 6;

  /// Chunk size
  uint64_t chunk_size = 1024;

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

/// Create result dataset
std::pair<HighFive::DataSet, HighFive::DataSet> create_result_datasets(
HighFive::Group& result_group,
const Parameters& simulation_parameters, const Simulation& simulation);

/// Write data frames
void write_data_frames(
HighFive::DataSet& dset_aggs, HighFive::DataSet& dset_configs,
uint64_t offset, uint64_t size, uint64_t N,
const arma::mat& buffer_aggs, const arma::mat& buffer_configs);

}  // namespace mch::app

#endif  // APP_APP_H_
