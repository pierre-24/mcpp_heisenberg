#ifndef APP_APP_H_
#define APP_APP_H_

#include <vector>
#include <string>
#include <utility>
#include <map>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include <toml++/toml.hpp>

namespace mch::app {

struct Parameters {
  /// Size of the supercell
  std::array<uint64_t, 3> supercell_size = {1, 1, 1};

  /// Magnetic sites
  std::vector<std::string> magnetic_sites;

  /// Pair definitions
  std::vector<mch::jpairdef_t> pair_defs;

  /// Spin values
  std::map<std::string, double> spin_values;

  /// Magnetic anisotropies
  std::map<std::string, double> magnetic_anisotropies;

  /// Initial configs
  std::map<std::string, ConfigType> initial_configs;

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

/// Runner for the simulation
class Runner {
 protected:
  mch::Geometry _geometry;
  mch::IsingHamiltonian _hamiltonian;
  mch::IsingMonteCarloRunner _runner;
  HighFive::File _h5_file;

  /// Write data frames in H5 datasets
  static void _write_data_frames(
     HighFive::DataSet& dset_aggs, HighFive::DataSet& dset_configs,
     uint64_t offset, uint64_t size, uint64_t N,
     const arma::mat& buffer_aggs, const arma::mat& buffer_configs);

 public:
  Runner() = delete;
  explicit Runner(const Parameters& parameters, const Geometry& initial_geometry, HighFive::File&& output_file);

  /// Run the simulation
  void run(const Parameters& parameters);
};

/// Set log level
void set_log_level(plog::Severity default_ = plog::Severity::info);

/// Read TOML
void read_toml_input(const std::string& input_file, Parameters& parameters);

}  // namespace mch::app

#endif  // APP_APP_H_
