#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include "app.h"

namespace mch::app {

void Parameters::update(toml::table& input) {
  // system
  auto supercell_node = input["supercell"];
  if (!!supercell_node) {
    if (!supercell_node.is_array()) {
      throw std::runtime_error("`supercell` must be an array");
    }

    auto& array = *(supercell_node.as_array());
    if (array.size() != 3) {
      throw std::runtime_error("`supercell` must be an array of size 3");
    }

    if (!array.at(0).is_integer() || !array.at(1).is_integer() || !array.at(1).is_integer()) {
      throw std::runtime_error("`supercell` must be an array of integers");
    }

    for (int iarr = 0; iarr < 3; ++iarr) {
      supercell_size[iarr] = array.at(iarr).as_integer()->get();

      if (supercell_size[iarr] < 1) {
        throw std::runtime_error("`supercell` must be an array of strictly positive integers");
      }
    }
  }

  auto ms_node = input["magnetic_sites"];
  if (!!ms_node) {
    if (!ms_node.is_array()) {
      throw std::runtime_error("`magnetic_sites` must be an array");
    }

    auto& array = *(ms_node.as_array());
    for (auto& ms : array) {
      magnetic_sites.push_back(ms.as_string()->get());
    }
  }

  auto js_node = input["pair_defs"];
  if (!!js_node) {
    if (!js_node.is_array()) {
      throw std::runtime_error("`pair_defs` must be an array");
    }

    auto& array = *(js_node.as_array());
    for (auto& md_node  : array) {
      if (!md_node.is_array()) {
        throw std::runtime_error("elements of `pair_defs` must be an array");
      }

      auto& subarray = *(md_node.as_array());

      if (subarray.size() != 4) {
        throw std::runtime_error("elements of `pair_defs` must contain 4 elements");
      }

      if (!subarray.at(0).is_string()
          || !subarray.at(1).is_string()
          || (!subarray.at(2).is_floating_point() && !subarray.at(2).is_integer())
          || (!subarray.at(3).is_floating_point() && !subarray.at(3).is_integer())) {
        throw std::runtime_error("elements of `pair_defs` must be of the form `['A', 'B', d, J]`");
      }

      pair_defs.push_back({
          subarray.at(0).as_string()->get(),
          subarray.at(1).as_string()->get(),
          subarray.at(2).as_floating_point()->get(),
          subarray.at(3).as_floating_point()->get()});
    }
  }

  auto sv_node = input["spin_values"];
  if (!!sv_node) {
    if (!sv_node.is_table()) {
      throw std::runtime_error("`spin_values` must be a table");
    }

    auto& table = *sv_node.as_table();

    table.for_each([this](const toml::key& key, auto&& val){
      if (!val.is_integer() && !val.is_floating_point()) {
        throw std::runtime_error("`spin_values` must contain floats as values");
      }

      spin_values[std::string(key.str())]= val.as_floating_point()->get();
    });
  }

  auto ma_node = input["magnetic_anisotropies"];
  if (!!ma_node) {
    if (!ma_node.is_table()) {
      throw std::runtime_error("`magnetic_anisotropies` must be a table");
    }

    auto& table = *ma_node.as_table();

    table.for_each([this](const toml::key& key, auto&& val){
      if (!val.is_integer() && !val.is_floating_point()) {
        throw std::runtime_error("`magnetic_anisotropies` must contain floats as values");
      }

      magnetic_anisotropies[std::string(key.str())]= val.as_floating_point()->get();
    });
  }

  auto sc_node = input["start_config"];
  if (!!sc_node) {
    std::string m = sc_node.value_or("random");

    if (m == "R" || m == "random") {
      start_config = StartConfig::Random;
    } else if (m == "F" || m == "ferri") {
      start_config = StartConfig::Ferri;
    } else if (m == "D" || m == "ferridown") {
      start_config = StartConfig::FerriDown;
    } else {
      throw std::runtime_error("`start_config` must be either `random` or `ferri` or `ferridown`");
    }
  }

  // simulation
  T = input["T"].value_or(T);
  H = input["H"].value_or(H);
  kB = input["kB"].value_or(kB);
  muB = input["muB"].value_or(muB);
  N = input["N"].value_or(N);

  // results
  save_interval = input["save_interval"].value_or(save_interval);
  deflate_level = input["deflate_level"].value_or(deflate_level);
  chunk_size = input["chunk_size"].value_or(chunk_size);

  if (deflate_level > 9) {
    throw std::runtime_error("`deflate_level` must be lower or equal to 9");
  }
}

void Parameters::print(std::ostream& stream) const {
  stream << "# system\n"
         << "supercell = [" << supercell_size[0] << ", " << supercell_size[1] << ", " << supercell_size[2] << "]\n";

  stream << "magnetic_sites = [";
  for (auto& ms : magnetic_sites) {
    stream << "'" << ms << "', ";
  }

  stream << "]\n";

  stream << "pair_defs = [\n  # A, B, d, J\n";
  for (auto& pair : pair_defs) {
    stream << "  ['" << pair.site_1 << "', '" << pair.site_2 << "', " << pair.distance << ", " << pair.J << "]\n";
  }

  stream << "]\n";

  stream << "spin_values = {\n";
  for (auto& pair : spin_values) {
    stream << "  " << pair.first << " = " << pair.second << ",\n";
  }

  stream << "}\n";

  stream << "magnetic_anisotropies = {\n";
  for (auto& pair : magnetic_anisotropies) {
    stream << "  " << pair.first << " = " << pair.second << ",\n";
  }

  stream << "}\n";

  stream << "start_config = '";

  if (start_config == Random) {
    stream << "random";
  } else if (start_config == Ferri) {
    stream << "ferri";
  } else if (start_config == FerriDown) {
    stream << "ferridown";
  }

  stream << "'\n";

  stream << "# simulation\n"
         << "kB = " << kB << "\n"
         << "muB = " << muB << "\n"
         << "T = " << T << "\n"
         << "H = " << H << "\n"
         << "N = " << N << "\n";

  stream << "# data frames\n"
         << "save_interval = " << save_interval << "\n"
         << "deflate_level = " << deflate_level << "\n"
         << "chunk_size = " << chunk_size << "\n";
}

void set_log_level(plog::Severity default_) {
  char* str_loglevel = getenv("MCH_LOG_LEVEL");
  if (str_loglevel != NULL) {
    auto loglevel = std::string(str_loglevel);
    if (loglevel == "none") {
      mch::init_logging(plog::Severity::none);
    } else if (loglevel == "fatal") {
      mch::init_logging(plog::Severity::fatal);
    } else if (loglevel == "error") {
      mch::init_logging(plog::Severity::error);
    } else if (loglevel == "warning") {
      mch::init_logging(plog::Severity::warning);
    } else if (loglevel == "info") {
      mch::init_logging(plog::Severity::info);
    } else if (loglevel == "debug") {
      mch::init_logging(plog::Severity::debug);
    } else if (loglevel == "verbose") {
      mch::init_logging(plog::Severity::verbose);
    } else {
      mch::init_logging(default_);
    }
  } else {
    mch::init_logging(default_);
  }
}

void read_toml_input(const std::string& input_file, Parameters& parameters) {
  if (!input_file.empty()) {
    std::cout << "*!> Using `" << input_file << "`\n";
    auto new_input_parameters = toml::parse_file(input_file);
    parameters.update(new_input_parameters);
  }

  std::cout << "*!> Input is now::\n```toml\n";
  parameters.print(std::cout);
  std::cout << "```\n";
}

Simulation prepare_simulation(const Parameters& parameters, const std::string& geometry_file) {
  // read geometry
  std::cout << "*!> Using `" << geometry_file << "`\n";

  auto geometry_fs = std::make_shared<std::ifstream>(geometry_file);

  if (!geometry_fs->is_open()) {
    throw std::runtime_error(fmt::format("Could not open geometry file `{}`", geometry_file));
  }

  auto geometry = mch::Geometry::from_poscar(geometry_fs);
  geometry_fs->close();

  LOGI << "Geometry of the unit cell is::\n```\n" << geometry.to_poscar() << "```\n";

  // Make supercell
  std::cout << "*!> Make system\n";
  auto supercell =
      geometry
          .filter_atoms(parameters.magnetic_sites)
          .to_supercell(parameters.supercell_size[0], parameters.supercell_size[1], parameters.supercell_size[2]);

  LOGI << "Geometry of the system is::\n```\n" << supercell.to_poscar() << "```\n";

  // Prepare hamiltonian
  std::cout << "*!> Make Hamiltonian\n";
  auto hamiltonian = mch::IsingHamiltonian::from_geometry(supercell, parameters.pair_defs);

  double magnetic_anisotropy = .0;

  for (auto& iondef : supercell.ions()) {
    if (parameters.magnetic_anisotropies.contains(iondef.first)) {
      double spin_value = 1.0;
      if (parameters.spin_values.contains(iondef.first)) {
        spin_value = parameters.spin_values.at(iondef.first);
      }

      magnetic_anisotropy += parameters.magnetic_anisotropies.at(iondef.first) * pow(spin_value, 2) * iondef.second;
    }
  }

  hamiltonian.set_magnetic_anisotropy(magnetic_anisotropy);

  // Make runner
  std::cout << "*!> Set initial config\n";

  auto spin_values = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::value(1.0));

  uint64_t nx = 0;
  for (auto& iondef : supercell.ions()) {
    if (parameters.spin_values.contains(iondef.first)) {
      spin_values.subvec(nx, nx + iondef.second - 1) *= parameters.spin_values.at(iondef.first);
    }

    nx += iondef.second;
  }

  if (parameters.start_config == StartConfig::Random) {
    spin_values = mch::RandomInitialConfig(spin_values).make();
  } else if (parameters.start_config == StartConfig::Ferri) {
    spin_values = mch::FerriInitialConfig(spin_values).make();
  } else if (parameters.start_config == StartConfig::FerriDown) {
    spin_values = mch::FerriDownInitalConfig(spin_values).make();
  }

  // Make runner
  std::cout << "*!> Make runner\n";

  return {
      .geometry = supercell,
      .hamiltonian = hamiltonian,
      .runner = mch::IsingMonteCarloRunner(hamiltonian, spin_values, parameters.kB)
  };
}

void save_simulation(HighFive::File& file, const Parameters& simulation_parameters, const Simulation& simulation) {
  // geometry itself
  auto geometry_group = file.createGroup("geometry");
  simulation.geometry.to_h5_group(geometry_group);

  // spin values
  auto spin_values = arma::vec(simulation.hamiltonian.number_of_magnetic_sites(), arma::fill::value(1.0));

  uint64_t nx = 0;
  for (auto& iondef : simulation.geometry.ions()) {
    if (simulation_parameters.spin_values.contains(iondef.first)) {
      spin_values.subvec(nx, nx + iondef.second - 1) *= fabs(simulation_parameters.spin_values.at(iondef.first));
    }

    nx += iondef.second;
  }

  geometry_group
      .createDataSet<double>("spin_values", HighFive::DataSpace({simulation.hamiltonian.number_of_magnetic_sites()}))
      .write_raw(spin_values.memptr());

  // magnetic anisotropies
  auto magnetic_anisotropies = arma::vec(simulation.hamiltonian.number_of_magnetic_sites(), arma::fill::value(.0));

  nx = 0;
  for (auto& iondef : simulation.geometry.ions()) {
    if (simulation_parameters.magnetic_anisotropies.contains(iondef.first)) {
      magnetic_anisotropies.subvec(nx, nx + iondef.second - 1).fill(
          simulation_parameters.magnetic_anisotropies.at(iondef.first));
    }

    nx += iondef.second;
  }

  geometry_group
      .createDataSet<double>(
          "magnetic_anisotropies", HighFive::DataSpace({simulation.hamiltonian.number_of_magnetic_sites()}))
      .write_raw(magnetic_anisotropies.memptr());

  // hamiltonian
  auto hamiltonian_group = file.createGroup("hamiltonian");
  simulation.hamiltonian.to_h5_group(hamiltonian_group);

  // infos
  std::array<double, 4> info = {
      simulation_parameters.kB, simulation_parameters.T, simulation_parameters.muB, simulation_parameters.H};
  hamiltonian_group.createDataSet("parameters", info).write(info);
}

std::pair<HighFive::DataSet, HighFive::DataSet> create_result_datasets(
HighFive::Group& result_group,
const Parameters& simulation_parameters, const Simulation& simulation) {
  LOGI << "Create result datasets";

  // aggs
  auto dset_aggs = result_group.createDataSet<double>(
      "aggregated_data", HighFive::DataSpace({simulation_parameters.N, 2}));

  // configs
  HighFive::DataSetCreateProps dapl_configs;

  dapl_configs.add(HighFive::Chunking(std::vector<hsize_t>{
      std::min(simulation_parameters.chunk_size, simulation_parameters.N),
      std::min(simulation_parameters.chunk_size, simulation.hamiltonian.number_of_magnetic_sites())}));

  dapl_configs.add(HighFive::Deflate{simulation_parameters.deflate_level});

  auto dset_configs = result_group.createDataSet<int8_t>(
      "configs",
      HighFive::DataSpace({simulation_parameters.N, simulation.hamiltonian.number_of_magnetic_sites()}),
      dapl_configs);

  return {dset_aggs, dset_configs};
}

void write_data_frames(
HighFive::DataSet& dset_aggs, HighFive::DataSet& dset_configs,
uint64_t offset, uint64_t size, uint64_t N,
const arma::mat& buffer_aggs, const arma::mat& buffer_configs) {
  LOGI << "Write frames [" << offset << "," << offset + size << ")";

  // write energies
  dset_aggs.select({offset, 0}, {size, 2})
      .write_raw(buffer_aggs.memptr());

  // write spins after transformation
  std::vector<std::vector<int8_t>> configs;
  buffer_configs.each_col([&configs, &N](auto& col) {
    std::vector<int8_t> config(N);
    std::transform(
        col.cbegin(), col.cend(), config.begin(), [](auto& val) { return (val < 0 ? -1 : 1); });

    configs.push_back(config);
  });

  dset_configs
      .select({offset, 0}, {size, N})
      .write(configs);
}

}  // namespace mch::app
