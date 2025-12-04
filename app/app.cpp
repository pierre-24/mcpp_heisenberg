#include <string>
#include <iostream>
#include <memory>

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

  // simulation
  T = input["T"].value_or(T);
  H = input["H"].value_or(H);
  N = input["N"].value_or(N);

  auto st_node = input["step_type"];
  if (!!st_node) {
    std::string m = st_node.value_or("sweep");

    if (m == "S" || m == "sweep") {
      step_type = Sweep;
    } else if (m == "C" || m == "cluster") {
      step_type = Cluster;
    } else {
      throw std::runtime_error("`step_type` must be either `sweep` or `cluster`");
    }
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

  stream << "# simulation\n"
         << "T = " << T << "\n"
         << "H = " << H << "\n"
         << "N = " << N << "\n"
         << "step_type = '" << (step_type == Sweep ? "sweep" : "cluster") << "'\n";
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

  std::cout << "*!> Geometry is::\n```\n" << geometry.to_poscar() << "```\n";

  // Make supercell
  std::cout << "*!> Make system\n";
  auto supercell =
      geometry
          .filter_atoms(parameters.magnetic_sites)
          .to_supercell(parameters.supercell_size[0], parameters.supercell_size[1], parameters.supercell_size[2]);

  std::cout << "*!> Geometry of the system is::\n```\n" << supercell.to_poscar() << "```\n";

  // Prepare hamiltonian
  std::cout << "*!> Make Hamiltonian\n";
  auto hamiltonian = mch::IsingHamiltonian::from_geometry(supercell, parameters.pair_defs);

  // Make runner
  std::cout << "*!> Make runner\n";
  return {
      .geometry = supercell,
      .hamiltonian = hamiltonian,
      .runner = mch::IsingMonteCarloRunner(hamiltonian)
  };
}

}  // namespace mch::app
