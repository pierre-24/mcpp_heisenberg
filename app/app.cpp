#include <string>
#include <iostream>
#include <memory>

#include "app.h"

namespace mch::app {

void Parameters::update(toml::table& input) {
  T = input["T"].value_or(T);
  H = input["H"].value_or(H);
  N = input["N"].value_or(N);
}

void Parameters::print(std::ostream& stream) const {
  stream << "T = " << T << "\n"
         << "H = " << H << "\n"
         << "N = " << N << "\n";
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
  std::cout << "*!> Make supercell\n";
  auto supercell = geometry.to_supercell(
      parameters.supercell_size[0], parameters.supercell_size[1], parameters.supercell_size[2]);

  // Prepare hamiltonian
  std::cout << "*!> Make Hamiltonian\n";
  auto hamiltonian = mch::IsingHamiltonian::from_geometry(supercell, parameters.magnetic_sites, parameters.pair_defs);

  // Make runner
  std::cout << "*!> Make runner\n";
  return {
      .geometry = supercell,
      .hamiltonian = hamiltonian,
      .runner = mch::IsingMonteCarloRunner(hamiltonian)
  };
}

}  // namespace mch::app
