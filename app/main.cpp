#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include "mcpp_heisenberg/app.h"

void print_title(std::ostream & stream, const std::string& title) {
  stream << "\n  _";
  for (uint64_t i = 0; i < title.size(); ++i) {
    stream << "_";
  }

  stream << "_\n | " << title << " |\n `-";
  for (uint64_t i = 0; i < title.size(); ++i) {
    stream << "-";
  }

  stream << "-´\n\n";
}

int main(int argc, char** argv) {
  std::cout << "*!> Welcome, this is "
            << PROJECT_NAME
            << " (" << PROJECT_VERSION << " @ " << PROJECT_BUILD_COMMIT << ", built " << PROJECT_BUILD_DATE << ")\n";

  // read CLI
  CLI::App app;

  std::string input_file;
  std::string geometry_file;
  std::string output_file = "MC.h5";

  app.add_option("-i,--input", input_file, "Input (TOML) file")
      ->check(CLI::ExistingFile);
  app.add_option("-g,--geometry", geometry_file, "Geometry (POSCAR) file")
      ->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", output_file, "Output (H5) file")
      ->capture_default_str();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  mch::elapsed::Chrono chrono_all;

  // Set the log level
  mch::app::set_log_level();

  print_title(std::cout, "Read inputs");

  // Read parameters
  mch::app::Parameters simulation_parameters;
  mch::app::read_toml_input(input_file, simulation_parameters);

  std::cout << "*!> Input is now::\n```toml\n";
  simulation_parameters.print(std::cout);
  std::cout << "```\n";

  // Read geometry
  std::cout << "!*> Using `" << geometry_file << "`\n";
  auto geometry_fs = std::make_shared<std::ifstream>(geometry_file);

  if (!geometry_fs->is_open()) {
    std::cerr << fmt::format("Could not open geometry file `{}`", geometry_file) << "\n";
    return EXIT_FAILURE;
  }

  auto initial_geometry = mch::Geometry::from_poscar(geometry_fs);
  geometry_fs->close();

  std::cout << "Geometry of the unit cell is::\n```\n" << initial_geometry.to_poscar() << "```\n";

  // Create H5 output
  std::cout << "!*> Using `" << output_file << "`\n";
  auto h5_file = HighFive::File(
      output_file, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  // Prepare simulation
  print_title(std::cout, "Prepare");
  auto simulation = mch::app::Runner(simulation_parameters, initial_geometry, std::move(h5_file));

  // run simulation
  print_title(std::cout, "Run MC");
  simulation.run(simulation_parameters);

  // The end ;)
  print_title(std::cout, "Epilog");
  std::cout << "*!> Done (took " << chrono_all.format() << ")\n";

  return EXIT_SUCCESS;
}
