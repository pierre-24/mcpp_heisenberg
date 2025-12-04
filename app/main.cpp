#include <string>
#include <iostream>
#include <memory>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mc.hpp>

#include "app.h"

int main(int argc, char** argv) {
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

  // Set the log level
  mch::app::set_log_level();

  // Read input
  mch::app::Parameters simulation_parameters;
  mch::app::read_toml_input(input_file, simulation_parameters);

  // Prepare simulation
  auto simulation = mch::app::prepare_simulation(simulation_parameters, geometry_file);

  // Run simulation
  double mean_energy = .0;
  double mean_magnetization = .0;

  for (uint64_t istep = 0; istep < simulation_parameters.N; ++istep) {
    simulation.runner.sweep(simulation_parameters.T, simulation_parameters.H);

    mean_energy += simulation.runner.energy();
    mean_magnetization += fabs(arma::sum(simulation.runner.spins()));
  }

  std::cout << "<E> = " << mean_energy / static_cast<double>(simulation_parameters.N * simulation.runner.N()) << ", "
            << "<|m|> = " << mean_magnetization / static_cast<double>(simulation_parameters.N * simulation.runner.N())
            << "\n";

  // Save
  {
    HighFive::File file(output_file, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    auto group = file.createGroup("results");
    simulation.runner.save(group);
  }

  return EXIT_SUCCESS;
}
