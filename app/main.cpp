#include <string>
#include <iostream>
#include <memory>
#include <algorithm>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include "app.h"

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

  // Read input
  mch::app::Parameters simulation_parameters;
  mch::app::read_toml_input(input_file, simulation_parameters);

  // Create H5
  HighFive::File h5_file(output_file, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  // Prepare simulation
  auto simulation = mch::app::prepare_simulation(simulation_parameters, geometry_file);

  mch::app::save_simulation(h5_file, simulation);

  print_title(std::cout, "Run MC");

  mch::elapsed::Chrono chrono_run;
  std::cout << "*!> Run for " << simulation_parameters.N << " steps, using "
            << (simulation_parameters.step_type == mch::app::Sweep ? "sweep" : "cluster") << " updates\n";

  // Create group and buffer
  auto result_group = h5_file.createGroup("results");
  auto pair = mch::app::create_result_datasets(result_group, simulation_parameters, simulation);
  auto dset_aggs = pair.first;
  auto dset_configs = pair.second;

  arma::mat buffer_aggs(2, simulation_parameters.save_interval);
  arma::mat buffer_configs(
      simulation.hamiltonian.number_of_magnetic_sites(), simulation_parameters.save_interval);

  // Run simulation
  uint64_t isave = 0;
  double mean_energy = .0;
  double mean_magnetization = .0;
  uint64_t offset_first_frame = 0;

  for (uint64_t istep = 0; istep < simulation_parameters.N; ++istep) {
    if (simulation_parameters.step_type == mch::app::Sweep) {
      simulation.runner.sweep(simulation_parameters.T, simulation_parameters.H);
    } else {
      simulation.runner.cluster_update(simulation_parameters.T, simulation_parameters.H);
    }

    // write buffer
    if (istep > 0 && istep % simulation_parameters.save_interval == 0) {
      mch::app::write_data_frames(
          dset_aggs, dset_configs,
          offset_first_frame, simulation_parameters.save_interval, simulation.hamiltonian.number_of_magnetic_sites(),
          buffer_aggs, buffer_configs);

      isave++;
      offset_first_frame = isave * simulation_parameters.save_interval;
    }

    buffer_aggs.col(istep % simulation_parameters.save_interval) = {
        simulation.runner.energy(),
        arma::sum(simulation.runner.spins())};

    buffer_configs.col(istep % simulation_parameters.save_interval) = simulation.runner.spins();

    mean_energy += simulation.runner.energy();
    mean_magnetization += fabs(arma::sum(simulation.runner.spins()));
  }

  // write last data frames
  uint64_t remaining_frames = simulation_parameters.N  - offset_first_frame;

  mch::app::write_data_frames(
      dset_aggs, dset_configs,
      offset_first_frame, remaining_frames, simulation.hamiltonian.number_of_magnetic_sites(),
      buffer_aggs, buffer_configs);

  std::cout << "*!> Done running! (took " << chrono_run.format() << ")\n";

  print_title(std::cout, "Epilog");

  // Statistics
  auto dN = static_cast<double>(simulation_parameters.N * simulation.hamiltonian.number_of_magnetic_sites());

  std::cout << "*!> Statistics over " << simulation_parameters.N << " steps "
            << "(N=" << simulation.hamiltonian.number_of_magnetic_sites() << ") :: "
            << "<E>/N = " << mean_energy / dN << ", "
            << "<|m|>/N = " << mean_magnetization / dN << "\n";

  std::cout << "*!> Done (took " << chrono_all.format() << ")\n";

  return EXIT_SUCCESS;
}
