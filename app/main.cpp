#include <string>
#include <iostream>
#include <memory>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include "app.h"

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

  // Set the log level
  mch::app::set_log_level();

  // Read input
  mch::app::Parameters simulation_parameters;
  mch::app::read_toml_input(input_file, simulation_parameters);

  // Create H5
  HighFive::File h5_file(output_file, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  // Prepare simulation
  auto simulation = mch::app::prepare_simulation(simulation_parameters, geometry_file);

  mch::app::save_simulation(h5_file, simulation_parameters, simulation);

  std::cout << "*!> Run for " << simulation_parameters.N << " steps\n";

  // Create group and buffer
  auto result_group = h5_file.createGroup("results");
  auto dset_energies = result_group.createDataSet<double>("energies", HighFive::DataSpace({simulation_parameters.N}));
  auto dset_configs = result_group.createDataSet<double>(
      "configs", HighFive::DataSpace({simulation_parameters.N, simulation.hamiltonian.number_of_magnetic_sites()}));

  arma::vec buffer_energies(simulation_parameters.save_interval);
  arma::mat buffer_configurations(
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
      LOGD << "write frames [" << offset_first_frame << "," << istep << ")";

      dset_energies
          .select({offset_first_frame}, {simulation_parameters.save_interval})
          .write_raw(buffer_energies.memptr());
      dset_configs
          .select(
              {offset_first_frame, 0},
              {simulation_parameters.save_interval, simulation.hamiltonian.number_of_magnetic_sites()})
          .write_raw(buffer_configurations.memptr());

      isave++;
      offset_first_frame = isave * simulation_parameters.save_interval;
    }

    buffer_energies.at(istep % simulation_parameters.save_interval) = simulation.runner.energy();
    buffer_configurations.col(istep % simulation_parameters.save_interval) = simulation.runner.spins();

    mean_energy += simulation.runner.energy();
    mean_magnetization += fabs(arma::sum(simulation.runner.spins()));
  }

  // write last buffer
  uint64_t remaining_frames = simulation_parameters.N  - offset_first_frame;
  LOGD << "write frames [" << offset_first_frame << "," << simulation_parameters.N << ")";
  dset_energies
      .select({offset_first_frame}, {remaining_frames})
      .write_raw(buffer_energies.memptr());
  dset_configs
      .select({offset_first_frame, 0}, {remaining_frames, simulation.hamiltonian.number_of_magnetic_sites()})
      .write_raw(buffer_configurations.memptr());

  std::cout << "*!> Done running!\n";

  // Statistics
  auto dN = static_cast<double>(simulation_parameters.N * simulation.hamiltonian.number_of_magnetic_sites());

  std::cout << "*> Statistics (N=" << simulation.hamiltonian.number_of_magnetic_sites() << ") :: "
            << "<E>/N = " << mean_energy / dN << ", "
            << "<|m|>/N = " << mean_magnetization / dN << "\n";

  // Save
  std::cout << "*!> Saving results\n";

  std::cout << "*!> Goodbye =)\n";

  return EXIT_SUCCESS;
}
