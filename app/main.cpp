#include <string>
#include <iostream>
#include <memory>
#include <vector>

#include <CLI/CLI.hpp>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>

#include "app.h"

/// Write data frames
void write_data_frames(
    HighFive::DataSet& dset_energies, HighFive::DataSet& dset_configs,
    uint64_t offset, uint64_t size, uint64_t N,
    const arma::mat& buffer_energies, const arma::mat& buffer_configurations) {
  // write energies
  dset_energies
      .select({offset}, {size})
      .write_raw(buffer_energies.memptr());

  // write spins
  std::vector<std::vector<int8_t>> configs;
  buffer_configurations.each_col([&configs, &N](auto& col) {
    std::vector<int8_t> config(N);
    std::transform(
        col.cbegin(), col.cend(), config.begin(), [](auto& val) { return (val < 0 ? -1 : 1); });

    configs.push_back(config);
  });

  dset_configs
      .select({offset, 0}, {size, N})
      .write(configs);
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

  // Set the log level
  mch::app::set_log_level();

  // Read input
  mch::app::Parameters simulation_parameters;
  mch::app::read_toml_input(input_file, simulation_parameters);

  // Create H5
  HighFive::File h5_file(output_file, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  // Prepare simulation
  auto simulation = mch::app::prepare_simulation(simulation_parameters, geometry_file);

  mch::app::save_simulation(h5_file, simulation);

  std::cout << "*!> Run for " << simulation_parameters.N << " steps\n";

  // Create group and buffer
  auto result_group = h5_file.createGroup("results");
  auto dset_energies = result_group.createDataSet<double>("energies", HighFive::DataSpace({simulation_parameters.N}));
  auto dset_configs = result_group.createDataSet<int8_t>(
      "configs", HighFive::DataSpace({simulation_parameters.N, simulation.hamiltonian.number_of_magnetic_sites()}));

  std::array<double, 2> info = {simulation_parameters.T, simulation_parameters.H};
  result_group.createDataSet("T&H", info).write(info);

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

      write_data_frames(dset_energies, dset_configs, offset_first_frame, simulation_parameters.save_interval,
                        simulation.hamiltonian.number_of_magnetic_sites(), buffer_energies, buffer_configurations);

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

  write_data_frames(dset_energies, dset_configs, offset_first_frame, remaining_frames,
                    simulation.hamiltonian.number_of_magnetic_sites(), buffer_energies, buffer_configurations);

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
