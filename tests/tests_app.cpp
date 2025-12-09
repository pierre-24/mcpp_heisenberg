#include <filesystem>
#include <cstdio>
#include <utility>
#include <string>
#include <vector>

#include <mcpp_heisenberg/mcpp_heisenberg.hpp>
#include "../app/app.h"

#include "tests.hpp"

class AppTestsSuite : public MCHTestsSuite {
 protected:
  mch::Geometry geometry;
  mch::app::Parameters parameters;
  std::filesystem::path temp_path;

  AppTestsSuite() {
    // Create temporary
    temp_path = std::filesystem::temp_directory_path() / std::tmpnam(nullptr);

    // Create a square lattice geometry
    arma::mat lattice = arma::eye(3, 3);
    lattice.at(0, 0) = 1;
    lattice.at(1, 1) = 1;
    lattice.at(2, 2) = 100;

    arma::mat positions = arma::mat(1, 3);

    geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions);

    // Tweak a few parameters
    parameters.T = 2.4;
    parameters.N = 10000;
    parameters.supercell_size = {10, 10, 1};
    parameters.magnetic_sites = {"H"};
    parameters.pair_defs = {{"H", "H", 1.0, 1.0}};
  }

  void TearDown() {
    EXPECT_TRUE(std::filesystem::remove(temp_path));
  }
};

/// Test if H5 contains what it should
TEST_F(AppTestsSuite, TestSaveSimple) {
  parameters.N = 200;
  parameters.supercell_size = {5, 5, 1};
  parameters.H = 0.1;
  parameters.initial_configs["H"] = mch::ConfigType::Random;

  // run & save
  {
    HighFive::File file(temp_path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    mch::app::Runner runner{parameters, geometry, std::move(file)};
    runner.run(parameters);
  }

  // read back and check
  {
    HighFive::File file(temp_path, HighFive::File::ReadOnly);

    // geometry
    auto geometry_group = file.getGroup("geometry");

    arma::mat re_lattice(3, 3);
    geometry_group.getDataSet("lattice_vectors").read_raw(re_lattice.memptr());
    arma::vec d(3);
    d = {
        parameters.supercell_size[0] * geometry.lattice().at(0, 0),
        parameters.supercell_size[1] * geometry.lattice().at(1, 1),
        parameters.supercell_size[2] * geometry.lattice().at(2, 2)};
    EXPECT_TRUE(arma::approx_equal(d, re_lattice.diag(), "abstol", 1e-5));

    std::vector<std::string> ion_types;
    geometry_group.getDataSet("ion_types").read(ion_types);
    EXPECT_EQ(ion_types, std::vector<std::string>{"H"});

    std::vector<uint64_t> ion_numbers;
    geometry_group.getDataSet("ion_numbers").read(ion_numbers);
    EXPECT_EQ(ion_numbers, std::vector<uint64_t>{25});

    std::vector<double> spin_values;
    geometry_group.getDataSet("spin_values").read(spin_values);
    EXPECT_EQ(spin_values, std::vector<double>(25, 1.0));

    std::vector<double> magnetic_anisotropies;
    geometry_group.getDataSet("magnetic_anisotropies").read(magnetic_anisotropies);
    EXPECT_EQ(magnetic_anisotropies, std::vector<double>(25, 0.0));

    // hamiltonian
    auto hamiltonian_group = file.getGroup("hamiltonian");

    std::vector<double> Js;
    hamiltonian_group.getDataSet("J").read(Js);
    EXPECT_EQ(Js, std::vector<double>(50, 1.0));

    std::vector<double> xparameters;
    hamiltonian_group.getDataSet("parameters").read(xparameters);
    std::vector<double> yparameters = {parameters.kB, parameters.T, parameters.muB, parameters.H};
    EXPECT_EQ(xparameters, yparameters);

    // results
    std::vector<std::array<uint64_t, 2>> pairs;
    hamiltonian_group.getDataSet("pairs").read(pairs);
    std::vector<mch::jpair_t> xpairs;
    for (uint64_t ipair = 0; ipair < pairs.size(); ++ipair) {
      xpairs.push_back({{pairs.at(ipair)[0], pairs.at(ipair)[1]}, Js.at(ipair)});
    }

    auto results_group = file.getGroup("results");

    arma::mat configs(25, parameters.N);
    results_group.getDataSet("configs").read_raw(configs.memptr());

    arma::mat aggs(2, parameters.N);
    results_group.getDataSet("aggregated_data").read_raw(aggs.memptr());

    auto hamiltonian = mch::IsingHamiltonian(25, xpairs);

    for (uint64_t istep = 0; istep < parameters.N; ++istep) {
      arma::vec config(configs.col(istep));

      EXPECT_EQ(aggs.at(1, istep), static_cast<double>(arma::sum(config)));
      EXPECT_NEAR(
          aggs.at(0, istep),
          hamiltonian.energy(config, parameters.muB * parameters.H), 1e-5);
    }
  }
}
