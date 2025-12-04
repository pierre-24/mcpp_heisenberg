#include <vector>
#include <filesystem>
#include <cstdio>

#include <mcpp_heisenberg/mc.hpp>

#include "tests.hpp"

class MCTestsSuite : public MCHTestsSuite {
 protected:
  mch::IsingHamiltonian square_hamiltonian;
  uint64_t N = 10;

  MCTestsSuite() {
    arma::mat lattice = arma::eye(3, 3);
    lattice.at(0, 0) = 2;
    lattice.at(1, 1) = 2;
    lattice.at(2, 2) = 100;

    arma::mat positions = arma::mat(1, 3);

    auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(N, N, 1);

    square_hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});
  }
};

/// Test the low temperature limit
TEST_F(MCTestsSuite, TestSquareLowTemp) {
  uint64_t MAX = 10000;
  uint64_t PROD_START = MAX / 10 * 2;

  auto runner = mch::IsingMonteCarloRunner(square_hamiltonian);
  double T = 0.1;

  arma::mat stats(MAX - PROD_START, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner.sweep(T);

    if (i > PROD_START) {
      stats.row(i - PROD_START) = {runner.energy(), fabs(arma::sum(runner.spins()))};
    }
  }

  EXPECT_NEAR(arma::mean(stats.col(0)) / static_cast<double>(N * N), -2, 1e-3);
  EXPECT_NEAR(fabs(arma::mean(stats.col(1))) / static_cast<double>(N * N), 1., 1e-3);
}

/// Test the high temperature limit
TEST_F(MCTestsSuite, TestSquareHighTemp) {
  uint64_t MAX = 10000;
  uint64_t PROD_START = MAX / 10 * 2;

  auto runner = mch::IsingMonteCarloRunner(square_hamiltonian);
  double T = 10.0;

  arma::mat stats(MAX - PROD_START, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner.sweep(T);

    if (i > PROD_START) {
      stats.row(i - PROD_START) = {runner.energy(), fabs(arma::sum(runner.spins()))};
    }
  }

  EXPECT_NEAR(fabs(arma::mean(stats.col(1))) / static_cast<double>(N * N), .0, 5e-1);
}

/// Test save
TEST_F(MCTestsSuite, TestSquareSave) {
  uint64_t MAX = 10;

  auto runner = mch::IsingMonteCarloRunner(square_hamiltonian);
  double T = 2.0;

  arma::mat stats(MAX, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner.sweep(T);
    stats.row(i) = {runner.energy(), arma::sum(runner.spins())};
  }

  // Save
  std::filesystem::path temp_path = std::filesystem::temp_directory_path() /= std::tmpnam(nullptr);
  {
    HighFive::File file(temp_path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

    auto group = file.createGroup("results");
    runner.save(group);
  }

  // Load
  {
    HighFive::File file(temp_path, HighFive::File::ReadOnly);
    auto group = file.getGroup("results");

    auto dset_energy = group.getDataSet("energies");
    arma::vec energies(MAX + 1);
    dset_energy.read_raw(energies.memptr());

    auto dset_configs = group.getDataSet("configs");
    arma::mat configs(square_hamiltonian.N(), MAX + 1);
    dset_configs.read_raw(configs.memptr());

    EXPECT_TRUE(arma::approx_equal(stats.col(0), energies.subvec(1, MAX), "abstol", 1e-5));
    EXPECT_TRUE(arma::approx_equal(arma::sum(configs.cols(1, MAX), 0), stats.col(1).t(), "abstol", 1e-5));
  }

  EXPECT_TRUE(std::filesystem::remove(temp_path));
}

/// Test square cluster update, in the low temperature limit
TEST_F(MCTestsSuite, TestSquareClusterUpdate) {
  uint64_t MAX = 10000;
  uint64_t PROD_START = MAX / 10 * 2;

  auto runner = mch::IsingMonteCarloRunner(square_hamiltonian);
  double T = 0.1;

  arma::mat stats(MAX - PROD_START, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner.cluster_update(T);

    if (i > PROD_START) {
      stats.row(i - PROD_START) = {runner.energy(), fabs(arma::sum(runner.spins()))};
    }
  }

  EXPECT_NEAR(arma::mean(stats.col(0)) / static_cast<double>(N * N), -2, 1e-3);  // <E>
  EXPECT_NEAR(fabs(arma::mean(stats.col(1))) / static_cast<double>(N * N), 1., 1e-3);  // <|m|>
}
