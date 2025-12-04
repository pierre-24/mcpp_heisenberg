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

  EXPECT_NEAR(arma::mean(stats.col(0)) / static_cast<double>(N * N), -2, 1e-3);  // <E>
  EXPECT_NEAR(arma::mean(stats.col(1)) / static_cast<double>(N * N), 1., 1e-3);  // <|m|>
}

/// Test square cluster update, in the low temperature limit
TEST_F(MCTestsSuite, TestSquareClusterUpdate) {
  uint64_t MAX = 10000;
  uint64_t PROD_START = MAX / 10 * 2;

  auto runner_sweep = mch::IsingMonteCarloRunner(square_hamiltonian);
  auto runner_cluster = mch::IsingMonteCarloRunner(square_hamiltonian);
  double T = 1.5;

  arma::mat stats_sweep(MAX - PROD_START, 2);
  arma::mat stats_cluster(MAX - PROD_START, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner_sweep.sweep(T);
    runner_cluster.cluster_update(T);

    if (i > PROD_START) {
      stats_sweep.row(i - PROD_START) = {runner_sweep.energy(), fabs(arma::sum(runner_sweep.spins()))};
      stats_cluster.row(i - PROD_START) = {runner_cluster.energy(), fabs(arma::sum(runner_cluster.spins()))};
    }
  }

  EXPECT_NEAR(arma::mean(stats_sweep.col(0)), arma::mean(stats_cluster.col(0)), 5e-1);  // <E>
  EXPECT_NEAR(arma::mean(stats_sweep.col(1)), arma::mean(stats_cluster.col(1)), 5e-1);  // <|m|>
}
