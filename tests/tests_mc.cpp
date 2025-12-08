#include <vector>
#include <filesystem>

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

  double T = 0.1;
  auto initial = mch::FerriInitialConfig(arma::vec(
      square_hamiltonian.number_of_magnetic_sites(), arma::fill::value(1.0))).make();
  auto runner = mch::IsingMonteCarloRunner(square_hamiltonian, initial);
  runner.reset_energy(T);

  arma::mat stats(MAX, 2);

  for (uint64_t i = 0; i < MAX; ++i) {
    runner.sweep(T);
    stats.row(i) = {runner.energy(), fabs(arma::sum(runner.spins()))};
  }

  LOGD << "last config is " << runner.spins();

  EXPECT_NEAR(arma::mean(stats.col(0)) / static_cast<double>(N * N), -2, 1e-3);  // <E>
  EXPECT_NEAR(arma::mean(stats.col(1)) / static_cast<double>(N * N), 1., 1e-3);  // <|m|>
}
