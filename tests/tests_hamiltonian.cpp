#include <vector>

#include <mcpp_heisenberg/hamiltonian.hpp>

#include "tests.hpp"

class HamiltonianTestsSuite : public MCHTestsSuite {
 protected:
  HamiltonianTestsSuite() = default;
};

/// Manual definition
TEST_F(HamiltonianTestsSuite, TestLinear) {
  auto hamiltonian = mch::Hamiltonian(4, {
      {{0, 1}, 1.},
      {{0, 3}, 1.},
      {{1, 2}, 1.},
      {{2, 3}, 1.},
  });

  auto spins = arma::uvec(hamiltonian.N(), arma::fill::ones);

  EXPECT_NEAR(hamiltonian.energy(spins), -1.0 * hamiltonian.N(), 1e-5);  // 2 interactions per site
}

/// Test chain (1D Ising)
TEST_F(HamiltonianTestsSuite, TestChain) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 100;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("chain", lattice, {{"H", 1}}, positions).to_supercell(4, 1, 1);

  auto hamiltonian = mch::Hamiltonian::from_geometry(geometry, {"H"}, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.N(), 4);

  auto spins = arma::uvec(hamiltonian.N(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -1.0 * hamiltonian.N(), 1e-5);
}

/// Test square lattice (2D Ising)
TEST_F(HamiltonianTestsSuite, TestSquare) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(4, 4, 1);

  auto hamiltonian = mch::Hamiltonian::from_geometry(geometry, {"H"}, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.N(), 16);

  auto spins = arma::uvec(hamiltonian.N(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -2.0 * hamiltonian.N(), 1e-5);  // 4 interactions per site
}

/// Test cube lattice (3D Ising), with selection of magnetic sites
TEST_F(HamiltonianTestsSuite, TestCubic) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 2;

  arma::mat positions = arma::mat(2, 3);
  positions.row(1) += 0.5;

  auto geometry = mch::Geometry("cube", lattice, {{"H", 1}, {"X", 1}}, positions).to_supercell(4, 4, 4);

  auto hamiltonian = mch::Hamiltonian::from_geometry(geometry, {"H"}, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.N(), 64);

  auto spins = arma::uvec(hamiltonian.N(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -3.0 * hamiltonian.N(), 1e-5);  // 6 interactions per site
}
