#include <vector>

#include <mcpp_heisenberg/hamiltonian.hpp>

#include "tests.hpp"

class HamiltonianTestsSuite : public MCHTestsSuite {
 protected:
  HamiltonianTestsSuite() = default;
};

/// Manual definition
TEST_F(HamiltonianTestsSuite, TestLinear) {
  auto hamiltonian = mch::IsingHamiltonian(4, {
      {{0, 1}, 1.},
      {{0, 3}, 1.},
      {{1, 2}, 1.},
      {{2, 3}, 1.},
  });

  auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::ones);

  EXPECT_NEAR(hamiltonian.energy(spins), -1.0 * hamiltonian.number_of_magnetic_sites(), 1e-5);
  // 2 interactions per site
}

/// Test chain (1D Ising)
TEST_F(HamiltonianTestsSuite, TestChain) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 100;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("chain", lattice, {{"H", 1}}, positions).to_supercell(4, 1, 1);

  auto hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.number_of_magnetic_sites(), 4);

  auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -1.0 * hamiltonian.number_of_magnetic_sites(), 1e-5);
}

/// Test square lattice (2D Ising)
TEST_F(HamiltonianTestsSuite, TestSquare) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(4, 4, 1);

  auto hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.number_of_magnetic_sites(), 16);

  auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -2.0 * hamiltonian.number_of_magnetic_sites(), 1e-5);
  // 4 interactions per site
}

/// Test cube lattice (3D Ising), with selection of magnetic sites
TEST_F(HamiltonianTestsSuite, TestCubic) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 2;

  arma::mat positions = arma::mat(2, 3);
  positions.row(1) += 0.5;

  auto geometry = mch::Geometry("cube", lattice, {{"H", 1}, {"X", 1}}, positions)
                      .filter_atoms({"H"})
                      .to_supercell(4, 4, 4);

  auto hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});

  EXPECT_EQ(hamiltonian.number_of_magnetic_sites(), 64);

  auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::ones);
  EXPECT_NEAR(hamiltonian.energy(spins), -3.0 * hamiltonian.number_of_magnetic_sites(), 1e-5);
  // 6 interactions per site
}

/// Test ΔE due to flippling one spin
TEST_F(HamiltonianTestsSuite, TestSquareDeltaE) {
  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(4, 4, 1);

  auto hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});

  for (int i = 0; i < 50; ++i) {
    auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::randu);
    spins.for_each([](double& e) { e = e > .5 ? 1.0 : -1.0; });

    double Eb = hamiltonian.energy(spins);
    double dE = hamiltonian.delta_energy(spins, 1);

    spins.at(1) *= -1;
    EXPECT_NEAR(dE, hamiltonian.energy(spins) - Eb, 1e-4);
  }
}

/// Test ΔE due to flippling one spin (with H)
TEST_F(HamiltonianTestsSuite, TestSquareDeltaEwithH) {
  double H = 1.0;

  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(4, 4, 1);

  auto hamiltonian = mch::IsingHamiltonian::from_geometry(geometry, {{"H", "H", 2.0, 1.0}});

  for (int i = 0; i < 50; ++i) {
    auto spins = arma::vec(hamiltonian.number_of_magnetic_sites(), arma::fill::randu);
    spins.for_each([](double& e) { e = e > .5 ? 1.0 : -1.0; });

    double Eb = hamiltonian.energy(spins, H);
    double dE = hamiltonian.delta_energy(spins, 1, H);

    spins.at(1) *= -1;
    EXPECT_NEAR(dE, hamiltonian.energy(spins, H) - Eb, 1e-4);
  }
}
