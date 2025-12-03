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

  EXPECT_NEAR(hamiltonian.energy(spins), -1.0 * hamiltonian.N(), 1e-5);
}
