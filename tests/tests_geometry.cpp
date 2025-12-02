#include <vector>
#include <memory>

#include <mcpp_heisenberg/geometry.hpp>

#include "tests.hpp"

class GeometryTestsSuite : public MCHTestsSuite {
 protected:
  GeometryTestsSuite() = default;
};

TEST_F(GeometryTestsSuite, TestCreate) {
  std::vector<mch::ion_type_t> ions = {mch::ion_type_t("Mg", 1), mch::ion_type_t("O", 1)};

  auto lattice = arma::mat(3, 3);
  lattice.diag() += 2.0;
  auto positions = arma::mat(2, 3);
  positions.row(1) += 0.5;

  auto geometry = mch::Geometry("test", lattice, ions, positions);

  EXPECT_NEAR(geometry.get_abc().at(0), 2.0, 1e-5);
  EXPECT_NEAR(geometry.get_alpha_beta_gamma().at(0), M_PI / 2, 1e-5);
}

constexpr char* MgO_frac = "MgO (fractional)\n"
    "     4.211\n"
    "       0.0000000000      0.5000000000      0.5000000000\n"
    "       0.5000000000      0.0000000000      0.5000000000\n"
    "       0.5000000000      0.5000000000      0.0000000000\n"
    "Mg O\n"
    "1 1\n"
    "Direct\n"
    "       0.0      0.0      0.0 Mg\n"
    "       0.5      0.5      0.5 O";

TEST_F(GeometryTestsSuite, TestReadPoscar) {
  auto ss = std::make_shared<std::stringstream>();

  (*ss) << MgO_frac;
  ss->seekg(0);

  auto geometry = mch::Geometry::from_poscar(ss);

  arma::vec abc = geometry.get_abc(), abg = geometry.get_alpha_beta_gamma();
  EXPECT_NEAR(abc.at(0), abc.at(1), 1e-5);
  EXPECT_NEAR(abc.at(1), abc.at(2), 1e-5);
  EXPECT_NEAR(abc.at(0), 2.97763, 1e-3);

  EXPECT_NEAR(abg.at(0), abg.at(1), 1e-5);
  EXPECT_NEAR(abg.at(1), abg.at(2), 1e-5);
  EXPECT_NEAR(abg.at(0) * 180 / M_PI, 60., 1e-5);
}

constexpr char* MgO_cart = "MgO (Cartesian)\n"
    "1.0\n"
    "        2.9776265621         0.0000000000         0.0000000000\n"
    "        1.4888132811         2.5787002458         0.0000000000\n"
    "        1.4888132811         0.8595667486         2.4312219072\n"
    "   Mg    O\n"
    "    1    1\n"
    "Selective dynamics\n"
    "Cartesian\n"
    "     0.000000000         0.000000000         0.000000000 T T T\n"
    "     2.977626562         1.719133497         1.215610954 F F F";


TEST_F(GeometryTestsSuite, TestReadPoscarCartesian) {
  auto ss = std::make_shared<std::stringstream>();

  (*ss) << MgO_cart;
  ss->seekg(0);

  auto geometry = mch::Geometry::from_poscar(ss);

  arma::vec abc = geometry.get_abc(), abg = geometry.get_alpha_beta_gamma();
  EXPECT_NEAR(abc.at(0), abc.at(1), 1e-5);
  EXPECT_NEAR(abc.at(1), abc.at(2), 1e-5);
  EXPECT_NEAR(abc.at(0), 2.97763, 1e-3);

  EXPECT_NEAR(abg.at(0), abg.at(1), 1e-5);
  EXPECT_NEAR(abg.at(1), abg.at(2), 1e-5);
  EXPECT_NEAR(abg.at(0) * 180 / M_PI, 60., 1e-5);
}
