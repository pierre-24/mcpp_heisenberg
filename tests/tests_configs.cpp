#include <map>
#include <string>
#include <vector>

#include <mcpp_heisenberg/configs.hpp>

#include "tests.hpp"

class ConfigsTestsSuite : public MCHTestsSuite {
 protected:
  arma::mat lattice = arma::mat(3, 3);
  arma::mat positions = arma::mat(2, 3);
  std::vector<mch::ion_type_t> ions = {mch::ion_type_t("Mg", 1), mch::ion_type_t("O", 1)};

  mch::Geometry geometry;

  ConfigsTestsSuite() {
    // BCC lattice
    lattice.diag() += 2.0;
    positions.row(1) += 0.5;

    geometry = mch::Geometry("BCC", lattice, ions, positions);
  }
};

TEST_F(ConfigsTestsSuite, TestAllUp) {
  auto test_geometry = geometry.to_supercell(1, 1, 4);

  std::map<std::string, double> spin_values;
  spin_values["Mg"] = 2.0;
  spin_values["O"] = .5;

  std::map<std::string, mch::ConfigType> initial_config;
  initial_config["Mg"] = mch::ConfigType::AllUp;
  initial_config["O"] = mch::ConfigType::AllUp;

  auto config = mch::make_config(test_geometry, spin_values, initial_config);

  arma::vec spins = {2., 2., 2., 2., .5, .5, .5, .5};
  EXPECT_TRUE(arma::approx_equal(config, spins, "abstol", 1e-5));
}

TEST_F(ConfigsTestsSuite, TestAllDown) {
  auto test_geometry = geometry.to_supercell(1, 1, 4);

  std::map<std::string, double> spin_values;
  spin_values["Mg"] = 2.0;
  spin_values["O"] = .5;

  std::map<std::string, mch::ConfigType> initial_config;
  initial_config["Mg"] = mch::ConfigType::AllDown;
  initial_config["O"] = mch::ConfigType::AllUp;

  auto config = mch::make_config(test_geometry, spin_values, initial_config);

  arma::vec spins = {-2., -2., -2., -2., .5, .5, .5, .5};
  EXPECT_TRUE(arma::approx_equal(config, spins, "abstol", 1e-5));
}

TEST_F(ConfigsTestsSuite, TestAntifFerroAZ) {
  auto test_geometry = geometry.to_supercell(1, 1, 4);

  std::map<std::string, double> spin_values;
  spin_values["Mg"] = 2.0;
  spin_values["O"] = .5;

  std::map<std::string, mch::ConfigType> initial_config;
  initial_config["Mg"] = mch::ConfigType::AntiFerroAZ;
  initial_config["O"] = mch::ConfigType::AllUp;

  auto config = mch::make_config(test_geometry, spin_values, initial_config);

  arma::vec spins = {2., -2., 2., -2., .5, .5, .5, .5};
  EXPECT_TRUE(arma::approx_equal(config, spins, "abstol", 1e-5));
}

TEST_F(ConfigsTestsSuite, TestAntifFerroAXY) {
  auto test_geometry = geometry.to_supercell(4, 2, 2);

  std::map<std::string, double> spin_values;
  spin_values["Mg"] = 2.0;
  spin_values["O"] = .5;

  std::map<std::string, mch::ConfigType> initial_config;
  initial_config["Mg"] = mch::ConfigType::AntiFerroAX;
  initial_config["O"] = mch::ConfigType::AntiFerroAY;

  auto config = mch::make_config(test_geometry, spin_values, initial_config);

  arma::vec spins = {
      2, 2, 2, 2,
      -2, -2, -2, -2,
      2, 2, 2, 2,
      -2, -2, -2, -2,
      .5, .5, -.5, -.5,
      .5, .5, -.5, -.5,
      .5, .5, -.5, -.5,
      .5, .5, -.5, -.5
  };

  EXPECT_TRUE(arma::approx_equal(config, spins, "abstol", 1e-5));
}
