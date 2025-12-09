#include <map>
#include <string>

#include <mcpp_heisenberg/configs.hpp>

namespace mch {

arma::vec make_config(
    const Geometry& geometry, const std::map<std::string, double>& spin_values,
    const std::map<std::string, ConfigType>& ions_configs) {
  std::random_device rd;
  auto rng = std::mt19937(rd());
  std::bernoulli_distribution dis(0.5);

  arma::vec config(geometry.number_of_atoms(), arma::fill::value(1.0));

  uint64_t nx = 0;
  for (auto& iondef : geometry.ions()) {
    // set initial value
    if (spin_values.contains(iondef.first)) {
      config.subvec(nx, nx + iondef.second - 1) *= spin_values.at(iondef.first);
    }

    // switch if any
    if (ions_configs.contains(iondef.first)) {
      auto config_type = ions_configs.at(iondef.first);

      if (config_type == ConfigType::AllDown) {
        config.subvec(nx, nx + iondef.second - 1) *= -1;
      } else if (config_type == ConfigType::AntiFerroAZ
          || config_type == ConfigType::AntiFerroAY
          || config_type == ConfigType::AntiFerroAX) {
        uint64_t axis = 0;
        switch (config_type) {
          case ConfigType::AntiFerroAX:
            axis = 0;
            break;
          case ConfigType::AntiFerroAY:
            axis = 1;
            break;
          case ConfigType::AntiFerroAZ:
            axis = 2;
            break;
          default:
            break;
        }

        arma::mat positions = geometry.positions().rows(nx, nx + iondef.second - 1);
        arma::uvec idx = arma::sort_index(positions.col(axis), "ascend");

        double c = positions.at(idx[0], axis);
        double spin = 1.0;
        for (auto& ipos : idx) {
          if (positions.at(ipos, axis) != c) {
            c = positions.at(ipos, axis);
            spin *= -1.0;
          }

          config.at(nx + ipos) *= spin;
        }
      } else if (config_type == ConfigType::Random) {
        config.subvec(nx, nx + iondef.second - 1).for_each([&dis, &rng](auto& val) {
          if (dis(rng)) {
            val *= -1;
          }
        });
      }
    }

    nx += iondef.second;
  }

  return config;
}

}  // namespace mch
