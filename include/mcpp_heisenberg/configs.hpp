#ifndef INCLUDE_MCPP_HEISENBERG_CONFIGS_HPP_
#define INCLUDE_MCPP_HEISENBERG_CONFIGS_HPP_

#include <map>
#include <string>

#include <mcpp_heisenberg/geometry.hpp>

namespace mch {

enum ConfigType {
  Random,
  AllUp,
  AllDown,
  AntiFerroAX,
  AntiFerroAY,
  AntiFerroAZ,
};

arma::vec make_config(
    const Geometry& geometry, const std::map<std::string, double>& spin_values,
    const std::map<std::string, ConfigType>& ions_configs);

}   // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_CONFIGS_HPP_
