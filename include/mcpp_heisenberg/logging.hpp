#ifndef  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_
#define  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_

#include <plog/Log.h>

namespace mch {

void init_logging(const plog::Severity& severity = plog::info);

}  // namespace mch

#endif  //  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_
