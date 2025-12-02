#ifndef INCLUDE_MCPP_HEISENBERG_ARMA_HPP_
#define INCLUDE_MCPP_HEISENBERG_ARMA_HPP_

#define ARMA_USE_OPENMP
#include <armadillo>

#include <sstream>
#include <string>

namespace mch::armautils {

/**
 * Check whether `a` and `b` are within `mult` $\times \varepsilon$ of each other.
 *
 * @tparam T floating point type
 * @param a number
 * @param b number
 * @param mult multiplier, default is 1
 * @return `true` if they are close
 */
template <typename T>
bool is_close(T a, T b, T mult = 1) {
    return (a-b) <= mult * arma::Datum<T>::eps;
}

}  // namespace mch::armautils

#endif  // INCLUDE_MCPP_HEISENBERG_ARMA_HPP_
