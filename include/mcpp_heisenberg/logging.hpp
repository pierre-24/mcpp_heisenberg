#ifndef  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_
#define  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_

#include <chrono>
#include <string>

#include <plog/Log.h>
#include <fmt/format.h>

namespace mch {

void init_logging(const plog::Severity& severity = plog::info);

namespace elapsed {

using clock_type = std::chrono::high_resolution_clock;
using float_sec = std::chrono::duration<double, std::chrono::seconds::period>;
using float_time_point = std::chrono::time_point<clock_type, float_sec>;

constexpr int one_minute = 60;
constexpr int one_hour = 60 * one_minute;
constexpr int one_day = 24 * one_hour;

class Chrono {
 private:
  float_time_point _start;

 public:
  Chrono(): _start{clock_type::now()} {}

  [[nodiscard]] std::string format() const {
    double elapsed = std::chrono::duration_cast<float_sec>(clock_type::now() - _start).count();

    int ndays = static_cast<int>(floor(elapsed / one_day));
    elapsed -= static_cast<float>(ndays * one_day);
    int nhours = static_cast<int>(floor(elapsed / one_hour));
    elapsed -= static_cast<float>(nhours * one_hour);
    int nmins = static_cast<int>(floor(elapsed / one_minute));
    elapsed -= static_cast<float>(nmins * one_minute);

    if (elapsed < 1) {
      return fmt::format("{:.2f} ms", elapsed * 1000);
    }

    if (nmins < 1) {
      return fmt::format("{:.2f} s", elapsed);
    }

    if (nhours < 1) {
      return fmt::format(
          "{} min{} and {:.1f} s", nmins, nmins < 2 ? "s": "", elapsed);
    }

    if (ndays < 1) {
      return fmt::format(
          "{} hour{} and {:.1f} mins", nhours, nhours < 2 ? "s": "",
          static_cast<float>(nmins) + (elapsed / one_minute));
    }

    return fmt::format(
        "{} day{}, {} hour{}, and {:.1f} mins",
        ndays, ndays < 2 ? "s": "", nhours, nhours < 2 ? "s": "",
        static_cast<float>(nmins) + (elapsed / one_minute));
  }
};

}  // namespace elapsed

}  // namespace mch

#endif  //  INCLUDE_MCPP_HEISENBERG_LOGGING_HPP_
