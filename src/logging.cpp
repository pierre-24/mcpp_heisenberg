#include <mcpp_heisenberg/logging.hpp>

#include <plog/Init.h>
#include <plog/Appenders/ConsoleAppender.h>

namespace mch {

class CustomFormatter {
 public:
  static plog::util::nstring header() { return ""; }
  static plog::util::nstring format(const plog::Record& record) {
    plog::util::nostringstream ss;

    if (record.getSeverity() <= plog::warning) {
      ss << PLOG_NSTR("!! ") << plog::severityToString(record.getSeverity()) << PLOG_NSTR(": ");
    }

    ss << record.getMessage() << PLOG_NSTR("\n");
    return ss.str();
  }
};

void init_logging(const plog::Severity& severity) {
  static plog::ConsoleAppender<CustomFormatter> consoleAppender;
  plog::init(severity, &consoleAppender);
}

}  // namespace mch
