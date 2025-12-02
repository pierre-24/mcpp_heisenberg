#include <iostream>

#include <plog/Init.h>
#include <plog/Appenders/ConsoleAppender.h>

#include "tests.hpp"

class CustomDebugFormatter {
 public:
  static plog::util::nstring header() { return ""; }
  static plog::util::nstring format(const plog::Record& record) {
    plog::util::nostringstream ss;
    ss << plog::severityToString(record.getSeverity())
       << PLOG_NSTR(" (") << record.getFile() << PLOG_NSTR(":") << record.getLine() << PLOG_NSTR("): ")
       << record.getMessage() << PLOG_NSTR("\n");
    return ss.str();
  }
};

void init_logging_debug(const plog::Severity& severity = plog::debug) {
  static plog::ConsoleAppender<CustomDebugFormatter> consoleAppender;
  plog::init(severity, &consoleAppender);
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  init_logging_debug();

  return RUN_ALL_TESTS();
}
