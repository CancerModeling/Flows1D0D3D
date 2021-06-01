////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_LOGGER_H
#define UTILS_LOGGER_H

#include "utils.hpp"

namespace macrocirculation {

/*! @brief Prints message to std::cout and also writes to the file */
class Logger {

public:
  Logger(const std::string &log_file, int rank,
         bool screen_out = true);

  ~Logger();

  // Output string message
  void log(const std::string &str, std::string tag,
           std::string severity);

  // All methods below overload above log function
  void log(const std::string &str, std::string tag) {
    log(str, tag, "info");
  };
  void log(const std::string &str) {
    log(str, "debug", "info");
  };
  void log(std::ostringstream &oss, std::string tag,
           std::string severity) {
    log(oss.str(), tag, severity);
    oss.str("");
    oss.clear();
  };
  void log(std::ostringstream &oss, std::string tag) {
    log(oss, tag, "info");
  };
  void log(std::ostringstream &oss) {
    log(oss, "debug", "info");
  };

  // overload function call
  void operator()(std::ostringstream &oss, std::string tag,
                  std::string severity) {
    log(oss, tag, severity);
  }
  void operator()(std::ostringstream &oss, std::string tag) {
    log(oss, tag, "info");
  }
  void operator()(std::ostringstream &oss) {
    log(oss, "debug", "info");
  }
  void operator()(const std::string &str, std::string tag,
                  std::string severity) {
    log(str, tag, severity);
  }
  void operator()(const std::string &str, std::string tag) {
    log(str, tag, "info");
  }
  void operator()(const std::string &str) {
    log(str, "debug", "info");
  }

private:
  int d_rank;
  bool d_screen_out;
  std::ofstream d_dbg_file;
};

} // namespace macrocirculation

#endif // UTILS_LOGGER_H
