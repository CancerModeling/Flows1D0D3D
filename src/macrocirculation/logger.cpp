////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "logger.hpp"
#include "aixlog.hpp"

namespace {

/*! @brief Formatted logging to libmesh out */
struct SinkLibmeshOut : public AixLog::SinkFormat {

  SinkLibmeshOut(int rank, AixLog::Severity severity, AixLog::Type type,
                 const std::string &format = "%Y-%m-%d %H-%M-%S.#ms [#severity] (#tag_func)")
      : d_rank(rank), AixLog::SinkFormat(severity, type, format) {}

  void log(const AixLog::Metadata &metadata,
           const std::string &message) override {
    if (d_rank == 0)
      do_log(std::cout, metadata, message);
  }

  int d_rank;
};

struct MySinkFile : public AixLog::SinkFormat {
  MySinkFile(int rank, AixLog::Severity severity, AixLog::Type type,
             const std::string &filename,
             const std::string &format =
               "%Y-%m-%d %H-%M-%S.#ms [#severity] (#tag_func)")
      : d_rank(rank), AixLog::SinkFormat(severity, type, format) {
    if (d_rank == 0)
      ofs.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
  }

  ~MySinkFile() override { if (d_rank == 0) ofs.close(); }

  void log(const AixLog::Metadata &metadata,
           const std::string &message) override {
    if (d_rank == 0)
      do_log(ofs, metadata, message);
  }

protected:
  mutable std::ofstream ofs;
  int d_rank;
};

std::string helper(int width, int s_width, std::string &tag) {
  int len = s_width + tag.length() + 5;
  if (width < len) {
    return tag;
  }

  int diff = width - len;
  int pad1 = diff / 2;
  int pad2 = diff - pad1;
  return std::string(pad1, ' ') + tag + std::string(pad2, ' ');
}
} // namespace

macrocirculation::Logger::Logger(const std::string &log_file, int rank,
                     bool screen_out)
    : d_rank(rank), d_screen_out(screen_out) {

  if (screen_out) {
    auto sink_cout = std::make_shared<SinkLibmeshOut>(d_rank, AixLog::Severity::trace,
                                                      AixLog::Type::normal);
    auto sink_file = std::make_shared<MySinkFile>(d_rank,
                                                  AixLog::Severity::trace, AixLog::Type::all, log_file + ".log");
    AixLog::Log::init({sink_cout, sink_file});
  } else {
    auto sink_cout = std::make_shared<SinkLibmeshOut>(d_rank, AixLog::Severity::notice,
                                                      AixLog::Type::normal);
    auto sink_file = std::make_shared<MySinkFile>(d_rank,
                                                  AixLog::Severity::trace, AixLog::Type::all, log_file + ".log");
    AixLog::Log::init({sink_cout, sink_file});
  }

}

macrocirculation::Logger::~Logger() {

  if (d_rank == 0 and d_dbg_file.is_open())
    d_dbg_file.close();
}

void macrocirculation::Logger::log(const std::string &str, std::string tag,
                       std::string severity) {

  if (d_rank != 0)
    return;

  if (severity == "info")
    LOG(INFO, helper(22, 4, tag)) << str;
  else if (severity == "notice")
    LOG(NOTICE, helper(22, 6, tag)) << str;
  else if (severity == "warning")
    LOG(WARNING, helper(22, 7, tag)) << str;
  else if (severity == "debug")
    LOG(DEBUG, helper(22, 5, tag)) << str;
}