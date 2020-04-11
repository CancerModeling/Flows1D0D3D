////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "logger.hpp"
#include "aixlog.hpp"
#include "utilIO.hpp"

using namespace AixLog;

namespace {

/**
 * @brief
 * Formatted logging to libmesh out
 */

struct SinkLibmeshOut : public SinkFormat {

  SinkLibmeshOut(Type type,
                 const std::string &format =
                     "%Y-%m-%d %H-%M-%S.#ms [#severity] (#tag_func)")
      : SinkFormat(severity, type, format) {}

  void log(const Metadata &metadata,
           const std::string &message) override {
    do_log(out, metadata, message);
  }
};

}

util::Logger::Logger(const std::string &log_file, Parallel::Communicator *comm,
         bool screen_out)
      : TimeStepLog(), d_comm_p(comm), d_screen_out(screen_out) {

    if (d_dbg_file.is_open())
      d_dbg_file.close();
    d_dbg_file = std::ofstream(log_file + "_dbg.log", std::ios_base::out);

    if (d_ts_file.is_open())
      d_ts_file.close();
    d_ts_file = std::ofstream(log_file + "_time.log", std::ios_base::out);
}

util::Logger::~Logger() {

  if (d_dbg_file.is_open())
    d_dbg_file.close();

  if (d_ts_file.is_open())
    d_ts_file.close();
}

void util::Logger::log(const std::string &str, bool screen_out) {

    if (screen_out)
      out << str;
    if (d_comm_p->rank() == 0)
      d_dbg_file << str;
}

std::string util::Logger::log_ts_base(const int i, const int ns) {

    std::ostringstream oss;

    auto spS = util::io::getSpaceS(ns);

    if (d_comm_p->rank() == 0) {

      // initial stuff
      if (i == 0) {
        oss << "# TimeStepLog\n";
        oss << "Setup_Time " << d_setup_time.time_diff() << "\n";

        oss << "#\n";
        oss << "SysNames\n";
        for (unsigned int j=0; j < d_n; j++)
          oss << d_sys_names[j] << " " << j << "\n";

        oss << "#\n";
        oss << "SysAssemblyTimes\n";
        for (unsigned int j=0; j < d_n; j++) {

          oss << d_sys_assembly_time[j].time_diff();

          if (j < d_n - 1)
            oss << ", ";
          else
            oss << "\n";
        }

        oss << "#\n";

        // add header for next step
        oss << "StepLog\n";
        oss << "T, ST, PST, NETT, NLI, PNLI, ";

        for (unsigned int j=0; j < d_n; j++) {

          oss << "S" + std::to_string(j);

          if (j < d_sys_solve_time[i].size() - 1)
            oss << ", ";
          else
            oss << "\n";
        }
      }

      oss << d_times[i] << ", "
          << d_solve_time[i].time_diff() << ", "
          << d_pres_solve_time[i].time_diff() << ", "
          << d_update_network_time[i].time_diff() << ", "
          << d_nonlinear_iters[i] << ", "
          << d_pres_nonlinear_iters[i] << ", ";

      for (unsigned int j=0; j < d_sys_solve_time[i].size(); j++) {

        oss << d_sys_solve_time[i][j].time_diff();

        if (j < d_sys_solve_time[i].size() - 1)
          oss << ", ";
        else
          oss << "\n";
      }
    }

    return oss.str();
}