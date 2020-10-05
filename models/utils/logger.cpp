////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "logger.hpp"
#include "aixlog.hpp"
#include "utilIO.hpp"

// using namespace AixLog;

namespace {

/**
 * @brief
 * Formatted logging to libmesh out
 */

struct SinkLibmeshOut : public AixLog::SinkFormat {

  SinkLibmeshOut(AixLog::Severity severity, AixLog::Type type,
                 const std::string &format = "%Y-%m-%d %H-%M-%S.#ms [#severity] (#tag_func)")
      : AixLog::SinkFormat(severity, type, format) {}

  //  SinkLibmeshOut(AixLog::Severity severity, AixLog::Type type,
  //                 const std::string &format = "%Y-%m-%d %H-%M-%S.#ms")
  //      : AixLog::SinkFormat(severity, type, format) {}

  void log(const AixLog::Metadata &metadata,
           const std::string &message) override {
    do_log(out, metadata, message);
  }
};

struct MySinkFile : public AixLog::SinkFormat {
  MySinkFile(AixLog::Severity severity, AixLog::Type type,
             const std::string &filename,
             const std::string &format =
               "%Y-%m-%d %H-%M-%S.#ms [#severity] (#tag_func)")
      : AixLog::SinkFormat(severity, type, format) {
    ofs.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
  }

  //  MySinkFile(AixLog::Severity severity, AixLog::Type type,
  //             const std::string &filename,
  //             const std::string &format =
  //             "%Y-%m-%d %H-%M-%S.#ms")
  //      : AixLog::SinkFormat(severity, type, format) {
  //    ofs.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
  //  }

  ~MySinkFile() override { ofs.close(); }

  void log(const AixLog::Metadata &metadata,
           const std::string &message) override {
    do_log(ofs, metadata, message);
  }

protected:
  mutable std::ofstream ofs;
};

std::string helper(int width, const std::string &severity, std::string &tag) {
  int len = severity.length() + tag.length() + 5;
  if (width < len) {
    return "[" + severity + "] (" + tag + ")";
  }

  int diff = width - len;
  int pad1 = diff / 2;
  int pad2 = diff - pad1;
  return "[" + severity + "]" + std::string(pad1, ' ') + "(" + tag + ")" +
         std::string(pad2, ' ');
}

std::string helper(int width, int s_width, std::string &tag) {
  int len = s_width + tag.length() + 5;
  if (width < len) {
    return tag;
  }

  int diff = width - len;
  return tag + std::string(diff, ' ');
}

std::string helper2(int width, int s_width, std::string &tag) {
  int len = s_width + tag.length() + 5;
  if (width < len) {
    return tag;
  }

  int diff = width - len;
  int pad1 = diff / 2;
  int pad2 = diff - pad1;
  return std::string(pad1, ' ') + tag + std::string(pad2, ' ');
}

template<class T>
inline T get_avg(const std::vector<T> &list) {

  if (list.size() == 0)
    return 0.;

  T avg = 0.;
  for (const auto &l : list)
    avg += l;
  avg = avg / (T(list.size()));

  return avg;
}

template<class T>
inline T get_total(const std::vector<T> &list) {

  if (list.size() == 0)
    return 0.;

  T avg = 0.;
  for (const auto &l : list)
    avg += l;
  return avg;
}

} // namespace

util::Logger::Logger(const std::string &log_file, Parallel::Communicator *comm,
                     bool screen_out)
    : TimeStepLog(), d_comm_p(comm), d_screen_out(screen_out) {

  if (screen_out) {
    auto sink_cout = std::make_shared<SinkLibmeshOut>(AixLog::Severity::trace,
                                                      AixLog::Type::normal);
    auto sink_file = std::make_shared<MySinkFile>(
      AixLog::Severity::trace, AixLog::Type::all, log_file + ".log");
    AixLog::Log::init({sink_cout, sink_file});
  } else {
    auto sink_cout = std::make_shared<SinkLibmeshOut>(AixLog::Severity::notice,
                                                      AixLog::Type::normal);
    auto sink_file = std::make_shared<MySinkFile>(
      AixLog::Severity::trace, AixLog::Type::all, log_file + ".log");
    AixLog::Log::init({sink_cout, sink_file});
  }

  if (d_ts_file.is_open())
    d_ts_file.close();
  d_ts_file = std::ofstream(log_file + "_ts.txt", std::ios_base::out);

  if (d_qoi_file.is_open())
    d_qoi_file.close();
  d_qoi_file = std::ofstream(log_file + "_qoi.txt", std::ios_base::out);
}

util::Logger::~Logger() {

  if (d_dbg_file.is_open())
    d_dbg_file.close();

  if (d_ts_file.is_open())
    d_ts_file.close();
}

void util::Logger::log(const std::string &str, std::string tag,
                       std::string severity) {

  if (severity == "info")
    LOG(INFO, helper2(22, 4, tag)) << str;
  else if (severity == "notice")
    LOG(NOTICE, helper2(22, 6, tag)) << str;
  else if (severity == "warning")
    LOG(WARNING, helper2(22, 7, tag)) << str;
  else if (severity == "debug")
    LOG(DEBUG, helper2(22, 5, tag)) << str;
}

std::string util::Logger::log_ts_base(const int i, const int ns) {

  std::ostringstream oss;

  auto spS = util::io::getSpaceS(ns);

  if (d_comm_p->rank() == 0) {

    // initial stuff
    if (i == 0) {
      oss << spS << "# TimeStepLog\n";
      oss << spS << "SetupTime " << d_setup_time.time_diff() << "\n";

      oss << spS << "#\n";
      oss << spS << "SysNames " << d_n << "\n";
      for (unsigned int j = 0; j < d_n; j++)
        oss << spS << d_sys_names[j] << " " << j << "\n";

      oss << spS << "#\n";
      oss << "SysAssemblyTimes\n";
      for (unsigned int j = 0; j < d_n; j++) {

        if (j == 0)
          oss << spS;

        oss << d_sys_assembly_time[j].time_diff();

        if (j < d_n - 1)
          oss << ", ";
        else
          oss << "\n";
      }

      oss << spS << "#\n";

      // add header for next step
      oss << spS << "StepLog " << 6 + d_n << "\n";
      oss << spS << "'T' 'ST' 'PST' 'NETT' 'NLI' 'PNLI' ";

      for (unsigned int j = 0; j < d_n; j++) {

        oss << "'S" + std::to_string(j) + "'";

        if (j < d_sys_solve_time[i].size() - 1)
          oss << " ";
        else
          oss << "\n";
      }
    }

    oss << spS << d_times[i] << " " << d_solve_time[i].time_diff() << " "
        << d_pres_solve_time[i].time_diff() << " "
        << d_update_network_time[i].time_diff() << " " << d_nonlinear_iters[i]
        << " " << d_pres_nonlinear_iters[i] << " ";

    for (unsigned int j = 0; j < d_sys_solve_time[i].size(); j++) {

      oss << d_sys_solve_time[i][j].time_diff();

      if (j < d_sys_solve_time[i].size() - 1)
        oss << " ";
      else
        oss << "\n";
    }
  }

  return oss.str();
}

std::string util::Logger::log_ts_base_final_avg(double sim_time, const int ns) {

  std::ostringstream oss;

  auto spS = util::io::getSpaceS(ns);

  if (d_comm_p->rank() == 0) {

    oss << spS << "#\n";

    // add header for next step
    oss << spS << "AvgStepLog " << 6 + d_n << "\n";
    oss << spS << "'T' 'ST' 'PST' 'NETT' 'NLI' 'PNLI' ";

    for (unsigned int j = 0; j < d_n; j++) {

      oss << "'S" + std::to_string(j) + "'";

      if (j < d_sys_solve_time[0].size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    // compute avg and write

    oss << spS << get_avg(d_times) << " "
        << get_avg(get_delta_t_vec(d_solve_time)) << " "
        << get_avg(get_delta_t_vec(d_pres_solve_time)) << " "
        << get_avg(get_delta_t_vec(d_update_network_time)) << " "
        << get_avg(d_nonlinear_iters) << " "
        << get_avg(d_pres_nonlinear_iters) << " ";

    for (unsigned int j = 0; j < d_sys_solve_time[0].size(); j++) {

      std::vector<float> delta_t;
      for (auto &i : d_sys_solve_time)
        delta_t.push_back(i[j].time_diff());

      oss << get_avg(delta_t);

      if (j < d_sys_solve_time[0].size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    oss << spS << "#\n";

    // add header for next step
    oss << spS << "TotalStepLog " << 6 + d_n << "\n";
    oss << spS << "'T' 'ST' 'PST' 'NETT' 'NLI' 'PNLI' ";

    for (unsigned int j = 0; j < d_n; j++) {

      oss << "'S" + std::to_string(j) + "'";

      if (j < d_sys_solve_time[0].size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    // compute avg and write

    oss << spS << get_total(d_times) << " "
        << get_total(get_delta_t_vec(d_solve_time)) << " "
        << get_total(get_delta_t_vec(d_pres_solve_time)) << " "
        << get_total(get_delta_t_vec(d_update_network_time)) << " "
        << get_total(d_nonlinear_iters) << " "
        << get_total(d_pres_nonlinear_iters) << " ";

    for (unsigned int j = 0; j < d_sys_solve_time[0].size(); j++) {

      std::vector<float> delta_t;
      for (auto &i : d_sys_solve_time)
        delta_t.push_back(i[j].time_diff());

      oss << get_total(delta_t);

      if (j < d_sys_solve_time[0].size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    // write total sim time
    oss << "#\n";
    oss << "FinalSimTime " << sim_time << "\n";
  }

  return oss.str();
}