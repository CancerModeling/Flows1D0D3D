////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_LOGGER_H
#define UTILS_LOGGER_H

#include "utilLibs.hpp"
#include "aixlog.hpp"

namespace {

float time_diff(std::chrono::steady_clock::time_point begin,
                       std::chrono::steady_clock::time_point end) {

  return std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                               begin).count();
}

int locate_in_set(const std::string &key, const std::vector<std::string> &set) {

  for (int i =0; i<set.size(); i++)
    if (set[i] == key)
      return i;

  return -1;
}
}

namespace util {

// easiy to use name to hold pair of time
struct TimePair {
  double d_dt;
  bool d_use_dt;
  steady_clock::time_point d_begin;
  steady_clock::time_point d_end;

  TimePair()
      : d_dt(0.), d_use_dt(false),
        d_begin(steady_clock::now()), d_end(steady_clock::now()) {}

  TimePair(steady_clock::time_point begin, steady_clock::time_point end)
      : d_dt(0.), d_use_dt(false),
        d_begin(begin), d_end(end) {}

  float time_diff() const {

    if (d_use_dt)
      return d_dt;
    else
      return std::chrono::duration_cast<std::chrono::microseconds>(d_end -
                                                                 d_begin)
        .count();
  }

  void add_dt(const float &dt) {
    if (!d_use_dt)
      d_use_dt = true;
    d_dt += dt;
  }
};

/*!
 * @brief Struct to store simulation information at each step
 */
struct TimeStepLog {

  /*! @brief Total systems */
  int d_n;

  /*! @brief Current step counter */
  int d_cur_step;

  /*! @brief Initial setup time */
  TimePair d_setup_time;

  /*! @brief  Each systems assembly time. Measured at the first step */
  std::vector<std::string> d_sys_names;

  /*! @brief  Each systems assembly time. Measured at the first step */
  std::vector<TimePair> d_sys_assembly_time;

  /*! @brief  Current time at this step */
  std::vector<float> d_times;

  /*! @brief  Total solve time per time step. Includes assembly, nonlinear
   * iterations, etc */
  std::vector<TimePair> d_solve_time;

  /*! @brief  Each systems total solve time per time step. Includes assembly,
   * nonlinear
   * iterations, etc */
  std::vector<std::vector<TimePair>> d_sys_solve_time;

  /*! @brief  Total pressure solve time per time step. Includes assembly,
   * nonlinear
   * iterations, etc */
  std::vector<TimePair> d_pres_solve_time;

  /*! @brief Time to update the network */
  std::vector<TimePair> d_update_network_time;

  /*! @brief  number of nonlinear iteration for the system */
  std::vector<int> d_nonlinear_iters;

  /*! @brief  number of nonlinear iteration for the pressure */
  std::vector<int> d_pres_nonlinear_iters;


  TimeStepLog() : d_n(0), d_cur_step(-1), d_setup_time(TimePair()) {}

  TimeStepLog(const std::vector<std::string> &sys_names)
      : d_n(sys_names.size()), d_cur_step(-1), d_setup_time(TimePair()),
        d_sys_names(sys_names),
        d_sys_assembly_time(std::vector<TimePair>(d_n, TimePair())) {}

  void init(const std::vector<std::string> &sys_names) {
    d_sys_names = sys_names;
    d_n = d_sys_names.size();
    d_sys_assembly_time = std::vector<TimePair>(d_n, TimePair());
  }

  void ready_new_step(int step) {

    // check
    if (d_cur_step >= step or (d_cur_step >= 0 and d_cur_step < step - 1)) {

      libmesh_error_msg("Check TimeStepLog::ready_new_step. It should be "
                        "called once per time step.\n "
                        "Current step in ts: " +
                        std::to_string(d_cur_step) +
                        " step in model: " + std::to_string(step));
    }

    d_times.push_back(0.);
    d_cur_step = int(d_times.size()) - 1;

    d_solve_time.emplace_back();
    d_sys_solve_time.emplace_back(d_n, TimePair());
    d_pres_solve_time.emplace_back();
    d_update_network_time.emplace_back();
    d_nonlinear_iters.push_back(0);
    d_pres_nonlinear_iters.push_back(0);
  }

  void add_sys_assembly_time(TimePair time_pair, const std::string &sys_name) {

    auto i = locate_in_set(sys_name, d_sys_names);
    if (i != -1)
      d_sys_assembly_time[i] = time_pair;
  }

  void add_cur_time(const float &time) {d_times[d_cur_step] = time;}

  void add_solve_time(TimePair time_pair) {d_solve_time[d_cur_step] = time_pair;}

  void add_sys_solve_time(TimePair time_pair, const std::string &sys_name) {

    auto i = locate_in_set(sys_name, d_sys_names);
    if (i != -1)
      d_sys_solve_time[d_cur_step][i] = time_pair;
  }

  void add_sys_solve_time(const float &dt, const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(dt);
  }

  void add_sys_solve_time(const std::vector<steady_clock::time_point>
      &sys_clocks, const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(
        time_diff(sys_clocks[n], steady_clock::now()));
  }

  void add_sys_solve_time(const steady_clock::time_point
                          &sys_clock, const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(
        time_diff(sys_clock, steady_clock::now()));
  }

  void add_pres_solve_time(TimePair time_pair) {d_pres_solve_time[d_cur_step]
                                                = time_pair;}

  void add_update_net_time(TimePair time_pair) {d_update_network_time[d_cur_step] =
                                                    time_pair;}

  void add_nonlin_iter(const int &i) {d_nonlinear_iters[d_cur_step] = i;}

  void add_pres_nonlin_iter(const int &i) {d_pres_nonlinear_iters[d_cur_step] = i;}
};

/*!
 * @brief Prints message to std::cout and also writes to the file
 */
class Logger {

public:

  Logger() : d_ts(), d_comm_p(nullptr), d_screen_out(true) {};

  Logger(const std::string &log_file, Parallel::Communicator *comm,
         bool screen_out = true)
      : d_ts(), d_comm_p(comm), d_screen_out(screen_out) {

    if (d_dbg_file.is_open())
      d_dbg_file.close();
    d_dbg_file = std::ofstream(log_file + "_dbg.log", std::ios_base::out);

    if (d_ts_file.is_open())
      d_ts_file.close();
    d_ts_file = std::ofstream(log_file + "_time.log", std::ios_base::out);
  };

  // overload function call
  void operator()(std::ostringstream &oss, bool screen_out) {
    log(oss, screen_out);
  }
  void operator()(std::ostringstream &oss) {
    log(oss, d_screen_out);
  }

  void operator()(const std::string &str, bool screen_out) {
    log(str, screen_out);
  }
  void operator()(const std::string &str) {
    log(str, d_screen_out);
  }

  void log(std::ostringstream &oss, bool screen_out) {

    if (screen_out)
      out << oss.str();
    if (d_comm_p->rank() == 0)
      d_dbg_file << oss.str();

    oss.str("");
    oss.clear();
  };

  void log(std::ostringstream &oss) {
    log(oss, d_screen_out);
  };

  void log(const std::string &str, bool screen_out) {

    if (screen_out)
      out << str;
    if (d_comm_p->rank() == 0)
      d_dbg_file << str;
  };

  void log(const std::string &str) {

    log(str, d_screen_out);
  };


  std::string log_ts_step(const int i) {

    std::ostringstream oss;

    if (d_comm_p->rank() == 0) {

      // initial stuff
      if (i == 0) {
        oss << "# TimeStepLog\n";
        oss << "Setup_Time " << d_ts.d_setup_time.time_diff() << "\n";

        oss << "#\n";
        oss << "SysNames\n";
        for (unsigned int j=0; j < d_ts.d_n; j++)
          oss << d_ts.d_sys_names[j] << " " << j << "\n";

        oss << "#\n";
        oss << "SysAssemblyTimes\n";
        for (unsigned int j=0; j < d_ts.d_n; j++) {

          oss << d_ts.d_sys_assembly_time[j].time_diff();

          if (j < d_ts.d_n - 1)
            oss << ", ";
          else
            oss << "\n";
        }

        oss << "#\n";

        // add header for next step
        oss << "StepLog\n";
        oss << "T, ST, PST, NETT, NLI, PNLI, ";

        for (unsigned int j=0; j < d_ts.d_n; j++) {

          oss << "S" + std::to_string(j);

          if (j < d_ts.d_sys_solve_time[i].size() - 1)
            oss << ", ";
          else
            oss << "\n";
        }
      }

      oss << d_ts.d_times[i] << ", "
          << d_ts.d_solve_time[i].time_diff() << ", "
          << d_ts.d_pres_solve_time[i].time_diff() << ", "
          << d_ts.d_update_network_time[i].time_diff() << ", "
          << d_ts.d_nonlinear_iters[i] << ", "
          << d_ts.d_pres_nonlinear_iters[i] << ", ";

      for (unsigned int j=0; j < d_ts.d_sys_solve_time[i].size(); j++) {

        oss << d_ts.d_sys_solve_time[i][j].time_diff();

        if (j < d_ts.d_sys_solve_time[i].size() - 1)
          oss << ", ";
        else
          oss << "\n";
      }

      d_ts_file << oss.str();
    }

    return oss.str();
  }

  void log_ts() {

    const auto i = d_ts.d_cur_step;
    std::string str = "\n\n___________\nTimeStepLog\n  " + log_ts_step(i);

    if (i > 0)
      log(str);
  }

  void log_ts_all() {
    for (int i=0; i<d_ts.d_solve_time.size(); i++)
      log_ts_step(i);
  }

  ~Logger() {

    if (d_dbg_file.is_open())
      d_dbg_file.close();

    if (d_ts_file.is_open())
      d_ts_file.close();
  }

  void init(const std::string &log_file, Parallel::Communicator *comm) {
    d_comm_p = comm;
    d_dbg_file.open(log_file, std::ios_base::out);
  }

public:
  TimeStepLog d_ts;

private:
  Parallel::Communicator *d_comm_p;
  bool d_screen_out;
  std::ofstream d_dbg_file;
  std::ofstream d_ts_file;
};

} // namespace util

#endif // UTILS_LOGGER_H
