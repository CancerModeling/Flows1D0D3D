////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_LOGGER_H
#define UTILS_LOGGER_H

#include "utilLibs.hpp"

namespace {

float time_diff(std::chrono::steady_clock::time_point begin,
                std::chrono::steady_clock::time_point end) {

  return std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                               begin)
    .count();
}

int locate_in_set(const std::string &key, const std::vector<std::string> &set) {

  for (int i = 0; i < set.size(); i++)
    if (set[i] == key)
      return i;

  return -1;
}
} // namespace

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
 * @brief class to store simulation information at each step
 */
class TimeStepLog {

public:
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

  void init_ts(const std::vector<std::string> &sys_names) {
    d_sys_names = sys_names;
    d_n = d_sys_names.size();
    d_sys_assembly_time = std::vector<TimePair>(d_n, TimePair());
  }

  void ready_new_step(int step, double cur_time = 0.) {

    // check
    if (d_cur_step >= step or (d_cur_step >= 0 and d_cur_step < step - 1)) {

      libmesh_error_msg("Check TimeStepLog::ready_new_step. It should be "
                        "called once per time step.\n "
                        "Current step in ts: " +
                        std::to_string(d_cur_step) +
                        " step in model: " + std::to_string(step));
    }

    d_times.push_back(cur_time);
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

  void add_cur_time(const float &time) { d_times[d_cur_step] = time; }

  void add_solve_time(TimePair time_pair) { d_solve_time[d_cur_step] = time_pair; }

  void add_sys_solve_time(TimePair time_pair, const std::string &sys_name) {

    auto i = locate_in_set(sys_name, d_sys_names);
    if (i != -1)
      d_sys_solve_time[d_cur_step][i] = time_pair;
  }

  void add_sys_solve_time(const float &dt, const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(dt);
  }

  void add_sys_solve_time(const std::vector<steady_clock::time_point>
                            &sys_clocks,
                          const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(
      time_diff(sys_clocks[n], steady_clock::now()));
  }

  void add_sys_solve_time(const steady_clock::time_point
                            &sys_clock,
                          const int &n) {

    d_sys_solve_time[d_cur_step][n].add_dt(
      time_diff(sys_clock, steady_clock::now()));
  }

  void add_pres_solve_time(TimePair time_pair) { d_pres_solve_time[d_cur_step] = time_pair; }

  void add_update_net_time(TimePair time_pair) { d_update_network_time[d_cur_step] =
                                                   time_pair; }

  void add_nonlin_iter(const int &i) { d_nonlinear_iters[d_cur_step] = i; }

  void add_pres_nonlin_iter(const int &i) { d_pres_nonlinear_iters[d_cur_step] = i; }

  std::vector<float> get_delta_t_vec(const std::vector<TimePair> &list) {

    auto delta_t = std::vector<float>(list.size(), 0.);
    for (size_t i = 0; i < list.size(); i++)
      delta_t[i] = list[i].time_diff();

    return delta_t;
  }
};

/*!
 * @brief Prints message to std::cout and also writes to the file
 */
class Logger : public TimeStepLog {

public:
  Logger(const std::string &log_file, Parallel::Communicator *comm,
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

  // output time step log
  std::string log_ts_base(const int i, const int ns);
  std::string log_ts_base_final_avg(double sim_time, const int ns);

  void log_ts() {

    const auto i = d_cur_step;

    // print to file without any format
    if (d_comm_p->rank() == 0)
      d_ts_file << log_ts_base(i, 0) << std::flush;

    // print to screen with format
    if (i > 0) {
      log(log_ts_base(i, 2), "TS log");
      log(" \n");
    }
  }

  void log_ts_final(double sim_time) {

    // print to file without any format
    if (d_comm_p->rank() == 0)
      d_ts_file << log_ts_base_final_avg(sim_time, 0) << std::flush;

    // print to screen with format
    log(log_ts_base_final_avg(sim_time, 2), "TS log");
  }

  void log_ts_all() {
    for (int i = 0; i < d_solve_time.size(); i++) {
      if (d_comm_p->rank() == 0)
        d_ts_file << log_ts_base(i, 0) << std::flush;
    }
  }

  void log_qoi_header(const double &time, const std::vector<std::string> &qoi_names) {

    std::ostringstream oss;
    oss << "time ";
    for (unsigned int i = 0; i < qoi_names.size(); i++) {
      oss << qoi_names[i];
      if (i < qoi_names.size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    // log to file
    if (d_comm_p->rank() == 0)
      d_qoi_file << oss.str() << std::flush;

    // log to screen
    std::string str = "\n  QoI log header\n  " + oss.str();
    log(str, "debug");

    //log_qoi(time, qoi);
  }

  void log_qoi(const double &time, const std::vector<double> &qoi) {

    std::ostringstream oss;
    oss << time << " ";
    for (unsigned int i = 0; i < qoi.size(); i++) {
      oss << qoi[i];
      if (i < qoi.size() - 1)
        oss << " ";
      else
        oss << "\n";
    }

    // log to file
    if (d_comm_p->rank() == 0)
      d_qoi_file << oss.str() << std::flush;

    // log to screen
    std::string str = "  " + oss.str();
    log(str, "QoI");
    log(" \n");
  }

private:
  Parallel::Communicator *d_comm_p;
  bool d_screen_out;
  std::ofstream d_dbg_file;
  std::ofstream d_ts_file;
  std::ofstream d_qoi_file;
};

} // namespace util

#endif // UTILS_LOGGER_H
