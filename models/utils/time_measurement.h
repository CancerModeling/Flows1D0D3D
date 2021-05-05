#ifndef TUMORMODELS_TIME_MEASUREMENT_H
#define TUMORMODELS_TIME_MEASUREMENT_H

#include <iostream>
#include <string>
#include <chrono>

namespace util {

class TimeMeasurement {
public:
  TimeMeasurement(const std::string &name)
      : d_name(name),
        d_begin(std::chrono::steady_clock::now()) {
    for (std::size_t i = 0; i < indent; i += 1)
      std::cout << " ";
    std::cout << "start " << d_name << std::endl;

    indent += 1;
  }

  ~TimeMeasurement() {
    indent -= 1;

    for (std::size_t i = 0; i < indent; i += 1)
      std::cout << " ";

    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<float> elapsed = end - d_begin;

    std::cout << "end " << d_name << " (took " << std::scientific << elapsed.count() << "s)" << std::endl;
  }

  TimeMeasurement(const TimeMeasurement &) = delete;
  void operator=(const TimeMeasurement &) = delete;
  TimeMeasurement(const TimeMeasurement &&) = delete;
  void operator=(const TimeMeasurement &&) = delete;

private:
  std::string d_name;
  std::chrono::steady_clock::time_point d_begin;

  static std::size_t indent;
};

std::size_t TimeMeasurement::indent = 0;

}

#endif //TUMORMODELS_TIME_MEASUREMENT_H
