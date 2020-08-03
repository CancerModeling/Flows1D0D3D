////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_RANDOM_DIST_H
#define UTILS_RANDOM_DIST_H

#include <random>
#include <boost/random.hpp>

// boost method
//typedef boost::mt19937 RandGen;
//typedef boost::lognormal_distribution<> LogNormalDist;
//typedef boost::uniform_real<> UniformDist;
//typedef boost::normal_distribution<> NormalDist;

typedef std::mt19937 RandGen;
typedef std::lognormal_distribution<> LogNormalDist;
typedef std::uniform_real_distribution<> UniformDist;
typedef std::normal_distribution<> NormalDist;

namespace util {

inline RandGen get_rgen(int seed) {
  if (seed < 0) {
    std::random_device rd;
    return RandGen(rd());
  } else {
    return RandGen(seed);
  }
}

class DistLogNormal {
public:
  DistLogNormal(double mean, double std, int seed = -1) : d_seed(seed), d_gen(get_rgen(seed)), d_dist(mean, std) {};

  double operator ()() {
    return d_dist(d_gen);
  }

  int d_seed;
  RandGen d_gen;
  LogNormalDist d_dist;
};

class DistNormal {
public:

  DistNormal(double mean, double std, int seed = -1) : d_seed(seed), d_gen(get_rgen(seed)), d_dist(mean, std) {}

  double operator ()() {
    return d_dist(d_gen);
  }

  int d_seed;
  RandGen d_gen;
  NormalDist d_dist;
};

class DistUniform {
public:

  DistUniform(double min, double max, int seed = -1) : d_seed(seed), d_gen(get_rgen(seed)), d_dist(min, max) {}

  double operator ()() {
    return d_dist(d_gen);
  }

  int d_seed;
  RandGen d_gen;
  UniformDist d_dist;
};

} // namespace util

#endif // UTILS_RANDOM_DIST_H
