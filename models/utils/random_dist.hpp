////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_RANDOM_DIST_H
#define UTILS_RANDOM_DIST_H

#include <boost/random.hpp>
#include <random>

#include "libmesh/parallel.h"

#include "utilIO.hpp"

// boost method
//typedef boost::mt19937 RandGen;
//typedef boost::lognormal_distribution<> LogNormalDist;
//typedef boost::uniform_real<> UniformDist;
//typedef boost::normal_distribution<> NormalDist;

typedef std::mt19937 RandGenerator;
typedef std::lognormal_distribution<> LogNormalDistribution;
typedef std::uniform_real_distribution<> UniformDistribution;
typedef std::normal_distribution<> NormalDistribution;

namespace util {

inline RandGenerator get_rd_gen(int &seed) {

  //return RandGenerator();

  if (seed < 0) {
    std::random_device rd;
    seed = rd();
  }

  return RandGenerator(seed);
}

inline std::default_random_engine get_rd_engine(int &seed) {

  //return std::default_random_engine();

  if (seed < 0) {
    std::random_device rd;
    seed = rd();
  }

  return std::default_random_engine(seed);
}

inline double transform_to_normal_dist(double mean, double std, double sample) {
  return std * sample + mean;
}

inline double transform_to_uniform_dist(double min, double max, double sample) {
  return min + sample * (max - min);
}


template<class T>
class DistributionSample {
public:
  DistributionSample(double arg1, double arg2, int seed = -1)
      : d_seed(seed), d_gen(get_rd_gen(seed)), d_dist(arg1, arg2){};

  DistributionSample() : d_seed(-1){};

  void init(double arg1, double arg2, int seed = -1) {

    d_seed = seed;
    d_gen = RandGenerator(get_rd_gen(seed));
    d_dist = T(arg1, arg2);
  }

  double operator()() {
    return d_dist(d_gen);
  }

  // not implemented for this
  void debug_out(const std::string &filename) {
    return;
  }

  int d_seed;
  RandGenerator d_gen;
  T d_dist;
};

template<class T>
class DistributionSampleParallel {
public:
  DistributionSampleParallel(double arg1, double arg2,
                             libMesh::Parallel::Communicator *comm,
                             int seed = -1,
                             unsigned int sample_size = 1000)
      : d_seed(seed), d_curSample(0), d_sampleSize(sample_size),
        d_gen(get_rd_gen(seed)), d_dist(arg1, arg2), d_procRank(0),
        d_procSize(0), d_comm_p(comm) {

    if (comm) {
      d_procRank = d_comm_p->rank();
      d_procSize = d_comm_p->size();
    }

    d_samples.reserve(d_sampleSize);
    if (d_procRank == 0)
      d_samples.resize(d_sampleSize);

    update();
  };

  DistributionSampleParallel()
      : d_seed(-1), d_curSample(0), d_sampleSize(0), d_procRank(0),
        d_procSize(0), d_comm_p(nullptr){};

  void init(double arg1, double arg2,
            libMesh::Parallel::Communicator *comm,
            int seed = -1,
            unsigned int sample_size = 1000) {

    d_seed = seed;
    d_curSample = 0;
    d_sampleSize = sample_size;
    d_gen = RandGenerator(get_rd_gen(seed));
    d_dist = T(arg1, arg2);
    d_comm_p = comm;

    if (comm) {
      d_procRank = d_comm_p->rank();
      d_procSize = d_comm_p->size();
    }

    d_samples.reserve(d_sampleSize);
    if (d_procRank == 0)
      d_samples.resize(d_sampleSize);

    update();
  }

  double operator()() {
    if (d_curSample < d_sampleSize) {
      return d_samples[d_curSample++];
    } else {
      update();
      return d_samples[d_curSample++];
    }
  }

  void update() {

    d_curSample = 0;
    if (d_procRank == 0) {
      for (unsigned int i = 0; i < d_sampleSize; i++)
        d_samples[i] = d_dist(d_gen);
    } else {
      d_samples.resize(0);
    }

    d_comm_p->allgather(d_samples);
    if (d_samples.size() != d_sampleSize)
      libmesh_error_msg("Error collecting samples");
  }

  void debug_out(const std::string &filename) {
    util::io::printFile(filename + "_" + std::to_string(d_procRank) + ".txt", d_samples, 0, " ");
  }

  int d_seed;
  unsigned int d_curSample;
  unsigned int d_sampleSize;
  RandGenerator d_gen;
  T d_dist;
  unsigned int d_procRank;
  unsigned int d_procSize;
  libMesh::Parallel::Communicator *d_comm_p;
  std::vector<double> d_samples;
};

} // namespace util

#endif // UTILS_RANDOM_DIST_H
