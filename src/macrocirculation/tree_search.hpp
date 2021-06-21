////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TREE_SEARCH_HPP
#define TREE_SEARCH_HPP

#include <vector>
#include "utils.hpp"
#include "tree_search_base.hpp"

namespace macrocirculation {

/*!
 * @brief A class for nearest neighbor search
 */
class BaseNSearch {

public:
  /*!
   * @brief Constructor
   */
  BaseNSearch(std::string name, size_t debug = 0)
    : d_debug(debug), d_tree_type(name) {}

  virtual double update_point_cloud(const std::vector<TreePt> &x,
                                  bool parallel = true) = 0;

  virtual double set_input_cloud() = 0;

  virtual size_t nearest_search(
    const TreePt &searchPoint, const size_t &num_elems,
    std::vector<size_t> &neighs,
    std::vector<double> &sqr_dist) = 0;

  virtual size_t radius_search(
    const TreePt &searchPoint, const double &search_r,
    std::vector<size_t> &neighs,
    std::vector<double> &sqr_dist) = 0;

public:
  /*@brief control the verbosity */
  size_t d_debug;

  /*@brief name of tree: nflann_kdtree */
  std::string d_tree_type;
};

/*!
 * @brief A class for nearest neighbor search using nanoflann library
 */
class NFlannSearchKd : public BaseNSearch {

public:
  /*!
   * @brief Constructor
   */
  NFlannSearchKd(const PointCloud &x, size_t debug = 0, double tree_resolution = 1.)
    : BaseNSearch("nflann_kdtree", debug), d_cloud(x), d_tree(3, d_cloud,
                                                              nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)) {
    d_params.sorted = false;
  }

  double set_input_cloud() override {
    auto t1 = std::chrono::steady_clock::now();
    d_tree.buildIndex();
    auto t2 = std::chrono::steady_clock::now();
    return time_diff(t1, t2);
  }

  double update_point_cloud(const std::vector<TreePt> &x,
                          bool parallel = true) override {
    return 0;
  }

  size_t nearest_search(
    const TreePt &searchPoint, const size_t &num_elems,
    std::vector<size_t> &neighs,
    std::vector<double> &sqr_dist) override {

    double query_pt[3] = {searchPoint(0), searchPoint(1), searchPoint(2)};

    neighs.resize(num_elems);
    sqr_dist.resize(num_elems);
    auto num_elems_act = d_tree.knnSearch(&query_pt[0], num_elems, &neighs[0], &sqr_dist[0]);
    neighs.resize(num_elems_act);
    sqr_dist.resize(num_elems_act);
    return num_elems_act;
  }

  size_t radius_search(
    const TreePt &searchPoint, const double &search_r,
    std::vector<size_t> &neighs,
    std::vector<double> &sqr_dist) override {

    double query_pt[3] = {searchPoint(0), searchPoint(1), searchPoint(2)};

    neighs.clear();
    sqr_dist.clear();
    TreeSearchRes resultSet(search_r * search_r, neighs, sqr_dist);
    return d_tree.radiusSearchCustomCallback(&query_pt[0], resultSet, d_params);
  }

public:
  /*! @brief coordinates of the points */
  PointCloudAdaptor d_cloud;

  /*! @brief Tree */
  NFlannKdTree d_tree;

  /*! @brief Tree search parameters */
  nanoflann::SearchParams d_params;
};

} // namespace macrocirculation

#endif //TREE_SEARCH_HPP
