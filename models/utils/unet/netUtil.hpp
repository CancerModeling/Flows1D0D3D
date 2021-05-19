////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNET_NETUTIL_H
#define UTIL_UNET_NETUTIL_H

#include "utils.hpp"

namespace util {

namespace unet {

void angle_correction(const Point &parent_d, Point &child_d,
                      const double &max_angle);

void compute_bifurcate_child_direction(const Point &parent_d,
                                       const Point &child_d, Point &child_d_2,
                                       const double &branch_angle);

void angle_correction_bifurcation(const Point &parent_d, Point &child_d,
                                  Point &child_d_2, const double &max_angle,
                                  const double &branch_angle);

std::vector<double> getElementCenter(int i, int j, int k, double h_3D);

std::vector<double> getCenterNeighbor(std::vector<double> center,
                                      std::vector<double> direction,
                                      double h_3D);

std::vector<std::vector<double>> defineDirections();

std::vector<double> getCenterFromIndex(int index, int N_3D, double h_3D);

bool isCenterInDomain(std::vector<double> center, double L_x);

int getElementIndex(std::vector<double> center, double h_3D, int N_3D);

std::vector<double> getCenterFace(std::vector<double> center,
                                  std::vector<double> direction, double h_3D);

double getDirichletValue(std::vector<double> center_face, double L_p,
                         double radius);

int getIndex(double value, double h_3D);

/*! @brief Returns the norm of a vector AND normalizes it. */
double normVector(std::vector<double> &vec);

/*! @brief Determines an arbitrary vector orthogonal to dir. */
std::vector<double> determineRotator(std::vector<double> dir);

/*! @brief Rotates the "rotator"-vector around the axis given by "dir" by theta radians, scales it with the given radius. */
std::vector<double> computeNodesOnCylinders(const std::vector<double>& dir,
                                            const std::vector<double>& rotator,
                                            const std::vector<double>& midpoint,
                                            double radius, double theta);

void updateWeightsAndIds(int N_s, int N_theta, int elementIndex,
                         std::vector<double> &weights,
                         std::vector<int> &id_3D_elements);

void determineWeightsAndIds(int N_s, int N_theta, int N_3D,
                            std::vector<double> coord,
                            std::vector<double> coord_neighbor, double radius,
                            double h_3D, double length_edge,
                            std::vector<double> &weights,
                            std::vector<int> &id_3D_elements);

void determineWeightsAndIds(int N_s, int N_theta, int N_3D,
                            std::vector<double> coord,
                            std::vector<double> coord_neighbor, double radius,
                            double h_3D, double length_edge,
                            std::vector<double> &weights,
                            std::vector<int> &id_3D_elements,
                            const int &integration_method,
                            const MeshBase &mesh,
                            bool check_elem_owner = false);

void determineWeightsAndIdsLineSource(int N_s, int N_theta, int N_3D,
                                      std::vector<double> coord,
                                      std::vector<double> coord_neighbor, double radius,
                                      double h_3D, double length_edge,
                                      std::vector<double> &weights,
                                      std::vector<int> &id_3D_elements,
                                      const int &integration_method,
                                      const MeshBase &mesh,
                                      bool check_elem_owner = false);

std::vector<int> getNeighboringElementIndices(int index, int N_3D, double h_3D,
                                              double L_x);

} // namespace unet

} // namespace util

#endif // UTIL_UNET_NETUTIL_H
