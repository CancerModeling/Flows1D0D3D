////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef TUMORMODELS_NETUTIL_H
#define TUMORMODELS_NETUTIL_H

#include "utilLibs.hpp"
#include "utils.hpp"

namespace netfc {

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

std::vector<std::vector<double>> defineDirectionsNeighboring();

std::vector<double> getCenterFromIndex( int index, int N_3D, double h_3D );

bool isCenterInDomain(std::vector<double> center, double L_x);

int getElementIndex(std::vector<double> center, double h_3D, int N_3D);

std::vector<double> getCenterFace(std::vector<double> center,
                                  std::vector<double> direction, double h_3D);

double getDirichletValue(std::vector<double> center_face, double L_p,
                         double radius);

int getIndex(double value, double h_3D);

double normVector(std::vector<double> &vec);

std::vector<double> determineRotator(std::vector<double> dir);

std::vector<double> computeNodesOnCylinders(std::vector<double> dir,
                                            std::vector<double> rotator,
                                            std::vector<double> midpoint,
                                            double radius, double theta);

void updateWeightsAndIds(int N_s, int N_theta, int elementIndex,
                         std::vector<double> &weights,
                         std::vector<int> &id_3D_elements);

void determineWeightsAndIds(int N_s, int N_theta, int N_3D, std::vector<double> coord, std::vector<double> coord_neighbor, 
                            double radius, double h_3D, double length_edge, std::vector<double> &weights, std::vector<int> &id_3D_elements);

std::vector<int> getNeighboringElementIndices( int index, int N_3D, double h_3D, double L_x );

}

#endif // TUMORMODELS_NETUTIL_H

