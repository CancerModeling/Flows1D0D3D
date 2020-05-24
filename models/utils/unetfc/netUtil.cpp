////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "netUtil.hpp"
#include "network.hpp"
#include <random>

void util::unetfc::angle_correction(const Point &parent_d, Point &child_d,
                      const double &max_angle) {

  auto child_angle = util::angle(child_d, parent_d);
  if (std::abs(child_angle) > max_angle) {

    // axis for rotation
    Point axis = parent_d.cross(child_d);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    child_d = util::rotate(parent_d, max_angle, axis);
  }

}

void util::unetfc::compute_bifurcate_child_direction(const Point &parent_d,
                                       const Point &child_d, Point &child_d_2,
                                       const double &branch_angle) {

  // get angle between child direction and parent direction and offset
  // it by branching angle
  double angle_1 = util::angle(parent_d, child_d) + branch_angle;

  // axis for rotation
  Point axis = parent_d.cross(child_d);
  axis = axis / axis.norm();

  // rotate parent_direction by -ve angle
  child_d_2 = util::rotate(parent_d, -angle_1, axis);

  out << "\n angle_1: " << angle_1 << ", axis: " << axis
      << ", child_d_2: " << child_d_2 << "\n";
}

void util::unetfc::angle_correction_bifurcation(const Point &parent_d, Point &child_d,
                                  Point &child_d_2, const double &max_angle,
                                  const double &branch_angle) {

  // check if angle of direction from parent direction is within
  // permissible range
  auto child_angle = util::angle(child_d_2, parent_d);
  if (std::abs(child_angle) > max_angle) {

    // axis for rotation
    Point axis = parent_d.cross(child_d_2);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    child_d_2 = util::rotate(parent_d, max_angle, axis);

  } else if (std::abs(child_angle) < branch_angle) {

    // axis for rotation
    Point axis = parent_d.cross(child_d_2);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    child_d_2 = util::rotate(parent_d, branch_angle, axis);
  }
}

std::vector<double> util::unetfc::getElementCenter(int i, int j, int k, double h_3D) {

  std::vector<double> center;

  center.push_back((double)i * h_3D + 0.5 * h_3D);

  center.push_back((double)j * h_3D + 0.5 * h_3D);

  center.push_back((double)k * h_3D + 0.5 * h_3D);

  return center;
}

std::vector<double> util::unetfc::getCenterNeighbor(std::vector<double> center,
                                      std::vector<double> direction,
                                      double h_3D) {

  std::vector<double> center_neighbor;

  for (int i = 0; i < 3; i++) {

    center_neighbor.push_back(center[i] + h_3D * direction[i]);
  }

  return center_neighbor;
}

std::vector<std::vector<double>> util::unetfc::defineDirections(){

                                 std::vector<std::vector<double>> directions;

				 for(int i = 0; i < 3; i++){

				     std::vector<double> direction_p;
				     std::vector<double> direction_m;

				     for(int j = 0; j < 3; j++){

				         if(i == j){

				            direction_p.push_back(1.0);
				            direction_m.push_back(-1.0);

				         } 
                                         else{

				            direction_p.push_back(0.0);
				            direction_m.push_back(0.0);
				         }

				     }

				     directions.push_back(direction_p);
				     directions.push_back(direction_m);

				}

				return directions;

}


std::vector<std::vector<double>> util::unetfc::defineDirectionsNeighboring(){

                                 std::vector<std::vector<double>> directions;

				 for(int i = 0; i < 3; i++){

				     std::vector<double> direction_p;
				     std::vector<double> direction_m;

				     for(int j = 0; j < 3; j++){

				         if(i == j){

				            direction_p.push_back(1.0);
				            direction_m.push_back(-1.0);

				         } 
                                         else{

				            direction_p.push_back(0.0);
				            direction_m.push_back(0.0);

				         }

				     }

				     directions.push_back(direction_p);
				     directions.push_back(direction_m);

				}

				for(int i = 0; i < 3; i++){

				     std::vector<double> direction_p;
				     std::vector<double> direction_m;

				     for(int j = 0; j < 3; j++){

				         if(i != j){

				            direction_p.push_back(1.0);
				            direction_m.push_back(-1.0);

				         } 
                                         else{

				            direction_p.push_back(0.0);
				            direction_m.push_back(0.0);

				         }

				     }

				     directions.push_back(direction_p);
				     directions.push_back(direction_m);

				}


				for(int i = 0; i < 3; i++){

				     std::vector<double> direction_p;
				     std::vector<double> direction_m;

				     for(int j = 0; j < 3; j++){

				         if(i == j){

				            direction_p.push_back(0.0);
				            direction_m.push_back(-1.0);

				         } 
                                         else{

				            direction_p.push_back(1.0);
				            direction_m.push_back(1.0);

				         }

				     }

				     directions.push_back(direction_p);
				     directions.push_back(direction_m);

				}				

				for(int i = 0; i < 3; i++){

				     std::vector<double> direction_p;
				     std::vector<double> direction_m;

				     for(int j = 0; j < 3; j++){

				         if(i != j){

				            direction_p.push_back(0.0);
				            direction_m.push_back(-1.0);

				         } 
                                         else{

				            direction_p.push_back(1.0);
				            direction_m.push_back(1.0);

				         }

				     }

				     directions.push_back(direction_p);
				     directions.push_back(direction_m);

				}

                                return directions;

}

bool util::unetfc::isCenterInDomain(std::vector<double> center, double L_x) {

  bool isInDomain = true;

  for (int i = 0; i < 3; i++) {

    if (center[i] < 0.0 || center[i] > L_x) {

      isInDomain = false;

      return isInDomain;
    }
  }

  return isInDomain;
}

int util::unetfc::getElementIndex(std::vector<double> center, double h_3D, int N_3D){

  int index = 0;

  index = index + std::floor(center[0] / h_3D);

  index = index + std::floor(center[1] / h_3D) * N_3D;

  index = index + std::floor(center[2] / h_3D) * N_3D * N_3D;

  // checking this error is costly so later we may want to remove this check
  if (index < 0) {
    std::ostringstream oss;
    oss << "index: " << index << ", loc: (" << util::io::printStr(center) <<
    ")";
    libmesh_error_msg("Index should not be negative. Check center location as"
                      " it may be outside the domain boundary.\n"
                      + oss.str());
  }

  return index;
}

std::vector<double> util::unetfc::getCenterFace(std::vector<double> center, std::vector<double> direction, double h_3D) {

  std::vector<double> center_face;

  for (int i = 0; i < 3; i++) {

    center_face.push_back(center[i] + 0.5 * h_3D * direction[i]);
  }

  return center_face;
}

double util::unetfc::getDirichletValue(std::vector<double> center_face, double L_p, double radius) {

  double dirichlet_value = 0.0;

  double dist = (center_face[0] - 0.5) * (center_face[0] - 0.5) +
                (center_face[1] - 0.5) * (center_face[1] - 0.5);

  dist = std::sqrt(dist);

  if (dist < radius) {

    dirichlet_value = (1.0 + center_face[2]) * L_p / (1.0 + L_p);

  } else {

    dirichlet_value = (1.0 + center_face[2]) * L_p / (1.0 + L_p) *
                      (1.0 - radius * std::log(dist / radius));
  }

  return dirichlet_value;
}

int util::unetfc::getIndex(double value, double h_3D) {

  int index = 0;

  while (index * h_3D < value - 1.0e-8) {

    index = index + 1;
  }

  return index;
}

double util::unetfc::normVector(std::vector<double> &vec) {

  double norm_vec = 0.0;

  for (int j = 0; j < 3; j++) {

    norm_vec = norm_vec + vec[j] * vec[j];
  }

  norm_vec = std::sqrt(norm_vec);

  for (int j = 0; j < 3; j++) {

    vec[j] = vec[j] / norm_vec;
  }

  return norm_vec;
}

std::vector<double> util::unetfc::determineRotator(std::vector<double> dir) {

  std::vector<double> rotator;

  if (std::abs(dir[0]) < 1.0e-13 && std::abs(dir[1]) < 1.0e-13) {

    rotator.push_back(1.0);

    rotator.push_back(0.0);

    rotator.push_back(0.0);

  } else {

    rotator.push_back(-dir[1]);

    rotator.push_back(dir[0]);

    rotator.push_back(0.0);
  }

  double length_rotator = normVector( rotator );   

  for(int j=0;j<3;j++){ 
                     
      rotator[ j ] = rotator[ j ]/length_rotator;

  }

  return rotator;

}

std::vector<double> util::unetfc::computeNodesOnCylinders(std::vector<double> dir,
                                            std::vector<double> rotator,
                                            std::vector<double> midpoint,
                                            double radius, double theta) {

  std::vector<double> cylinder_node(3);

  double c = std::cos(theta);

  double s = std::sin(theta);

  gmm::row_matrix<std::vector<double>> R(3, 3);

  R(0, 0) = dir[0] * dir[0] * (1.0 - c) + c;
  R(1, 1) = dir[1] * dir[1] * (1.0 - c) + c;
  R(2, 2) = dir[2] * dir[2] * (1.0 - c) + c;
  R(0, 1) = dir[0] * dir[1] * (1.0 - c) - dir[2] * s;
  R(0, 2) = dir[0] * dir[2] * (1.0 - c) + dir[1] * s;
  R(1, 0) = dir[0] * dir[1] * (1.0 - c) + dir[2] * s;
  R(2, 0) = dir[0] * dir[2] * (1.0 - c) - dir[1] * s;
  R(1, 2) = dir[1] * dir[2] * (1.0 - c) - dir[0] * s;
  R(2, 1) = dir[1] * dir[2] * (1.0 - c) + dir[0] * s;

  gmm::mult(R, rotator, cylinder_node);

  gmm::scale(cylinder_node, radius);

  gmm::add(midpoint, cylinder_node);

  return cylinder_node;
}

void util::unetfc::updateWeightsAndIds(int N_s, int N_theta, int elementIndex,
                                 std::vector<double> &weights,
                                 std::vector<int> &id_3D_elements) {

  bool elementfound = false;

  int indexfound = 0;

  int numberOfElements = id_3D_elements.size();

  for (int i = 0; i < numberOfElements; i++) {

    if (elementIndex == id_3D_elements[i]) {

      indexfound = i;

      elementfound = true;
      break;
    }
  }

  if (elementfound == true) {

    weights[indexfound] = weights[indexfound] + 1.0 / (N_s * N_theta);

  } else {

    id_3D_elements.push_back(elementIndex);
    weights.push_back(1.0 / (N_s * N_theta));
  }
}

void util::unetfc::determineWeightsAndIds(int N_s, int N_theta, int N_3D, std::vector<double> coord, std::vector<double> coord_neighbor,
                                    double radius, double h_3D, double& length_edge, std::vector<double> &weights, std::vector<int> &id_3D_elements){

     std::vector<double> direction, rotator;

     for(int j=0;j<3;j++){
  
         direction.push_back( coord_neighbor[ j ] - coord[ j ] );

     }

     length_edge = normVector( direction );

     rotator = determineRotator( direction );

     double length_rotator = normVector( rotator );

     for(int i_s=1;i_s<N_s-1;i_s++){

         std::vector<double> midpoint(3);

         double theta = 0.0;

         for(int j=0;j<3;j++){

             midpoint[ j ] = coord[ j ] + ( ( (double) i_s/(double) N_s ) * (0.5*length_edge) ) * direction[ j ];

         }

         for(int i_theta=0;i_theta<N_theta;i_theta++){

             theta = ((double) i_theta)/((double) N_theta)*2.0*M_PI;

             std::vector<double> cylinder_node = computeNodesOnCylinders( direction,  rotator, midpoint,  radius, theta );

             int elementIndex = getElementIndex( cylinder_node, h_3D, N_3D );

             // Compute weights and element ids
             updateWeightsAndIds( N_s-2, N_theta, elementIndex, weights, id_3D_elements );

         }

     }

}

std::vector<double> util::unetfc::getCenterFromIndex( int index, int N_3D, double h_3D ){

                    std::vector<double> center = std::vector<double>(3,0.0);

                    int k = index/( N_3D*N_3D );

                    index = index - ( k*N_3D*N_3D );

                    int j = index/N_3D;

                    int i = index - ( j*N_3D );

                    center[ 0 ] = (double)i * h_3D + 0.5 * h_3D;
                    center[ 1 ] = (double)j * h_3D + 0.5 * h_3D;
                    center[ 2 ] = (double)k * h_3D + 0.5 * h_3D; 

                    return center;

}

std::vector<int> util::unetfc::getNeighboringElementIndices( int index, int N_3D, double h_3D, double L_x ){

                 std::vector<int> indicesNeighbors;

                 std::vector<int> directions;

                 directions.push_back( -2 );
                 directions.push_back(  2 );

                 directions.push_back( -1 );
                 directions.push_back(  1 );

                 directions.push_back( -N_3D );
                 directions.push_back(  N_3D );

                 directions.push_back( -2*N_3D );
                 directions.push_back(  2*N_3D );

                 directions.push_back(  N_3D-1 );
                 directions.push_back(  N_3D+1 );

                 directions.push_back(  2*N_3D-1 );
                 directions.push_back(  2*N_3D+1 );

                 directions.push_back( -2*N_3D-1 );
                 directions.push_back(  2*N_3D+1 );

                 directions.push_back( -N_3D-1 );
                 directions.push_back( -N_3D+1 );

                 directions.push_back( -1+( N_3D*N_3D ) );
                 directions.push_back(  1+( N_3D*N_3D ) );

                 directions.push_back( -1+( 2*N_3D*N_3D ) );
                 directions.push_back(  1+( 2*N_3D*N_3D ) );

                 directions.push_back( -N_3D+( N_3D*N_3D ) );
                 directions.push_back(  N_3D+( N_3D*N_3D ) );

                 directions.push_back( -N_3D+( 2*N_3D*N_3D ) );
                 directions.push_back(  N_3D+( 2*N_3D*N_3D ) );

                 directions.push_back(  N_3D-1+( N_3D*N_3D ) );
                 directions.push_back(  N_3D+1+( N_3D*N_3D ) );

                 directions.push_back(  N_3D-1+2*( N_3D*N_3D ) );
                 directions.push_back(  N_3D+1+2*( N_3D*N_3D ) );

                 directions.push_back( -N_3D-1+( N_3D*N_3D ) );
                 directions.push_back( -N_3D+1+( N_3D*N_3D ) );

                 directions.push_back( -N_3D-1+( 2*N_3D*N_3D ) );
                 directions.push_back( -N_3D+1+( 2*N_3D*N_3D ) );

                 directions.push_back( -1-( N_3D*N_3D ) );
                 directions.push_back(  1-( N_3D*N_3D ) );

                 directions.push_back( -1-( 2*N_3D*N_3D ) );
                 directions.push_back(  1-( 2*N_3D*N_3D ) );

                 directions.push_back( -N_3D-( N_3D*N_3D ) );
                 directions.push_back(  N_3D-( N_3D*N_3D ) );

                 directions.push_back(  N_3D-1-( N_3D*N_3D ) );
                 directions.push_back(  N_3D+1-( N_3D*N_3D ) );

                 directions.push_back( -N_3D-1-( N_3D*N_3D ) );
                 directions.push_back( -N_3D+1-( N_3D*N_3D ) );

                 for(int i=0;i<directions.size();i++){

                     int index_neighbor = index+directions[ i ];

                     std::vector<double> center = getCenterFromIndex( index_neighbor, N_3D, h_3D );

                     bool isInDomain = true;

                     for(int j=0;j<3;j++){

                         if( center[ i ]<0.0 || center[ i ]>L_x ){

                             isInDomain = false;

                             break;
                             
                         }

                     }

                     indicesNeighbors.push_back( index_neighbor );

                 }

                 return indicesNeighbors;

}
