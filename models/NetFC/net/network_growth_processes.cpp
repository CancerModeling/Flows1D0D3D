////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "netUtil.hpp"
#include "network.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

void netfc::Network::updateNetwork() {

  std::cout << " " << std::endl;
  std::cout << "Update the network" << std::endl;

  int numberOfNodes = VGM.getNumberOfNodes();

  std::cout << " " << std::endl;
  std::cout << "Number of nodes: " << numberOfNodes << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Mark nodes for apical growth " << std::endl;
  markApicalGrowth();

  std::cout << " " << std::endl;
  std::cout << "Process apical growth" << std::endl;
  processApicalGrowth();

  numberOfNodes = VGM.getNumberOfNodes();

  std::cout << " " << std::endl;
  std::cout << "Number of nodes after growing the network: " << numberOfNodes << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Mark edges for sprouting " << std::endl;
  markSproutingGrowth();

  std::cout << " " << std::endl;
  std::cout << "Process sprouting growth" << std::endl;
  processSproutingGrowth();

  numberOfNodes = VGM.getNumberOfNodes();

  std::cout << " " << std::endl;
  std::cout << "Number of nodes after growing the network: " << numberOfNodes << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Remove terminal vessels" << std::endl;
  removeRedundantTerminalVessels();

  std::cout << " " << std::endl;
  std::cout << "Link terminal vessels" << std::endl;
  linkTerminalVessels();

  std::cout << " " << std::endl;
  std::cout << "Adapt radius" << std::endl;
  adaptRadius();

  numberOfNodes = VGM.getNumberOfNodes();

  std::cout << " " << std::endl;
  std::cout << "Number of nodes after linking terminal vessels to the network: " << numberOfNodes << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Rescale the 1D matrices and vectors" << std::endl;

  Ac_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  b_c = std::vector<double>(numberOfNodes, 0.0);

  A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
  b_flow_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
  b_nut_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  phi_sigma = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);
  phi_sigma_old = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  P_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  for (int i = 0; i < N_tot_3D; i++) {

    phi_sigma[i] = phi_sigma_3D[i];
    phi_sigma_old[i] = phi_sigma_3D[i];
    P_3D1D[i] = P_3D[i];
  }

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    phi_sigma[N_tot_3D + indexOfNode] = pointer->c_v;
    phi_sigma_old[N_tot_3D + indexOfNode] = pointer->c_v;
    P_3D1D[N_tot_3D + indexOfNode] = pointer->p_v;

    pointer = pointer->global_successor;
  }
}

void netfc::Network::linkTerminalVessels() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    if (pointer->neighbors.size() == 1) {

      std::vector<double> coord = pointer->coord;

      std::shared_ptr<VGNode> pointer_1 = VGM.getHead();

      double radius = pointer->radii[0];

      int index = pointer->index;

      std::vector<double> dir_term_vessel = std::vector<double>(3, 0.0);

      for (int i = 0; i < 3; i++) {

        dir_term_vessel[i] = coord[i] - pointer->neighbors[0]->coord[i];
      }

      double length_dir = gmm::vect_norm2(dir_term_vessel);

      std::vector<double> normal_plane = std::vector<double>(3, 0.0);

      for (int i = 0; i < 3; i++) {

        normal_plane[i] = dir_term_vessel[i] / length_dir;
      }

      double p_v = pointer->p_v;

      int numberOfNeighbors = pointer->neighbors.size();

      while (pointer_1) {

        std::vector<double> coord_1 = pointer_1->coord;

        int index_1 = pointer_1->index;

        double dist = util::dist_between_points(coord, coord_1);

        std::vector<double> diff = std::vector<double>(3, 0.0);

        for (int i = 0; i < 3; i++) {

          diff[i] = coord_1[i] - coord[i];
        }

        double dist_plane = 0.0;

        for (int i = 0; i < 3; i++) {

          dist_plane += normal_plane[i] * diff[i];
        }

        double pv_1 = pointer_1->p_v;

        int numberOfNeighbors_1 = pointer_1->neighbors.size();

        if (dist_plane > 0.05 && index != index_1 && dist < h_3D) {

          if (numberOfNeighbors_1 == 1 && pv_1 < 0.95 * p_v) {

            std::cout << " " << std::endl;
            std::cout << "dist: " << dist << "\n";
            std::cout << "index: " << index << "\n";
            std::cout << "index_1: " << index_1 << "\n";
            std::cout << "Update pointer"
                      << "\n";

            pointer->p_boundary = 0.0;
            pointer->c_boundary = 0.0;
            pointer->typeOfVGNode = TypeOfNode::InnerNode;
            pointer->apicalGrowth = false;
            pointer->radii.push_back(radius);
            pointer->radii_initial.push_back(radius);
            pointer->L_p.push_back(input.d_tissue_flow_L_p);
            pointer->L_s.push_back(input.d_tissue_nut_L_s);
            pointer->edge_touched.push_back(true);
            pointer->sprouting_edge.push_back(false);
            pointer->neighbors.push_back(pointer_1);

            double length = util::dist_between_points(pointer->coord, pointer_1->coord);
            double p_node = pointer->p_v;
            double p_neighbor = pointer_1->p_v;
            double delta_p = p_neighbor - p_node;
            double tau_w_ini = (radius * std::abs(delta_p)) / (2.0 * length);

            pointer->tau_w_initial.push_back(tau_w_ini);

            std::cout << "Update pointer 1"
                      << "\n";

            pointer_1->p_boundary = 0.0;
            pointer_1->c_boundary = 0.0;
            pointer_1->typeOfVGNode = TypeOfNode::InnerNode;
            pointer_1->apicalGrowth = false;
            pointer_1->radii.push_back(radius);
            pointer_1->radii_initial.push_back(radius);
            pointer_1->L_p.push_back(input.d_tissue_flow_L_p);
            pointer_1->L_s.push_back(input.d_tissue_nut_L_s);
            pointer_1->edge_touched.push_back(true);
            pointer_1->sprouting_edge.push_back(false);
            pointer_1->neighbors.push_back(pointer);
            pointer_1->tau_w_initial.push_back(tau_w_ini);

          } else if (numberOfNeighbors_1 > 1) {

            std::cout << " " << std::endl;
            std::cout << "dist: " << dist << "\n";
            std::cout << "index: " << index << "\n";
            std::cout << "index_1: " << index_1 << "\n";
            std::cout << "Update pointer"
                      << "\n";

            pointer->p_boundary = 0.0;
            pointer->c_boundary = 0.0;
            pointer->typeOfVGNode = TypeOfNode::InnerNode;
            pointer->apicalGrowth = false;
            pointer->radii.push_back(radius);
            pointer->radii_initial.push_back(radius);
            pointer->L_p.push_back(input.d_tissue_flow_L_p);
            pointer->L_s.push_back(input.d_tissue_nut_L_s);
            pointer->edge_touched.push_back(true);
            pointer->sprouting_edge.push_back(false);
            pointer->neighbors.push_back(pointer_1);

            double length = util::dist_between_points(pointer->coord, pointer_1->coord);
            double p_node = pointer->p_v;
            double p_neighbor = pointer_1->p_v;
            double delta_p = p_neighbor - p_node;
            double tau_w_ini = (radius * std::abs(delta_p)) / (2.0 * length);

            pointer->tau_w_initial.push_back(tau_w_ini);

            std::cout << "Update pointer 1"
                      << "\n";

            pointer_1->p_boundary = 0.0;
            pointer_1->c_boundary = 0.0;
            pointer_1->typeOfVGNode = TypeOfNode::InnerNode;
            pointer_1->apicalGrowth = false;
            pointer_1->radii.push_back(radius);
            pointer_1->radii_initial.push_back(radius);
            pointer_1->L_p.push_back(input.d_tissue_flow_L_p);
            pointer_1->L_s.push_back(input.d_tissue_nut_L_s);
            pointer_1->edge_touched.push_back(true);
            pointer_1->sprouting_edge.push_back(false);
            pointer_1->neighbors.push_back(pointer);
            pointer_1->tau_w_initial.push_back(tau_w_ini);
          }

          break;
        }

        pointer_1 = pointer_1->global_successor;
      }
    }

    pointer = pointer->global_successor;
  }

  //Reset values
  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    pointer->apicalGrowth = false;

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->sprouting_edge[i] = false;
      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }
}

void netfc::Network::markApicalGrowth() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::cout << " " << std::endl;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1) {

      const auto &coord = pointer->coord;

      // if node is near the boundary, we do not process Omega = (0,L)^3
      if (0.0001 < coord[0] && coord[0] < L_x - 0.0001 && 0.0001 < coord[1] && coord[1] < L_x - 0.0001 && 0.0001 < coord[2] && coord[2] < L_x - 0.0001) {

        int index = getElementIndex(coord, h_3D, N_3D);

        double taf_node = phi_TAF[index];

        std::cout << "taf_node: " << taf_node << std::endl;

        if (taf_node > input.d_network_update_taf_threshold) {

          std::cout << "index: " << pointer->index << std::endl;

          pointer->apicalGrowth = true;
        }
      }
    }

    pointer = pointer->global_successor;
  }
}

void netfc::Network::processApicalGrowth() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  std::cout << "L_x: " << L_x << "\n";

  // Initialize random objects
  std::lognormal_distribution<> log_normal_distribution(input.d_log_normal_mean, input.d_log_normal_std_dev);
  std::random_device rd;
  std::mt19937 generator(rd());

  int numberOfNodes_old = VGM.getNumberOfNodes();

  // mark node for growth based on a certain criterion
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::cout << "Number of nodes before: " << VGM.getNumberOfNodes() << "\n";

  while (pointer) {

    if (pointer->apicalGrowth) {

      std::cout << " "
                << "\n";
      std::cout << "Processing node: " << pointer->index << "\n";

      std::vector<double> coord = pointer->coord;
      std::cout << "Compute direction based on TAF"
                << "\n";

      int element_index = getElementIndex(coord, h_3D, N_3D);
      std::cout << "element_index: " << element_index << "\n";

      std::vector<int> indicesNeighbors = getNeighboringElementIndices(element_index, N_3D, h_3D, L_x);
      std::vector<double> TAF_neighbors;

      double TAF_point = phi_TAF[element_index];

      double TAF_max = 0.0;

      double TAF_max_2 = 0.0;

      double TAF_min = 0.0;

      std::vector<double> new_point_1 = std::vector<double>(3, 0.0);
      std::vector<double> diff = std::vector<double>(3, 0.0);
      std::vector<double> dir_term_vessel = std::vector<double>(3, 0.0);
      std::vector<double> normal_plane = std::vector<double>(3, 0.0);

      for (int i = 0; i < 3; i++) {

        dir_term_vessel[i] = coord[i] - pointer->neighbors[0]->coord[i];
      }

      double length_dir = gmm::vect_norm2(dir_term_vessel);

      double prod_coord, dist_new_point = 0.0;

      for (int i = 0; i < 3; i++) {

        normal_plane[i] = dir_term_vessel[i] / length_dir;
      }

      // lognormal distribution
      double log_dist = log_normal_distribution(generator);
      double radius_p = pointer->radii[0];

      // get length
      double length = log_dist * radius_p;

      if (radius_p < 8.5e-3) {

        radius_p = 8.5e-3;
      }

      if (length > 3.0 * h_3D) {

        length = 3.0 * h_3D;
        std::cout << "Length is shortened !!!" << std::endl;
      }

      std::vector<double> rotator = determineRotator(normal_plane);

      double length_rotator = normVector(rotator);

      std::vector<double> midpoint(3, 0.0);
      std::vector<double> max_vec(3, 0.0);
      std::vector<double> max_vec_2(3, 0.0);
      std::vector<double> min_vec(3, 0.0);

      double theta = 0.0;

      for (int j = 0; j < 3; j++) {

        midpoint[j] = coord[j] + (length * normal_plane[j]);
        rotator[j] = rotator[j] / length_rotator;
      }

      int N_theta = 35;

      int N_r = 35;

      for (int i_theta = 0; i_theta < N_theta; i_theta++) {

        theta = ((double) i_theta) / ((double) N_theta) * 2.0 * M_PI;

        for (int j_r = 1; j_r < N_r; j_r++) {

          double r = ((double) j_r) / ((double) N_r) * length * std::tan(70.0 / 180.0 * M_PI);

          std::vector<double> cylinder_node = computeNodesOnCylinders(normal_plane, rotator, midpoint, r, theta);

          if (isCenterInDomain(cylinder_node, L_x) && length_rotator > 0.0) {

            int index_cone = getElementIndex(cylinder_node, h_3D, N_3D);

            double TAF = phi_TAF[index_cone];

            if (TAF_max < 1.0e-16) {

              TAF_max = TAF;
            }

            if (TAF > TAF_max - 1.0e-8) {

              TAF_max = TAF;

              max_vec = cylinder_node;

            } else if (TAF_max - 1.0e-8 > TAF && TAF > TAF_max_2 - 1.0e-8) {

              TAF_max_2 = TAF;

              max_vec_2 = cylinder_node;
            }

            if (TAF_min - 1.0e-8 > TAF) {

              TAF_min = TAF;

              min_vec = cylinder_node;
            }

          } else {

            //std::cout << "Node not considered!" << std::endl;
          }
        }
      }

      std::cout << " " << std::endl;
      std::vector<double> direction(3, 0.0);

      if (TAF_point < TAF_max) {

        for (int i = 0; i < 3; i++) {

          direction[i] = (max_vec[i] - coord[i]) + (0.5 * normal_plane[i]);
          //( 0.5*max_vec[ i ] + 0.5*max_vec_2[ i ] - coord[ i ] ) + ( 0.5 * normal_plane[ i ] );
        }

      } else {

        std::cout << "Find nearest node!" << std::endl;

        std::vector<double> near_node(3, 0.0);

        near_node = findNearNetworkNode(coord, normal_plane);

        std::cout << "near_node: " << near_node << std::endl;

        for (int i = 0; i < 3; i++) {

          direction[i] = (0.5 * max_vec[i] + 0.5 * max_vec_2[i] - coord[i]) + (1.0 * normal_plane[i]);
        }
      }

      double length_d = gmm::vect_norm2(direction);

      if (length_d > 0.0) {

        for (int i = 0; i < 3; i++) {

          direction[i] = direction[i] / length_d;
        }

        for (int i = 0; i < 3; i++) {

          new_point_1[i] = coord[i] + length * direction[i];
        }
      }

      std::vector<double> new_point_link = std::vector<double>(3, 0.0);

      bool isIntersecting = testIntersection(coord, new_point_1, radius_p, pointer);

      std::cout << "coord: " << coord << "\n";
      std::cout << "rotator: " << rotator << "\n";
      std::cout << "max_vec: " << max_vec << "\n";
      std::cout << "normal_plane: " << normal_plane << "\n";
      std::cout << "length_d: " << length_d << "\n";
      std::cout << "dir_term_vessel: " << dir_term_vessel << "\n";
      std::cout << "direction: " << direction << "\n";
      std::cout << "new_point: " << new_point_1 << "\n";
      std::cout << "TAF_point: " << TAF_point << "\n";

      double global_max_TAF = gmm::vect_norminf(phi_TAF);
      std::cout << "global_max_TAF: " << global_max_TAF << "\n";
      std::cout << "TAF_max: " << TAF_max << "\n";
      std::cout << "TAF_min: " << TAF_min << "\n";

      // check if we bifurcate at this node
      bool bifurcate = false;

      double prob = 0.5 + 0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) / std::sqrt(2.0 * input.d_log_normal_std_dev * input.d_log_normal_std_dev));
      std::cout << "prob: " << prob << "\n";

      if (prob > 0.935) {

        bifurcate = true;
      }

      if (!bifurcate && length_d > 0.0) {

        if (!isIntersecting) {

          createASingleNode(new_point_1, radius_p, pointer);
        }

      } else if (bifurcate) {

        std::cout << "Create bifuraction"
                  << "\n";

        double gamma = input.d_net_radius_exponent_gamma;
        double R_c = std::pow(2.0, -1.0 / gamma) * radius_p;

        if (R_c > radius_p) {

          R_c = radius_p;
        }

        // create normal distribution function
        std::normal_distribution<> normal_distribution(R_c, R_c / 35.0);

        double radius_b1 = normal_distribution(generator);
        double radius_b2 = normal_distribution(generator);

        if (radius_b1 < 8.5e-3) {

          radius_b1 = 8.5e-3;
        }

        if (radius_b2 < 8.5e-3) {

          radius_b2 = 8.5e-3;
        }

        std::vector<double> new_point_2 = std::vector<double>(3, 0.0);

        double branch_angle_1 = 0.0;
        double branch_angle_2 = 0.0;

        double angle_arg_1 = ((radius_p * radius_p * radius_p * radius_p) + (radius_b1 * radius_b1 * radius_b1 * radius_b1) - (radius_b2 * radius_b2 * radius_b2 * radius_b2)) /
                             (2.0 * radius_p * radius_p * radius_b1 * radius_b1);

        double angle_arg_2 = ((radius_p * radius_p * radius_p * radius_p) + (radius_b2 * radius_b2 * radius_b2 * radius_b2) - (radius_b1 * radius_b1 * radius_b1 * radius_b1)) /
                             (2.0 * radius_p * radius_p * radius_b2 * radius_b2);

        if (std::abs(angle_arg_1) < 1.0 && std::abs(angle_arg_2) < 1.0) {

          branch_angle_1 = std::acos(angle_arg_1);
          branch_angle_2 = std::acos(angle_arg_2);

          std::cout << "radius_p: " << radius_p << "\n";
          std::cout << "radius_b1: " << radius_b1 << "\n";
          std::cout << "radius_b2: " << radius_b2 << "\n";
          std::cout << "branch_angle_1: " << branch_angle_1 * 180.0 / M_PI << "\n";
          std::cout << "branch_angle_2: " << branch_angle_2 * 180.0 / M_PI << "\n";

          double branch_angle = branch_angle_1 + branch_angle_2;

          std::cout << "branch_angle: " << branch_angle * 180.0 / M_PI << "\n";

          if (branch_angle * 180.0 / M_PI < 160.0 && branch_angle * 180.0 / M_PI > 40.0) {

            std::vector<double> rotation_axis = util::cross_prod(direction, dir_term_vessel);
            std::vector<double> diff_2 = util::rotate(normal_plane, branch_angle_2, rotation_axis);
            std::vector<double> diff_1 = util::rotate(normal_plane, -branch_angle_1, rotation_axis);

            if (util::dist_between_points(diff_1, direction) < util::dist_between_points(diff_2, direction)) {

              for (int i = 0; i < 3; i++) {

                diff_1[i] = 0.5 * diff_1[i] + 0.5 * direction[i];
              }

            } else {

              for (int i = 0; i < 3; i++) {

                diff_2[i] = 0.5 * diff_2[i] + 0.5 * direction[i];
              }
            }

            double length_diff_1 = gmm::vect_norm2(diff_1);
            double length_diff_2 = gmm::vect_norm2(diff_2);

            // lognormal distribution
            double log_dist = log_normal_distribution(generator);

            // get length
            double length_1 = log_dist * radius_b1;
            double length_2 = log_dist * radius_b2;

            if (length_1 > 3.0 * h_3D) {

              length_1 = 3.0 * h_3D;
            }

            if (length_2 > 3.0 * h_3D) {

              length_2 = 3.0 * h_3D;
            }

            if (length_diff_2 > 0.0 && length_diff_1 > 0.0) {

              for (int i = 0; i < 3; i++) {

                new_point_1[i] = coord[i] + (length_1 * diff_1[i] / length_diff_1);
                new_point_2[i] = coord[i] + (length_2 * diff_2[i] / length_diff_2);
              }
            }

            std::cout << "rotation_axis: " << rotation_axis << "\n";
            std::cout << "branch_angle: " << branch_angle << "\n";
            std::cout << "length_1: " << length_1 << "\n";
            std::cout << "length_2: " << length_2 << "\n";
            std::cout << "new_point_1: " << new_point_1 << "\n";
            std::cout << "new_point_2: " << new_point_2 << "\n";

            if (gmm::vect_norm2(direction) > 0.0 && length_diff_2 > 0.0 && length_diff_1 > 0.0) {

              bool isIntersecting_1 = testIntersection(coord, new_point_1, radius_p, pointer);
              bool isIntersecting_2 = testIntersection(coord, new_point_2, radius_p, pointer);

              if (!isIntersecting_1) {

                createASingleNode(new_point_1, radius_b1, pointer);
              }

              if (!isIntersecting_2) {

                createASingleNode(new_point_2, radius_b2, pointer);
              }
            }
          }
        }
      }
    }

    pointer = pointer->global_successor;
  }

  // reset the boolean values
  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    pointer->apicalGrowth = false;

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->sprouting_edge[i] = false;
      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }
}


void netfc::Network::createASingleNode(std::vector<double> new_point, double radius, std::shared_ptr<VGNode> &pointer) {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  bool isInside = isCenterInDomain(new_point, L_x);

  std::cout << "Create new node"
            << "\n";

  if (isInside) {

    VGNode new_node;

    new_node.index = VGM.getNumberOfNodes();
    new_node.coord = new_point;
    new_node.p_boundary = 0.95 * pointer->p_v;
    new_node.p_v = 0.95 * pointer->p_v;
    new_node.c_boundary = input.d_in_nutrient;
    new_node.c_v = 0.0;
    new_node.typeOfVGNode = TypeOfNode::DirichletNode;
    new_node.apicalGrowth = false;
    new_node.radii.push_back(radius);
    new_node.radii_initial.push_back(radius);
    new_node.L_p.push_back(input.d_tissue_flow_L_p);
    new_node.L_s.push_back(input.d_tissue_nut_L_s);
    new_node.edge_touched.push_back(true);
    new_node.sprouting_edge.push_back(false);
    new_node.neighbors.push_back(pointer);
    new_node.notUpdated = 0;
    new_node.typeOfVGNode = TypeOfNode::NeumannNode;

    double length = util::dist_between_points(pointer->coord, new_point);
    double p_node = pointer->p_v;
    double p_neighbor = 0.95 * pointer->p_v;
    double delta_p = p_neighbor - p_node;
    double tau_w_ini = (radius * std::abs(delta_p)) / (2.0 * length);

    new_node.tau_w_initial.push_back(tau_w_ini);

    auto sp_newNode = std::make_shared<VGNode>(new_node);
    std::cout << "New index: " << new_node.index << "\n";
    std::cout << "Neighbor index: " << pointer->index << "\n";

    std::cout << "Update old node"
              << "\n";

    pointer->p_boundary = 0.0;
    pointer->c_boundary = 0.0;
    pointer->typeOfVGNode = TypeOfNode::InnerNode;
    pointer->apicalGrowth = false;
    pointer->radii.push_back(radius);
    pointer->radii_initial.push_back(radius);
    pointer->L_p.push_back(input.d_tissue_flow_L_p);
    pointer->L_s.push_back(input.d_tissue_nut_L_s);
    pointer->edge_touched.push_back(true);
    pointer->sprouting_edge.push_back(false);
    pointer->neighbors.push_back(sp_newNode);
    pointer->notUpdated = 0;
    pointer->tau_w_initial.push_back(tau_w_ini);
    std::cout << "Attach new node as pointer"
              << "\n";

    VGM.attachPointerToNode(sp_newNode);
  }
}


bool netfc::Network::testCollision(std::vector<double> point) {

  bool isColliding = false;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::cout << "Test Collision"
            << "\n";

  while (pointer) {

    std::vector<double> coord = pointer->coord;

    std::vector<double> diff = std::vector<double>(3, 0.0);

    for (int i = 0; i < 3; i++) {

      diff[i] = coord[i] - point[i];
    }

    double dist = gmm::vect_norm2(diff);

    if (dist < h_3D) {

      //std::cout << "Node not inserted, dist: " << dist << "\n";

      isColliding = true;

      break;
    }

    pointer = pointer->global_successor;
  }

  return isColliding;
}

void netfc::Network::removeRedundantTerminalVessels() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1) {

      const auto &coord = pointer->coord;

      // if node is near the boundary, we do not process Omega = (0,L)^3
      if (0.0001 < coord[0] && coord[0] < L_x - 0.0001 && 0.0001 < coord[1] && coord[1] < L_x - 0.0001 && 0.0001 < coord[2] && coord[2] < L_x - 0.0001) {

        int updateNumber = pointer->notUpdated;

        pointer->notUpdated = updateNumber + 1;

        std::cout << "pointer->notUpdated: " << pointer->notUpdated << std::endl;
      }
    }

    pointer = pointer->global_successor;
  }

  std::cout << " " << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1 && pointer->notUpdated > 2) {

      int index = pointer->index;

      std::cout << "Remove node with index: " << index << std::endl;
      std::cout << "coord: " << pointer->coord << std::endl;

      int numberOfNeighbors_neighbor = pointer->neighbors[0]->neighbors.size();

      int index_neighbor = pointer->neighbors[0]->index;

      std::cout << "numberOfNeighbors: " << pointer->neighbors.size() << std::endl;
      std::cout << "numberOfNeighbors_neighbor: " << numberOfNeighbors_neighbor << std::endl;
      std::cout << "index_neighbor: " << index_neighbor << std::endl;

      std::vector<std::shared_ptr<VGNode>> new_neighbors;

      std::vector<bool> new_edge_touched;
      std::vector<bool> new_sprouting_edge;
      std::vector<double> new_radii;
      std::vector<double> new_L_p;
      std::vector<double> new_L_s;
      std::vector<double> new_tau_w_initial;

      for (int i = 0; i < numberOfNeighbors_neighbor; i++) {

        std::cout << "pointer->neighbors[ 0 ]->neighbors[ i ]->index: " << pointer->neighbors[0]->neighbors[i]->index << std::endl;

        if (index != pointer->neighbors[0]->neighbors[i]->index) {

          new_neighbors.push_back(pointer->neighbors[0]->neighbors[i]);
          new_edge_touched.push_back(pointer->neighbors[0]->edge_touched[i]);
          new_sprouting_edge.push_back(pointer->neighbors[0]->sprouting_edge[i]);
          new_radii.push_back(pointer->neighbors[0]->radii[i]);
          new_L_p.push_back(pointer->neighbors[0]->L_p[i]);
          new_L_s.push_back(pointer->neighbors[0]->L_s[i]);
          new_tau_w_initial.push_back(pointer->neighbors[0]->tau_w_initial[i]);
        }
      }

      pointer->neighbors[0]->neighbors = new_neighbors;
      pointer->neighbors[0]->edge_touched = new_edge_touched;
      pointer->neighbors[0]->sprouting_edge = new_sprouting_edge;
      pointer->neighbors[0]->radii = new_radii;
      pointer->neighbors[0]->radii_initial = new_radii;
      pointer->neighbors[0]->tau_w_initial = new_tau_w_initial;
      pointer->neighbors[0]->L_p = new_L_p;
      pointer->neighbors[0]->L_s = new_L_s;
      pointer->neighbors[0]->notUpdated = 0;

      if (numberOfNeighbors_neighbor == 2) {

        pointer->neighbors[0]->typeOfVGNode = TypeOfNode::NeumannNode;
        pointer->neighbors[0]->p_boundary = 0.0;
        pointer->neighbors[0]->c_boundary = 1.0;

      } else {

        pointer->neighbors[0]->typeOfVGNode = TypeOfNode::InnerNode;
        pointer->neighbors[0]->p_boundary = 0.0;
        pointer->neighbors[0]->c_boundary = 0.0;
      }

      std::shared_ptr<VGNode> old_pointer = pointer;

      if (pointer->global_predecessor) {

        pointer->global_predecessor->global_successor = pointer->global_successor;

      } else {

        pointer->global_successor->global_predecessor = NULL;
      }

      if (pointer->global_successor) {

        pointer->global_successor->global_predecessor = pointer->global_predecessor;

      } else {

        pointer->global_predecessor->global_successor = NULL;
      }

      pointer = pointer->global_successor;

      old_pointer.reset();

    } else {

      pointer = pointer->global_successor;
    }
  }

  VGM.determineNumberOfNodes();
  int numberOfNodes = VGM.getNumberOfNodes();

  std::cout << " " << std::endl;
  std::cout << "Number of nodes after removing redundant nodes: " << numberOfNodes << std::endl;

  //renumber nodes
  std::cout << "Renumber nodes" << std::endl;

  pointer = VGM.getHead();

  int counter = 0;

  while (pointer) {

    pointer->index = counter;

    counter = counter + 1;

    const auto &coord = pointer->coord;

    if (0.0001 < coord[0] && coord[0] < L_x - 0.0001 && 0.0001 < coord[1] && coord[1] < L_x - 0.0001 && 0.0001 < coord[2] && coord[2] < L_x - 0.0001) {

      if (pointer->typeOfVGNode == TypeOfNode::InnerNode && pointer->neighbors.size() == 1) {

        pointer->typeOfVGNode = TypeOfNode::NeumannNode;
      }
    }

    pointer = pointer->global_successor;
  }
}


void netfc::Network::markSproutingGrowth() {

  const auto &input = d_model_p->get_input_deck();

  std::lognormal_distribution<> log_normal_distribution(input.d_log_normal_mean, input.d_log_normal_std_dev);
  std::random_device rd;
  std::mt19937 generator(rd());

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    for (int i = 0; i < numberOfEdges; i++) {

      if (pointer->edge_touched[i] == false) {

        double radius = pointer->radii[i];

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        int local_index = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        double sproutingProbability = 0.0;

        std::vector<double> diff = std::vector<double>(3, 0.0);
        std::vector<double> mid_point = std::vector<double>(3, 0.0);

        for (int j = 0; j < 3; j++) {

          diff[j] = coord_neighbor[j] - coord[j];
          mid_point[j] = 0.5 * (coord_neighbor[j] + coord[j]);
        }

        int element_index = getElementIndex(mid_point, h_3D, N_3D);

        double TAF = phi_TAF[element_index];

        double TAF_th = input.d_network_update_taf_threshold;

        double length = gmm::vect_norm2(diff);

        double log_dist = log_normal_distribution(generator);

        sproutingProbability = 0.5 + 0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) / std::sqrt(2.0 * input.d_log_normal_std_dev * input.d_log_normal_std_dev));

        double global_max_TAF = gmm::vect_norminf(phi_TAF);

        if (sproutingProbability > 0.9 && TAF > TAF_th) {

          pointer->neighbors[i]->sprouting_edge[local_index] = true;
          pointer->sprouting_edge[i] = true;
        }

        pointer->neighbors[i]->edge_touched[local_index] = true;
        pointer->edge_touched[i] = true;
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {

      std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }
}


void netfc::Network::processSproutingGrowth() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double gamma = input.d_net_radius_exponent_gamma;

  std::lognormal_distribution<> log_normal_distribution(input.d_log_normal_mean, input.d_log_normal_std_dev);
  std::random_device rd;
  std::mt19937 generator_log(rd());

  std::default_random_engine generator;

  double L_x = input.d_domain_params[1];

  std::cout << " " << std::endl;

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    for (int i = 0; i < numberOfEdges; i++) {

      if (pointer->sprouting_edge[i] == true && pointer->edge_touched[i] == false) {

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        std::vector<double> mid_point = std::vector<double>(3, 0.0);

        double radius = pointer->radii[i];

        double PSI = 1.05;

        double radius_prime = std::pow((std::pow(PSI, gamma) - 1.0), 1.0 / gamma) * radius;

        double radius_min = 8.5e-3;

        double p_v_neighbor = pointer->neighbors[i]->p_v;

        std::uniform_real_distribution<double> distribution_uniform(radius_min, radius_prime);

        double radius_new = radius_min;

        if (radius_prime > radius_min) {

          radius_new = distribution_uniform(generator);
        }

        for (int j = 0; j < 3; j++) {

          mid_point[j] = 0.5 * (coord_neighbor[j] + coord[j]);
        }

        std::cout << "Compute direction of growth" << std::endl;

        int element_index = getElementIndex(mid_point, h_3D, N_3D);

        std::vector<int> indicesNeighbors = getNeighboringElementIndices(element_index, N_3D, h_3D, L_x);
        std::vector<double> TAF_neighbors;

        for (int j = 0; j < numberOfEdges; j++) {

          TAF_neighbors.push_back(phi_TAF[indicesNeighbors[j]]);
        }

        int max_index = 0;
        int max_index_2 = 0;

        double TAF_max = 0.0;
        double TAF_max_2 = 0.0;

        for (int j = 0; j < indicesNeighbors.size(); j++) {

          double TAF = TAF_neighbors[j];

          if (TAF > TAF_max - 1.0e-8) {

            max_index = indicesNeighbors[j];

            TAF_max = TAF;

          } else if (TAF_max - 1.0e-8 > TAF && TAF > TAF_max_2 - 1.0e-8) {

            max_index_2 = indicesNeighbors[j];

            TAF_max_2 = TAF;
          }
        }

        std::vector<double> dir_new_vessel = std::vector<double>(3, 0.0);
        std::vector<double> dir_vessel = std::vector<double>(3, 0.0);
        std::vector<double> new_point = std::vector<double>(3, 0.0);

        std::vector<double> max_center = getCenterFromIndex(max_index, N_3D, h_3D);
        std::vector<double> max_center_2 = getCenterFromIndex(max_index_2, N_3D, h_3D);

        for (int j = 0; j < 3; j++) {

          dir_new_vessel[j] = max_center[j] - mid_point[j];
          dir_vessel[j] = coord_neighbor[j] - coord[j];
        }

        double norm_dir_new_vessel = gmm::vect_norm2(dir_new_vessel);
        double norm_dir_vessel = gmm::vect_norm2(dir_vessel);

        for (int j = 0; j < 3; j++) {

          dir_new_vessel[j] = dir_new_vessel[j] / norm_dir_new_vessel;
          dir_vessel[j] = dir_vessel[j] / norm_dir_vessel;
        }

        double prod_angle = 0.0;

        for (int j = 0; j < 3; j++) {

          prod_angle += dir_new_vessel[j] * dir_vessel[j];
        }

        double angle = std::acos(prod_angle);

        // lognormal distribution
        double log_dist = log_normal_distribution(generator_log);

        // get length
        double length_new = log_dist * radius_new;

        if (length_new > 3.0 * h_3D) {

          length_new = 3.0 * h_3D;
        }

        if (length_new < 2.0 * radius) {

          length_new = 2.0 * radius;
        }

        for (int j = 0; j < 3; j++) {

          new_point[j] = mid_point[j] + length_new * dir_new_vessel[j];
        }


        bool isColliding = testCollision(new_point);

        if (angle * 180.0 / M_PI > 20.0 && angle * 180.0 / M_PI < 90.0 && norm_dir_new_vessel > 0.0 && !isColliding) {

          std::cout << "Create new node_1" << std::endl;
          VGNode new_node_1;

          new_node_1.index = VGM.getNumberOfNodes();
          new_node_1.coord = mid_point;
          new_node_1.p_boundary = 0.0;
          new_node_1.p_v = pointer->p_v;
          new_node_1.c_boundary = input.d_in_nutrient;
          new_node_1.c_v = pointer->c_v;
          new_node_1.typeOfVGNode = TypeOfNode::InnerNode;
          new_node_1.apicalGrowth = false;

          new_node_1.radii.push_back(radius);
          new_node_1.radii.push_back(radius);
          new_node_1.radii.push_back(radius_new);

          new_node_1.radii_initial.push_back(radius);
          new_node_1.radii_initial.push_back(radius);
          new_node_1.radii_initial.push_back(radius_new);

          new_node_1.L_p.push_back(input.d_tissue_flow_L_p);
          new_node_1.L_p.push_back(input.d_tissue_flow_L_p);
          new_node_1.L_p.push_back(input.d_tissue_flow_L_p);

          new_node_1.L_s.push_back(input.d_tissue_nut_L_s);
          new_node_1.L_s.push_back(input.d_tissue_nut_L_s);
          new_node_1.L_s.push_back(input.d_tissue_nut_L_s);

          new_node_1.edge_touched.push_back(true);
          new_node_1.edge_touched.push_back(true);
          new_node_1.edge_touched.push_back(true);

          new_node_1.sprouting_edge.push_back(false);
          new_node_1.sprouting_edge.push_back(false);
          new_node_1.sprouting_edge.push_back(false);

          new_node_1.neighbors.push_back(pointer);
          new_node_1.neighbors.push_back(pointer->neighbors[i]);
          new_node_1.notUpdated = 0;

          double length_old = util::dist_between_points(mid_point, coord_neighbor);
          double length_new = util::dist_between_points(mid_point, coord_neighbor);

          double delta_p_old = p_v_neighbor - pointer->p_v;
          double delta_p_new = 0.95 * pointer->p_v - pointer->p_v;

          double tau_w_ini_old = (radius * std::abs(delta_p_old)) / (2.0 * length_old);
          double tau_w_ini_new = (radius * std::abs(delta_p_new)) / (2.0 * length_new);

          new_node_1.tau_w_initial.push_back(tau_w_ini_old);
          new_node_1.tau_w_initial.push_back(tau_w_ini_old);
          new_node_1.tau_w_initial.push_back(tau_w_ini_new);

          auto sp_newNode_1 = std::make_shared<VGNode>(new_node_1);

          std::cout << "Create new node_2" << std::endl;
          VGNode new_node_2;

          new_node_2.index = VGM.getNumberOfNodes() + 1;
          new_node_2.coord = new_point;
          new_node_2.p_boundary = 0.95 * pointer->p_v;
          new_node_2.p_v = 0.95 * pointer->p_v;
          new_node_2.c_boundary = input.d_in_nutrient;
          new_node_2.c_v = pointer->c_v;
          new_node_2.typeOfVGNode = TypeOfNode::DirichletNode;
          new_node_2.apicalGrowth = false;
          new_node_2.radii.push_back(radius_new);
          new_node_2.radii_initial.push_back(radius_new);
          new_node_2.L_p.push_back(input.d_tissue_flow_L_p);
          new_node_2.L_s.push_back(input.d_tissue_nut_L_s);
          new_node_2.edge_touched.push_back(true);
          new_node_2.sprouting_edge.push_back(false);
          new_node_2.typeOfVGNode = TypeOfNode::NeumannNode;
          new_node_2.notUpdated = 0;
          new_node_2.tau_w_initial.push_back(tau_w_ini_new);
          new_node_2.neighbors.push_back(sp_newNode_1);

          auto sp_newNode_2 = std::make_shared<VGNode>(new_node_2);

          new_node_1.neighbors.push_back(sp_newNode_2);

          std::cout << "Update connectivity"
                    << "\n";

          pointer->replacePointerWithIndex(i, sp_newNode_1);

          int index_to_be_replaced = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

          pointer->neighbors[i]->replacePointerWithIndex(index_to_be_replaced, sp_newNode_1);
          pointer->markEdgeLocalIndex(i);
          pointer->sprouting_edge[i] = false;
          pointer->neighbors[i]->markEdgeLocalIndex(index_to_be_replaced);
          pointer->neighbors[i]->sprouting_edge[index_to_be_replaced] = false;

          std::cout << "Attach new node_1 as pointer"
                    << "\n";
          VGM.attachPointerToNode(sp_newNode_1);

          std::cout << "Attach new node_2 as pointer"
                    << "\n";
          VGM.attachPointerToNode(sp_newNode_2);

        } else {

          pointer->edge_touched[i] = true;
          pointer->sprouting_edge[i] = false;

          int index_neighbor = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);
          pointer->neighbors[i]->markEdgeLocalIndex(index_neighbor);
        }

      } else {

        pointer->edge_touched[i] = true;
        pointer->sprouting_edge[i] = false;

        int index_neighbor = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);
        pointer->neighbors[i]->markEdgeLocalIndex(index_neighbor);
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {

      pointer->edge_touched[i] = false;
      pointer->sprouting_edge[i] = false;
    }

    pointer = pointer->global_successor;
  }
}

std::vector<double> netfc::Network::findNearNetworkNode(std::vector<double> coord, std::vector<double> normal_plane) {

  std::vector<double> coord_near_node(3, 0.0);

  std::vector<double> coord_min_node(3, 0.0);

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double dist_min = 0.0;

  int numberOfNodes = 0;

  while (pointer) {

    coord_near_node = pointer->coord;

    double dist = util::dist_between_points(coord, coord_near_node);

    std::vector<double> diff = std::vector<double>(3, 0.0);

    for (int i = 0; i < 3; i++) {

      diff[i] = coord_near_node[i] - coord[i];
    }

    double dist_plane = 0.0;

    for (int i = 0; i < 3; i++) {

      dist_plane += normal_plane[i] * diff[i];
    }

    if (dist_plane > 0.0 && dist_plane < 0.4) {

      dist_min = dist_plane;

      numberOfNodes = numberOfNodes + 1;

      for (int i = 0; i < 3; i++) {

        coord_min_node[i] = coord_min_node[i] + coord_near_node[i];
      }
    }

    pointer = pointer->global_successor;
  }

  std::cout << "numberOfNodes: " << numberOfNodes << std::endl;

  if (numberOfNodes > 0) {

    for (int i = 0; i < 3; i++) {

      coord_min_node[i] = coord_min_node[i] / ((double) numberOfNodes);
    }

  } else {

    coord_min_node = normal_plane;
  }

  return coord_min_node;
}

bool netfc::Network::testIntersection(std::vector<double> point_1, std::vector<double> point_2, double radius, std::shared_ptr<VGNode> &pointer_test) {

  const auto &input = d_model_p->get_input_deck();

  bool isIntersecting = false;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::cout << "Test intersection " << std::endl;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    int index = pointer->index;

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (pointer->edge_touched[i] == false && !isIntersecting) {

        std::vector<double> point_p1 = pointer->coord;

        std::vector<double> point_p2 = pointer->neighbors[i]->coord;

        double radius_p = pointer->radii[i];

        if (util::dist_between_points(point_1, point_p1) < 1.0e-11 || util::dist_between_points(point_2, point_p1) < 1.0e-11 ||
            util::dist_between_points(point_1, point_p2) < 1.0e-11 || util::dist_between_points(point_2, point_p2) < 1.0e-11) {

          isIntersecting = false;

        } else {

          double length_p = util::dist_between_points(point_p1, point_p2);
          double length_test = util::dist_between_points(point_1, point_2);

          std::vector<double> dir_p(3, 0.0);
          std::vector<double> dir_test(3, 0.0);

          for (int j = 0; j < 3; j++) {

            dir_p[j] = point_p2[j] - point_p1[j];
            dir_test[j] = point_2[j] - point_1[j];
          }

          for (int j = 1; j < 10; j++) {

            std::vector<double> new_point_p(3, 0.0);
            std::vector<double> new_point_test(3, 0.0);

            for (int k = 0; k < 3; k++) {

              new_point_p[k] = point_p1[k] + (length_p * (double) j / 10.0 * dir_p[k]);
              new_point_test[k] = point_1[k] + (length_test * (double) j / 10.0 * dir_test[k]);
            }

            if (util::dist_between_points(new_point_p, new_point_test) < radius + radius_p) {

              isIntersecting = true;

              std::cout << "Vessel is intersecting" << std::endl;

              break;
            }
          }
        }
      }

      pointer->edge_touched[i] = true;
      pointer->neighbors[i]->markEdge(index);
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {

      pointer->edge_touched[i] = false;
      pointer->sprouting_edge[i] = false;
    }

    pointer = pointer->global_successor;
  }

  return isIntersecting;
}


void netfc::Network::adaptRadius() {

  const auto &input = d_model_p->get_input_deck();

  double dt = d_model_p->d_dt;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double tau_w_avg = 0.0;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (pointer->edge_touched[i] == false) {

        int local_index = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        double radius = pointer->radii[i];

        std::vector<double> coord_node = pointer->coord;
        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord_node, coord_neighbor);

        double p_node = pointer->p_v;
        double p_neighbor = pointer->neighbors[i]->p_v;
        double delta_p = p_neighbor - p_node;

        double tau_w = (radius * std::abs(delta_p)) / (2.0 * length);

        if (tau_w > tau_w_avg) {

          tau_w_avg = tau_w;
        }

        pointer->neighbors[i]->edge_touched[local_index] = true;
        pointer->edge_touched[i] = true;
      }
    }

    pointer = pointer->global_successor;
  }

  int numberOfNodes = VGM.getNumberOfNodes();

  tau_w_avg = tau_w_avg / (double) numberOfNodes;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {

      pointer->edge_touched[i] = false;
      pointer->sprouting_edge[i] = false;
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<bool> remove_edge(false, numberOfNeighbors);
    std::vector<int> edgesToBeRemoved;

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (pointer->edge_touched[i] == false) {

        int local_index = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        double radius = pointer->radii[i];

        double radius_initial = pointer->radii_initial[i];

        std::vector<double> coord_node = pointer->coord;
        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord_node, coord_neighbor);

        double p_node = pointer->p_v;
        double p_neighbor = pointer->neighbors[i]->p_v;
        double delta_p = p_neighbor - p_node;

        double c_node = pointer->c_v;
        double c_neighbor = pointer->neighbors[i]->c_v;

        double tau_w_ini = pointer->tau_w_initial[i];

        double tau_w = (radius * std::abs(delta_p)) / (2.0 * length);

        double A_r = 2.0 * (tau_w_ini - tau_w) - 20.0 * (0.5 * (c_node + c_neighbor));

        double R_initial = (8.0 * length) / (M_PI * radius_initial * radius_initial * radius_initial * radius_initial);

        double R_U = 1.3 * R_initial;
        double R_L = 0.7 * R_initial;

        double C_initial = -std::log((R_initial - R_L) / (R_U - R_initial));

        double R_new = (R_L + R_U * std::exp(A_r - C_initial)) / (1.0 + std::exp(A_r - C_initial));

        double radius_new = std::pow((8.0 * length) / (R_new * M_PI), 0.25);

        if (radius_new < 8.5e-3) {

          radius_new = 8.5e-3;

          std::cout << "Remove vessel with index: " << i << std::endl;

          edgesToBeRemoved.push_back(i);

          remove_edge[i] = true;

          int numberOfNeighbors_Neighbor = pointer->neighbors.size();
        }

        if (radius < 0.02) {

          std::cout << " " << std::endl;
          std::cout << "tau_avg: " << tau_w_avg << std::endl;
          std::cout << "tau_w: " << tau_w << std::endl;
          std::cout << "2.0*( tau_w_ini - tau_w ): " << 2.0 * (tau_w_ini - tau_w) << std::endl;
          std::cout << "-20.0*( 0.5*(c_node+c_neighbor) ): " << -20.0 * (0.5 * (c_node + c_neighbor)) << std::endl;
          std::cout << "A_r: " << A_r << std::endl;
          std::cout << "R_initial: " << R_initial << std::endl;
          std::cout << "R_U: " << R_U << std::endl;
          std::cout << "R_L: " << R_L << std::endl;
          std::cout << "R_new: " << R_new << std::endl;
          std::cout << "radius_initial: " << radius_initial << std::endl;
          std::cout << "radius: " << radius << std::endl;
          std::cout << "radius_new: " << radius_new << std::endl;
          std::cout << "coord_node: " << coord_node << std::endl;

          pointer->neighbors[i]->radii[local_index] = radius_new;
          pointer->radii[i] = radius_new;
        }

        pointer->neighbors[i]->edge_touched[local_index] = true;
        pointer->edge_touched[i] = true;
      }
    }

    // Remove redundant vessels
    for (int i = 0; i < numberOfNeighbors; i++) {

      if (remove_edge[i] == true) {

        int index_neighbor = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        pointer->neighbors[i]->removeComponent(index_neighbor);
      }
    }

    pointer->removeComponents(edgesToBeRemoved);

    pointer = pointer->global_successor;
  }

  // Remove nodes without neighbor
  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    if (numberOfEdges == 0) {
    }

    pointer = pointer->global_successor;
  }

  // Renumber nodes
  pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {

      pointer->edge_touched[i] = false;
      pointer->sprouting_edge[i] = false;
    }

    pointer = pointer->global_successor;
  }
}
