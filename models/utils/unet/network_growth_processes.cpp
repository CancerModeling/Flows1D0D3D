////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "modelUtil.hpp"
#include "netUtil.hpp"
#include "network.hpp"
#include "random_dist.hpp"

void util::unet::Network::updateNetwork(BaseAssembly &taf_sys,
                                        BaseAssembly &grad_taf_sys) {

  d_model_p->d_log("Update the network \n", "net update");

  int numberOfNodesOld = VGM.getNumberOfNodes();

  if (d_update_number % d_update_interval == 0) {

    // update TAF vector from libmesh taf system
    if (taf_sys.d_sys_name != "TAF")
      libmesh_error_msg("Must pass TAF system to update network.");

    // get TAF at element centroid
    taf_sys.localize_solution_with_elem_id_numbering_non_const_elem(phi_TAF_3D,
                                                                    {0}, false);

    if (d_procRank == 0) {

      d_model_p->d_log("Number of nodes: " + std::to_string(numberOfNodesOld) +
                           " \n",
                       "net update");

      // util::get_elem_sol(taf_sys, phi_TAF);

      d_model_p->d_log("Mark nodes for apical growth \n", "net update");
      markApicalGrowth();

      d_model_p->d_log("Process apical growth \n", "net update");
      processApicalGrowth();
      check_vessel_length();

      auto numberOfNodes = VGM.getNumberOfNodes();

      d_model_p->d_log("Number of nodes after growing the network: " +
                           std::to_string(numberOfNodes) + " \n",
                       "net update");

      d_model_p->d_log("Mark edges for sprouting \n", "net update");
      markSproutingGrowth();

      d_model_p->d_log("Process sprouting growth \n", "net update");
      processSproutingGrowth();
      // check_vessel_length();

      numberOfNodes = VGM.getNumberOfNodes();

      d_model_p->d_log("Number of nodes after growing the network: " +
                           std::to_string(numberOfNodes) + " \n",
                       "net update");
      d_model_p->d_log("Remove redundant vessels \n", "net update");
      removeRedundantTerminalVessels();
      // check_vessel_length();

      numberOfNodes = VGM.getNumberOfNodes();

      d_model_p->d_log("Link terminal vessels \n", "net update");
      linkTerminalVessels();
      // check_vessel_length();

      numberOfNodes = VGM.getNumberOfNodes();

      d_model_p->d_log(
          "Number of nodes after linking terminal vessels to the network: " +
              std::to_string(numberOfNodes) + " \n",
          "net update");

      d_model_p->d_log("Remove short vessels \n", "net update");
      auto pointer = VGM.getHead();

      while (pointer) {

        int numberOfNeighbors = pointer->neighbors.size();

        for (int i = 0; i < numberOfNeighbors; i++) {

          double length = util::dist_between_points(
              pointer->coord, pointer->neighbors[i]->coord);

          if (util::definitelyLessThan(length, 1.e-8)) {

            pointer->removeComponent(i);
            d_model_p->d_log("neighbor removed length=0!!! \n", "net update");
          }
        }

        pointer = pointer->global_successor;
      } // loop for zero length

      check_vessel_length();

      d_model_p->d_log("Remove isolated nodes \n", "net update");

      pointer = VGM.getHead();

      while (pointer) {

        int numberOfEdges = pointer->neighbors.size();

        if (numberOfEdges == 0) {

          d_model_p->d_log("Remove node \n", "net update");

          std::shared_ptr<VGNode> old_pointer = pointer;

          if (pointer->global_predecessor) {

            pointer->global_predecessor->global_successor =
                pointer->global_successor;

          } else {

            pointer->global_successor->global_predecessor = NULL;
          }

          if (pointer->global_successor) {

            pointer->global_successor->global_predecessor =
                pointer->global_predecessor;

          } else {

            pointer->global_predecessor->global_successor = NULL;
          }

          pointer = pointer->global_successor;

          old_pointer.reset();
        }

        pointer = pointer->global_successor;
      } // loop for remove node

      check_vessel_length();

      d_model_p->d_log("Reset nodes \n", "net update");

      pointer = VGM.getHead();

      while (pointer) {

        int numberOfEdges = pointer->neighbors.size();

        for (int i = 0; i < numberOfEdges; i++) {

          pointer->edge_touched[i] = false;
          pointer->sprouting_edge[i] = false;
        }

        pointer = pointer->global_successor;
      } // loop for reset node

      VGM.determineNumberOfNodes();

      numberOfNodes = VGM.getNumberOfNodes();

      d_model_p->d_log("Rescale the 1D matrices and vectors \n", "net update");
      if (!d_coupled_solver and numberOfNodesOld != numberOfNodes) {
        A_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes,
                                                       numberOfNodes);
        b.resize(numberOfNodes);
        P_v.resize(numberOfNodes);

        C_v.resize(numberOfNodes);
        C_v_old.resize(numberOfNodes);

        // Ac_VGM =
        //    gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes,
        //    numberOfNodes);
        b_c.resize(numberOfNodes);
      } // update matrix and vector

      if (d_coupled_solver and numberOfNodesOld != numberOfNodes) {

        A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
            N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
        b_flow_3D1D.resize(N_tot_3D + numberOfNodes);

        // A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
        //    N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
        b_nut_3D1D.resize(N_tot_3D + numberOfNodes);

        // resize function does not change the value of existing elements
        phi_sigma.resize(N_tot_3D + numberOfNodes);
        phi_sigma_old.resize(N_tot_3D + numberOfNodes);
        P_3D1D.resize(N_tot_3D + numberOfNodes);

        for (int i = 0; i < N_tot_3D; i++) {

          phi_sigma[i] = phi_sigma_3D[i];
          phi_sigma_old[i] = phi_sigma_3D[i];
          P_3D1D[i] = P_3D[i];
        }
      } // update matrix and vector

      pointer = VGM.getHead();
      while (pointer) {

        int indexOfNode = pointer->index;

        if (d_coupled_solver) {
          phi_sigma[N_tot_3D + indexOfNode] = pointer->c_v;
          phi_sigma_old[N_tot_3D + indexOfNode] = pointer->c_v;
          P_3D1D[N_tot_3D + indexOfNode] = pointer->p_v;
        } else {
          C_v[indexOfNode] = pointer->c_v;
          C_v_old[indexOfNode] = pointer->c_v;
          P_v[indexOfNode] = pointer->p_v;
        }

        pointer = pointer->global_successor;
      } // loop for update solution in network

      // compute element and weights (iterative method may require this data)
      const auto &input = d_model_p->get_input_deck();
      if (input.d_compute_elem_weights and input.d_model_name != "NetFCFVFE")
        compute_elem_weights();
      if (numberOfNodes != numberOfNodesOld) {
        d_has_network_changed = true;
        d_model_p->d_log("Added " +
                             std::to_string(numberOfNodes - numberOfNodesOld) +
                             " vertices to the network \n",
                         "net update");
      }
    } // if zero processor
  }   // if network update step

  // adopt radius
  if ((d_update_number + 1) % d_update_interval == 0) {

    if (d_procRank == 0) {
      d_model_p->d_log("Adapt radius \n", "net update");
      adaptRadius();
    }
  } // loop for adopt radius

  // mark nodes too close to boundary as dirichlet
  if (d_update_number % d_update_interval == 0 or
      (d_update_number + 1) % d_update_interval == 0) {
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while( pointer ) {

         int numberOfNeighbors = pointer->neighbors.size();

         if( numberOfNeighbors == 1 ) {

             const auto &coord = pointer->coord;

             // if node is near the boundary, we do not process Omega = (0,L)^3 
             if( 0.00001>coord[0] || coord[0]>L_x-0.00001 || 0.00001>coord[1] || coord[1]>L_x-0.00001 || 0.00001>coord[2] || coord[2]>L_x-0.00001 ){

                   pointer->typeOfVGNode = TypeOfNode::DirichletNode;

             }

         }

         pointer = pointer->global_successor;

  }

  pointer = VGM.getHead();
}

  // compute and communicate updated network
  if (d_update_number % d_update_interval == 0 or
      (d_update_number + 1) % d_update_interval == 0) {
    check_vessel_length();
    prepare_and_communicate_network();
  }
}

void util::unet::Network::linkTerminalVessels() {

  const auto &input = d_model_p->get_input_deck();

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

        if (dist_plane > 0.05 && index != index_1 && dist < h_3D &&
            dist > 0.0 && length_dir > 0.0) {

          if (numberOfNeighbors_1 < 3) {

            oss << " " << std::endl;
            oss << "dist: " << dist << "\n";
            oss << "index: " << index << "\n";
            oss << "index_1: " << index_1 << "\n";
            oss << "Update pointer"
                << "\n";
            d_model_p->d_log(oss, "net update");

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

            double length =
                util::dist_between_points(pointer->coord, pointer_1->coord);
            double p_node = pointer->p_v;
            double p_neighbor = pointer_1->p_v;
            double delta_p = p_neighbor - p_node;
            double tau_w_ini = (radius * std::abs(delta_p)) / (2.0 * length);

            pointer->tau_w_initial.push_back(tau_w_ini);

            oss << "Update pointer 1 \n";
            d_model_p->d_log(oss, "net update");

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

  // Reset values
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

void util::unet::Network::markApicalGrowth() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1) {

      const auto &coord = pointer->coord;

      // if node is near the boundary, we do not process Omega = (0,L)^3
      if (0.00001 < coord[0] && coord[0] < L_x - 0.00001 &&
          0.00001 < coord[1] && coord[1] < L_x - 0.00001 &&
          0.00001 < coord[2] && coord[2] < L_x - 0.00001) {

        int index = getElementIndex(coord, h_3D, N_3D);

        double taf_node = phi_TAF_3D[index];

        // std::cout << "taf_node: " << taf_node << std::endl;

        if (taf_node > input.d_network_update_taf_threshold) {

          // std::cout << "index: " << pointer->index << std::endl;

          pointer->apicalGrowth = true;
        }
      } else { // fix Dirichlet boundaries

        oss << "index: " << pointer->index << std::endl;
        oss << "coord: " << pointer->coord << std::endl;
        oss << " " << std::endl;
        d_model_p->d_log(oss, "net update");

        pointer->typeOfVGNode = TypeOfNode::DirichletNode;
      }
    }

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::processApicalGrowth() {

  const auto &input = d_model_p->get_input_deck();

  // mark node for growth based on a certain criterion
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  oss << "Number of nodes before: " << VGM.getNumberOfNodes() << "\n";
  d_model_p->d_log(oss, "net update");

  int counter = 0;

  while (pointer) {

    if (pointer->apicalGrowth) {

      oss << " \n";
      oss << "Processing node: " << pointer->index << "\n";

      std::vector<double> coord = pointer->coord;
      oss << "Compute direction based on TAF\n";

      int element_index = getElementIndex(coord, h_3D, N_3D);
      oss << "element_index: " << element_index << "\n";
      oss << "counter: " << counter << "\n";
      d_model_p->d_log(oss, "net update");

      std::vector<int> indicesNeighbors =
          getNeighboringElementIndices(element_index, N_3D, h_3D, L_x);
      std::vector<double> TAF_neighbors;

      double TAF_point = phi_TAF_3D[element_index];

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

      // TAG: Random
      auto log_dist = d_logNormalDist();

      double radius_p = pointer->radii[0];

      // get length
      double length = log_dist * radius_p;

      if (radius_p < input.d_min_radius) {

        radius_p = input.d_min_radius;
      }

      if (length > 3.0 * h_3D) {

        length = 3.0 * h_3D;
        d_model_p->d_log("Length is shortened !!!", "net update");
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

        theta = ((double)i_theta) / ((double)N_theta) * 2.0 * M_PI;

        for (int j_r = 1; j_r < N_r; j_r++) {

          double r = ((double)j_r) / ((double)N_r) * length *
                     std::tan(70.0 / 180.0 * M_PI);

          std::vector<double> cylinder_node = computeNodesOnCylinders(
              normal_plane, rotator, midpoint, r, theta);

          if (isCenterInDomain(cylinder_node, L_x) &&
              util::definitelyGreaterThan(length_rotator, 0.)) {

            int index_cone = getElementIndex(cylinder_node, h_3D, N_3D);

            double TAF = phi_TAF_3D[index_cone];

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
          }
        } // radial loop
      }   // tangential loop

      std::vector<double> direction(3, 0.0);

      if (TAF_point < TAF_max) {

        for (int i = 0; i < 3; i++) {

          direction[i] = (max_vec[i] - coord[i]) + (0.65 * normal_plane[i]);
        }

      } else {

        oss << "Find nearest node!" << std::endl;

        std::vector<double> near_node(3, 0.0);

        near_node = findNearNetworkNode(coord, normal_plane);

        oss << "near_node: " << near_node << std::endl;
        d_model_p->d_log(oss, "net update");

        for (int i = 0; i < 3; i++) {

          direction[i] = (0.5 * max_vec[i] + 0.5 * max_vec_2[i] - coord[i]) +
                         (1.0 * normal_plane[i]);
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

      bool isIntersecting =
          testIntersection(coord, new_point_1, radius_p, pointer);

      oss << "coord: " << coord << "\n";
      oss << "rotator: " << rotator << "\n";
      oss << "max_vec: " << max_vec << "\n";
      oss << "normal_plane: " << normal_plane << "\n";
      oss << "length_d: " << length_d << "\n";
      oss << "dir_term_vessel: " << dir_term_vessel << "\n";
      oss << "direction: " << direction << "\n";
      oss << "length: " << length << "\n";
      oss << "new_point: " << new_point_1 << "\n";
      oss << "TAF_point: " << TAF_point << "\n";
      d_model_p->d_log(oss, "net update");

      double global_max_TAF = gmm::vect_norminf(phi_TAF_3D);
      /*
                     std::cout << "global_max_TAF: " << global_max_TAF << "\n";
                     std::cout << "TAF_max: " << TAF_max << "\n";
                     std::cout << "TAF_min: " << TAF_min << "\n";
      */
      // check if we bifurcate at this node
      bool bifurcate = false;

      double prob =
          0.5 + 0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) /
                               std::sqrt(2.0 * input.d_log_normal_std_dev *
                                         input.d_log_normal_std_dev));
      /*
                     std::cout << "prob: " << prob << "\n";
                     std::cout << "input.d_network_bifurcate_prob: " <<
         input.d_network_bifurcate_prob << "\n";
      */
      if (prob > input.d_network_bifurcate_prob) {

        bifurcate = true;
      }

      if (!bifurcate && length_d > 0.0 && length > 0.0) {

        if (!isIntersecting) {

          createASingleNode(new_point_1, radius_p, pointer);
          counter++;
        }

      } else if (bifurcate && radius_p > input.d_min_radius) {

        // std::cout << "Create bifuraction" << "\n";

        double gamma = input.d_net_radius_exponent_gamma;
        double R_c = std::pow(2.0, -1.0 / gamma) * radius_p;

        if (R_c > radius_p) {

          R_c = radius_p;
        }

        // TAG: Random
        double radius_b1 =
            util::transform_to_normal_dist(R_c, R_c / 35.0, d_normalDist());
        double radius_b2 =
            util::transform_to_normal_dist(R_c, R_c / 35.0, d_normalDist());

        if (radius_b1 < input.d_min_radius) {

          radius_b1 = input.d_min_radius;
        }

        if (radius_b2 < input.d_min_radius) {

          radius_b2 = input.d_min_radius;
        }

        std::vector<double> new_point_2 = std::vector<double>(3, 0.0);

        double branch_angle_1 = 0.0;
        double branch_angle_2 = 0.0;

        double angle_arg_1 =
            ((radius_p * radius_p * radius_p * radius_p) +
             (radius_b1 * radius_b1 * radius_b1 * radius_b1) -
             (radius_b2 * radius_b2 * radius_b2 * radius_b2)) /
            (2.0 * radius_p * radius_p * radius_b1 * radius_b1);

        double angle_arg_2 =
            ((radius_p * radius_p * radius_p * radius_p) +
             (radius_b2 * radius_b2 * radius_b2 * radius_b2) -
             (radius_b1 * radius_b1 * radius_b1 * radius_b1)) /
            (2.0 * radius_p * radius_p * radius_b2 * radius_b2);

        if (std::abs(angle_arg_1) < 1.0 && std::abs(angle_arg_2) < 1.0) {

          branch_angle_1 = std::acos(angle_arg_1);
          branch_angle_2 = std::acos(angle_arg_2);

          oss << "radius_p: " << radius_p << "\n";
          oss << "radius_b1: " << radius_b1 << "\n";
          oss << "radius_b2: " << radius_b2 << "\n";
          oss << "branch_angle_1: " << branch_angle_1 * 180.0 / M_PI << "\n";
          oss << "branch_angle_2: " << branch_angle_2 * 180.0 / M_PI << "\n";
          d_model_p->d_log(oss, "net update");

          double branch_angle = branch_angle_1 + branch_angle_2;

          if (branch_angle * 180.0 / M_PI < 160.0 &&
              branch_angle * 180.0 / M_PI > 40.0) {

            std::vector<double> rotation_axis =
                util::cross_prod(direction, dir_term_vessel);
            std::vector<double> diff_2 =
                util::rotate(normal_plane, branch_angle_2, rotation_axis);
            std::vector<double> diff_1 =
                util::rotate(normal_plane, -branch_angle_1, rotation_axis);

            if (util::dist_between_points(diff_1, direction) <
                util::dist_between_points(diff_2, direction)) {

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

            // TAG: Random
            log_dist = d_logNormalDist();

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

                new_point_1[i] =
                    coord[i] + (length_1 * diff_1[i] / length_diff_1);
                new_point_2[i] =
                    coord[i] + (length_2 * diff_2[i] / length_diff_2);
              }
            }

            oss << "rotation_axis: " << rotation_axis << "\n";
            oss << "branch_angle: " << branch_angle << "\n";
            oss << "length_1: " << length_1 << "\n";
            oss << "length_2: " << length_2 << "\n";
            oss << "new_point_1: " << new_point_1 << "\n";
            oss << "new_point_2: " << new_point_2 << "\n";
            d_model_p->d_log(oss, "net update");

            if (gmm::vect_norm2(direction) > 0.0 && length_diff_2 > 0.0 &&
                length_diff_1 > 0.0) {

              bool isIntersecting_1 =
                  testIntersection(coord, new_point_1, radius_p, pointer);
              bool isIntersecting_2 =
                  testIntersection(coord, new_point_2, radius_p, pointer);

              if (!isIntersecting_1) {

                createASingleNode(new_point_1, radius_b1, pointer);
              }

              if (!isIntersecting_2) {

                createASingleNode(new_point_2, radius_b2, pointer);
              }
            }
          }
        }
      } // if bifurcate
    }   // if apical

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

bool util::unet::Network::testIntersection(
    std::vector<double> point_1, std::vector<double> point_2, double radius,
    std::shared_ptr<VGNode> &pointer_test) {

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

        if (util::dist_between_points(point_1, point_p1) < 1.0e-11 ||
            util::dist_between_points(point_2, point_p1) < 1.0e-11 ||
            util::dist_between_points(point_1, point_p2) < 1.0e-11 ||
            util::dist_between_points(point_2, point_p2) < 1.0e-11) {

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

          for (int j = 1; j < 20; j++) {

            std::vector<double> new_point_p(3, 0.0);
            std::vector<double> new_point_test(3, 0.0);

            for (int k = 0; k < 3; k++) {

              new_point_p[k] =
                  point_p1[k] + (length_p * (double)j / 20.0 * dir_p[k]);
              new_point_test[k] =
                  point_1[k] + (length_test * (double)j / 20.0 * dir_test[k]);
            }

            if (util::dist_between_points(new_point_p, new_point_test) <
                radius + radius_p) {

              isIntersecting = true;

              oss << "Vessel is intersecting" << std::endl;
              d_model_p->d_log(oss, "net update");

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

void util::unet::Network::createALinkingNode(std::vector<double> new_point,
                                             double radius,
                                             std::shared_ptr<VGNode> &pointer) {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  bool isInside = isCenterInDomain(new_point, L_x);

  // std::cout << "Create new linking node"  << "\n";

  if (isInside) {

    VGNode new_node;

    new_node.index = VGM.getNumberOfNodes();
    new_node.coord = new_point;
    new_node.p_boundary = 0.95 * pointer->p_v;
    new_node.p_v = 0.95 * pointer->p_v;
    new_node.c_boundary = input.d_in_nutrient;
    new_node.c_v = pointer->c_v; // 0.0;
    new_node.typeOfVGNode = TypeOfNode::NeumannNode;
    new_node.apicalGrowth = false;
    new_node.radii.push_back(radius);
    new_node.radii_initial.push_back(radius);
    new_node.L_p.push_back(input.d_tissue_flow_L_p);
    new_node.L_s.push_back(input.d_tissue_nut_L_s);
    new_node.edge_touched.push_back(true);
    new_node.sprouting_edge.push_back(false);
    new_node.neighbors.push_back(pointer);
    new_node.notUpdated = 0;

    auto sp_newNode = std::make_shared<VGNode>(new_node);
    /*
        std::cout << "New index: " << new_node.index << "\n";
        std::cout << "Neighbor index: " << pointer->index << "\n";

        std::cout << "Update old node" << "\n";
    */

    pointer->p_boundary = 0.0;
    pointer->c_boundary = 0.0;
    pointer->typeOfVGNode = TypeOfNode::InnerNode;
    pointer->apicalGrowth = false;
    pointer->radii.push_back(radius);
    new_node.radii_initial.push_back(radius);
    pointer->L_p.push_back(input.d_tissue_flow_L_p);
    pointer->L_s.push_back(input.d_tissue_nut_L_s);
    pointer->edge_touched.push_back(true);
    pointer->sprouting_edge.push_back(false);
    pointer->neighbors.push_back(sp_newNode);
    pointer->notUpdated = 0;
    //   std::cout << "Attach new node as pointer" << "\n";

    VGM.attachPointerToNode(sp_newNode);
  }
}

void util::unet::Network::createASingleNode(std::vector<double> new_point,
                                            double radius,
                                            std::shared_ptr<VGNode> &pointer) {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  bool isInside = isCenterInDomain(new_point, L_x);

  //     std::cout << "Create new node" << "\n";

  if (isInside) {

    VGNode new_node;

    new_node.index = VGM.getNumberOfNodes();
    new_node.coord = new_point;
    new_node.p_boundary = 0.95 * pointer->p_v;
    new_node.p_v = 0.95 * pointer->p_v;
    new_node.c_boundary = 0.0;
    new_node.c_v = pointer->c_v; // 0.0;
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
    // new_node.typeOfVGNode = TypeOfNode::DirichletNode;

    double length = util::dist_between_points(pointer->coord, new_point);
    double p_node = pointer->p_v;
    double p_neighbor = 0.95 * pointer->p_v;
    double delta_p = p_neighbor - p_node;
    double tau_w_ini = (radius * std::abs(delta_p)) / (2.0 * length);

    new_node.tau_w_initial.push_back(tau_w_ini);

    auto sp_newNode = std::make_shared<VGNode>(new_node);
    /*
             std::cout << "New index: " << new_node.index << "\n";
             std::cout << "Neighbor index: " << pointer->index << "\n";

             std::cout << "Update old node" << "\n";
    */
    pointer->p_boundary = 0.0;
    pointer->c_boundary = 0.0;
    // pointer->c_v = 0.0;
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
    //	 std::cout << "Attach new node as pointer" << "\n";

    VGM.attachPointerToNode(sp_newNode);
  }
}

bool util::unet::Network::testCollision(std::vector<double> point) {

  bool isColliding = false;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  //  std::cout << "Test Collision" << "\n";

  while (pointer) {

    std::vector<double> coord = pointer->coord;
    std::vector<double> diff = std::vector<double>(3, 0.0);

    for (int i = 0; i < 3; i++) {

      diff[i] = coord[i] - point[i];
    }

    double dist = gmm::vect_norm2(diff);

    if (dist < h_3D) {

      // std::cout << "Node not inserted, dist: " << dist << "\n";

      isColliding = true;

      break;
    }

    pointer = pointer->global_successor;
  }

  return isColliding;
}

void util::unet::Network::removeRedundantTerminalVessels() {

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1) {

      const auto &coord = pointer->coord;

      // if node is near the boundary, we do not process Omega = (0,L)^3
      if (0.0001 < coord[0] && coord[0] < L_x - 0.0001 && 0.0001 < coord[1] &&
          coord[1] < L_x - 0.0001 && 0.0001 < coord[2] &&
          coord[2] < L_x - 0.0001) {

        int updateNumber = pointer->notUpdated;

        pointer->notUpdated = updateNumber + 1;

        oss << "pointer->notUpdated: " << pointer->notUpdated << std::endl;
        d_model_p->d_log(oss, "net update");
      }
    }

    pointer = pointer->global_successor;
  }

  d_model_p->d_log(" \n", "net update");

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1 && pointer->notUpdated > 2) {

      int index = pointer->index;

      oss << "Remove node with index: " << index << std::endl;
      oss << "coord: " << pointer->coord << std::endl;

      int numberOfNeighbors_neighbor = pointer->neighbors[0]->neighbors.size();

      int index_neighbor = pointer->neighbors[0]->index;

      oss << "numberOfNeighbors: " << pointer->neighbors.size() << std::endl;
      oss << "numberOfNeighbors_neighbor: " << numberOfNeighbors_neighbor
          << std::endl;
      oss << "index_neighbor: " << index_neighbor << std::endl;
      d_model_p->d_log(oss, "net update");

      std::vector<std::shared_ptr<VGNode>> new_neighbors;

      std::vector<bool> new_edge_touched;
      std::vector<bool> new_sprouting_edge;
      std::vector<double> new_radii;
      std::vector<double> new_L_p;
      std::vector<double> new_L_s;
      std::vector<double> new_tau_w_initial;

      for (int i = 0; i < numberOfNeighbors_neighbor; i++) {

        if (index != pointer->neighbors[0]->neighbors[i]->index) {

          new_neighbors.push_back(pointer->neighbors[0]->neighbors[i]);
          new_edge_touched.push_back(pointer->neighbors[0]->edge_touched[i]);
          new_sprouting_edge.push_back(
              pointer->neighbors[0]->sprouting_edge[i]);
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

        pointer->global_predecessor->global_successor =
            pointer->global_successor;

      } else {

        pointer->global_successor->global_predecessor = NULL;
      }

      if (pointer->global_successor) {

        pointer->global_successor->global_predecessor =
            pointer->global_predecessor;

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

  oss << " " << std::endl;
  oss << "Number of nodes after removing redundant nodes: " << numberOfNodes
      << std::endl;
  d_model_p->d_log(oss, "net update");

  // renumber nodes
  // std::cout << "Renumber nodes" << std::endl;

  pointer = VGM.getHead();

  int counter = 0;

  while (pointer) {

    pointer->index = counter;

    counter = counter + 1;

    const auto &coord = pointer->coord;

    if (0.0001 < coord[0] && coord[0] < L_x - 0.0001 && 0.0001 < coord[1] &&
        coord[1] < L_x - 0.0001 && 0.0001 < coord[2] &&
        coord[2] < L_x - 0.0001) {

      if (pointer->typeOfVGNode == TypeOfNode::InnerNode &&
          pointer->neighbors.size() == 1) {

        pointer->typeOfVGNode = TypeOfNode::NeumannNode;
      }
    }

    pointer = pointer->global_successor;
  }

  // Reset values
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

void util::unet::Network::markSproutingGrowth() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    for (int i = 0; i < numberOfEdges; i++) {

      int numberOfEdges_neighbor = pointer->neighbors[i]->neighbors.size();

      if (pointer->edge_touched[i] == false && numberOfEdges > 1 &&
          numberOfEdges_neighbor > 1) {

        double radius = pointer->radii[i];

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        int local_index =
            pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        double sproutingProbability = 0.0;

        std::vector<double> diff = std::vector<double>(3, 0.0);
        std::vector<double> mid_point = std::vector<double>(3, 0.0);

        for (int j = 0; j < 3; j++) {

          diff[j] = coord_neighbor[j] - coord[j];
          mid_point[j] = 0.5 * (coord_neighbor[j] + coord[j]);
        }

        int element_index = getElementIndex(mid_point, h_3D, N_3D);

        double TAF = phi_TAF_3D[element_index];

        double TAF_th = input.d_network_update_taf_threshold;

        double length = gmm::vect_norm2(diff);

        // TAG: Random
        double log_dist = d_logNormalDist();

        sproutingProbability =
            0.5 +
            0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) /
                           std::sqrt(2.0 * input.d_log_normal_std_dev *
                                     input.d_log_normal_std_dev));

        // double global_max_TAF = gmm::vect_norminf(phi_TAF_3D);
        /*
                            std::cout << "TAF: " << TAF << std::endl;
                            std::cout << "TAF_th: " << TAF_th << std::endl;
                            std::cout << "sproutingProbability: " <<
           sproutingProbability << std::endl; std::cout << "coord: " << coord <<
           std::endl;
        */

        oss << " " << std::endl;
        oss << "radius_initial: " << pointer->radii_initial[i] << std::endl;
        oss << "radius: " << radius << std::endl;
        oss << "coord: " << coord << std::endl;
        d_model_p->d_log(oss, "net update");

        if (sproutingProbability > input.d_sprouting_prob && TAF > TAF_th) {

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

void util::unet::Network::processSproutingGrowth() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double gamma = input.d_net_radius_exponent_gamma;

  // std::cout << " " << std::endl;

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    for (int i = 0; i < numberOfEdges; i++) {

      if (pointer->sprouting_edge[i] == true &&
          pointer->edge_touched[i] == false && numberOfEdges > 1) {

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        std::vector<double> mid_point = std::vector<double>(3, 0.0);

        double radius = pointer->radii[i];

        double PSI = 1.05;

        double radius_prime =
            std::pow((std::pow(PSI, gamma) - 1.0), 1.0 / gamma) * radius;

        double radius_min = input.d_min_radius;

        double p_v_neighbor = pointer->neighbors[i]->p_v;

        double radius_new = 1.25 * radius_min;

        if (1.25 * radius_min < radius_prime) {

          // TAG: Random
          radius_new = util::transform_to_uniform_dist(
              1.25 * radius_min, radius_prime, d_uniformDist());
        }

        if (radius_new > 0.035) {

          radius_new = 0.035;
        }

        std::vector<double> dir_vessel = std::vector<double>(3, 0.0);
        std::vector<double> max_vec = std::vector<double>(3, 0.0);

        for (int j = 0; j < 3; j++) {

          mid_point[j] = 0.5 * (coord_neighbor[j] + coord[j]);
          dir_vessel[j] = coord_neighbor[j] - coord[j];
        }

        // std::cout<< "mid_point: " << mid_point << std::endl;

        double length_vessel = normVector(dir_vessel);

        std::vector<double> normed_dir_vessel = std::vector<double>(3, 0.0);

        for (int j = 0; j < 3; j++) {

          normed_dir_vessel[j] = dir_vessel[j] / length_vessel;
        }

        std::vector<double> rotator = determineRotator(dir_vessel);

        double length_rotator = normVector(rotator);

        int N_s = 10;
        int N_theta = 35;

        std::vector<double> vessel_point(3, 0.0);
        std::vector<double> max_vessel_point(3, 0.0);

        double theta = 0.0;

        double TAF_vessel_surf = 0.0;
        double Max_TAF_vessel_surf = 0.0;

        double counter = 0.0;

        for (int i_s = 0; i_s < N_s + 1; i_s++) {

          for (int j = 0; j < 3; j++) {

            vessel_point[j] =
                (double)i_s / (double)N_s * (0.75 - 0.25) * length_vessel *
                    normed_dir_vessel[j] +
                (0.25 * length_vessel * normed_dir_vessel[j] + coord[j]);
          }

          double TAF_vessel_surf = 0.0;

          for (int i_theta = 0; i_theta < N_theta; i_theta++) {

            theta = ((double)i_theta) / ((double)N_theta) * 2.0 * M_PI;

            std::vector<double> cylinder_node = computeNodesOnCylinders(
                normed_dir_vessel, rotator, vessel_point, radius, theta);

            if (isCenterInDomain(cylinder_node, L_x) &&
                util::definitelyGreaterThan(length_rotator, 0.)) {

              int index = getElementIndex(cylinder_node, h_3D, N_3D);

              double TAF = phi_TAF_3D[index];

              TAF_vessel_surf = TAF_vessel_surf + TAF;

              counter = counter + 1.0;
            }
          }

          if (counter > 0.0) {

            TAF_vessel_surf = TAF_vessel_surf / counter;
          }

          if (TAF_vessel_surf > Max_TAF_vessel_surf) {

            Max_TAF_vessel_surf = TAF_vessel_surf;

            max_vessel_point = vessel_point;
          }

          if (Max_TAF_vessel_surf < 1.0e-16) {

            Max_TAF_vessel_surf = TAF_vessel_surf;

            max_vessel_point = vessel_point;
          }
        }

        // std::cout<< "max_vessel_point: " << max_vessel_point << std::endl;

        theta = 0.0;

        double TAF_max = 0.0;

        for (int i_theta = 0; i_theta < N_theta; i_theta++) {

          theta = ((double)i_theta) / ((double)N_theta) * 2.0 * M_PI;

          std::vector<double> cylinder_node =
              computeNodesOnCylinders(normed_dir_vessel, rotator,
                                      max_vessel_point, 2.0 * radius, theta);

          if (isCenterInDomain(cylinder_node, L_x) &&
              util::definitelyGreaterThan(length_rotator, 0.)) {

            int index_cone = getElementIndex(cylinder_node, h_3D, N_3D);

            double TAF = phi_TAF_3D[index_cone];

            if (TAF_max < 1.0e-16) {

              TAF_max = TAF;
            }

            if (TAF > TAF_max - 1.0e-8) {

              TAF_max = TAF;

              max_vec = cylinder_node;
            }
          }
        }

        std::vector<double> dir_new_vessel = std::vector<double>(3, 0.0);
        std::vector<double> new_point = std::vector<double>(3, 0.0);

        for (int j = 0; j < 3; j++) {

          dir_new_vessel[j] = max_vec[j] - max_vessel_point[j];
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
        double angle_deg = angle * 180. / M_PI;

        // TAG: Random
        double log_dist = d_logNormalDist();

        // get length
        double length_new = log_dist * radius_new;

        if (length_new < 2.0 * radius) {

          length_new = 2.0 * radius;
        }

        if (length_new > 3.0 * h_3D) {

          length_new = 3.0 * h_3D;
        }

        for (int j = 0; j < 3; j++) {

          new_point[j] = mid_point[j] + length_new * dir_new_vessel[j];
        }

        if (util::definitelyGreaterThan(
                util::dist_between_points(new_point, mid_point), 3.0 * h_3D))
          libmesh_error_msg(
              "Error new vessel length " +
              std::to_string(util::dist_between_points(new_point, mid_point)) +
              " is larger than permissible value " + std::to_string(3. * h_3D));
        if (util::definitelyLessThan(
                util::dist_between_points(new_point, mid_point), 1.e-8))
          libmesh_error_msg(
              "Error new vessel length " +
              std::to_string(util::dist_between_points(new_point, mid_point)) +
              " is too small");

        bool isColliding = testCollision(new_point);

        if (angle * 180.0 / M_PI > 10.0 && angle * 180.0 / M_PI < 170.0 &&
            length_vessel > 0.0 && !isColliding && length_vessel > 0.13) {

          d_model_p->d_log("New vessel with length = " +
                               std::to_string(util::dist_between_points(
                                   new_point, mid_point)) +
                               "\n",
                           "net update");
          d_model_p->d_log("Create new node_1\n", "net update");
          VGNode new_node_1;

          new_node_1.index = VGM.getNumberOfNodes();
          new_node_1.coord = mid_point;
          new_node_1.p_boundary = 0.0;
          new_node_1.p_v = pointer->p_v;
          new_node_1.c_boundary = 0.0;   // input.d_in_nutrient;
          new_node_1.c_v = pointer->c_v; // 0.0;
          new_node_1.typeOfVGNode = TypeOfNode::InnerNode;
          new_node_1.apicalGrowth = false;
          new_node_1.radii_initial.push_back(radius);
          new_node_1.radii_initial.push_back(radius);
          new_node_1.radii.push_back(radius);
          new_node_1.radii.push_back(radius);
          /*
                                  if( std::pow( ( std::pow( radius,gamma
             )-std::pow( radius_new,gamma ) ), 1.0/gamma )>radius ){

                                      new_node_1.radii.push_back( std::pow( (
             std::pow( radius,gamma )-std::pow( radius_new,gamma ) ), 1.0/gamma
             ) ); new_node_1.radii.push_back( radius );

                                      new_node_1.radii_initial.push_back(
             std::pow( ( std::pow( radius,gamma )-std::pow( radius_new,gamma )
             ), 1.0/gamma ) ); new_node_1.radii_initial.push_back( radius );

                                  }
                                  else{

                                      new_node_1.radii.push_back( radius );
                                      new_node_1.radii.push_back( std::pow( (
             std::pow( radius,gamma )-std::pow( radius_new,gamma ) ), 1.0/gamma
             ) );

                                      new_node_1.radii_initial.push_back( radius
             ); new_node_1.radii_initial.push_back( std::pow( ( std::pow(
             radius,gamma )-std::pow( radius_new,gamma ) ), 1.0/gamma ) );

                                  }
          */
          new_node_1.radii.push_back(radius_new);

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

          // TODO Both lengths are same. Bug here??
          double length_old =
              util::dist_between_points(mid_point, coord_neighbor);
          double length_new =
              util::dist_between_points(mid_point, coord_neighbor);

          double delta_p_old = p_v_neighbor - pointer->p_v;
          double delta_p_new = 0.95 * pointer->p_v - pointer->p_v;

          double tau_w_ini_old =
              (radius * std::abs(delta_p_old)) / (2.0 * length_old);
          double tau_w_ini_new =
              (radius * std::abs(delta_p_new)) / (2.0 * length_new);

          new_node_1.tau_w_initial.push_back(tau_w_ini_old);
          new_node_1.tau_w_initial.push_back(tau_w_ini_old);
          new_node_1.tau_w_initial.push_back(tau_w_ini_new);

          auto sp_newNode_1 = std::make_shared<VGNode>(new_node_1);

          d_model_p->d_log("Create new node_2\n", "net update");
          VGNode new_node_2;

          new_node_2.index = VGM.getNumberOfNodes() + 1;
          new_node_2.coord = new_point;
          new_node_2.p_boundary = 0.95 * pointer->p_v;
          new_node_2.p_v = 0.95 * pointer->p_v;
          new_node_2.c_boundary = 0.0;   // pointer->c_v;
          new_node_2.c_v = pointer->c_v; // 0.0;
          new_node_2.apicalGrowth = false;
          new_node_2.radii.push_back(radius_new);
          new_node_2.radii_initial.push_back(radius_new);
          new_node_2.L_p.push_back(input.d_tissue_flow_L_p);
          new_node_2.L_s.push_back(input.d_tissue_nut_L_s);
          new_node_2.edge_touched.push_back(true);
          new_node_2.sprouting_edge.push_back(false);
          new_node_2.typeOfVGNode = TypeOfNode::NeumannNode;
          // new_node_2.typeOfVGNode = TypeOfNode::DirichletNode;
          new_node_2.notUpdated = 0;
          new_node_2.tau_w_initial.push_back(tau_w_ini_new);
          new_node_2.neighbors.push_back(sp_newNode_1);

          auto sp_newNode_2 = std::make_shared<VGNode>(new_node_2);

          new_node_1.neighbors.push_back(sp_newNode_2);

          oss << "Update connectivity"
              << "\n";
          d_model_p->d_log(oss, "net update");

          pointer->replacePointerWithIndex(i, sp_newNode_1);

          int index_to_be_replaced =
              pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

          pointer->neighbors[i]->replacePointerWithIndex(index_to_be_replaced,
                                                         sp_newNode_1);

          pointer->sprouting_edge[i] = false;
          pointer->edge_touched[i] = true;

          pointer->neighbors[i]->sprouting_edge[index_to_be_replaced] = false;
          pointer->neighbors[i]->edge_touched[index_to_be_replaced] = true;

          oss << "Attach new node_1 as pointer"
              << "\n";
          VGM.attachPointerToNode(sp_newNode_1);

          oss << "Attach new node_2 as pointer"
              << "\n";
          d_model_p->d_log(oss, "net update");
          VGM.attachPointerToNode(sp_newNode_2);

        } else {

          int localIndex =
              pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

          pointer->edge_touched[i] = true;
          pointer->sprouting_edge[i] = false;

          pointer->neighbors[i]->edge_touched[localIndex] = true;
          pointer->neighbors[i]->sprouting_edge[localIndex] = false;
        }

      } else {

        int localIndex =
            pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        pointer->edge_touched[i] = true;
        pointer->sprouting_edge[i] = false;

        pointer->neighbors[i]->edge_touched[localIndex] = true;
        pointer->neighbors[i]->sprouting_edge[localIndex] = false;
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

std::vector<double>
util::unet::Network::findNearNetworkNode(std::vector<double> coord,
                                         std::vector<double> normal_plane) {

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

    if (dist_plane > 0.0) {

      dist_min = dist_plane;

      numberOfNodes = numberOfNodes + 1;

      for (int i = 0; i < 3; i++) {

        coord_min_node[i] = coord_min_node[i] + coord_near_node[i];
      }
    }

    pointer = pointer->global_successor;
  }

  if (numberOfNodes > 0) {

    for (int i = 0; i < 3; i++) {

      coord_min_node[i] = coord_min_node[i] / ((double)numberOfNodes);
    }
  }

  return coord_min_node;
}

void util::unet::Network::adaptRadius() {

  const auto &input = d_model_p->get_input_deck();

  double dt = d_model_p->d_dt;

  double L_x = input.d_domain_params[1];

  std::shared_ptr<VGNode> pointer = VGM.getHead();

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

    std::vector<bool> remove_edge(numberOfNeighbors, false);
    std::vector<int> edgesToBeRemoved;

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (pointer->edge_touched[i] == false) {

        int local_index =
            pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

        double radius = pointer->radii[i];

        double radius_initial = pointer->radii_initial[i];

        std::vector<double> coord_node = pointer->coord;
        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord_node, coord_neighbor);

        double p_node = pointer->p_v;
        double p_neighbor = pointer->neighbors[i]->p_v;
        double delta_p = p_neighbor - p_node;

        double tau_w = (radius * std::abs(delta_p)) / (2.0 * length);

        double k_WSS = 0.05;

        double S_WSS = 0.0;

        if (1.0 + tau_w > 0.0) {

          S_WSS = std::log(1.0 + tau_w);
        }

        double S_tot = k_WSS * S_WSS;
        double delta_r = dt * radius * S_tot;

        double radius_new = radius + delta_r;

        int numberOfNeighbors_Neighbor =
            pointer->neighbors[i]->neighbors.size();

        if ((radius < input.d_min_radius && numberOfNeighbors < 2) ||
            (numberOfNeighbors_Neighbor == 1 && numberOfNeighbors == 1)) {

          radius_new = input.d_min_radius;

          oss << "Remove vessel with index: " << i << std::endl;
          d_model_p->d_log(oss, "net update");

          edgesToBeRemoved.push_back(i);

          remove_edge[i] = true;
        }

        if (radius_new < 1.2 * radius && radius_new > 0.9 * radius) {

          oss << " " << std::endl;
          oss << "tau_w: " << tau_w << std::endl;
          oss << "S_WSS: " << S_WSS << std::endl;
          oss << "S_tot: " << S_tot << std::endl;
          oss << "radius_initial: " << radius_initial << std::endl;
          oss << "radius: " << radius << std::endl;
          oss << "radius_new: " << radius_new << std::endl;
          oss << "coord_node: " << coord_node << std::endl;
          d_model_p->d_log(oss, "net update");

          // pointer->neighbors[ i ]->radii[ local_index ] = radius_new;
          pointer->radii[i] = radius_new;
        }

        // pointer->neighbors[ i ]->edge_touched[ local_index ] = true;
        pointer->edge_touched[i] = true;
      }
    }
    /*
                // Remove redundant vessels
                for(int i=0;i<numberOfNeighbors;i++){

                    if( remove_edge[ i ] == true ){

                        int index_neighbor = pointer->neighbors[ i
       ]->getLocalIndexOfNeighbor( pointer );

                        pointer->neighbors[ i ]->removeComponent( index_neighbor
       );

                    }

                }

                if( edgesToBeRemoved.size()>0 ){

                    pointer->removeComponents( edgesToBeRemoved );

                }
    */
    pointer = pointer->global_successor;
  }

  // std::cout << " " << std::endl;
  /*
       // Remove nodes without neighbor
       pointer = VGM.getHead();

       while( pointer ){

              int numberOfEdges = pointer->neighbors.size();

              if( numberOfEdges == 0 ){

                  std::cout << "Remove node" << std::endl;

                  std::shared_ptr<VGNode> old_pointer = pointer;

                  if( pointer->global_predecessor ){

                      pointer->global_predecessor->global_successor =
     pointer->global_successor;

                  }
                  else{

                      pointer->global_successor->global_predecessor = NULL;

                  }

                  if( pointer->global_successor ){

                      pointer->global_successor->global_predecessor =
     pointer->global_predecessor;

                  }
                  else{

                      pointer->global_predecessor->global_successor = NULL;

                  }

                  pointer = pointer->global_successor;

                  old_pointer.reset();

              }
              else if( numberOfEdges == 1 ){

                  std::vector<double> coord = pointer->coord;

                  // if node is near the boundary, we do not process Omega =
     (0,L)^3 if( 0.0001<coord[0] && coord[0]<L_x-0.0001 && 0.0001<coord[1] &&
     coord[1]<L_x-0.0001 && 0.0001<coord[2] && coord[2]<L_x-0.0001 ){

                      pointer->typeOfVGNode = TypeOfNode::NeumannNode;

                  }

              }
              else{

                  pointer->typeOfVGNode = TypeOfNode::InnerNode;

              }

              pointer = pointer->global_successor;

       }

       // Renumber nodes
       std::cout << " " << std::endl;
       std::cout << "Renumber nodes" << std::endl;

       pointer = VGM.getHead();

       int counter = 0;

       while( pointer ){

              pointer->index = counter;

              counter = counter + 1;

              pointer = pointer->global_successor;

       }
  */
  VGM.determineNumberOfNodes();

  // Reset
  oss << " " << std::endl;
  oss << "Reset nodes" << std::endl;
  d_model_p->d_log(oss, "net update");

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
