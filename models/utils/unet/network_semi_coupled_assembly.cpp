////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "netUtil.hpp"
#include "network.hpp"

void util::unet::Network::assembleVGMSystemForPressure(BaseAssembly &pres_sys) {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  // factor to enhance condition of matrix
  const double factor_p = input.d_assembly_factor_p_t;

  int numberOfNodes = VGM.getNumberOfNodes();
  if (A_VGM.nrows() != numberOfNodes)
    A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  for (unsigned int i = 0; i < A_VGM.nrows(); i++)
      A_VGM[i].clear();

  if (b.size() != numberOfNodes)
    b.resize(numberOfNodes);

  for (unsigned int i = 0; i < b.size(); i++)
    b[i] = 0.;

  int indexOfNode = 0;
  int numberOfNeighbors = 0;
  int indexNeighbor = 0;
  std::vector<double> coord;
  std::vector<double> coord_neighbor;
  double radius = 0.;
  double length = 0.;
  double L_p = 0.;
  double K_1D = 0.;
  double surface_area = 0.;
  std::vector<double> weights;
  std::vector<int> id_3D_elements;

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    numberOfNeighbors = pointer->neighbors.size();

    // loop over segments and compute 1D and 1D-3D coupling
    for (int i = 0; i < numberOfNeighbors; i++) {

      radius = pointer->radii[i];
      indexNeighbor = pointer->neighbors[i]->index;
      coord_neighbor = pointer->neighbors[i]->coord;
      length = util::dist_between_points(coord, coord_neighbor);
      L_p = pointer->L_p[i];
      K_1D = getK1D(0.5 * (coord[2] + coord_neighbor[2]), L_p, radius);

      // Surface area of cylinder
      surface_area = 2.0 * M_PI * (0.5 * length) * radius;
      determineWeightsAndIds(input.d_num_points_length,
                             input.d_num_points_angle, N_3D, coord,
                             coord_neighbor, radius, h_3D, 0.5 * length,
                             weights, id_3D_elements, mesh, false);

      if (numberOfNeighbors == 1 and
          pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        A_VGM(indexOfNode, indexOfNode) += factor_p * 1.0;
        b[indexOfNode] += factor_p * pointer->p_boundary;

        // we do not consider 1D-3D coupling for this case
        continue;
      }

      // diffusion term
      A_VGM(indexOfNode, indexOfNode) += factor_p * K_1D / length;
      A_VGM(indexOfNode, indexNeighbor) += -factor_p * K_1D / length;

      // Add coupling entry to 1D1D matrix
      A_VGM(indexOfNode, indexOfNode) += factor_p * L_p * surface_area;

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      int numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0)
          continue;

        b[indexOfNode] += factor_p * L_p * surface_area * weights[j] *
                          P_3D[id_3D_elements[j]];
      } // loop over 3D elements
    }   // loop over neighbor segments

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::assembleVGMSystemForNutrient(BaseAssembly &pres_sys,
                                                       BaseAssembly &nut_sys) {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  // factor to enhance condition of matrix
  const double factor_c = input.d_assembly_factor_c_t;

  int numberOfNodes = VGM.getNumberOfNodes();
  if (A_VGM.nrows() != numberOfNodes)
    A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  for (unsigned int i = 0; i < A_VGM.nrows(); i++)
      A_VGM[i].clear();

  if (b_c.size() != numberOfNodes)
    b_c.resize(numberOfNodes);

  for (unsigned int i = 0; i < b_c.size(); i++)
    b_c[i] = 0.;

  std::vector<std::vector<double>> directions;
  directions = defineDirections();
  double dt = d_model_p->d_dt;

  int indexOfNode = 0;
  int numberOfNeighbors = 0;
  int indexNeighbor = 0;
  std::vector<double> coord;
  std::vector<double> coord_neighbor;
  double radius = 0.;
  double length = 0.;
  double L_s = 0.;
  double L_p = 0.;
  double surface_area = 0.;
  std::vector<double> weights;
  std::vector<int> id_3D_elements;
  int numberOfElements = 0;
  double p_v = 0.;
  double p_neighbor = 0.;
  double v_interface = 0.;
  double p_t = 0.0;
  double c_t = 0.0;
  double phi_sigma_boundary = 0.0;

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    p_v = pointer->p_v;
    numberOfNeighbors = pointer->neighbors.size();

    // loop over segments and compute 1D and 1D-3D coupling
    for (int i = 0; i < numberOfNeighbors; i++) {

      radius = pointer->radii[i];
      indexNeighbor = pointer->neighbors[i]->index;
      coord_neighbor = pointer->neighbors[i]->coord;
      length = util::dist_between_points(coord, coord_neighbor);
      p_neighbor = pointer->neighbors[i]->p_v;
      v_interface =
          -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v);
      L_s = pointer->L_s[i];
      L_p = pointer->L_p[i];

      // Surface area of cylinder
      surface_area = 2.0 * M_PI * (0.5 * length) * radius;
      determineWeightsAndIds(input.d_num_points_length,
                             input.d_num_points_angle, N_3D, coord,
                             coord_neighbor, radius, h_3D, 0.5 * length,
                             weights, id_3D_elements, mesh, false);

      if (numberOfNeighbors == 1) {

        if (v_interface > 0.0 &&
            pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

          // mass matrix
          A_VGM(indexOfNode, indexOfNode) += factor_c * length;

          // old time step term
          b_c[indexOfNode] += factor_c * length * C_v_old[indexOfNode];

          if (p_v < input.d_identify_vein_pres)
            phi_sigma_boundary = input.d_in_nutrient_vein;
          else
            phi_sigma_boundary = input.d_in_nutrient;

          // diffusion and advection
          A_VGM(indexOfNode, indexNeighbor) += -factor_c * dt * D_v / length;
          A_VGM(indexOfNode, indexOfNode) +=
              factor_c * (dt * v_interface + 2.0 * dt * D_v / length);

          // advection, diffusion flux due to boundary condition
          b_c[indexOfNode] +=
              factor_c *
              (dt * v_interface - dt * D_v / length) * phi_sigma_boundary;

        } else {

          // mass matrix
          A_VGM(indexOfNode, indexOfNode) += factor_c * length;

          // old time step term
          b_c[indexOfNode] += factor_c * length * C_v_old[indexOfNode];

          // diffusion and advection
          A_VGM(indexOfNode, indexNeighbor) +=
              factor_c * (dt * v_interface - dt * D_v / length);
          A_VGM(indexOfNode, indexOfNode) +=
              -factor_c * dt * v_interface + dt * D_v / length;

          // 1D part of the coupling
          A_VGM(indexOfNode, indexOfNode) +=
              factor_c * dt * L_s * surface_area;
        }
      } // if numberOfNeighbors = 1
      else {

        // mass matrix
        A_VGM(indexOfNode, indexOfNode) += factor_c * 0.5 * length;

        // old time step term
        b_c[indexOfNode] += factor_c * 0.5 * length * C_v_old[indexOfNode];

        // advection term
        if (v_interface > 0.0)
          A_VGM(indexOfNode, indexOfNode) += factor_c * dt * v_interface;
        else
          A_VGM(indexOfNode, indexNeighbor) += factor_c * dt * v_interface;

        // diffusion term
        A_VGM(indexOfNode, indexOfNode) += factor_c * dt * D_v / length;
        A_VGM(indexOfNode, indexNeighbor) += -factor_c * dt * D_v / length;

        // 1D part of the coupling
        A_VGM(indexOfNode, indexOfNode) +=
            factor_c * dt * L_s * surface_area;
      } // if numberOfNeighbors > 1

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0)
          continue;

        c_t = phi_sigma_3D[id_3D_elements[j]];
        b_c[indexOfNode] +=
            factor_c * dt * L_s * surface_area * weights[j] * c_t;

        // osmotic reflection term
        p_t = P_3D[id_3D_elements[j]];

        if (p_v - p_t > 0.0) {

          // 1D equation
          // -2pi R (p_v - p_t) phi_v term in right hand side of 1D equation
          A_VGM(indexOfNode, indexOfNode) +=
              factor_c * dt * (1. - osmotic_sigma) * L_p * surface_area *
              weights[j] * (p_v - p_t);
        } else {

          // 1D equation
          // -2pi R (p_v - p_t) phi_sigma term in right hand side of 1D
          // equation
          b_c[indexOfNode] += -factor_c * dt * (1. - osmotic_sigma) * L_p *
                              surface_area * weights[j] * (p_v - p_t) * c_t;
        }
      } // loop over 3D elements
    }   // loop over neighbor segments

    pointer = pointer->global_successor;
  }
}