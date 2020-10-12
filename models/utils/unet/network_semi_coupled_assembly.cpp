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
  unsigned int assembly_cases;

  double row_sum = 0.;

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    numberOfNeighbors = pointer->neighbors.size();

    // find cases
    assembly_cases = d_vertexBdFlag[indexOfNode];

    std::string assembly_cases_str =
      get_assembly_cases_pres_str(pointer, input.d_identify_vein_pres);

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
                             weights, id_3D_elements,
                             input.d_coupling_3d1d_integration_method, mesh, false);

      // case specific implementation
      if (assembly_cases & UNET_PRES_BDRY_DIRIC) {

        libmesh_assert_equal_to_msg(assembly_cases_str, "boundary_dirichlet",
                                    "Error assembly case " + std::to_string(assembly_cases) + " does not match expected case " + "boundary_dirichlet");

        A_VGM(indexOfNode, indexOfNode) += factor_p * 1.0;
        b[indexOfNode] += factor_p * pointer->p_boundary;
      } else {

        // diffusion term
        A_VGM(indexOfNode, indexOfNode) += factor_p * K_1D / length;
        A_VGM(indexOfNode, indexNeighbor) += -factor_p * K_1D / length;

        // Add coupling entry to 1D1D matrix
        A_VGM(indexOfNode, indexOfNode) += factor_p * L_p * surface_area;
      }

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      int numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0 or assembly_cases & UNET_PRES_BDRY_DIRIC)
          continue;

        b[indexOfNode] += factor_p * L_p * surface_area * weights[j] *
                          P_3D[id_3D_elements[j]];
      } // loop over 3D elements
    }   // loop over neighbor segments

    row_sum = A_VGM(indexOfNode, indexOfNode);
    for (int i = 0; i < numberOfNeighbors; i++) {
      row_sum += A_VGM(indexOfNode, pointer->neighbors[i]->index);
    }

    if (row_sum < 0.)
      libmesh_warning("Network node " + std::to_string(indexOfNode) + " is not diagonally dominated. Sum of row = " + std::to_string(row_sum));

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
  unsigned int assembly_cases;
  bool dirichlet_fixed = false;
  double row_sum = 0.;

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    p_v = pointer->p_v;
    numberOfNeighbors = pointer->neighbors.size();

    // find cases
    //assembly_cases = get_assembly_cases_nut(pointer, input.d_identify_vein_pres);
    assembly_cases = d_vertexBdFlag[indexOfNode];

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
                             weights, id_3D_elements,
                             input.d_coupling_3d1d_integration_method, mesh, false);

      dirichlet_fixed = false;
      // case specific implementation
      if (assembly_cases & UNET_NUT_BDRY_ARTERY_INLET) {

        A_VGM(indexOfNode, indexOfNode) += factor_c * 1.0;
        b_c[indexOfNode] += factor_c * input.d_in_nutrient;

        dirichlet_fixed = true;

      } else if (assembly_cases & UNET_NUT_BDRY_VEIN_INLET) {

        A_VGM(indexOfNode, indexOfNode) += factor_c * 1.0;
        b_c[indexOfNode] += factor_c * input.d_in_nutrient_vein;

        dirichlet_fixed = true;

      } else if (assembly_cases & UNET_NUT_BDRY_ARTERY_OUTLET or
                 assembly_cases & UNET_NUT_BDRY_VEIN_OUTLET) {

        // advection
        A_VGM(indexOfNode, indexOfNode) += -factor_c * dt * v_interface;
        A_VGM(indexOfNode, indexNeighbor) += +factor_c * dt * v_interface;

      } else if (assembly_cases & UNET_NUT_BDRY_INNER_OUTLET or
                 assembly_cases & UNET_NUT_BDRY_INNER_INLET) {

        // advection
        if (v_interface > -1.0e-8)
          A_VGM(indexOfNode, indexOfNode) += factor_c * dt * v_interface;
        else
          A_VGM(indexOfNode, indexNeighbor) += factor_c * dt * v_interface;

      } else if (assembly_cases & UNET_NUT_INNER) {

        // advection term
        if (v_interface > -1.0e-8)
          A_VGM(indexOfNode, indexOfNode) += factor_c * dt * v_interface;
        else
          A_VGM(indexOfNode, indexNeighbor) += factor_c * dt * v_interface;
      }

      // common entries to various cases
      if (dirichlet_fixed == false) {

        // mass matrix
        A_VGM(indexOfNode, indexOfNode) += factor_c * 0.5 * length;

        // diffusion
        A_VGM(indexOfNode, indexOfNode) += factor_c * dt * D_v / length;
        A_VGM(indexOfNode, indexNeighbor) += -factor_c * dt * D_v / length;

        // from previous time step
        b_c[indexOfNode] += factor_c * 0.5 * length * C_v_old[indexOfNode];

        // 1D part of the coupling (Check this term for v > 0 and Dirichlet node)
        A_VGM(indexOfNode, indexOfNode) += factor_c * dt * L_s * surface_area;
      }

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0 or dirichlet_fixed)
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

    row_sum = A_VGM(indexOfNode, indexOfNode);
    for (int i = 0; i < numberOfNeighbors; i++) {
      row_sum += A_VGM(indexOfNode, pointer->neighbors[i]->index);
    }

    if (row_sum < 0.)
      libmesh_warning("Network node " + std::to_string(indexOfNode) + " is not diagonally dominated. Sum of row = " + std::to_string(row_sum));

    pointer = pointer->global_successor;
  }
}