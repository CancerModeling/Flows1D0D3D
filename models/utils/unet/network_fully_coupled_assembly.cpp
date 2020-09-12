////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "netUtil.hpp"
#include "network.hpp"

void util::unet::Network::assemble3D1DSystemForPressure(BaseAssembly &nut_sys,
                                                        BaseAssembly &tum_sys) {

  // get some new 3D systems
  auto &pro = d_model_p->get_assembly("Prolific");
  auto &hyp = d_model_p->get_assembly("Hypoxic");
  auto &ecm = d_model_p->get_assembly("ECM");
  auto &nut = d_model_p->get_assembly("Nutrient");

  const auto &input = d_model_p->get_input_deck();
  const auto &mesh = d_model_p->get_mesh();

  int numberOfNodes = VGM.getNumberOfNodes();

  if (A_flow_3D1D.nrows() != N_tot_3D + numberOfNodes)
    A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
        N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);

  for (int i = 0; i < A_flow_3D1D.nrows(); i++)
    A_flow_3D1D[i].clear();

  if (b_flow_3D1D.size() != N_tot_3D + numberOfNodes)
    b_flow_3D1D.resize(N_tot_3D + numberOfNodes);

  for (int i = 0; i < b_flow_3D1D.size(); i++)
    b_flow_3D1D[i] = 0.0;

  std::vector<std::vector<double>> directions;
  directions = defineDirections();
  double vol_elem = h_3D * h_3D * h_3D;
  double area_face = h_3D * h_3D;
  double area_face_by_h = h_3D;

  // declare local variables so that we reuse them
  int index = 0;
  int index_neighbor = 0;
  auto center = std::vector<double>(3, 0.);
  auto center_neighbor = std::vector<double>(3, 0.);
  bool isInnerFace = false;

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
  bool dirichlet_fixed = false;

  Real nut_cur = 0.;
  Real ecm_cur = 0.;
  Real nut_proj = 0.;
  Real ecm_proj = 0.;
  Real chem_pro_cur = 0.;
  Real chem_hyp_cur = 0.;

  Gradient pro_grad = 0.;
  Gradient hyp_grad = 0.;

  Gradient Sp = 0.;

  double row_sum = 0.;

  // assemble 3D
  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        index = i + j * N_3D + k * N_3D * N_3D;

        // Get element center
        center = getElementCenter(i, j, k, h_3D);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          center_neighbor = getCenterNeighbor(center, directions[face], h_3D);
          isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            // diffusion term
            A_flow_3D1D(index, index) += K_3D * area_face_by_h;
            A_flow_3D1D(index, index_neighbor) += -K_3D * area_face_by_h;
          } else {

            if (scenario == "test_single_vessel") {

              std::vector<double> center_face =
                  getCenterFace(center, directions[face], h_3D);

              double dirichlet_value = getDirichletValue(
                  center_face, VGM.getHead()->L_p[0], VGM.getHead()->radii[0]);

              A_flow_3D1D(index, index) += 2.0 * K_3D * area_face_by_h;
              b_flow_3D1D[index] +=
                  2.0 * K_3D * area_face_by_h * dirichlet_value;
            }
          }
        } // loop over faces

        // Libmesh element
        auto elem = mesh.elem_ptr(index);
        pro.init_dof(elem);
        hyp.init_dof(elem);
        nut.init_dof(elem);
        ecm.init_dof(elem);

        // get finite-volume quantities
        nut_cur = nut.get_current_sol(0);
        nut_proj = util::project_concentration(nut_cur);

        // loop over sides of the element
        for (auto side : elem->side_index_range()) {

          if (elem->neighbor_ptr(side) != nullptr) {

            // div(chem_tum * grad(tum)) term
            // requires integration over face of an element
            pro.d_fe_face->reinit(elem, side);

            // loop over quadrature points
            for (unsigned int qp = 0; qp < pro.d_qrule_face.n_points();
                 qp++) {

              chem_pro_cur = 0.; chem_hyp_cur = 0.;
              pro_grad = 0.; hyp_grad = 0.;
              ecm_cur = 0.; ecm_proj = 0.;
              for (unsigned int l = 0; l < tum_sys.d_phi_face.size(); l++) {

                chem_pro_cur +=
                    pro.d_phi_face[l][qp] * pro.get_current_sol_var(l, 1);

                pro_grad.add_scaled(pro.d_dphi_face[l][qp],
                                    pro.get_current_sol_var(l, 0));

                chem_hyp_cur +=
                    pro.d_phi_face[l][qp] * hyp.get_current_sol_var(l, 1);

                hyp_grad.add_scaled(pro.d_dphi_face[l][qp],
                                    hyp.get_current_sol_var(l, 0));

                ecm_cur +=
                    pro.d_phi_face[l][qp] * ecm.get_current_sol(l);
              }

              ecm_proj = util::project_concentration(ecm_cur);

              if (input.d_assembly_method == 1) {
                Sp = (chem_pro_cur + input.d_chi_c * nut_cur +
                    input.d_chi_h * ecm_cur) *
                     pro_grad +
                     (chem_hyp_cur + input.d_chi_c * nut_cur +
                         input.d_chi_h * ecm_cur) *
                     hyp_grad;
              } else {
                Sp = (chem_pro_cur + input.d_chi_c * nut_proj +
                    input.d_chi_h * ecm_proj) *
                     pro_grad +
                     (chem_hyp_cur + input.d_chi_c * nut_proj +
                         input.d_chi_h * ecm_proj) *
                     hyp_grad;
              }

              // add to force
              b_flow_3D1D[index] += -pro.d_JxW_face[qp] *
                                    input.d_tissue_flow_coeff * Sp *
                                    pro.d_qface_normals[qp];
            } // loop over quadrature points on face

          } // elem neighbor is not null
        }   // loop over faces
      }     // z-loop
    }       // y-loop
  }         // x-loop

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    numberOfNeighbors = pointer->neighbors.size();

    // find cases
    assembly_cases = d_vertexBdFlag[indexOfNode];

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
      dirichlet_fixed = false;
      if (assembly_cases & UNET_PRES_BDRY_DIRIC) {

        if ((!input.d_outlet_apply_neumann and !input.d_inlet_apply_neumann) or
            !pointer->is_initial_node) {

          dirichlet_fixed = true;
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += 1.0;
          b_flow_3D1D[N_tot_3D + indexOfNode] += pointer->p_boundary;
        }
        else {

          // Experiment
          if (pointer->p_v < input.d_identify_vein_pres) {

            // vein's end node
            if (input.d_outlet_apply_neumann) {

              std::cout << "P1D: " << P_3D1D[N_tot_3D + indexOfNode] << "\n";
              if (P_3D1D[N_tot_3D + indexOfNode] < -1.e-8 or
                  P_3D1D[N_tot_3D + indexOfNode] > 1.e10)
                libmesh_error_msg("Error: Pressure should not be negative or "
                                  "very large\n");

              // diffusion term
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
                  K_1D / length;
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
                  -K_1D / length;

              // nonhomogeneous neumann flux
              b_flow_3D1D[N_tot_3D + indexOfNode] +=
                  M_PI * radius * radius * input.d_outlet_neumann_val;

            } else {
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
                  1.0;
              b_flow_3D1D[N_tot_3D + indexOfNode] += pointer->p_boundary;

              dirichlet_fixed = true;
            }
          } // vein's end point
          else {

            // artery's end node
            if (input.d_inlet_apply_neumann) {

              std::cout << "P1D: " << P_3D1D[N_tot_3D + indexOfNode] << "\n";
              if (P_3D1D[N_tot_3D + indexOfNode] < -1.e-8 or
                  P_3D1D[N_tot_3D + indexOfNode] > 1.e10)
                libmesh_error_msg("Error: Pressure should not be negative or "
                                  "very large\n");

              // diffusion term
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
                  K_1D / length;
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
                  -K_1D / length;

              // nonhomogeneous neumann flux
              b_flow_3D1D[N_tot_3D + indexOfNode] +=
                  M_PI * radius * radius * input.d_inlet_neumann_val;

            } else {
              A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
                  1.0;
              b_flow_3D1D[N_tot_3D + indexOfNode] += pointer->p_boundary;

              dirichlet_fixed = true;
            }
          } // artery's end point
        } // apply non-homogeneous Neumann bc

      } else {

        // diffusion term
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            K_1D / length;
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
            -K_1D / length;

        // Add coupling entry to 1D1D matrix
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            L_p * surface_area;
      }

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      int numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0 or assembly_cases & UNET_PRES_BDRY_DIRIC)
          continue;

        // A_3D1D
        A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) +=
            -L_p * surface_area * weights[j];

        // A_3D3D
        A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) +=
            L_p * surface_area * weights[j];

        // A_1D3D
        A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +=
            -L_p * surface_area * weights[j];

      } // loop over 3D elements
    }   // loop over neighbor segments

    row_sum = A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode);
    for (int i = 0; i < numberOfNeighbors; i++) {
      row_sum += A_flow_3D1D(N_tot_3D + indexOfNode,
                             N_tot_3D + pointer->neighbors[i]->index);
    }

    if (row_sum < 0.)
      libmesh_warning("Network node " + std::to_string(indexOfNode) +
                      " is not diagonally dominated. Sum of row = " +
                      std::to_string(row_sum));

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::assemble3D1DSystemForNutrients(
    BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

  // get some new 3D systems
  auto &pro = d_model_p->get_assembly("Prolific");
  auto &hyp = d_model_p->get_assembly("Hypoxic");
  auto &nec = d_model_p->get_assembly("Necrotic");
  auto &ecm = d_model_p->get_assembly("ECM");
  auto &nut = d_model_p->get_assembly("Nutrient");
  auto &taf = d_model_p->get_assembly("TAF");
  auto &mde = d_model_p->get_assembly("MDE");

  const auto &input = d_model_p->get_input_deck();
  const auto &mesh = d_model_p->get_mesh();

  int numberOfNodes = VGM.getNumberOfNodes();
  if (A_flow_3D1D.nrows() != N_tot_3D + numberOfNodes)
    A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
        N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);

  for (int i = 0; i < A_flow_3D1D.nrows(); i++)
    A_flow_3D1D[i].clear();

  if (b_nut_3D1D.size() != N_tot_3D + numberOfNodes)
    b_nut_3D1D.resize(N_tot_3D + numberOfNodes);

  for (int i = 0; i < b_nut_3D1D.size(); i++)
    b_nut_3D1D[i] = 0.0;

  // Store current and old solution
  Real pro_cur = 0.;
  Real hyp_cur = 0.;
  Real nec_cur = 0.;
  Real ecm_cur = 0.;
  Real mde_cur = 0.;

  Real pro_proj = 0.;
  Real hyp_proj = 0.;
  Real nec_proj = 0.;
  Real ecm_proj = 0.;
  Real mde_proj = 0.;

  Real pres_cur = 0.;
  Real nut_old = 0.;
  Real ecm_old = 0.;
  Real chem_pro_old = 0.;
  Real chem_hyp_old = 0.;

  Real nut_old_proj = 0.;
  Real ecm_old_proj = 0.;

  Gradient pro_grad = 0.;
  Gradient hyp_grad = 0.;
  Gradient nec_grad = 0.;

  Gradient pro_old_grad = 0.;
  Gradient hyp_old_grad = 0.;

  Gradient Sp_old = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  std::vector<std::vector<double>> directions;
  directions = defineDirections();
  double vol_elem = h_3D * h_3D * h_3D;
  double area_face = h_3D * h_3D;
  double area_face_by_h = h_3D;
  double dt = d_model_p->d_dt;

  // declare local variables so that we reuse them
  int index = 0;
  int index_neighbor = 0;
  auto center = std::vector<double>(3, 0.);
  auto center_neighbor = std::vector<double>(3, 0.);
  bool isInnerFace = false;
  Number chem_tum_cur = 0.;
  Gradient tum_grad = 0.;
  double v = 0.;
  Number v_mu = 0.;

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
  double phi_sigma_boundary = 0.0;
  unsigned int assembly_cases;
  bool dirichlet_fixed = false;

  double row_sum = 0.;

  // assemble 3D
  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        index = i + j * N_3D + k * N_3D * N_3D;

        A_flow_3D1D(index, index) += vol_elem;
        b_nut_3D1D[index] += vol_elem * phi_sigma_old[index];

        // Get element center
        center = getElementCenter(i, j, k, h_3D);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          center_neighbor = getCenterNeighbor(center, directions[face], h_3D);
          isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            v = -K_3D * (P_3D1D[index_neighbor] - P_3D1D[index]) / h_3D;

            // upwinding
            if (v > 0.0)
              A_flow_3D1D(index, index) += dt * area_face * v;
            else
              A_flow_3D1D(index, index_neighbor) += dt * area_face * v;

            // diffusion term
            A_flow_3D1D(index, index) += dt * D_v_3D * area_face_by_h;
            A_flow_3D1D(index, index_neighbor) += -dt * D_v_3D * area_face_by_h;
          }
        } // loop over faces

        // Libmesh element
        auto elem = mesh.elem_ptr(index);

        // loop over sides of the element
        for (auto side : elem->side_index_range()) {

          if (elem->neighbor_ptr(side) != nullptr) {

            // grad(tum) term
            // requires integration over face of an element
            pro.d_fe_face->reinit(elem, side);

            // loop over quadrature points
            for (unsigned int qp = 0; qp < pro.d_qrule_face.n_points();
                 qp++) {

              chem_pro_old = 0.; chem_hyp_old = 0.;
              pro_grad = 0.; hyp_grad = 0.; nec_grad = 0.;
              pro_old_grad = 0.; hyp_old_grad = 0.;
              ecm_old = 0.; ecm_old_proj = 0.;
              for (unsigned int l = 0; l < pro.d_phi_face.size(); l++) {

                chem_pro_old +=
                    pro.d_phi_face[l][qp] * pro.get_old_sol_var(l, 1);

                pro_grad.add_scaled(pro.d_dphi_face[l][qp],
                                    pro.get_current_sol_var(l, 0));

                pro_old_grad.add_scaled(pro.d_dphi_face[l][qp],
                                        pro.get_old_sol_var(l, 0));

                chem_hyp_old +=
                    pro.d_phi_face[l][qp] * hyp.get_old_sol_var(l, 1);

                hyp_grad.add_scaled(pro.d_dphi_face[l][qp],
                                    hyp.get_current_sol_var(l, 0));

                hyp_old_grad.add_scaled(pro.d_dphi_face[l][qp],
                                        hyp.get_old_sol_var(l, 0));

                nec_grad.add_scaled(pro.d_dphi_face[l][qp],
                                    nec.get_current_sol(l));

                ecm_old +=
                    pro.d_phi_face[l][qp] * ecm.get_old_sol(l);
              }

              ecm_old_proj = util::project_concentration(ecm_old);

              if (input.d_assembly_method == 1) {
                Sp_old = (chem_pro_old + input.d_chi_c * nut_old +
                    input.d_chi_h * ecm_old) *
                         pro_old_grad +
                         (chem_hyp_old + input.d_chi_c * nut_old +
                             input.d_chi_h * ecm_old) *
                         hyp_old_grad;
              } else {
                Sp_old = (chem_pro_old + input.d_chi_c * nut_old_proj +
                    input.d_chi_h * ecm_old_proj) *
                         pro_old_grad +
                         (chem_hyp_old + input.d_chi_c * nut_old_proj +
                             input.d_chi_h * ecm_old_proj) *
                         hyp_old_grad;
              }

              // chemotactic term
              b_nut_3D1D[index] += -pro.d_JxW_face[qp] * dt * input.d_chi_c *
                                   (pro_grad + hyp_grad + nec_grad) *
                                   pro.d_qface_normals[qp];

              // advection term
              v_mu = pro.d_JxW_face[qp] * dt * K_3D * Sp_old *
                     pro.d_qface_normals[qp];

              // goes to the dof of element (not the neighbor)
              A_flow_3D1D(index, index) += v_mu;
            } // loop over quadrature points on face

          } // elem neighbor is not null
        }   // loop over faces

        //
        // compute source terms
        //
        nut.init_dof(elem);
        pro.init_dof(elem);
        hyp.init_dof(elem);
        nec.init_dof(elem);
        taf.init_dof(elem);
        ecm.init_dof(elem);
        mde.init_dof(elem);

        // init fe and element matrix and vector
        hyp.init_fe(elem);

        // loop over quadrature points
        for (unsigned int qp = 0; qp < hyp.d_qrule.n_points(); qp++) {
          // Computing solution
          pro_cur = 0.;
          hyp_cur = 0.;
          nec_cur = 0.;
          ecm_cur = 0.;
          mde_cur = 0.;

          for (unsigned int l = 0; l < hyp.d_phi.size(); l++) {

            pro_cur += hyp.d_phi[l][qp] * pro.get_current_sol_var(l, 0);
            hyp_cur += hyp.d_phi[l][qp] * hyp.get_current_sol_var(l, 0);
            nec_cur += hyp.d_phi[l][qp] * nec.get_current_sol(l);
            ecm_cur += hyp.d_phi[l][qp] * ecm.get_current_sol(l);
            mde_cur += hyp.d_phi[l][qp] * mde.get_current_sol(l);
          }

          if (input.d_assembly_method == 1) {

            compute_rhs = hyp.d_JxW[qp] * dt *
                          (input.d_lambda_A * (pro_cur + hyp_cur) +
                           input.d_lambda_ECM_D * ecm_cur * mde_cur);

            compute_mat =
                hyp.d_JxW[qp] * dt *
                (input.d_lambda_P * pro_cur +
                 input.d_lambda_Ph * hyp_cur +
                 input.d_lambda_ECM_P * (1. - ecm_cur) *
                     util::heaviside(ecm_cur - input.d_bar_phi_ECM_P));

          } else {

            mde_proj = util::project_concentration(mde_cur);
            ecm_proj = util::project_concentration(ecm_cur);
            pro_proj = util::project_concentration(pro_cur);
            hyp_proj = util::project_concentration(hyp_cur);
            nec_proj = util::project_concentration(nec_cur);

            compute_rhs = hyp.d_JxW[qp] * dt *
                          (input.d_lambda_A * (pro_proj + hyp_proj) +
                           input.d_lambda_ECM_D * ecm_proj * mde_proj);

            compute_mat =
                hyp.d_JxW[qp] * dt *
                (input.d_lambda_P * pro_proj +
                 input.d_lambda_Ph * hyp_proj +
                 input.d_lambda_ECM_P * (1. - ecm_proj) *
                 util::heaviside(ecm_proj - input.d_bar_phi_ECM_P));
          }

          // add rhs
          b_nut_3D1D[index] += compute_rhs;

          // add matrix
          A_flow_3D1D(index, index) += compute_mat;

        } // loop over quad points
      }   // z-loop
    }     // y-loop
  }       // x-loop

  // assemble 1D and 1D-3D coupling
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    indexOfNode = pointer->index;
    coord = pointer->coord;
    p_v = pointer->p_v;
    numberOfNeighbors = pointer->neighbors.size();

    // find cases
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

        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += 1.0;
        b_nut_3D1D[N_tot_3D + indexOfNode] += input.d_in_nutrient;

        dirichlet_fixed = true;

      } else if (assembly_cases & UNET_NUT_BDRY_VEIN_INLET) {

        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += 1.0;
        b_nut_3D1D[N_tot_3D + indexOfNode] += input.d_in_nutrient_vein;

        dirichlet_fixed = true;

      } else if (assembly_cases & UNET_NUT_BDRY_ARTERY_OUTLET or
                 assembly_cases & UNET_NUT_BDRY_VEIN_OUTLET) {

        // advection
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            -dt * v_interface;
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
            +dt * v_interface;

      } else if (assembly_cases & UNET_NUT_BDRY_INNER_OUTLET or
                 assembly_cases & UNET_NUT_BDRY_INNER_INLET) {

        // advection
        if (v_interface > -1.0e-8)
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
              dt * v_interface;
        else
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
              dt * v_interface;

      } else if (assembly_cases & UNET_NUT_INNER) {

        // advection term
        if (v_interface > -1.0e-8)
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
              dt * v_interface;
        else
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
              dt * v_interface;
      }

      // common entries to various cases
      if (dirichlet_fixed == false) {

        // mass matrix
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            0.5 * length;

        // diffusion
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            dt * D_v / length;
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
            -dt * D_v / length;

        // from previous time step
        b_nut_3D1D[N_tot_3D + indexOfNode] +=
            0.5 * length * phi_sigma_old[N_tot_3D + indexOfNode];

        // 1D part of the coupling (Check this term for v > 0 and Dirichlet
        // node)
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            dt * L_s * surface_area;
      }

      // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
      numberOfElements = id_3D_elements.size();

      for (int j = 0; j < numberOfElements; j++) {

        if (id_3D_elements[j] < 0 or dirichlet_fixed)
          continue;

        // A_3D1D
        A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) +=
            -dt * L_s * surface_area * weights[j];

        // A_3D3D
        A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) +=
            dt * L_s * surface_area * weights[j];

        // A_1D3D
        A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +=
            -dt * L_s * surface_area * weights[j];

        // osmotic reflection term
        p_t = P_3D[id_3D_elements[j]];

        if (p_v - p_t > 0.0) {

          // 1D equation
          // -2pi R (p_v - p_t) phi_v term in right hand side of 1D equation
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
              dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] *
              (p_v - p_t);

          // A_3D3D
          A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) +=
              -dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] *
              (p_v - p_t);

        } else {

          // 3D equation
          // 2pi R (p_v - p_t) phi_sigma term in right hand side of 3D
          // equation
          A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) +=
              -dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] *
              (p_v - p_t);

          A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +=
              dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] *
              (p_v - p_t);
        }
      } // loop over 3D elements
    }   // loop over neighbor segments

    row_sum = A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode);
    for (int i = 0; i < numberOfNeighbors; i++) {
      row_sum += A_flow_3D1D(N_tot_3D + indexOfNode,
                             N_tot_3D + pointer->neighbors[i]->index);
    }

    if (row_sum < 0.)
      libmesh_warning("Network node " + std::to_string(indexOfNode) +
                      " is not diagonally dominated. Sum of row = " +
                      std::to_string(row_sum));

    pointer = pointer->global_successor;
  }
}
