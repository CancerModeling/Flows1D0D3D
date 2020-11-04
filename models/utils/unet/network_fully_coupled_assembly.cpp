////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "netUtil.hpp"
#include "network.hpp"

void util::unet::Network::assemble3D1DSystemForPressure(BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

  // get some new 3D systems
  auto &pro = d_model_p->get_assembly("Prolific");
  auto &hyp = d_model_p->get_assembly("Hypoxic");
  auto &ecm = d_model_p->get_assembly("ECM");
  auto &nut = d_model_p->get_assembly("Nutrient");

  const auto &input = d_model_p->get_input_deck();
  const auto &mesh = d_model_p->get_mesh();

  int numberOfNodes = VGM.getNumberOfNodes();

  if (A_flow_3D1D.nrows() != N_tot_3D + numberOfNodes)
    libmesh_error_msg("A_flow_3D1D error.");

  for (int i = 0; i < A_flow_3D1D.nrows(); i++) {

    A_flow_3D1D[i].clear();
  }

  for (int i = 0; i < b_flow_3D1D.size(); i++) {

    b_flow_3D1D[i] = 0.0;
  }

  std::vector<std::vector<double>> directions;
  directions = defineDirections();

  Real nut_cur = 0.;
  Real ecm_cur = 0.;
  Real nut_proj = 0.;
  Real ecm_proj = 0.;
  Real chem_pro_cur = 0.;
  Real chem_hyp_cur = 0.;

  Gradient pro_grad = 0.;
  Gradient hyp_grad = 0.;

  Gradient Sp = 0.;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        b_flow_3D1D[index] = 0.0;

        // Get element center
        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          std::vector<double> center_neighbor =
            getCenterNeighbor(center, directions[face], h_3D);

          bool isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            int index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            // diffusion term
            A_flow_3D1D(index, index) =
              A_flow_3D1D(index, index) + K_3D * h_3D * h_3D / h_3D;

            A_flow_3D1D(index, index_neighbor) = -K_3D * h_3D * h_3D / h_3D;

          } else {

            if (scenario == "test_single_vessel") {

              std::vector<double> center_face =
                getCenterFace(center, directions[face], h_3D);

              double dirichlet_value = getDirichletValue(
                center_face, VGM.getHead()->L_p[0], VGM.getHead()->radii[0]);

              A_flow_3D1D(index, index) =
                A_flow_3D1D(index, index) + 2.0 * K_3D * h_3D * h_3D / h_3D;

              b_flow_3D1D[index] =
                b_flow_3D1D[index] + (2.0 * K_3D * h_3D * h_3D * dirichlet_value) / h_3D;
            }
          }
        }

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

              chem_pro_cur = 0.;
              chem_hyp_cur = 0.;
              pro_grad = 0.;
              hyp_grad = 0.;
              ecm_cur = 0.;
              ecm_proj = 0.;
              for (unsigned int l = 0; l < pro.d_phi_face.size(); l++) {

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
      }
    }
  }

  // Assemble 1D Matrix, 1D right hand side and coupling matrices (pressure)
  // std::cout << " " << std::endl;
  // std::cout << "Assemble 1D Matrix and Coupling matrices (pressure)"
  //          << std::endl;

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    //std::cout << "\n----------------------------------------------\n";
    //std::cout << "indexOfNode: " << indexOfNode << std::endl;

    int numberOfNeighbors = pointer->neighbors.size();
    std::vector<double> coord = pointer->coord;

    if (numberOfNeighbors == 1) {

      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = 1.0;
        b_flow_3D1D[N_tot_3D + indexOfNode] = pointer->p_boundary;

      } else if (pointer->typeOfVGNode == TypeOfNode::NeumannNode) {

        double radius = pointer->radii[0];

        //pointer->printInformationOfNode();
        //std::cout << "radius: " << radius << std::endl;
        //std::cout << "numberOfNodes: " << VGM.getNumberOfNodes() << std::endl;
        //std::cout << "A_flow_3D1D: " << A_flow_3D1D.size() << std::endl;

        // Assemble coupling terms (nutrients)
        std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;
        int N_s = input.d_num_points_length;
        int N_theta = input.d_num_points_angle;
        int indexNeighbor = pointer->neighbors[0]->index;
        //std::cout << "indexNeighbor: " << indexNeighbor << std::endl;

        std::vector<double> weights;
        std::vector<int> id_3D_elements;

        // Permeabilty vessel wall for nutrients
        double L_p = pointer->L_p[0];
        double length = util::dist_between_points(coord, coord_neighbor);
        double K_1D = getK1D(0.5 * (coord[2] + coord_neighbor[2]), L_p, radius);
        //std::cout << "N_tot_3D+indexOfNode: " << N_tot_3D+indexOfNode << std::endl;
        //std::cout << "gmm::mat_nrows(A_flow_3D1D): " << gmm::mat_nrows(A_flow_3D1D) << std::endl;

        // diffusion term
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + K_1D / length;
        //std::cout << "A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode): " << A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) << std::endl;

        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = -K_1D / length;

        double length_edge = 0.5 * length;

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

        // Add coupling entry to 1D1D matrix
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + L_p * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();
        //std::cout << "numberOfElements: " << numberOfElements << std::endl;

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {
            // A_3D1D
            A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
              A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
              L_p * surface_area * weights[j];

            // A_3D3D
            A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) =
              A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) +
              L_p * surface_area * weights[j];

            // A_1D3D
            A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
              A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) -
              L_p * surface_area * weights[j];
          }

        } // loop over 3D elements

      } // if neumann node

    } // if number of neighbors = 1
    else {

      for (int i = 0; i < numberOfNeighbors; i++) {

        // Discretisation of the differential operator
        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord, coord_neighbor);
        double L_p = pointer->L_p[i];
        double K_1D = getK1D(0.5 * (coord[2] + coord_neighbor[2]), L_p, radius);

        // diffusion term
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + K_1D / length;
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = -K_1D / length;

        // Coupling terms
        int N_s = input.d_num_points_length;
        int N_theta = input.d_num_points_angle;

        std::vector<double> weights;
        std::vector<int> id_3D_elements;

        double length_edge = 0.5 * length;

        if (pointer->neighbors[i]->neighbors.size() == 1 && pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

          length_edge = length;
        }

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor,
                               radius, h_3D, length_edge, weights,
                               id_3D_elements);


        // Add coupling entry to 1D1D matrix
        A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_flow_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + L_p * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {
            // A_3D1D
            A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
              A_flow_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) - L_p * surface_area * weights[j];

            // A_3D3D
            A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) =
              A_flow_3D1D(id_3D_elements[j], id_3D_elements[j]) + L_p * surface_area * weights[j];

            // A_1D3D
            A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
              A_flow_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) - L_p * surface_area * weights[j];
          }
        } // loop over 3D elements

      } // loop over neighbor segments

      b_flow_3D1D[N_tot_3D + indexOfNode] = 0.0;
    } // if number of neighbors > 1

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::assemble3D1DSystemForNutrients(BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

  const auto &input = d_model_p->get_input_deck();
  const auto &mesh = d_model_p->get_mesh();

  // get some new 3D systems
  auto &pro = d_model_p->get_assembly("Prolific");
  auto &hyp = d_model_p->get_assembly("Hypoxic");
  auto &nec = d_model_p->get_assembly("Necrotic");
  auto &ecm = d_model_p->get_assembly("ECM");
  auto &nut = d_model_p->get_assembly("Nutrient");
  auto &taf = d_model_p->get_assembly("TAF");
  auto &mde = d_model_p->get_assembly("MDE");


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

  Gradient pro_old_grad = 0.;
  Gradient hyp_old_grad = 0.;

  Gradient Sp_old = 0.;

  Real compute_rhs = 0.;
  Real compute_mat = 0.;

  Number v_mu = 0.;

  // 3D-1D coupled flow problem on a cube
  // std::cout << " " << std::endl;
  //  std::cout << "3D-1D coupled nutrient transport problem on a cube \Omega =
  //  (0,"
  //            << L_x << ")^3" << std::endl;

  // Number of Elements (3D) in each space direction
  //  std::cout << " " << std::endl;
  //  std::cout << "Number of Elements (3D) in each space direction: " << N_3D
  //            << std::endl;
  //  std::cout << "Total number of Elements in 3D: " << N_tot_3D << std::endl;

  // Mesh size
  //  std::cout << " " << std::endl;
  //  std::cout << "Mesh size h_3D: " << h_3D << std::endl;

  double vol_elem = h_3D * h_3D * h_3D;
  // std::cout << "vol_elem: " << vol_elem << std::endl;

  double area_face = h_3D * h_3D;
  // std::cout << "area_face: " << area_face << std::endl;

  // Assemble 3D Matrix and right hand side (pressure)
  // std::cout << " " << std::endl;
  double K_3D = input.d_tissue_flow_coeff;
  //  std::cout << "Assemble 3D Matrix and right hand side (nutrients)"
  //            << std::endl;
  //  std::cout << "K_3D: " << K_3D << std::endl;
  //  std::cout << "D_v_3D: " << D_v_3D << std::endl;

  int numberOfNodes = VGM.getNumberOfNodes();

  // std::cout << "numberOfNodes: " << N_tot_3D + numberOfNodes << std::endl;

  A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D + numberOfNodes,
                                                      N_tot_3D + numberOfNodes);

  for (int i = 0; i < A_nut_3D1D.nrows(); i++) {

    A_nut_3D1D[i].clear();
  }

  for (int i = 0; i < b_nut_3D1D.size(); i++) {

    b_nut_3D1D[i] = 0.0;
  }

  std::vector<std::vector<double>> directions;

  directions = defineDirections();

  double dt = d_model_p->d_dt;
  // std::cout << "dt: " << dt << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        A_nut_3D1D(index, index) += vol_elem;

        b_nut_3D1D[index] = vol_elem * phi_sigma_old[index];

        // Get element center
        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          std::vector<double> center_neighbor =
            getCenterNeighbor(center, directions[face], h_3D);

          bool isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            int index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            double v = -K_3D * (P_3D1D[index_neighbor] - P_3D1D[index]) / h_3D;

            // upwinding
            if (v > 0.0) {

              A_nut_3D1D(index, index) += dt * area_face * v;

            } else {

              A_nut_3D1D(index, index_neighbor) += dt * area_face * v;
            }

            // diffusion term
            A_nut_3D1D(index, index) += dt * D_v_3D * area_face / h_3D;

            A_nut_3D1D(index, index_neighbor) += -dt * D_v_3D * area_face / h_3D;
          }
        }

        // Libmesh element
        auto elem = mesh.elem_ptr(index);
        pro.init_dof(elem);
        hyp.init_dof(elem);
        nec.init_dof(elem);
        ecm.init_dof(elem);
        nut.init_dof(elem);
        taf.init_dof(elem);
        mde.init_dof(elem);

        // loop over sides of the element
        for (auto side : elem->side_index_range()) {

          if (elem->neighbor_ptr(side) != nullptr) {

            // grad(tum) term
            // requires integration over face of an element
            pro.d_fe_face->reinit(elem, side);

            // loop over quadrature points
            for (unsigned int qp = 0; qp < pro.d_qrule_face.n_points();
                 qp++) {

              chem_pro_old = 0.;
              chem_hyp_old = 0.;
              pro_grad = 0.;
              hyp_grad = 0.;
              pro_old_grad = 0.;
              hyp_old_grad = 0.;
              ecm_old = 0.;
              ecm_old_proj = 0.;
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
                                   (pro_grad + hyp_grad) *
                                   pro.d_qface_normals[qp];

              // advection term
              v_mu = pro.d_JxW_face[qp] * dt * K_3D * Sp_old *
                     pro.d_qface_normals[qp];

              // goes to the dof of element (not the neighbor)
              A_nut_3D1D(index, index) += v_mu;
            } // loop over quadrature points on face

          } // elem neighbor is not null
        }   // loop over faces

        //
        // compute source terms
        //

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
          A_nut_3D1D(index, index) += compute_mat;

        } // loop over quad points

      } // z - loop
    }
  }

  // Assemble 1D Matrix, 1D right hand side and coupling matrices (nutrients)
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v = pointer->p_v;

    if (numberOfNeighbors == 1) {

      double radius = pointer->radii[0];

      int indexNeighbor = pointer->neighbors[0]->index;

      std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;

      double length = util::dist_between_points(coord, coord_neighbor);

      double p_neighbor = pointer->neighbors[0]->p_v;

      double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v);

      double volume = length / 2.0 * radius * radius * M_PI;

      if (v_interface > 0.0 && pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        double phi_sigma_boundary = 0.0;

        if (p_v < input.d_identify_vein_pres) {

          phi_sigma_boundary = input.d_in_nutrient_vein;

        } else {

          phi_sigma_boundary = input.d_in_nutrient;
        }

        // dirichlet bc
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = 1.0;
        b_nut_3D1D[N_tot_3D + indexOfNode] = phi_sigma_boundary;

      } // if dirichlet node
      else if (v_interface < 0.0 && pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        double length_edge = 0.5 * length;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = length_edge;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = dt * v_interface;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) - dt * v_interface;

        b_nut_3D1D[N_tot_3D + indexOfNode] = length_edge * phi_sigma_old[N_tot_3D + indexOfNode];

        // Assemble coupling terms (nutrients)
        int N_s = input.d_num_points_length;
        int N_theta = input.d_num_points_angle;

        std::vector<double> weights;

        std::vector<int> id_3D_elements;

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor,
                               radius, h_3D, length_edge, weights, id_3D_elements);

        // Permeabilty vessel wall for nutrients
        double L_s = pointer->L_s[0];
        double L_p = pointer->L_p[0];

        // 1D part of the coupling
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + dt * L_s * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();

        double p_t = 0.0;

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {

            // A_3D1D
            A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
              dt * L_s * surface_area * weights[j];

            // A_3D3D
            A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) +
              dt * L_s * surface_area * weights[j];

            // A_1D3D
            A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) -
              dt * L_s * surface_area * weights[j];

            // osmotic reflection term
            p_t = P_3D1D[id_3D_elements[j]];

            if (p_v - p_t > 0.0) {

              // 1D equation
              // -2pi R (p_v - p_t) phi_v term in right hand side of 1D equation
              A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
                A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              // A_3D3D
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
                A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);


            } else {

              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
                A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
                A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);
            }
          }

        } // loop over 3D elements

      } // if outflow dirichlet node
      else {

        double length_edge = 0.5 * length;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = length_edge;

        if (v_interface > -1.0e-8) {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += dt * v_interface;

        } else {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) += dt * v_interface;
        }

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = -dt * D_v / length;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) = A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + dt * D_v / length;

        b_nut_3D1D[N_tot_3D + indexOfNode] = length_edge * phi_sigma_old[N_tot_3D + indexOfNode];

        // Assemble coupling terms (nutrients)

        int N_s = input.d_num_points_length;
        int N_theta = input.d_num_points_angle;

        std::vector<double> weights;

        std::vector<int> id_3D_elements;

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor,
                               radius, h_3D, length_edge, weights, id_3D_elements);

        // Permeabilty vessel wall for nutrients
        double L_s = pointer->L_s[0];
        double L_p = pointer->L_p[0];

        // 1D part of the coupling
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) + dt * L_s * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();

        double p_t = 0.0;

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {

            // A_3D1D
            A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
              dt * L_s * surface_area * weights[j];

            // A_3D3D
            A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) +
              dt * L_s * surface_area * weights[j];

            // A_1D3D
            A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) -
              dt * L_s * surface_area * weights[j];

            // osmotic reflection term
            p_t = P_3D1D[id_3D_elements[j]];

            if (p_v - p_t > 0.0) {

              // -2pi R (p_v - p_t) phi_v term in right hand side of 1D equation
              A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
                A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              // A_3D3D
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
                A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);


            } else {

              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
                A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
                A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);
            }
          }

        } // loop over 3D elements

      } // if neumann node

    } // if number of neighbors = 1
    else {

      for (int i = 0; i < numberOfNeighbors; i++) {

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = util::dist_between_points(coord, coord_neighbor);

        double length_edge = 0.5 * length;

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        double p_neighbor = pointer->neighbors[i]->p_v;

        double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v);

        if (pointer->neighbors[i]->neighbors.size() == 1 && pointer->neighbors[i]->typeOfVGNode == TypeOfNode::DirichletNode && v_interface > 0.0) {
          length_edge = length;
        }

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = 0.0;

        // mass matrix
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += length_edge;

        // advection term
        if (v_interface > -1.0e-8) {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += dt * v_interface;

        } else {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) += dt * v_interface;
        }

        // diffusion term
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += dt * D_v / length;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) -= dt * D_v / length;

        // from previous time step
        b_nut_3D1D[N_tot_3D + indexOfNode] += length_edge * phi_sigma_old[N_tot_3D + indexOfNode];

        // Assemble coupling terms (nutrients)

        int N_s = input.d_num_points_length;

        int N_theta = input.d_num_points_angle;

        std::vector<double> weights;

        std::vector<int> id_3D_elements;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor,
                               radius, h_3D, length_edge, weights,
                               id_3D_elements);

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * length_edge * radius;

        // Permeabilty vessel wall for nutrients
        double L_s = pointer->L_s[i];
        double L_p = pointer->L_p[i];

        // 1D part of the coupling
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
          dt * L_s * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();

        double p_t = 0.0;

        for (int j = 0; j < numberOfElements; j++) {

          if (id_3D_elements[j] > -1) {

            // A_3D1D
            A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
              dt * L_s * surface_area * weights[j];

            // A_3D3D
            A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) +
              dt * L_s * surface_area * weights[j];

            // A_1D3D
            A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) -
              dt * L_s * surface_area * weights[j];

            // osmotic reflection term
            p_t = P_3D1D[id_3D_elements[j]];

            if (p_v - p_t > 0.0) {

              // 1D equation
              // -2pi R (p_v - p_t) phi_v term in right hand side of 1D equation
              A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
                A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              // A_3D3D
              A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) =
                A_nut_3D1D(id_3D_elements[j], N_tot_3D + indexOfNode) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);


            } else {

              // 3D equation
              // 2pi R (p_v - p_t) phi_sigma term in right hand side of 3D equation
              A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
                A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) -
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);

              A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) =
                A_nut_3D1D(N_tot_3D + indexOfNode, id_3D_elements[j]) +
                dt * (1. - osmotic_sigma) * L_p * surface_area * weights[j] * (p_v - p_t);
            }
          }

        } // loop over 3D elements

      } // loop over neighbor segments

    } // if number of neighbors > 1

    pointer = pointer->global_successor;
  }
}
