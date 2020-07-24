////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "modelUtil.hpp"
#include "network_data_structure.cpp"
#include "network_growth_processes.cpp"
#include "network_fully_coupled_assembly.cpp"
#include "network_semi_coupled_assembly.cpp"

void util::unet::Network::create_initial_network() {

  const auto &input = d_model_p->get_input_deck();
  auto comm = d_model_p->get_comm();
  d_is_network_changed = true;
  d_coupled_solver =
      d_model_p->d_name == "NetFCFVFE" or d_model_p->d_name == "NetFV";

  // equation system
  std::vector<std::vector<double>> vertices;
  std::vector<double> pressures;
  std::vector<double> radii;
  std::vector<std::vector<unsigned int>> elements;

  scenario = input.d_scenario;

  // std::cout << " " << std::endl;
  oss << "Scenario: " << scenario << std::endl;

  readData(vertices, pressures, radii, elements);

  transferDataToVGM(vertices, pressures, radii, elements);

  int numberOfNodes = VGM.getNumberOfNodes();

  int refinementLevel = input.d_network_init_refinement;

  for (int i = 0; i < refinementLevel; i++) {

    refine1DMesh();
  }

  numberOfNodes = VGM.getNumberOfNodes();

  // std::cout << " " << std::endl;
  oss << "Number of nodes in network: " << numberOfNodes << std::endl;
  d_model_p->d_log(oss, "debug");

  // get some fixed parameters
  N_3D = input.d_num_elems;
  N_tot_3D = N_3D * N_3D * N_3D;
  L_x = input.d_domain_params[1];
  h_3D = L_x / (double)N_3D;
  mu = input.d_init_vessel_mu;
  D_v = input.d_D_sigma_v;
  D_v_3D = input.d_D_sigma;
  D_TAF = input.d_D_TAF;
  osmotic_sigma = input.d_osmotic_sigma;
  K_3D = input.d_tissue_flow_coeff;

  // allocate space for solution of 3D species
  P_3D = std::vector<double>(N_tot_3D, 0.0);
  phi_sigma_3D = std::vector<double>(N_tot_3D, 0.0);
  phi_TAF_3D = std::vector<double>(N_tot_3D, 0.0);

  // initialize matrix and vector
  if (!d_coupled_solver) {
    // 1D pressure: matrix, rhs, and solution
    A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
    b = std::vector<double>(numberOfNodes, 0.);
    P_v = std::vector<double>(numberOfNodes, 0.);

    // 1D nutrient: matrix, rhs, and solution
    Ac_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
    b_c = std::vector<double>(numberOfNodes, 0.0);

    C_v = std::vector<double>(numberOfNodes, 0.0);
    C_v_old = std::vector<double>(numberOfNodes, 0.0);
  } else {
    // 3D1D flow problem: matrix, rhs and solution
    A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
        N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
    b_flow_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

    P_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

    A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
        N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
    b_nut_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

    phi_sigma = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);
    phi_sigma_old = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

    for (int i = 0; i < N_tot_3D; i++) {

      phi_sigma_3D[i] = input.d_nut_ic_value;
      phi_sigma_old[i] = input.d_nut_ic_value;
      phi_sigma[i] = input.d_nut_ic_value;
    }
  }

  // initialize nutrient as one in artery
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    if (pointer->radii[0] < input.d_identify_artery_radius) {

      if (d_coupled_solver) {
        phi_sigma_old[N_tot_3D + indexOfNode] = 1.0;
        phi_sigma[N_tot_3D + indexOfNode] = 1.0;
      } else {

        C_v[indexOfNode] = 1.;
        C_v_old[indexOfNode] = 1.;
      }
    } else {

      if (d_coupled_solver) {
        phi_sigma_old[N_tot_3D + indexOfNode] = 0.;
        phi_sigma[N_tot_3D + indexOfNode] = 0.;
      } else {

        C_v[indexOfNode] = 0.;
        C_v_old[indexOfNode] = 0.;
      }
    }

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::solve3D1DFlowProblem(BaseAssembly &pres_sys,
                                               BaseAssembly &tum_sys) {

  if (pres_sys.d_sys_name != "Pressure" or tum_sys.d_sys_name != "Tumor")
    libmesh_error_msg("Must pass Pressure and Tumor system to solve 3D-1D "
                      "pressure");

  const auto &input = d_model_p->get_input_deck();

  assemble3D1DSystemForPressure(pres_sys, tum_sys);

  // Solve linear system of equations
  // std::cout << " " << std::endl;
  // std::cout << "Solve linear system of equations (pressure)" << std::endl;

  // gmm::iteration iter(10E-11, 2);

  size_t restart = 500;

  if (d_update_number == 1) {

    P_3D1D = b_flow_3D1D;
  }

  gmm::iteration iter(1.0E-10, 2);

  // gmm::identity_matrix PR;

  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_flow_3D1D, 50,
                                                               1e-8);

  // gmm::ilutp_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_flow_3D1D,
  // 50, 1e-4);

  // gmm::ildlt_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_flow_3D1D);

  // gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_flow_3D1D);

  gmm::gmres(A_flow_3D1D, P_3D1D, b_flow_3D1D, PR, restart, iter);

  // gmm::bicgstab(A_flow_3D1D, P_3D1D, b_flow_3D1D, PR, iter);

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->p_v = P_3D1D[N_tot_3D + indexOfNode];
    pointer = pointer->global_successor;
  }

  for (int i = 0; i < N_tot_3D; i++) {

    P_3D[i] = P_3D1D[i];
  }

  // copy the 3D pressure to libmesh pressure system
  pres_sys.set_elem_sol(P_3D);
}

void util::unet::Network::solve3D1DNutrientProblem(BaseAssembly &nut_sys,
                                                   BaseAssembly &tum_sys) {

  if (nut_sys.d_sys_name != "Nutrient" or tum_sys.d_sys_name != "Tumor")
    libmesh_error_msg("Must pass Nutrient and Tumor system to solve 3D-1D "
                      "nutrient");

  const auto &input = d_model_p->get_input_deck();
  const auto timeStep = d_model_p->d_step;

  assemble3D1DSystemForNutrients(nut_sys, tum_sys);

  // if this is first call inside nonlinear loop, we guess current
  // concentration as old concentration

  if (d_model_p->d_nonlinear_step == 0)
    phi_sigma = phi_sigma_old;
  if (d_model_p->d_step == 1)
    phi_sigma = b_nut_3D1D;

  size_t restart = 150;

  gmm::iteration iter(1.0E-10);

  // gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_nut_3D1D,
  // 50, 1e-4);

  // gmm::identity_matrix PR;

  // gmm::ilutp_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_nut_3D1D,
  // 70, 1e-8);

  gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_nut_3D1D);

  gmm::gmres(A_nut_3D1D, phi_sigma, b_nut_3D1D, PR, restart, iter);

  // gmm::bicgstab(A_nut_3D1D, phi_sigma, b_nut_3D1D, PR, iter);

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->c_v = phi_sigma[N_tot_3D + indexOfNode];

    if (phi_sigma[N_tot_3D + indexOfNode] > 1.0) {

      pointer->c_v = 1.0;
      phi_sigma[N_tot_3D + indexOfNode] = 1.0;
    }
    /*
        std::cout << "index: " << pointer->index << " c_v: " << pointer->c_v
                  << " p_v: " << pointer->p_v << " coord: " << pointer->coord
                  << std::endl;
    */
    pointer = pointer->global_successor;
  }

  // do not modify old with current concentration as this solver could be
  // called inside nonlinear loop at given time step
  // Rather update the old with new where this solver is called at the end of
  // nonlinear loop phi_sigma_old = phi_sigma;

  // Extract nutrient concentrations
  // std::cout << " " << std::endl;
  // std::cout << "Extract nutrient concentrations" << std::endl;

  for (int i = 0; i < N_tot_3D; i++) {

    phi_sigma_3D[i] = phi_sigma[i];
  }

  // copy the 3D nutrient to libmesh nutrient system
  nut_sys.set_elem_sol(phi_sigma_3D);
}

void util::unet::Network::solveVGMforPressure(BaseAssembly &pres_sys) {

  if (pres_sys.d_sys_name != "Pressure")
    libmesh_error_msg("Must pass Pressure system to solve 1D pressure");

  // gather pressure solution in all processors
  pres_sys.localize_solution_with_elem_id_numbering_const_elem(P_3D, {0},
                                                               false);

  assembleVGMSystemForPressure(pres_sys);

  gmm::iteration iter(10E-18);

  //  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_VGM, 50,
  //  1e-5);

  gmm::identity_matrix P;

  P_v = b;

  gmm::bicgstab(A_VGM, P_v, b, P, iter);

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->p_v = P_v[indexOfNode];

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::solveVGMforNutrient(BaseAssembly &pres_sys,
                                              BaseAssembly &nut_sys) {

  if (pres_sys.d_sys_name != "Pressure" or nut_sys.d_sys_name != "Nutrient")
    libmesh_error_msg("Must pass Pressure and Nutrient system to solve 1D "
                      "nutrient");

  // we do not update 3D pressure assuming that pressure system is solved
  // before solving other systems and therefore 3D pressure in P_3D is
  // already updated

  // gather nutrient solution in all processors
  nut_sys.localize_solution_with_elem_id_numbering_const_elem(phi_sigma_3D, {0},
                                                              false);

  assembleVGMSystemForNutrient(pres_sys, nut_sys);

  // if this is first call inside nonlinear loop, we guess current
  // concentration as old concentration
  if (d_model_p->d_nonlinear_step == 0)
    C_v = C_v_old;
  if (d_model_p->d_step == 1)
    C_v = b_c;

  // get preconditioner
  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(Ac_VGM, 50, 1e-6);

  // solve
  gmm::iteration iter(5.0e-11);
  gmm::bicgstab(Ac_VGM, C_v, b_c, P, iter);

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->c_v = C_v[indexOfNode];

    pointer = pointer->global_successor;
  }

  // do not modify old with current concentration as this solver could be
  // called inside nonlinear loop at given time step
  // Rather update the old with new where this solver is called at the end of
  // nonlinear loop
  // C_v_old = C_v;
}

double util::unet::Network::getDirichletValue(std::vector<double> center_face,
                                              double L_p, double radius) {

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

double util::unet::Network::getK1D(double s, double L_p, double radius) {

  if (scenario == "test_single_vessel") {

    return L_p / (1.0 + L_p) * (1.0 + s + 0.5 * s * s);

  } else {

    return (radius * radius * radius * radius * M_PI) / (8.0 * mu);
  }
}
