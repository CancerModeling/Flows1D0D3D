////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "modelUtil.hpp"
#include "network_data_structure.cpp"
#include "network_fully_coupled_assembly.cpp"
#include "network_growth_processes.cpp"
#include "network_semi_coupled_assembly.cpp"

void util::unet::Network::unmark_network_nodes() {
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    pointer->node_marked = false;

    pointer = pointer->global_successor;
  } // loop over vertices
}

/// Simple depth first search in the graph
void depth_first_search(std::shared_ptr<util::unet::VGNode> pointer) {
  pointer->node_marked = true;

  for (auto &neighborPtr : pointer->neighbors) {
    if (!neighborPtr->node_marked) {
      depth_first_search(neighborPtr);
    }
  }
}

void util::unet::Network::mark_nodes_connected_with_initial_nodes() {
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    if (pointer->is_initial_node && !pointer->node_marked) {
      depth_first_search(pointer);
    }

    pointer = pointer->global_successor;
  } // loop over vertices
}

void util::unet::Network::add_lengths_and_volumes_of_unmarked_network(double &total_length, double &total_volume) {
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    for (std::size_t idx = 0; idx < pointer->neighbors.size(); idx += 1) {
      const auto neighbor = pointer->neighbors[idx];
      // we only count the edge if either the node or its neighbor is not marked:
      bool one_is_unmarked = !pointer->node_marked || !neighbor->node_marked;
      // to avoid counting every edge twice, we only count edges,
      // where the index of our node is smaller than the neighbor index.
      if (one_is_unmarked && pointer->index < neighbor->index) {
        const auto length = util::dist_between_points(pointer->coord, neighbor->coord);
        const auto r = pointer->radii[idx];
        total_length += length;
        total_volume += length * r * r * M_PI;
      }
    }

    pointer = pointer->global_successor;
  } // loop over vertices
}

void util::unet::Network::get_length_and_volume_of_network(double &total_length, double &total_volume) {
  total_length = 0;
  total_volume = 0;

  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    for (std::size_t idx = 0; idx < pointer->neighbors.size(); idx += 1) {
      const auto neighbor = pointer->neighbors[idx];
      // to avoid counting every edge twice, we only count edges,
      // where the index of our node is smaller than the neighbor index.
      if (pointer->index < neighbor->index) {
        const auto length = util::dist_between_points(pointer->coord, neighbor->coord);
        const auto r = pointer->radii[idx];
        total_length += length;
        total_volume += length * r * r * M_PI;
      }
    }

    pointer = pointer->global_successor;
  } // loop over vertices
}

int util::unet::Network::get_number_of_bifurcations() {
  int num_bifurcations = 0;

  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {

    if (pointer->neighbors.size() > 2)
      num_bifurcations += 1;

    pointer = pointer->global_successor;
  } // loop over vertices

  return num_bifurcations;
}


void util::unet::Network::delete_unmarked_nodes() {
  // make sure we have at least one node
  if (VGM.isEmpty())
    return;

  // 1. remove the nodes from the neighbor list
  {
    std::shared_ptr<VGNode> pointer = VGM.getHead();

    auto notMarkedFunctional = [](const VGNode &node) -> bool { return !node.node_marked; };

    while (pointer) {
      // remove all the marked nodes from the neighbors list
      // note we have to change several lists at once
      pointer->remove(notMarkedFunctional);

      pointer = pointer->global_successor;
    } // loop over vertices
  }

  // 2. remove the nodes from the global list
  {
    std::shared_ptr<VGNode> pointer = VGM.getHead();

    while (pointer) {

      // we back up the successor, since the VGM.remove method will set it to a nullpointer
      auto successor = pointer->global_successor;

      if (!pointer->node_marked)
        VGM.remove(pointer);

      pointer = successor;
    } // loop over vertices
  }
}

void util::unet::Network::reenumerate_dofs() {
  int next_index = 0;
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    pointer->index = next_index;

    pointer = pointer->global_successor;
    next_index += 1;
  } // loop over vertices
}

void util::unet::Network::delete_unconnected_nodes() {
  unmark_network_nodes();
  mark_nodes_connected_with_initial_nodes();
  add_lengths_and_volumes_of_unmarked_network(total_removed_length, total_removed_volume);
  delete_unmarked_nodes();
  reenumerate_dofs();
}

void util::unet::Network::resize_matrices_direct_solver() {
  // if we use the coupled solver, we do not need to update the matrices of the direct solver
  if (d_coupled_solver)
    return;

  const auto numberOfNodes = VGM.getNumberOfNodes();
  const auto numberOfNodesOld = b.size();

  // we rescale if the number of dof has changed:
  d_model_p->d_log("Rescale the 1D matrices and vectors \n", "net update");
  if (numberOfNodesOld != numberOfNodes) {
    A_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes,
                                                   numberOfNodes);
    b.resize(numberOfNodes);
    P_v.resize(numberOfNodes);

    C_v.resize(numberOfNodes);
    C_v_old.resize(numberOfNodes);

    b_c.resize(numberOfNodes);
  } // update matrix and vector
}

void util::unet::Network::copy_network_to_vectors() {
  auto pointer = VGM.getHead();
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
}

void util::unet::Network::resize_matrices_coupled_solver() {
  if (d_coupled_solver) {
    auto numberOfNodes = VGM.getNumberOfNodes();

    A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
      N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
    b_flow_3D1D.resize(N_tot_3D + numberOfNodes);

    A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
      N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
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
}

void util::unet::Network::set_bc_of_added_vessels_to_neumann() {
  auto pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    if (numberOfNeighbors == 1 && !pointer->is_given) {

      const auto &coord = pointer->coord;

      std::cout << "Set inner node to boundary node" << std::endl;

      pointer->typeOfVGNode = TypeOfNode::NeumannNode;
    }

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::reset_edge_flags() {
  auto pointer = VGM.getHead();

  while (pointer) {

    int numberOfEdges = pointer->neighbors.size();

    for (int i = 0; i < numberOfEdges; i++) {
      pointer->edge_touched[i] = false;
      pointer->sprouting_edge[i] = false;
    }

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::delete_old_sprouters() {
  // all none initial nodes between two initial nodes are removed
  std::shared_ptr<VGNode> pointer = VGM.getHead();
  while (pointer) {
    auto successor = pointer->global_successor;
    if (pointer->neighbors.size() == 2 && pointer->is_sprouting_node) {
      std::cout << "is sprouting " << pointer->index << std::endl;
      // remove node from its neighbors
      for (auto n : pointer->neighbors)
        n->remove([=](auto &p) { return &p == pointer.get(); });
      VGM.remove(pointer);
      // connect neighbors to each other
      {
        pointer->neighbors[0]->neighbors.push_back(pointer->neighbors[1]);
        pointer->neighbors[0]->radii.push_back(0.5 * (pointer->radii[0] + pointer->radii[1]));
        pointer->neighbors[0]->L_p.push_back(0.5 * (pointer->L_p[0] + pointer->L_p[1]));
        pointer->neighbors[0]->L_s.push_back(0.5 * (pointer->L_s[0] + pointer->L_s[1]));
        pointer->neighbors[0]->radii_initial.push_back(0.5 * (pointer->radii_initial[0] + pointer->radii_initial[1]));
        pointer->neighbors[0]->tau_w_initial.push_back(0.5 * (pointer->tau_w_initial[0] + pointer->tau_w_initial[1]));
        pointer->neighbors[0]->edge_touched.push_back(pointer->edge_touched[1] || pointer->edge_touched[0]);
        pointer->neighbors[0]->sprouting_edge.push_back(pointer->sprouting_edge[1] || pointer->sprouting_edge[0]);

        pointer->neighbors[1]->neighbors.push_back(pointer->neighbors[0]);
        pointer->neighbors[1]->radii.push_back(0.5 * (pointer->radii[0] + pointer->radii[1]));
        pointer->neighbors[1]->L_p.push_back(0.5 * (pointer->L_p[0] + pointer->L_p[1]));
        pointer->neighbors[1]->L_s.push_back(0.5 * (pointer->L_s[0] + pointer->L_s[1]));
        pointer->neighbors[1]->radii_initial.push_back(0.5 * (pointer->radii_initial[0] + pointer->radii_initial[1]));
        pointer->neighbors[1]->tau_w_initial.push_back(0.5 * (pointer->tau_w_initial[0] + pointer->tau_w_initial[1]));
        pointer->neighbors[1]->edge_touched.push_back(pointer->edge_touched[1] || pointer->edge_touched[0]);
        pointer->neighbors[1]->sprouting_edge.push_back(pointer->sprouting_edge[1] || pointer->sprouting_edge[0]);
      }
    }
    pointer = successor;
  } // loop over vertices
  reenumerate_dofs();
}

void util::unet::Network::create_initial_network() {

  //
  // Initialize the network
  //
  // Assembly and solve will only happen at processor zero so we create
  // some fields only on zero processor
  //
  // Seperate the variables from what is needed in processor zero and
  // on other processors
  //
  //
  // Only on processor zero:
  //  - VGM data
  //  - Assembly matrices
  //  - Right hand side
  //

  auto &input = d_model_p->get_input_deck();
  d_coupled_solver = d_model_p->d_name == "NetFCFVFE" or input.d_coupled_1d3d;

  if (d_procSize > 1 and d_coupled_solver)
    libmesh_error_msg(
      "Fully 1D-3D coupled solver only works in serial execution");

  scenario = input.d_scenario;
  oss << "Scenario: " << scenario << std::endl;

  // read file and create initial network
  // Do this only on first processor
  if (d_procRank == 0) {
    std::vector<std::vector<double>> vertices;
    std::vector<double> pressures;
    std::vector<double> radii;
    std::vector<std::vector<unsigned int>> elements;

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
  }

  // create data of vertices and segments in processor zero
  // and communicate with other processors
  prepare_and_communicate_network();

  // get minimum length of vessel segments and set it as minimum length
  // for sprouting
  if (d_procRank == 0) {
    double min_vessel_length = 100. * input.d_domain_params[1];
    std::shared_ptr<VGNode> pointer = VGM.getHead();
    while (pointer) {
      for (int i = 0; i < pointer->neighbors.size(); i++) {
        double length = util::dist_between_points(pointer->coord, pointer->neighbors[i]->coord);
        if (min_vessel_length > length)
          min_vessel_length = length;
      }
      pointer = pointer->global_successor;
    } // loop over vertices

    // if user did not specify min length, set now
    if (input.d_min_length_for_sprouting < 1.e-12)
      input.d_min_length_for_sprouting = 0.9 * min_vessel_length;
  }

  // Initialize common data
  d_has_network_changed = true;
  d_update_interval = input.d_network_update_interval;
  N_3D = input.d_num_elems;
  N_tot_3D = N_3D * N_3D * N_3D;
  L_x = input.d_domain_params[1];
  h_3D = L_x / (double) N_3D;
  mu = input.d_init_vessel_mu;
  D_v = input.d_D_sigma_v;
  D_v_3D = input.d_D_sigma;
  D_TAF = input.d_D_TAF;
  osmotic_sigma = input.d_osmotic_sigma;
  K_3D = input.d_tissue_flow_coeff;
  L_p = input.d_tissue_flow_L_p;
  L_s = input.d_tissue_nut_L_s;

  // allocate space for solution of 3D species
  if (d_procRank == 0) {
    P_3D = std::vector<double>(N_tot_3D, 0.0);
    phi_sigma_3D = std::vector<double>(N_tot_3D, 0.0);
    phi_TAF = std::vector<double>(N_tot_3D, 0.0);
  }

  // initialize matrix and vector
  if (!d_coupled_solver) {

    // exclusive to processor zero
    if (d_procRank == 0) {

      // 1D pressure: matrix, rhs, and solution
      A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(d_numVertices, d_numVertices);
      b = std::vector<double>(d_numVertices, 0.);

      // 1D nutrient: matrix, rhs, and solution
      b_c = std::vector<double>(d_numVertices, 0.0);
      C_v_old = std::vector<double>(d_numVertices, 0.0);
    }

    // common to all processors
    P_v = std::vector<double>(d_numVertices, 0.);
    C_v = std::vector<double>(d_numVertices, 0.0);

  } else {

    // 3D1D flow problem: matrix, rhs and solution
    A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
      N_tot_3D + d_numVertices, N_tot_3D + d_numVertices);
    b_flow_3D1D = std::vector<double>(N_tot_3D + d_numVertices, 0.0);

    P_3D1D = std::vector<double>(N_tot_3D + d_numVertices, 0.0);

    A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D + d_numVertices,
                                                        N_tot_3D + d_numVertices);
    b_nut_3D1D = std::vector<double>(N_tot_3D + d_numVertices, 0.0);

    phi_sigma = std::vector<double>(N_tot_3D + d_numVertices, 0.0);
    phi_sigma_old = std::vector<double>(N_tot_3D + d_numVertices, 0.0);

    for (int i = 0; i < N_tot_3D; i++) {

      phi_sigma_3D[i] = input.d_nut_ic_value;
      phi_sigma_old[i] = input.d_nut_ic_value;
      phi_sigma[i] = input.d_nut_ic_value;
    }

    for (int i = 0; i < d_numVertices; i++) {
      phi_sigma_old[N_tot_3D + i] = 0.0;
      phi_sigma[N_tot_3D + i] = 0.0;
    }
  }

  // initialize nutrient as one in artery
  // Only on processor zero
  //if (d_procRank == 0) {
  if (scenario == "two_vessels" and d_procRank == 0) {

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
    } // loop over vertices
  }   // if processor zero
}

void util::unet::Network::solve3D1DFlowProblem(BaseAssembly &pres_sys,
                                               BaseAssembly &tum_sys) {

  if (pres_sys.d_sys_name != "Pressure" or tum_sys.d_sys_name != "Tumor")
    libmesh_error_msg("Must pass Pressure and Tumor system to solve 3D-1D "
                      "pressure");

  const auto &input = d_model_p->get_input_deck();

  std::cout << " " << std::endl;
  std::cout << "Assemble system of equations (pressure)" << std::endl;

  assemble3D1DSystemForPressure(pres_sys, tum_sys);

  // Solve linear system of equations
  std::cout << " " << std::endl;
  std::cout << "Solve linear system of equations (pressure)" << std::endl;

  // gmm::iteration iter(10E-11, 2);

  size_t restart = 500;

  if (d_update_number == 1) {

    P_3D1D = b_flow_3D1D;
  }

  gmm::iteration iter(1.0E-10);

  gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_flow_3D1D);

  gmm::gmres(A_flow_3D1D, P_3D1D, b_flow_3D1D, PR, restart, iter);

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

  // also update the boundary flag
  update_and_communicate_bdry_flag();
}

void util::unet::Network::solve3D1DNutrientProblem(BaseAssembly &nut_sys,
                                                   BaseAssembly &tum_sys) {
  std::cout << "called" << std::endl;

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

  gmm::iteration iter(1.0E-12);

  gmm::ilu_precond<gmm::row_matrix<gmm::wsvector<double>>> PR(A_nut_3D1D);

  gmm::gmres(A_nut_3D1D, phi_sigma, b_nut_3D1D, PR, restart, iter);
  //gmm::lu_solve(A_nut_3D1D, phi_sigma, b_nut_3D1D);

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->c_v = phi_sigma[N_tot_3D + indexOfNode];

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

  std::cout << "pressure step 1" << std::endl;
  // solve only on processor zero and then communicate to all other processors
  if (d_procRank == 0) {

    assembleVGMSystemForPressure(pres_sys);

    std::cout << "pressure step assembly" << std::endl;

    gmm::iteration iter(10E-18);

    gmm::identity_matrix P;

    P_v = b;

    gmm::bicgstab(A_VGM, P_v, b, P, iter);

    std::cout << "pressure step solve" << std::endl;

    std::shared_ptr<VGNode> pointer = VGM.getHead();

    while (pointer) {

      int indexOfNode = pointer->index;

      pointer->p_v = P_v[indexOfNode];

      pointer = pointer->global_successor;
    }

  } // solve on processor zero

  // communicate solution to all processors
  if (d_procSize > 1) {

    // resize P_v in non zero processors before communication
    if (d_procRank > 0)
      P_v.resize(0);

    d_comm_p->allgather(P_v);
    if (P_v.size() != d_numVertices)
      libmesh_error_msg("Size of P_v " + std::to_string(P_v.size()) +
                        " after allgather does not match expected size " +
                        std::to_string(d_numVertices));
  }

  // also update the boundary flag
  update_and_communicate_bdry_flag();
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

  // solve only on processor zero and then communicate to all other processors
  if (d_procRank == 0) {

    assembleVGMSystemForNutrient(pres_sys, nut_sys);

    // if this is first call inside nonlinear loop, we guess current
    // concentration as old concentration
    if (d_model_p->d_nonlinear_step == 0)
      C_v = C_v_old;
    if (d_model_p->d_step == 1)
      C_v = b_c;

    // get preconditioner
    gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_VGM, 50,
                                                                1e-6);

    // solve
    gmm::iteration iter(5.0e-11);
    gmm::bicgstab(A_VGM, C_v, b_c, P, iter);

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

  } // solve on processor zero

  // communicate solution to all processors
  if (d_procSize > 1) {

    // resize C_v in non zero processors before communication
    if (d_procRank > 0)
      C_v.resize(0);

    d_comm_p->allgather(C_v);
    if (C_v.size() != d_numVertices)
      libmesh_error_msg("Size of C_v " + std::to_string(C_v.size()) +
                        " after allgather does not match expected size " +
                        std::to_string(d_numVertices));
  }
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

std::vector<double> util::unet::Network::compute_qoi() {

  if (d_procRank > 0)
    return {};

  // For Tobias analysis
  // TODO remove in final merge
  {
    // report 1D data
    auto pointer = VGM.getHead();

    int numberOfBifurcations = 0;

    int numberOfVessels = 0;

    double total_length = 0.0;

    double total_volume = 0.0;

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      if (numberOfNeighbors > 2) {

        numberOfBifurcations++;
      }

      for (int i = 0; i < numberOfNeighbors; i++) {

        if (pointer->edge_touched[i] == false) {

          double length = util::dist_between_points(pointer->coord, pointer->neighbors[i]->coord);

          total_length = total_length + length;

          total_volume = total_volume + length * pointer->radii[i] * pointer->radii[i] * M_PI;

          numberOfVessels++;

          int localIndex = pointer->neighbors[i]->getLocalIndexOfNeighbor(pointer);

          pointer->edge_touched[i] = true;

          pointer->neighbors[i]->edge_touched[localIndex] = true;
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

    std::string path_number_bifurcations = "two_vessels_number_bifurcations.txt";

    std::fstream file_number_bifurcations;
    file_number_bifurcations.open(path_number_bifurcations, std::ios::out | std::ios::app);
    file_number_bifurcations << numberOfBifurcations << std::endl;

    std::string path_total_length = "two_vessels_total_length.txt";

    std::fstream file_total_length;
    file_total_length.open(path_total_length, std::ios::out | std::ios::app);
    file_total_length << total_length << std::endl;

    std::string path_number_of_vessels = "two_vessels_number_of_vessels.txt";

    std::fstream file_number_of_vessels;
    file_number_of_vessels.open(path_number_of_vessels, std::ios::out | std::ios::app);
    file_number_of_vessels << numberOfVessels << std::endl;

    std::string path_total_volume = "two_vessels_total_volume.txt";

    std::fstream file_total_volume;
    file_total_volume.open(path_total_volume, std::ios::out | std::ios::app);
    file_total_volume << total_volume << std::endl;
  }

  // 0 - r_v_mean             1 - r_v_std
  // 2 - l_v_mean             3 - l_v_std         4 - l_v_total
  // 5 - vessel_vol           6 - vessel_density
  auto qoi = std::vector<double>(7, 0.);

  const auto &input = d_model_p->get_input_deck();

  // compute the mean and also store the array for calculation of std dev
  std::vector<double> r_v;
  std::vector<double> l_v;
  double domain_vol = std::pow(input.d_domain_params[1], 3);

  for (unsigned int i = 0; i < d_numSegments; i++) {

    auto node1 = d_segments[2 * i + 0];
    auto node2 = d_segments[2 * i + 1];
    auto radius = d_segmentData[i];

    std::vector<double> coord1;
    std::vector<double> coord2;
    for (unsigned int j = 0; j < 3; j++) {
      coord1.push_back(d_vertices[3 * node1 + j]);
      coord2.push_back(d_vertices[3 * node2 + j]);
    }

    double length = util::dist_between_points(coord1, coord2);

    // add to data
    r_v.push_back(radius);
    l_v.push_back(length);
    qoi[0] += radius;
    qoi[2] += length;
    qoi[4] += length;
    qoi[5] += M_PI * radius * radius * length;
  }

  // compute mean and std dev
  qoi[0] = qoi[0] / d_numSegments;
  qoi[2] = qoi[2] / d_numSegments;
  qoi[1] = util::get_std_dev(r_v, qoi[0]);
  qoi[3] = util::get_std_dev(r_v, qoi[2]);
  qoi[6] = qoi[5] / domain_vol;

  return qoi;
}

unsigned int util::unet::Network::get_assembly_cases_pres(
  const std::shared_ptr<VGNode> &pointer,
  const double &identify_vein_pres) const {

  // find various cases
  if (pointer->neighbors.size() == 1) {
    if (pointer->typeOfVGNode == TypeOfNode::DirichletNode)
      // return "boundary_dirichlet";
      return UNET_PRES_BDRY_DIRIC;
    else
      // return "boundary_inner";
      return UNET_PRES_BDRY_INNER;
  } // neighbor == 1
  else {
    // return "inner";
    return UNET_PRES_INNER;
  }
}

unsigned int util::unet::Network::get_assembly_cases_nut(
  const std::shared_ptr<VGNode> &pointer,
  const double &identify_vein_pres) const {

  // find various cases
  if (pointer->neighbors.size() == 1) {
    if (pointer->p_v > pointer->neighbors[0]->p_v) {
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {
        if (pointer->p_v >= identify_vein_pres)
          // return "boundary_artery_inlet";
          return UNET_NUT_BDRY_ARTERY_INLET;
        else
          // return "boundary_vein_inlet";
          return UNET_NUT_BDRY_VEIN_INLET;
      } // inlet dirichlet
      else {
        // return "boundary_inner_inlet";
        return UNET_NUT_BDRY_INNER_INLET;
      } // not dirichlet
    }   // v > 0
    else {
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        if (pointer->p_v >= identify_vein_pres)
          // return "boundary_artery_inlet";
          return UNET_NUT_BDRY_ARTERY_OUTLET;
        else
          // return "boundary_vein_inlet";
          return UNET_NUT_BDRY_VEIN_OUTLET;
      } else
        return UNET_NUT_BDRY_INNER_OUTLET;
    } // v < 0
  }   // neighbor == 1
  else {
    // return "inner";
    return UNET_NUT_INNER;
  }
}

std::string util::unet::Network::get_assembly_cases_pres_str(
  const std::shared_ptr<VGNode> &pointer,
  const double &identify_vein_pres) const {

  // find various cases
  if (pointer->neighbors.size() == 1) {
    if (pointer->typeOfVGNode == TypeOfNode::DirichletNode)
      return "boundary_dirichlet";
    // return UNET_PRES_BDRY_DIRIC;
    else
      return "boundary_inner";
    // return UNET_PRES_BDRY_INNER;
  } // neighbor == 1
  else {
    return "inner";
    // return UNET_PRES_INNER;
  }
}

std::string util::unet::Network::get_assembly_cases_nut_str(
  const std::shared_ptr<VGNode> &pointer,
  const double &identify_vein_pres) const {

  // find various cases
  if (pointer->neighbors.size() == 1) {
    if (pointer->p_v > pointer->neighbors[0]->p_v) {
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {
        if (pointer->p_v >= identify_vein_pres)
          return "boundary_artery_inlet";
        //return UNET_NUT_BDRY_ARTERY_INLET;
        else
          return "boundary_vein_inlet";
        //return UNET_NUT_BDRY_VEIN_INLET;
      } // inlet dirichlet
      else {
        return "boundary_inner_inlet";
        //return UNET_NUT_BDRY_INNER_INLET;
      } // not dirichlet
    }   // v > 0
    else {
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {

        if (pointer->p_v >= identify_vein_pres)
          return "boundary_artery_outlet";
        //return UNET_NUT_BDRY_ARTERY_OUTLET;
        else
          return "boundary_vein_outlet";
        //return UNET_NUT_BDRY_VEIN_OUTLET;
      } else
        return "boundary_inner_outlet";
      //return UNET_NUT_BDRY_INNER_OUTLET;
    } // v < 0
  }   // neighbor == 1
  else {
    return "inner";
    //return UNET_NUT_INNER;
  }
}

void util::unet::Network::prepare_and_communicate_network() {

  //  if (d_coupled_solver)
  //    return;

  // On processor zero, prepare the network data
  if (d_procRank == 0) {

    const auto &input = d_model_p->get_input_deck();

    //// vertex data
    d_numSegments = 0;
    d_numVertices = VGM.getNumberOfNodes();
    d_vertices.resize(d_numVertices * 3);
    d_vertexBdFlag.resize(d_numVertices);
    d_segments.resize(0);
    d_segmentData.resize(0);

    // loop over nodes
    auto pointer = VGM.getHead();
    while (pointer) {

      unsigned int i_start = 3 * pointer->index;
      for (int i = 0; i < 3; i++)
        d_vertices[i_start + i] = pointer->coord[i];

      // get boundary flag
      unsigned int bdry_flag = UNET_FREE_MASK;

      // for pressure and nutrient
      bdry_flag |= get_assembly_cases_pres(pointer, input.d_identify_vein_pres);
      bdry_flag |= get_assembly_cases_nut(pointer, input.d_identify_vein_pres);

      d_vertexBdFlag[pointer->index] = bdry_flag;

      for (int i = 0; i < pointer->neighbors.size(); i++) {
        if (!pointer->edge_touched[i]) {
          // add segment
          d_segments.push_back(pointer->index);
          d_segments.push_back(pointer->neighbors[i]->index);

          // add radius
          d_segmentData.push_back(pointer->radii[i]);

          // check for length
          double length = util::dist_between_points(pointer->coord,
                                                    pointer->neighbors[i]->coord);
          if (util::definitelyGreaterThan(length, 0.5))
            libmesh_error_msg("Error unusually long vessel detected. Length = " +
                              std::to_string(length) +
                              " is greater than tolerance 0.5");

          if (util::definitelyLessThan(length, 1.e-8))
            libmesh_error_msg("Error unusually short vessel detected. Length = " +
                              std::to_string(length));

          d_numSegments += 1;
          pointer->edge_touched[i] = true;
          pointer->neighbors[i]->markEdge(pointer->index);
        }
      }

      pointer = pointer->global_successor;
    }

    // reset the marked edges
    pointer = VGM.getHead();
    while (pointer) {

      for (int i = 0; i < pointer->neighbors.size(); i++)
        pointer->edge_touched[i] = false;

      pointer = pointer->global_successor;
    }
  } else {
    d_numVertices = 0;
    d_numSegments = 0;
    d_vertices.resize(0);
    d_vertexBdFlag.resize(0);
    d_segments.resize(0);
    d_segmentData.resize(0);
  }

  // communicate only if more than one processor
  if (d_procSize > 1) {
    // get the number of vertices and segments
    d_comm_p->max(d_numVertices);
    d_comm_p->max(d_numSegments);

    // get the remaining datas
    d_comm_p->allgather(d_vertices);
    if (d_vertices.size() != 3 * d_numVertices)
      libmesh_error_msg(
        "Error size of vertices = " + std::to_string(d_vertices.size()) +
        " does not match expected size = " +
        std::to_string(3 * d_numVertices));

    d_comm_p->allgather(d_vertexBdFlag);
    if (d_vertexBdFlag.size() != d_numVertices)
      libmesh_error_msg(
        "Error size of vertexBdFlag = " +
        std::to_string(d_vertexBdFlag.size()) +
        " does not match expected size = " + std::to_string(d_numVertices));

    d_comm_p->allgather(d_segments);
    if (d_segments.size() != 2 * d_numSegments)
      libmesh_error_msg("Error size of segments connectivity = " +
                        std::to_string(d_segments.size()) +
                        " does not match expected size = " +
                        std::to_string(2 * d_numSegments));

    d_comm_p->allgather(d_segmentData);
    if (d_segmentData.size() != d_numSegments)
      libmesh_error_msg(
        "Error size of segmentData = " +
        std::to_string(d_segmentData.size()) +
        " does not match expected size = " + std::to_string(d_numSegments));
  }
}

void util::unet::Network::update_and_communicate_bdry_flag() {

  //  if (d_coupled_solver)
  //    return;

  // On processor zero, prepare the network data
  if (d_procRank == 0) {

    const auto &input = d_model_p->get_input_deck();

    //// vertex data
    // d_vertexBdFlag.resize(d_numVertices);

    // loop over nodes
    auto pointer = VGM.getHead();
    while (pointer) {

      // get boundary flag
      unsigned int bdry_flag = UNET_FREE_MASK;

      // for pressure and nutrient
      bdry_flag |= get_assembly_cases_pres(pointer, input.d_identify_vein_pres);
      bdry_flag |= get_assembly_cases_nut(pointer, input.d_identify_vein_pres);

      d_vertexBdFlag[pointer->index] = bdry_flag;
      pointer = pointer->global_successor;
    }
  } else {
    d_numVertices = 0;
    d_vertexBdFlag.resize(0);
  }

  // communicate only if more than one processor
  if (d_procSize > 1) {
    // get the number of vertices and segments
    d_comm_p->max(d_numVertices);

    d_comm_p->allgather(d_vertexBdFlag);
    if (d_vertexBdFlag.size() != d_numVertices)
      libmesh_error_msg(
        "Error size of vertexBdFlag = " +
        std::to_string(d_vertexBdFlag.size()) +
        " does not match expected size = " + std::to_string(d_numVertices));
  }
}

void util::unet::Network::check_vessel_length() {

  // On processor zero, prepare the network data
  if (d_procRank == 0) {

    // loop over nodes
    auto pointer = VGM.getHead();
    while (pointer) {

      for (int i = 0; i < pointer->neighbors.size(); i++) {
        // get length of vessel
        auto length = util::dist_between_points(pointer->coord,
                                                pointer->neighbors[i]->coord);
        if (util::definitelyGreaterThan(length, 0.5))
          libmesh_error_msg("Error unusually long vessel detected. Length = " +
                            std::to_string(length) +
                            " is greater than tolerance 0.5");

        if (util::definitelyLessThan(length, 1.e-8))
          libmesh_error_msg("Error unusually short vessel detected. Length = " +
                            std::to_string(length));
      }

      pointer = pointer->global_successor;
    }
  }
}
