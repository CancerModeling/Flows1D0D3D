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

void util::unet::Network::create_initial_network() {

  const auto &input = d_model_p->get_input_deck();
  d_is_network_changed = true;

  // equation system
  std::vector<std::vector<double>> vertices;
  std::vector<double> pressures;
  std::vector<double> radii;
  std::vector<std::vector<unsigned int>> elements;

  scenario = input.d_scenario;

  //std::cout << " " << std::endl;
  oss << "Scenario: " << scenario << std::endl;

  readData(vertices, pressures, radii, elements);

  transferDataToVGM(vertices, pressures, radii, elements);

  int numberOfNodes = VGM.getNumberOfNodes();

  //std::cout << " " << std::endl;
  //oss << "Number of nodes: " << numberOfNodes << std::endl;

  // refine mesh
  //std::cout << " " << std::endl;
  //oss << "Refine mesh" << std::endl;

  int refinementLevel = input.d_network_init_refinement;

  for(int i = 0; i < refinementLevel; i++) {

      refine1DMesh();

  }

  numberOfNodes = VGM.getNumberOfNodes();

  //std::cout << " " << std::endl;
  oss << "Number of nodes: " << numberOfNodes << std::endl;
  d_model_p->d_delayed_msg += "Network Info \n";
  d_model_p->d_delayed_msg += oss.str() + " \n";
  oss.str("");
  oss.clear();

  // compute element and weights
  if (input.d_compute_elem_weights and input.d_model_name != "NetFCFVFE")
    compute_elem_weights();

  // Print data
  //std::cout << " " << std::endl;
  //std::cout << "Print data" << std::endl;
  // printDataVGM();

  // initialize matrix and vector
  // 1D pressure: matrix, rhs, and solution
  A_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  b = std::vector<double>(numberOfNodes, 0.);
  P_v = std::vector<double>(numberOfNodes, 0.);

  // 1D nutrient: matrix, rhs, and solution
  Ac_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  b_c = std::vector<double>(numberOfNodes, 0.0);

  C_v = std::vector<double>(numberOfNodes, 0.0);
  C_v_old = std::vector<double>(numberOfNodes, 0.0);

  // 3D1D flow problem: matrix, rhs and solution
  N_3D = input.d_num_elems;
  N_tot_3D = N_3D * N_3D * N_3D;

  double L_x = input.d_domain_params[1];
  h_3D = L_x / (double)N_3D;

  A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(
      N_tot_3D + numberOfNodes, N_tot_3D + numberOfNodes);
  b_flow_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  P_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);
  P_3D = std::vector<double>(N_tot_3D, 0.0);

  A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D + numberOfNodes,
                                                      N_tot_3D + numberOfNodes);
  b_nut_3D1D = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);

  phi_sigma = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);
  phi_sigma_old = std::vector<double>(N_tot_3D + numberOfNodes, 0.0);
  phi_sigma_3D = std::vector<double>(N_tot_3D, 0.0);

  for(int i = 0; i < N_tot_3D; i++) {

      phi_sigma_3D[i]  = input.d_nut_ic_value;
      phi_sigma_old[i] = input.d_nut_ic_value;
      phi_sigma[i] = input.d_nut_ic_value;

  }

  for(int i = 0; i < numberOfNodes; i++) {

      phi_sigma_old[N_tot_3D+i] = 0.0;
      phi_sigma[N_tot_3D+i] = 0.0;

  }

  phi_TAF = std::vector<double>(N_tot_3D, 0.0);
  phi_TAF_old = std::vector<double>(N_tot_3D, 0.0);

  mu = input.d_init_vessel_mu;

  D_v = input.d_D_sigma_v;

  D_v_3D = input.d_D_sigma;

  D_TAF = input.d_D_TAF;

  osmotic_sigma = input.d_osmotic_sigma;
}

void util::unet::Network::solve3D1DFlowProblem(BaseAssembly &pres_sys, BaseAssembly &tum_sys) {

  if (pres_sys.d_sys_name != "Pressure" or tum_sys.d_sys_name != "Tumor")
    libmesh_error_msg("Must pass Pressure and Tumor system to solve 3D-1D "
                      "pressure");

  const auto &input = d_model_p->get_input_deck();
  double L_x = input.d_domain_params[1];
  const auto timeStep = d_model_p->d_step;

  assemble3D1DSystemForPressure(pres_sys, tum_sys);

  // Solve linear system of equations
  //std::cout << " " << std::endl;
  //std::cout << "Solve linear system of equations (pressure)" << std::endl;

  //gmm::iteration iter(10E-11, 2);
  gmm::iteration iter(10E-11);

  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> Pr(A_flow_3D1D, 50,
                                                               1e-5);

  P_3D1D = b_flow_3D1D;

  gmm::bicgstab(A_flow_3D1D, P_3D1D, b_flow_3D1D, Pr, iter);

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->p_v = P_3D1D[N_tot_3D + indexOfNode];
    pointer = pointer->global_successor;
  }

  // Extract pressures
  //std::cout << " " << std::endl;
  //std::cout << "Extract pressures" << std::endl;

  for (int i = 0; i < N_tot_3D; i++) {

    P_3D[i] = P_3D1D[i];
  }

  // Compute velocity field
  /*
  std::cout << "Compute velocity field" << std::endl;
  std::vector<std::vector<double>> V_3D;

  // Space directions
  std::vector<std::vector<double>> directions;
  directions = defineDirections();

  // Goes K divided by mu in Darcy's law
  double K_3D = input.d_tissue_flow_coeff;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        std::vector<double> vel_element(3);

        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        // Iterate over the interfaces
        for (int face = 0; face < 6; face++) {

          std::vector<double> center_neighbor =
              getCenterNeighbor(center, directions[face], h_3D);

          bool isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            int index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            for (int l = 0; l < 3; l++) {

              vel_element[l] += -K_3D * (P_3D[index_neighbor] - P_3D[index]) /
                                h_3D * directions[face][l];
            }
          }
        }

        V_3D.push_back(vel_element);
      }
    }
  }

  if (timeStep % 2 == 0) {

    std::cout << " " << std::endl;
    std::cout << "Plot solutions" << std::endl;
    writeDataToVTK3D_Pressure(P_3D, V_3D, N_3D, h_3D, timeStep);
  }
  */

  // copy the 3D pressure to libmesh pressure system
  util::set_elem_sol(pres_sys, P_3D);
}

void util::unet::Network::solve3D1DNutrientProblem(BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

     if (nut_sys.d_sys_name != "Nutrient" or tum_sys.d_sys_name != "Tumor")
         libmesh_error_msg("Must pass Nutrient and Tumor system to solve 3D-1D " "nutrient");

     const auto &input = d_model_p->get_input_deck();
     const auto timeStep = d_model_p->d_step;

     // Solver
     gmm::iteration iter(4.0e-11);

     assemble3D1DSystemForNutrients(nut_sys, tum_sys);

     // if this is first call inside nonlinear loop, we guess current
     // concentration as old concentration
 
     if (d_model_p->d_nonlinear_step == 0)
         phi_sigma = phi_sigma_old;
     if (d_model_p->d_step == 1)
         phi_sigma = b_nut_3D1D;

     gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_nut_3D1D, 60, 1e-8);
     gmm::bicgstab(A_nut_3D1D, phi_sigma, b_nut_3D1D, P, iter);

     auto pointer = VGM.getHead();
   
     while (pointer) {

            int indexOfNode = pointer->index;

		pointer->c_v = phi_sigma[N_tot_3D + indexOfNode];
		// std::cout << "index: " << pointer->index << " c_v: " << pointer->c_v << "
		// p_v: " << pointer->p_v << " coord: " << pointer->coord << std::endl;
		pointer = pointer->global_successor;
     }

     // do not modify old with current concentration as this solver could be
     // called inside nonlinear loop at given time step
     // Rather update the old with new where this solver is called at the end of nonlinear loop
     // phi_sigma_old = phi_sigma;

     // Extract nutrient concentrations
     //std::cout << " " << std::endl;
     //std::cout << "Extract nutrient concentrations" << std::endl;

     for (int i = 0; i < N_tot_3D; i++) {

          phi_sigma_3D[i] = phi_sigma[i];
     }

     if (timeStep % input.d_dt_output_interval == 0) {

         std::cout << " " << std::endl;
         std::cout << "Plot solutions" << std::endl;
         writeDataToVTK3D_Nutrients(phi_sigma_3D, N_3D, h_3D, timeStep);
        // writeDataToVTKTimeStep_VGM(timeStep);
     }

     // copy the 3D nutrient to libmesh nutrient system
     util::set_elem_sol(nut_sys, phi_sigma_3D);

}

void util::unet::Network::solveVGMforPressure(BaseAssembly &pres_sys) {

  if (pres_sys.d_sys_name != "Pressure")
    libmesh_error_msg("Must pass Pressure system to solve 1D pressure");

  // std::cout << "Assemble pressure matrix and right hand side" << std::endl;
  assembleVGMSystemForPressure(pres_sys);

  gmm::iteration iter(10E-18);

  //  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_VGM, 50,
  //  1e-5);

  gmm::identity_matrix P;

  P_v = b;

  gmm::bicgstab(A_VGM, P_v, b, P, iter);

  if (P_v.size() < 20) {
    oss << "        P_v = (" << util::io::printStr(P_v) << ")" << std::endl;
    d_model_p->d_delayed_msg = oss.str();
    util::clear_oss(oss);
  }

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->p_v = P_v[indexOfNode];

    pointer = pointer->global_successor;
  }

  //  std::cout << "Write vtk" << std::endl;
  //  writeDataToVTK_VGM();
}

void util::unet::Network::solveVGMforNutrient(BaseAssembly &pres_sys,
                                                BaseAssembly &nut_sys) {

  if (pres_sys.d_sys_name != "Pressure" or nut_sys.d_sys_name != "Nutrient")
    libmesh_error_msg("Must pass Pressure and Nutrient system to solve 1D "
                      "nutrient");

  // std::cout << " " << std::endl;
  // std::cout << "Assemble nutrient matrix and right hand side" << std::endl;
  assembleVGMSystemForNutrient(pres_sys, nut_sys);

  // if this is first call inside nonlinear loop, we guess current
  // concentration as old concentration
  if (d_model_p->d_nonlinear_step == 0)
    C_v = C_v_old;
  if (d_model_p->d_step == 1)
    C_v = b_c;

  // get preconditioner
  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(Ac_VGM, 50, 1e-6);

  // debug
  if (false) {
    out << "nonlinear iter: " << d_model_p->d_nonlinear_step << "\n";
    out << "Ac_VGM rows: " << Ac_VGM.nrows()
        << " Ac_VGM cols: " << Ac_VGM.ncols() << " C_v size: " << C_v.size()
        << " b_c size: " << b_c.size() << "\n";

    out << Ac_VGM << "\n\n";
    //    out << C_v << "\n\n";
    //    out << b_c << "\n\n";
  }

  // solve
  gmm::iteration iter(5.0e-11);
  gmm::bicgstab(Ac_VGM, C_v, b_c, P, iter);

  // std::cout << C_v << std::endl;

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

void util::unet::Network::assemble3D1DSystemForPressure(BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

  const auto &input = d_model_p->get_input_deck();
  double L_x = input.d_domain_params[1];

  // 3D-1D coupled flow problem on a unit cube
  //std::cout << " " << std::endl;
  //std::cout << "3D-1D coupled flow problem on a cube \Omega = (0," << L_x
  //          << ")^3" << std::endl;

  // Number of Elements (3D) in each space direction
  //std::cout << " " << std::endl;
  //std::cout << "Number of Elements (3D) in each space direction: " << N_3D
  //          << std::endl;
  //std::cout << "Total number of Elements in 3D: " << N_tot_3D << std::endl;

  // Mesh size
  //std::cout << " " << std::endl;
  //std::cout << "Mesh size h_3D: " << h_3D << std::endl;

  // Assemble 3D Matrix and right hand side (pressure)
  //std::cout << " " << std::endl;
  // It is K divided by mu
  double K_3D = input.d_tissue_flow_coeff;
  //std::cout << "Assemble 3D Matrix and right hand side (pressure)" <<
  //std::endl;
  //std::cout << "K_3D: " << K_3D << std::endl;

  int numberOfNodes = VGM.getNumberOfNodes();

  //std::cout << "numberOfUnknowns: " << N_tot_3D + numberOfNodes << std::endl;

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

            A_flow_3D1D(index, index) = A_flow_3D1D(index, index) + K_3D * h_3D * h_3D / h_3D;

            A_flow_3D1D(index, index_neighbor) = -K_3D * h_3D * h_3D / h_3D;

          } 
          else {

            if (scenario == "test_single_vessel") {

              std::vector<double> center_face =
                  getCenterFace(center, directions[face], h_3D);

              double dirichlet_value = getDirichletValue(
                  center_face, VGM.getHead()->L_p[0], VGM.getHead()->radii[0]);

              A_flow_3D1D(index, index) = A_flow_3D1D(index, index) + 2.0 * K_3D * h_3D * h_3D / h_3D;

              b_flow_3D1D[index] = b_flow_3D1D[index] + (2.0 * K_3D * h_3D * h_3D * dirichlet_value) / h_3D;

            }

          }

        }

      }

    }

  }

  // Assemble 1D Matrix, 1D right hand side and coupling matrices (pressure)
  //std::cout << " " << std::endl;
  //std::cout << "Assemble 1D Matrix and Coupling matrices (pressure)"
  //          << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    // std::cout << "\n----------------------------------------------\n";
    // std::cout << "indexOfNode: " << indexOfNode << std::endl;

    int numberOfNeighbors = pointer->neighbors.size();

    //std::cout << "numberOfNeighbors: " << numberOfNeighbors << std::endl;

    std::vector<double> coord = pointer->coord;

    //std::cout << "coord: " << coord << std::endl;

    if (numberOfNeighbors == 1) {

	if( pointer->typeOfVGNode == TypeOfNode::DirichletNode ){

	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = 1.0;
	    b_flow_3D1D[ N_tot_3D+indexOfNode ] = pointer->p_boundary;

	    double radius = pointer->radii[ 0 ];

	    // Assemble coupling terms (nutrients)
	    std::vector<double> coord_neighbor = pointer->neighbors[ 0 ]->coord;

	    int N_s = input.d_num_points_length;
	    int N_theta = input.d_num_points_angle;

	    std::vector<double> weights;
	    std::vector<int> id_3D_elements;

	    double length = util::dist_between_points( coord, coord_neighbor );
	    double length_edge = 0.5 * length;
	  
	    // Surface area of cylinder
	    double surface_area = 2.0*M_PI*length_edge*radius;

	    determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

	    // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
	    int numberOfElements = id_3D_elements.size();

	    // Permeabilty vessel wall for nutrients
	    double L_p = pointer->L_p[ 0 ];

	    for(int j=0;j<numberOfElements;j++){

		// A_3D1D
		A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - L_p*surface_area*weights[ j ];

		// A_3D3D
		A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) = A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) + L_p*surface_area*weights[ j ];

	    }

	}
	else if( pointer->typeOfVGNode == TypeOfNode::NeumannNode ){

	    double radius = pointer->radii[ 0 ];

	    // Assemble coupling terms (nutrients)
	    std::vector<double> coord_neighbor = pointer->neighbors[ 0 ]->coord;

	    int N_s = input.d_num_points_length;
	    int N_theta = input.d_num_points_angle;
	    int indexNeighbor = pointer->neighbors[ 0 ]->index;

	    std::vector<double> weights;
	    std::vector<int> id_3D_elements;

	    // Permeabilty vessel wall for nutrients
	    double L_p = pointer->L_p[ 0 ];
	    double length = util::dist_between_points( coord, coord_neighbor );
	    double K_1D = getK1D( 0.5*( coord[2]+coord_neighbor[2] ), L_p, radius );

	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + K_1D/length;
	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexNeighbor) = - K_1D/length;

	    double length_edge = 0.5 * length;
	  
	    // Surface area of cylinder
	    double surface_area = 2.0*M_PI*length_edge*radius;

	    determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

	    // Add coupling entry to 1D1D matrix
	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + L_p*surface_area;

	    // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
	    int numberOfElements = id_3D_elements.size();

	    for(int j=0;j<numberOfElements;j++){

		if( id_3D_elements[ j ]>-1 ){

		    // A_3D1D
		    A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - L_p*surface_area*weights[ j ];

		    // A_3D3D
		    A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  = A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  + L_p*surface_area*weights[ j ];

		    // A_1D3D
		    A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) = A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) - L_p*surface_area*weights[ j ];

		}

	    }

	}

    } 
    else {

      for (int i = 0; i < numberOfNeighbors; i++) {


	    // Discretisation of the differential operator
	    double radius = pointer->radii[ i ];

	    int indexNeighbor = pointer->neighbors[ i ]->index;

	    std::vector<double> coord_neighbor = pointer->neighbors[ i ]->coord;

	    double length = util::dist_between_points( coord, coord_neighbor );
	    double L_p = pointer->L_p[ i ];
	    double K_1D = getK1D( 0.5*( coord[2]+coord_neighbor[2] ), L_p, radius );

	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + K_1D/length;
	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexNeighbor) = - K_1D/length;

	    // Coupling terms   
	    int N_s = input.d_num_points_length;
	    int N_theta = input.d_num_points_angle;

	    std::vector<double> weights;
	    std::vector<int> id_3D_elements;

	    double length_edge = 0.5 * length;

	    // Surface area of cylinder
	    double surface_area = 2.0*M_PI*length_edge*radius;

	    determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

	    // Add coupling entry to 1D1D matrix
	    A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + L_p*surface_area;

	    // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
	    int numberOfElements = id_3D_elements.size();

	    for(int j=0;j<numberOfElements;j++){

		if( id_3D_elements[ j ]>-1 ){

		     // A_3D1D
		     A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - L_p*surface_area*weights[ j ];

		     // A_3D3D
		     A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  = A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  + L_p*surface_area*weights[ j ];

		     // A_1D3D
		     A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) = A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) - L_p*surface_area*weights[ j ];

		}

	    }

      }
              
      b_flow_3D1D[ N_tot_3D+indexOfNode ] = 0.0;

    }

    pointer = pointer->global_successor;

  }

}

void util::unet::Network::assemble3D1DSystemForNutrients(BaseAssembly &nut_sys, BaseAssembly &tum_sys) {

  const auto &input = d_model_p->get_input_deck();
  const auto &mesh = d_model_p->get_mesh();
  double L_x = input.d_domain_params[1];

  // 3D-1D coupled flow problem on a cube
  //std::cout << " " << std::endl;
  //  std::cout << "3D-1D coupled nutrient transport problem on a cube \Omega = (0,"
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
  //std::cout << "vol_elem: " << vol_elem << std::endl;

  double area_face = h_3D * h_3D;
  //std::cout << "area_face: " << area_face << std::endl;

  // Assemble 3D Matrix and right hand side (pressure)
  //std::cout << " " << std::endl;
  double K_3D = input.d_tissue_flow_coeff;
  //  std::cout << "Assemble 3D Matrix and right hand side (nutrients)"
  //            << std::endl;
  //  std::cout << "K_3D: " << K_3D << std::endl;
  //  std::cout << "D_v_3D: " << D_v_3D << std::endl;

  int numberOfNodes = VGM.getNumberOfNodes();

  //std::cout << "numberOfNodes: " << N_tot_3D + numberOfNodes << std::endl;

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
  //std::cout << "dt: " << dt << std::endl;

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

          std::vector<double> center_neighbor = getCenterNeighbor(center, directions[face], h_3D);

          bool isInnerFace = isCenterInDomain(center_neighbor, L_x);

          if (isInnerFace) {

            int index_neighbor = getElementIndex(center_neighbor, h_3D, N_3D);

            double v = -K_3D * (P_3D1D[index_neighbor] - P_3D1D[index]) / h_3D;

            if (v > 0.0) {

              A_nut_3D1D(index, index) += dt * area_face * v;

            } else {

              A_nut_3D1D(index, index_neighbor) += dt * area_face * v;
            }

            A_nut_3D1D(index, index) += dt * D_v_3D * area_face / h_3D;

            A_nut_3D1D(index, index_neighbor) = -dt * D_v_3D * area_face / h_3D;

          }

        }

        // Libmesh element
        auto elem = mesh.elem_ptr(i);

        // loop over sides of the element
        for (auto side : elem->side_index_range()) {

          if (elem->neighbor_ptr(side) != nullptr) {

            const Elem *neighbor = elem->neighbor_ptr(side);

            // div(chem_tum * grad(tum)) term
            // requires integration over face of an element
            tum_sys.d_fe_face->reinit(elem, side);

            // loop over quadrature points
            for (unsigned int qp = 0; qp < tum_sys.d_qrule_face.n_points(); qp++) {

              Real chem_tum_cur = 0.;
              Gradient tum_grad = 0.;
              for (unsigned int l = 0; l < tum_sys.d_phi_face.size(); l++) {

                chem_tum_cur +=
                    tum_sys.d_phi_face[l][qp] * tum_sys.get_current_sol_var(l, 1);

                tum_grad.add_scaled(tum_sys.d_dphi_face[l][qp],
                                    tum_sys.get_current_sol_var(l, 0));
              }

              // add to force
              b_nut_3D1D[index] += -tum_sys.d_JxW_face[qp] * input.d_tissue_flow_coeff *
                       chem_tum_cur * tum_grad * tum_sys.d_qface_normals[qp];
            } // loop over quadrature points on face

          } // elem neighbor is not null
        }   // loop over faces

      }

    }

  }

  // Assemble 1D Matrix, 1D right hand side and coupling matrices (nutrients)
  //  std::cout << " " << std::endl;
  //  std::cout << "Assemble 1D Matrix and Coupling matrices (nutrients)"
  //            << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v = pointer->p_v;

    double c_v_k = pointer->c_v;

    if (numberOfNeighbors == 1) {

	double radius = pointer->radii[0];

	int indexNeighbor = pointer->neighbors[0]->index;

	std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;

	double length = 0.0;

	for(int j = 0; j < 3; j++){

	    length += (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
	}

	length = std::sqrt(length);

	double p_neighbor = pointer->neighbors[0]->p_v;

	double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v);

	double volume = length / 2.0 * radius * radius * M_PI;

	if( v_interface > 0.0 && pointer->typeOfVGNode == TypeOfNode::DirichletNode ){

	    A_nut_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = 1.0;

	    double phi_sigma_boundary = 0.0;

	    if(p_v < input.d_identify_vein_pres){

	       phi_sigma_boundary = input.d_in_nutrient_vein;

	    }
	    else{

	       phi_sigma_boundary = input.d_in_nutrient;

	    }

	    b_nut_3D1D[N_tot_3D+indexOfNode] = phi_sigma_boundary;

	    // Assemble coupling terms (nutrients)
	    std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;

	    int N_s = input.d_num_points_length;
	    int N_theta = input.d_num_points_angle;

	    std::vector<double> weights;

	    std::vector<int> id_3D_elements;

	    double length_edge = 0.5 * length;
	  
	    // Surface area of cylinder
	    double surface_area = 2.0*M_PI*length_edge*radius;

	    determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

	    // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
	    int numberOfElements = id_3D_elements.size();

	    // Permeabilty vessel wall for nutrients
	    double L_s = pointer->L_s[ 0 ];

	    double p_t = 0.0;

	    for(int j=0;j<numberOfElements;j++){

		// A_3D1D
		A_nut_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_nut_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - dt*L_s*surface_area*weights[ j ];

		// A_3D3D
		A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  = A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  + dt*L_s*surface_area*weights[ j ];

	    }

	} 
	else{

	   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = length;

	   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexNeighbor) = dt * v_interface - dt * D_v / length;

	   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) - dt * v_interface + dt * D_v / length;

	   b_nut_3D1D[N_tot_3D+indexOfNode] = length * phi_sigma_old[N_tot_3D+indexOfNode];

	   // Assemble coupling terms (nutrients)

	   int N_s = input.d_num_points_length; 
	   int N_theta = input.d_num_points_angle;

	   std::vector<double> weights;

	   std::vector<int> id_3D_elements;

	   double length_edge = 0.5 * length;
	  
	   // Surface area of cylinder
	   double surface_area = 2.0*M_PI*length_edge*radius;

	   determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

	   // Permeabilty vessel wall for nutrients
	   double L_s = pointer->L_s[ 0 ];

	   // 1D part of the coupling
	   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) + dt * L_s * surface_area;

	   // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
	   int numberOfElements = id_3D_elements.size();

	   double p_t = 0.0;

	   for(int j=0;j<numberOfElements;j++){

	       // A_3D1D
	       A_nut_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_nut_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - dt*L_s*surface_area*weights[ j ];

	       // A_3D3D
	       A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  = A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  + dt*L_s*surface_area*weights[ j ];

	       // A_1D3D
	       A_nut_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) = A_nut_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) - dt*L_s*surface_area*weights[ j ];

	       p_t = P_3D1D[ id_3D_elements[ j ] ];

	       if( p_v-p_t > 0.0 ){

		   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) + 
					                                         dt * (1. - osmotic_sigma) * pointer->L_p[ 0 ] * surface_area * weights[ j ] * (p_v - p_t);

	       }
	       else{

		   A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) = A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) - 
					                                         dt * (1. - osmotic_sigma) * pointer->L_p[ 0 ] * surface_area * weights[ j ] * (p_v - p_t);    

	       }

	   }

	}

    } 
    else {

      for (int i = 0; i < numberOfNeighbors; i++) {

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) = 0.0;

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = 0.0;

        for (int j = 0; j < 3; j++) {

          length +=
              (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
        }

        length = std::sqrt(length);

        double p_neighbor = pointer->neighbors[i]->p_v;

        double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) *
                             (p_neighbor - p_v);

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) += length;

        if (v_interface > 0.0) {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
              dt * v_interface;

        } else {

          A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) +=
              dt * v_interface;
        }

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +=
            dt * D_v / length;

        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexNeighbor) -=
            dt * D_v / length;

        b_nut_3D1D[N_tot_3D + indexOfNode] +=
            length * phi_sigma_old[N_tot_3D + indexOfNode];

        // Assemble coupling terms (nutrients)

        int N_s = input.d_num_points_length;

        int N_theta = input.d_num_points_angle;

        std::vector<double> weights;

        std::vector<int> id_3D_elements;

        double length_edge = 0.0;

        determineWeightsAndIds(N_s, N_theta, N_3D, coord, coord_neighbor,
                               radius, h_3D, length_edge, weights,
                               id_3D_elements);

        // Surface area of cylinder
        double surface_area = 2.0 * M_PI * 0.5 * length_edge * radius;

        // Permeabilty vessel wall for nutrients
        double L_s = pointer->L_s[i];

        // 1D part of the coupling
        A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
            A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
            dt * L_s * surface_area;

        // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
        int numberOfElements = id_3D_elements.size();

        double p_t = 0.0;

        for (int j = 0; j < numberOfElements; j++) {

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

          p_t = P_3D1D[id_3D_elements[j]];

          if (p_v - p_t > 0.0) {

            A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) =
                A_nut_3D1D(N_tot_3D + indexOfNode, N_tot_3D + indexOfNode) +
                dt * (1. - osmotic_sigma) * pointer->L_p[i] * surface_area *
                    weights[j] * (p_v - p_t);

          } else {

            A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) =
                A_nut_3D1D(id_3D_elements[j], id_3D_elements[j]) -
                dt * (1. - osmotic_sigma) * pointer->L_p[i] * surface_area *
                    weights[j] * (p_v - p_t);
          }

        }

      }

    }

    pointer = pointer->global_successor;

  }

}

void util::unet::Network::assembleVGMSystemForPressure(
    BaseAssembly &pres_sys) {

  const auto &mesh = d_model_p->get_mesh();

  const auto &input = d_model_p->get_input_deck();

  // factor to enhance condition of matrix
  const double factor_p = input.d_assembly_factor_p_t;

  int numberOfNodes = VGM.getNumberOfNodes();

  // reinitialize data
  //  A_VGM =
  //      gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  if (A_VGM.nrows() != numberOfNodes) {

    A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  } else {

    // A_VGM.clear_mat();

    for (unsigned int i = 0; i < A_VGM.nrows(); i++)
      A_VGM[i].clear();
  }

  if (b.size() != numberOfNodes) {

    b.resize(numberOfNodes);
  }

  for (unsigned int i = 0; i < b.size(); i++) {

    b[i] = 0.;
  }

  if (P_v.size() != numberOfNodes) {

    P_v.resize(numberOfNodes);
  }

  for (unsigned int i = 0; i < P_v.size(); i++) {

    P_v[i] = 0.;
  }

  std::vector<ElemWeights> J_b_points;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v_k = pointer->p_v;

    // std::cout << "p_v_k: " << p_v_k << std::endl;

    if (numberOfNeighbors == 1) {

      A_VGM(indexOfNode, indexOfNode) = 1.0;

      b[indexOfNode] = pointer->p_boundary;

    } else {

      // get element data at points on cylinder surface
      if (input.d_compute_elem_weights)
        J_b_points = pointer->J_b_points;
      else
        J_b_points = compute_elem_weights_at_node(pointer);

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = J_b_points[i];

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = 0.0;

        for (int j = 0; j < 3; j++) {

          length +=
              (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
        }

        length = std::sqrt(length);

        // multiplying both sides by factor
        double t_seg =
            factor_p * M_PI * std::pow(radius, 4) / (8.0 * length * mu);

        A_VGM(indexOfNode, indexOfNode) += t_seg;
        A_VGM(indexOfNode, indexNeighbor) -= t_seg;

        // implicit for p_v in source
        A_VGM(indexOfNode, indexOfNode) +=
            factor_p * input.d_coupling_method_theta * pointer->L_p[i] *
            J_b_data.half_cyl_surf;

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = mesh.elem_ptr(e_id);
          pres_sys.init_dof(elem);
          auto p_t_k = pres_sys.get_current_sol(0);

          // explicit for p_t in source
          b[indexOfNode] -=
              factor_p * pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
              ((1. - input.d_coupling_method_theta) * p_v_k - p_t_k);
        }

      } // loop over neighbor segments

      // b[indexOfNode] = 0.0;
    }

    //    std::cout << "p_v_k: " << pointer->p_v << " s: " << pointer->coord[2]
    //              << ", b: " << b[indexOfNode]
    //              << ", L_p: " << util::io::printStr(pointer->L_p) <<
    //              std::endl;

    pointer = pointer->global_successor;
  }

  // std::cout << "A_VGM: " << A_VGM << std::endl;
}

void util::unet::Network::assembleVGMSystemForNutrient(
    BaseAssembly &pres_sys, BaseAssembly &nut_sys) {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  // factor to enhance condition of matrix
  const double factor_c = input.d_assembly_factor_c_t;
  double coupling_theta = input.d_coupling_method_theta;
  if (input.d_decouple_nutrients)
    coupling_theta = 0.;

  int numberOfNodes = VGM.getNumberOfNodes();

  // reinitialize data
  //  Ac_VGM =
  //      gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  if (Ac_VGM.nrows() != numberOfNodes) {

    Ac_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  } else {

    // Ac_VGM.clear_mat();

    for (unsigned int i = 0; i < Ac_VGM.nrows(); i++)
      Ac_VGM[i].clear();
  }

  if (b_c.size() != numberOfNodes) {

    b_c.resize(numberOfNodes);
  }

  for (unsigned int i = 0; i < b_c.size(); i++) {

    b_c[i] = 0.;
  }

  if (C_v.size() != numberOfNodes) {

    C_v.resize(numberOfNodes);
  }

  //  for (unsigned int i=0; i<C_v.size(); i++) {
  //
  //    C_v[i] = 0.;
  //  }

  std::vector<ElemWeights> J_b_points;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double dt = d_model_p->d_dt;

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v_k = pointer->p_v;

    double c_v_k = pointer->c_v;

    double c_couple = c_v_k;
    if (input.d_decouple_nutrients)
      c_couple = C_v_old[indexOfNode];

    // init
    Ac_VGM(indexOfNode, indexOfNode) = 0.;

    if (numberOfNeighbors == 1) {

      double radius = pointer->radii[0];

      int indexNeighbor = pointer->neighbors[0]->index;

      std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;

      double length = 0.0;

      for (int j = 0; j < 3; j++) {

        length +=
            (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
      }

      length = std::sqrt(length);

      double p_neighbor = pointer->neighbors[0]->p_v;

      double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) *
                           (p_neighbor - p_v_k);

      double volume = length / 2.0 * radius * radius * M_PI;

      // init
      Ac_VGM(indexOfNode, indexNeighbor) = 0.;

      if (v_interface > 0.0) {
        // inlet

        // if artery
        if (p_v_k >= input.d_identify_vein_pres) {

          // Dirichlet on inlet
          Ac_VGM(indexOfNode, indexOfNode) += factor_c * 1.0;
          b_c[indexOfNode] += factor_c * input.d_in_nutrient;

        }
        else {

          // Dirchlet on inlet
          bool dirichlet_on_vein = false;
          if (dirichlet_on_vein) {
            Ac_VGM(indexOfNode, indexOfNode) += factor_c * 1.0;
            b_c[indexOfNode] += factor_c * input.d_in_nutrient_vein;
          } else {

            // mass matrix
            Ac_VGM(indexOfNode, indexOfNode) += factor_c * length;

            // diffusion
            Ac_VGM(indexOfNode, indexOfNode) += factor_c * dt * D_v / length;
            Ac_VGM(indexOfNode, indexNeighbor) += - factor_c * dt * D_v /
                                                  length;

            // advection
            Ac_VGM(indexOfNode, indexOfNode) += factor_c * dt * v_interface;
            Ac_VGM(indexOfNode, indexNeighbor) -= factor_c * dt * v_interface;

            // old time step
            b_c[indexOfNode] += factor_c * length * C_v_old[indexOfNode];
          }
        }

      } else {
        // outlet

        // mass matrix
        Ac_VGM(indexOfNode, indexOfNode) += factor_c * length;

        // diffusion
        Ac_VGM(indexOfNode, indexOfNode) += factor_c * dt * D_v / length;
        Ac_VGM(indexOfNode, indexNeighbor) += - factor_c * dt * D_v / length;

        // advection
        Ac_VGM(indexOfNode, indexOfNode) -= factor_c * dt * v_interface;
        Ac_VGM(indexOfNode, indexNeighbor) += factor_c * dt * v_interface;

        // old time step
        b_c[indexOfNode] += factor_c * length * C_v_old[indexOfNode];
      }

    } else {

      // get element data at points on cylinder surface
      if (input.d_compute_elem_weights)
        J_b_points = pointer->J_b_points;
      else
        J_b_points = compute_elem_weights_at_node(pointer);

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = J_b_points[i];

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

        // init
        Ac_VGM(indexOfNode, indexNeighbor) = 0.0;

        std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

        double length = 0.0;

        for (int j = 0; j < 3; j++) {

          length +=
              (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
        }

        length = std::sqrt(length);

        double p_neighbor = pointer->neighbors[i]->p_v;

        double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) *
                             (p_neighbor - p_v_k);

        Ac_VGM(indexOfNode, indexOfNode) += factor_c * length;

        if (v_interface > 0.0) {

          Ac_VGM(indexOfNode, indexOfNode) += factor_c * dt * v_interface;

        } else {

          Ac_VGM(indexOfNode, indexNeighbor) += factor_c * dt * v_interface;
        }

        Ac_VGM(indexOfNode, indexOfNode) += factor_c * dt * D_v / length;

        Ac_VGM(indexOfNode, indexNeighbor) -= factor_c * dt * D_v / length;

        b_c[indexOfNode] += factor_c * length * C_v_old[indexOfNode];

        // coupling between 3d and 1d nutrient
        {
          // implicit part of the coupling
          Ac_VGM(indexOfNode, indexOfNode) +=
              factor_c * dt * coupling_theta * pointer->L_s[i] *
              J_b_data.half_cyl_surf;

          // compute explicit part of the coupling
          for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

            auto e_id = J_b_data.elem_id[e];
            auto e_w = J_b_data.elem_weight[e];

            // get 3d pressure
            const auto *elem = mesh.elem_ptr(e_id);
            pres_sys.init_dof(elem);
            Real p_t_k = pres_sys.get_current_sol(0);

            // get 3d nutrient
            nut_sys.init_dof(elem);
            Real c_t_k = nut_sys.get_current_sol(0);

            // explicit
            double source =
                dt * pointer->L_s[i] * J_b_data.half_cyl_surf * e_w *
                ((1. - coupling_theta) * c_couple - c_t_k);
            b_c[indexOfNode] -= factor_c * source;

            //            if (pointer->p_v <
            //                    input.d_identify_vein_pres &&
            //                i == 0 && e < 2)
            //              out << "index: " << indexOfNode << ", neighbor: " <<
            //              indexNeighbor
            //                  << ", source: " << source << ", pressure: " <<
            //                  p_v_k
            //                  << ", radius: " << radius << ", c_t: " << c_t_k
            //                  << ", L_s: " << pointer->L_s[i] << "\n";

            // term due to pressure difference
            double c_transport = 0.;
            if (p_v_k - p_t_k >= 0.)
              c_transport = c_couple;
            else
              c_transport = c_t_k;
            b_c[indexOfNode] -= factor_c * dt * (1. - osmotic_sigma) *
                                pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                                (p_v_k - p_t_k) * c_transport;
          }
        } // coupling
      }
    }

//    if (indexOfNode == 14 or indexOfNode == 15 or indexOfNode == 21) {
//      out << "\ni: " << indexOfNode << ", bi: " << b_c[indexOfNode]
//          << ", Aii: " << Ac_VGM(indexOfNode, indexOfNode);
//      for (int i = 0; i < numberOfNeighbors; i++) {
//        int indexNeighbor = pointer->neighbors[i]->index;
//        out << ", Ai" << indexNeighbor << ": " << Ac_VGM(indexOfNode,
//                                                         indexNeighbor) ;
//      }
//      out << "\n\n";
//    }

    pointer = pointer->global_successor;
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

  double K_1D = 0.0;

  if (scenario == "test_single_vessel") {

    K_1D = L_p / (1.0 + L_p) * (1.0 + s + 0.5 * s * s);

  } else {

    K_1D = (radius * radius * radius * radius * M_PI) / (8.0 * mu);
  }

  return K_1D;
}
