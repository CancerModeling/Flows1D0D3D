////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "network_data_structure.cpp"
#include "netUtil.hpp"
#include "../model.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>


void netfcfvfe::Network::create_initial_network() {

     const auto &input = d_model_p->get_input_deck();
     d_is_network_changed = true;

     // equation system
     std::vector< std::vector<double> > vertices;
     std::vector< double > pressures;
     std::vector< double > radii;
     std::vector< std::vector<unsigned int> > elements;

     scenario = input.d_scenario;

     std::cout << " " << std::endl;
     std::cout << "Scenario: " << scenario << std::endl;

     readData(vertices, pressures, radii, elements);

     transferDataToVGM(vertices, pressures, radii, elements);

     int numberOfNodes = VGM.getNumberOfNodes();

     std::cout << " " << std::endl;
     std::cout << "Number of nodes: " << numberOfNodes << std::endl;

     // refine mesh
     std::cout << " " << std::endl;
     std::cout << "Refine mesh" << std::endl;

     int refinementLevel = input.d_network_init_refinement;

     for(int i=0;i<refinementLevel;i++){

         refine1DMesh();

     }

     numberOfNodes = VGM.getNumberOfNodes();

     std::cout << " " << std::endl;
     std::cout << "Number of nodes: " << numberOfNodes << std::endl;

     // Print data
     std::cout << " " << std::endl;
     std::cout << "Print data" << std::endl;
     //printDataVGM();

     // initialize matrix and vector
     // 1D nutrient: matrix, rhs, and solution
     Ac_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
     b_c = std::vector<double>(numberOfNodes, 0.0);

     C_v     = std::vector<double>(numberOfNodes, 0.0);
     C_v_old = std::vector<double>(numberOfNodes, 0.0);

     // 3D1D flow problem: matrix, rhs and solution
     N_3D = input.d_num_elems;
     N_tot_3D = N_3D * N_3D * N_3D;

     double L_x = input.d_domain_params[1];
     h_3D = L_x/(double) N_3D;

     A_flow_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D+numberOfNodes, N_tot_3D+numberOfNodes);
     b_flow_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     P_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);
     P_3D = std::vector<double>(N_tot_3D, 0.0);

     A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D+numberOfNodes, N_tot_3D+numberOfNodes);
     b_nut_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     A_TAF_3D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D, N_tot_3D);
     b_TAF_3D = std::vector<double>(N_tot_3D, 0.0);

     phi_sigma = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);
     phi_sigma_old = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);
     phi_sigma_3D = std::vector<double>(N_tot_3D, 0.0);

     phi_P = std::vector<double>(N_tot_3D, 0.0);
     phi_P_old = std::vector<double>(N_tot_3D, 0.0);

     phi_TAF = std::vector<double>(N_tot_3D, 0.0);
     phi_TAF_old = std::vector<double>(N_tot_3D, 0.0);

     sol_Vec_Phi_P = std::vector<double>(2*N_tot_3D, 0.0);
     sol_Vec_Phi_P_old = std::vector<double>(2*N_tot_3D, 0.0);

     mu = input.d_init_vessel_mu;

     D_v = input.d_D_sigma_v;

     D_TAF = input.d_D_TAF;

     osmotic_sigma = input.d_osmotic_sigma;

}

void netfcfvfe::Network::writeDataToVTKTimeStep_VGM(int timeStep) {

     std::string path = scenario;
     path += "_1D_";
     path += std::to_string(timeStep);
     path.append(".vtk");

     std::fstream filevtk;
     filevtk.open(path, std::ios::out);
     filevtk << "# vtk DataFile Version 2.0" << std::endl;
     filevtk << "Network Nutrient Transport" << std::endl;
     filevtk << "ASCII" << std::endl;
     filevtk << "DATASET POLYDATA" << std::endl;

     filevtk << "POINTS " << VGM.getNumberOfNodes() << " float" << std::endl;

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     while( pointer ) {

            std::vector<double> coord;

            coord = pointer->coord;

            filevtk << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;

            pointer = pointer->global_successor;

     }

     int numberOfNodes = 0;

     pointer = VGM.getHead();

     while( pointer ) {

            int numberOfNeighbors = pointer->neighbors.size();

            for (int i = 0; i < numberOfNeighbors; i++) {

                 if (!pointer->edge_touched[i]) {

                     numberOfNodes = numberOfNodes + 1;

                     pointer->edge_touched[i] = true;

                     pointer->neighbors[i]->markEdge(pointer->index);

                 }

            }

            pointer = pointer->global_successor;

     }

     pointer = VGM.getHead();

     while( pointer ){

            int numberOfNeighbors = pointer->neighbors.size();

            for(int i = 0; i < numberOfNeighbors; i++) {

                pointer->edge_touched[i] = false;
   
            }

            pointer = pointer->global_successor;

     }

     int polygons = 3 * numberOfNodes;

     filevtk << " " << std::endl;
     filevtk << "LINES " << numberOfNodes << " " << polygons << std::endl;

     pointer = VGM.getHead();

     while( pointer ){

            int numberOfNeighbors = pointer->neighbors.size();

            for(int i = 0; i < numberOfNeighbors; i++){

                if (!pointer->edge_touched[i]) {

                    filevtk << "2 " << pointer->index << " " << pointer->neighbors[i]->index  << std::endl;

                    pointer->edge_touched[i] = true;

                    pointer->neighbors[i]->markEdge(pointer->index);

                }
            }

            pointer = pointer->global_successor;

     }

     pointer = VGM.getHead();

     while( pointer ){

            int numberOfNeighbors = pointer->neighbors.size();

            for(int i = 0; i < numberOfNeighbors; i++){

                pointer->edge_touched[i] = false;

            }

            pointer = pointer->global_successor;

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_DATA " << numberOfNodes << std::endl;
     filevtk << "SCALARS radii float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     pointer = VGM.getHead();

     while( pointer ) {

            int numberOfNeighbors = pointer->neighbors.size();

            for(int i = 0; i < numberOfNeighbors; i++){

                if(!pointer->edge_touched[i]) {

                   filevtk << pointer->radii[i] << std::endl;

                   pointer->edge_touched[i] = true;

                   pointer->neighbors[i]->markEdge(pointer->index);

                }
            }

            pointer = pointer->global_successor;
     }

     pointer = VGM.getHead();

     while( pointer ){

            int numberOfNeighbors = pointer->neighbors.size();

            for(int i = 0; i < numberOfNeighbors; i++) {

                pointer->edge_touched[i] = false;

            }

            pointer = pointer->global_successor;

     }

     filevtk << " " << std::endl;
     filevtk << "POINT_DATA " << VGM.getNumberOfNodes() << std::endl;
     filevtk << "SCALARS Pressure_(1D) float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     pointer = VGM.getHead();

     while( pointer ) {

            filevtk << pointer->p_v << std::endl;

            pointer = pointer->global_successor;

     }

     filevtk << "SCALARS Phi_Sigma_(1D) float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     pointer = VGM.getHead();

     while( pointer ) {

            filevtk << pointer->c_v << std::endl;

            pointer = pointer->global_successor;

     }

}

void netfcfvfe::Network::assemble3D1DSystemForPressure(){

     const auto &input = d_model_p->get_input_deck();

     double L_x = input.d_domain_params[1];

     // 3D-1D coupled flow problem on a unit cube 
     std::cout << " " << std::endl;
     std::cout << "3D-1D coupled flow problem on a cube \Omega = (0," << L_x << ")^3" << std::endl;

     // Number of Elements (3D) in each space direction
     std::cout << " " << std::endl;
     std::cout << "Number of Elements (3D) in each space direction: " << N_3D << std::endl;
     std::cout << "Total number of Elements in 3D: " << N_tot_3D << std::endl;

     // Mesh size
     std::cout << " " << std::endl;
     std::cout << "Mesh size h_3D: " << h_3D << std::endl;

     // Assemble 3D Matrix and right hand side (pressure)
     std::cout << " " << std::endl;
     double K_3D = input.d_tissue_flow_K;
     std::cout << "Assemble 3D Matrix and right hand side (pressure)" << std::endl;
     std::cout << "K_3D: " << K_3D << std::endl;

     int numberOfNodes = VGM.getNumberOfNodes();

     std::cout << "numberOfUnknowns: " << N_tot_3D+numberOfNodes << std::endl;

     for(int i = 0; i < A_flow_3D1D.nrows(); i++){
 
         A_flow_3D1D[i].clear();

     }

     for(int i = 0; i < b_flow_3D1D.size(); i++) {

         b_flow_3D1D[i] = 0.0;

     } 

     std::vector<  std::vector< double > > directions;
     directions = defineDirections();

     for(int i=0;i<N_3D;i++){ // x-loop

         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop

                 int index = i + j*N_3D + k*N_3D*N_3D;

                 b_flow_3D1D[ index ] = 0.0;

                 // Get element center
                 std::vector< double > center = getElementCenter( i, j, k, h_3D );

                 // Iterate over the interfaces
                 for(int face=0;face<6;face++){

                     std::vector< double > center_neighbor = getCenterNeighbor( center, directions[ face ], h_3D );
    
                     bool isInnerFace = isCenterInDomain( center_neighbor, L_x );

                     if( isInnerFace ){

                       //  std::cout << "index: " << index << std::endl;

                         int index_neighbor = getElementIndex( center_neighbor, h_3D, N_3D );
                         
                         A_flow_3D1D( index, index ) = A_flow_3D1D( index, index ) + K_3D*h_3D*h_3D/h_3D;

                       //  std::cout << "index_neighbor: " << index_neighbor << std::endl;

                         A_flow_3D1D( index, index_neighbor ) = -K_3D*h_3D*h_3D/h_3D;

                         // std::cout << "A_flow_3D1D( index, index_neighbor ): " << A_flow_3D1D( index, index_neighbor ) << std::endl;

                     }
                     else{

                         if( scenario == "test_single_vessel" ){

                             std::vector< double > center_face = getCenterFace( center, directions[ face ], h_3D );

                             double dirichlet_value = getDirichletValue( center_face, VGM.getHead()->L_p[ 0 ], VGM.getHead()->radii[ 0 ] );
                    
                             A_flow_3D1D( index, index ) = A_flow_3D1D( index, index ) + 2.0*K_3D*h_3D*h_3D/h_3D;

                             b_flow_3D1D[ index ] = b_flow_3D1D[ index ] + (2.0*K_3D*h_3D*h_3D*dirichlet_value)/h_3D;

                         }

                     }

                     //std::cout << "A_flow_3D1D( index, index ): " << A_flow_3D1D( index, index ) << std::endl;

                 }

             }

         }

     }

     // Assemble 1D Matrix, 1D right hand side and coupling matrices (pressure)
     std::cout<< " " << std::endl;
     std::cout << "Assemble 1D Matrix and Coupling matrices (pressure)" << std::endl;

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     while( pointer ){

            int indexOfNode = pointer->index;

            //std::cout << "indexOfNode: " << indexOfNode << std::endl;

            int numberOfNeighbors = pointer->neighbors.size();

            //std::cout << "numberOfNeighbors: " << numberOfNeighbors << std::endl;

            std::vector<double> coord = pointer->coord;

            //std::cout << "coord: " << coord << std::endl;

            if( numberOfNeighbors == 1 ){

                A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = 1.0;

                b_flow_3D1D[ N_tot_3D+indexOfNode ] = pointer->p_boundary;
                
                int indexNeighbor = pointer->neighbors[ 0 ]->index;

                //std::cout << "indexNeighbor: " << indexNeighbor << std::endl;

            }
            else{

                for( int i=0; i<numberOfNeighbors; i++ ){

                     // Discretisation of the differential operator

                     double radius = pointer->radii[ i ];

                     //std::cout << "radius: " << radius << std::endl;

                     int indexNeighbor = pointer->neighbors[ i ]->index;

                     //std::cout << "indexNeighbor: " << indexNeighbor << std::endl;

                     std::vector<double> coord_neighbor = pointer->neighbors[ i ]->coord;

                     double length = 0.0;

                     for(int j=0;j<3;j++){

                         length += (coord[j]-coord_neighbor[j])*(coord[j]-coord_neighbor[j]);

                     }

                     length = std::sqrt( length );

                     double L_p = pointer->L_p[ i ];

                     double K_1D = getK1D( 0.5*( coord[2]+coord_neighbor[2] ), L_p, radius );

                     A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + K_1D/length;

                     A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexNeighbor) = - K_1D/length;
   
                     // Coupling terms
   
                     int N_s = input.d_num_points_length;

                     int N_theta = input.d_num_points_angle;

                     std::vector<double> weights;

                     std::vector<int> id_3D_elements;

                     double length_edge = 0.0;

                     determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);

                     // Surface area of cylinder
                     double surface_area = 2.0*M_PI*0.5*length_edge*radius;

                     // Add coupling entry to 1D1D matrix
                     A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = A_flow_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) + L_p*surface_area;

                     // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
                     int numberOfElements = id_3D_elements.size();

                     for(int j=0;j<numberOfElements;j++){

                         // A_3D1D
                         A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) = A_flow_3D1D(id_3D_elements[ j ],N_tot_3D+indexOfNode) - L_p*surface_area*weights[ j ];

                         // A_3D3D
                         A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  = A_flow_3D1D(id_3D_elements[ j ],id_3D_elements[ j ])  + L_p*surface_area*weights[ j ];

                         // A_1D3D
                         A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) = A_flow_3D1D(N_tot_3D+indexOfNode,id_3D_elements[ j ]) - L_p*surface_area*weights[ j ];

                     }
        
                }

                b_flow_3D1D[ N_tot_3D+indexOfNode ] = 0.0;

            }

            pointer = pointer->global_successor;

     }

}

void netfcfvfe::Network::assemble3D1DSystemForNutrients(){

     const auto &input = d_model_p->get_input_deck();

     double L_x = input.d_domain_params[1];

     // 3D-1D coupled flow problem on a cube 
     std::cout << " " << std::endl;
     std::cout << "3D-1D coupled nutrient transport problem on a cube \Omega = (0," << L_x << ")^3" << std::endl;

     // Number of Elements (3D) in each space direction
     std::cout << " " << std::endl;
     std::cout << "Number of Elements (3D) in each space direction: " << N_3D << std::endl;
     std::cout << "Total number of Elements in 3D: " << N_tot_3D << std::endl;

     // Mesh size
     std::cout << " " << std::endl;
     std::cout << "Mesh size h_3D: " << h_3D << std::endl;

     double vol_elem = h_3D * h_3D * h_3D;
     std::cout << "vol_elem: " << vol_elem << std::endl;

     double area_face = h_3D * h_3D;
     std::cout << "area_face: " << area_face << std::endl;

     // Assemble 3D Matrix and right hand side (pressure)
     std::cout << " " << std::endl;
     double K_3D = input.d_tissue_flow_K;
     std::cout << "Assemble 3D Matrix and right hand side (nutrients)" << std::endl;
     std::cout << "K_3D: " << K_3D << std::endl;
     std::cout << "D_v: " << D_v << std::endl;

     int numberOfNodes = VGM.getNumberOfNodes();

     std::cout << "numberOfNodes: " << N_tot_3D+numberOfNodes << std::endl;

     A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D+numberOfNodes, N_tot_3D+numberOfNodes);

     for(int i = 0; i < A_nut_3D1D.nrows(); i++){
 
         A_nut_3D1D[i].clear();

     }

     for(int i = 0; i < b_nut_3D1D.size(); i++) {

         b_nut_3D1D[i] = 0.0;

     }

     std::vector<  std::vector< double > > directions;

     directions = defineDirections();

     double dt = d_model_p->d_dt;
     std::cout << "dt: " << dt << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop

         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop

                 int index = i + j*N_3D + k*N_3D*N_3D;

                 A_nut_3D1D( index, index ) += vol_elem;

                 b_nut_3D1D[ index ] = vol_elem * phi_sigma_old[ index ];

                 // Get element center
                 std::vector< double > center = getElementCenter( i, j, k, h_3D );

                 // Iterate over the interfaces
                 for(int face=0;face<6;face++){

                     std::vector< double > center_neighbor = getCenterNeighbor( center, directions[ face ], h_3D );
    
                     bool isInnerFace = isCenterInDomain( center_neighbor, L_x );

                     if( isInnerFace ){

                         int index_neighbor = getElementIndex( center_neighbor, h_3D, N_3D );

                         double v = -K_3D*( P_3D1D[ index_neighbor ] - P_3D1D[ index ] )/h_3D;

                         if( v > 0.0 ){

                             A_nut_3D1D( index, index ) += dt * area_face * v;

                         }
                         else{

                             A_nut_3D1D( index, index_neighbor ) += dt * area_face * v;

                         }

                         A_nut_3D1D( index, index ) += dt * D_v * area_face/h_3D;

                         A_nut_3D1D( index, index_neighbor ) = - dt * D_v * area_face/h_3D;

                     }

                 }

             }

         }

     }

     // Assemble 1D Matrix, 1D right hand side and coupling matrices (nutrients)
     std::cout<< " " << std::endl;
     std::cout << "Assemble 1D Matrix and Coupling matrices (nutrients)" << std::endl;

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     while( pointer ){

            int indexOfNode = pointer->index;

            int numberOfNeighbors = pointer->neighbors.size();

            std::vector<double> coord = pointer->coord;

            double p_v = pointer->p_v;

            double c_v_k = pointer->c_v;

            if( numberOfNeighbors == 1 ){

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

                if(v_interface > 0.0){

                   A_nut_3D1D(N_tot_3D+indexOfNode,N_tot_3D+indexOfNode) = 1.0;

                   if(p_v < input.d_identify_vein_pres){

                      b_nut_3D1D[N_tot_3D+indexOfNode] = input.d_in_nutrient_vein;

                   }
                   else{

                      b_nut_3D1D[N_tot_3D+indexOfNode] = input.d_in_nutrient;

                   }

                } 
                else{

                   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = length;

                   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexNeighbor) = dt * v_interface - dt * D_v / length;

                   A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) = A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) - dt * v_interface + dt * D_v / length;

                   b_nut_3D1D[N_tot_3D+indexOfNode] = length * phi_sigma_old[N_tot_3D+indexOfNode];

                }

            }
            else{

                for(int i = 0; i < numberOfNeighbors; i++) {

                    double radius = pointer->radii[i];

                    int indexNeighbor = pointer->neighbors[i]->index;

                    A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexNeighbor) = 0.0;

                    std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

                    double length = 0.0;

                    for(int j=0; j<3; j++) {

                        length += (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);

                    }

                    length = std::sqrt(length);

                    double p_neighbor = pointer->neighbors[i]->p_v;

                    double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v);

                    A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) += length;

                    if(v_interface > 0.0){

                       A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) += dt * v_interface;

                    } 
                    else{

                       A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexNeighbor) += dt * v_interface;

                    }

                    A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexOfNode) += dt * D_v / length;
 
                    A_nut_3D1D(N_tot_3D+indexOfNode, N_tot_3D+indexNeighbor) -= dt * D_v / length;

                    b_nut_3D1D[N_tot_3D+indexOfNode] += length * phi_sigma_old[N_tot_3D+indexOfNode];

                    // Assemble coupling terms (nutrients)
   
                    int N_s = input.d_num_points_length;

                    int N_theta = input.d_num_points_angle;

                    std::vector<double> weights;

                    std::vector<int> id_3D_elements;

                    double length_edge = 0.0;

                    determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);
                  
                    // Surface area of cylinder
                    double surface_area = 2.0*M_PI*0.5*length_edge*radius;

                    // Permeabilty vessel wall for nutrients
                    double L_s = pointer->L_s[i];

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
                                                                                         dt * (1. - osmotic_sigma) * pointer->L_p[i] * surface_area * weights[ j ] * (p_v - p_t); 

                         }
		         else{

                             A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) = A_nut_3D1D(id_3D_elements[ j ],id_3D_elements[ j ]) - 
                                                                                      dt * (1. - osmotic_sigma) * pointer->L_p[i] * surface_area * weights[ j ] * (p_v - p_t);    

                         }

                    }

                }

            }

            pointer = pointer->global_successor;

     }

}

void netfcfvfe::Network::assembleVGMSystemForNutrients(){

     const auto &tum_sys = d_model_p->get_system();
     const auto &mesh = d_model_p->get_mesh();

     // 3d pressure
     const auto &pres = tum_sys.get_system<TransientLinearImplicitSystem>("Pressure");
     const DofMap &pres_map = pres.get_dof_map();
     const unsigned int v_pres = pres.variable_number("pressure");

     // 3d nutrient
     const auto &nut = tum_sys.get_system<TransientLinearImplicitSystem>("Nutrient");
     const DofMap &nut_map = nut.get_dof_map();
     const unsigned int v_nut = nut.variable_number("nutrient");

     const auto &input = d_model_p->get_input_deck();

     int numberOfNodes = VGM.getNumberOfNodes();

     if(Ac_VGM.nrows() != numberOfNodes){

        A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

     } 
     else{

        for(unsigned int i = 0; i < Ac_VGM.nrows(); i++){

            Ac_VGM[i].clear();

        }

     }

     if(b_c.size() != numberOfNodes){

        b_c.resize(numberOfNodes);

     }

     for(unsigned int i = 0; i < b_c.size(); i++){

         b_c[i] = 0.0;

     }

     if(C_v.size() != numberOfNodes) {

        C_v.resize(numberOfNodes);

     }

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     double dt = d_model_p->d_dt;

     while(pointer){

           int indexOfNode = pointer->index;

           int numberOfNeighbors = pointer->neighbors.size();

           std::vector<double> coord = pointer->coord;

           double p_v_k = pointer->p_v;

           double c_v_k = pointer->c_v;

           if(numberOfNeighbors == 1){

              double radius = pointer->radii[0];

              int indexNeighbor = pointer->neighbors[0]->index;

              std::vector<double> coord_neighbor = pointer->neighbors[0]->coord;

              double length = 0.0;

              for(int j = 0; j < 3; j++){

                  length += (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
              }

              length = std::sqrt(length);

              double p_neighbor = pointer->neighbors[0]->p_v;

              double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v_k);

              double volume = length / 2.0 * radius * radius * M_PI;

              if(v_interface > 0.0){

                 Ac_VGM(indexOfNode, indexOfNode) = 1.0;

                 if(p_v_k < input.d_identify_vein_pres){

                    b_c[indexOfNode] = input.d_in_nutrient_vein;

                 }
                 else{

                    b_c[indexOfNode] = input.d_in_nutrient;

                 }

              } 
              else{

                 Ac_VGM(indexOfNode, indexOfNode) = length;

                 Ac_VGM(indexOfNode, indexNeighbor) = dt * v_interface - dt * D_v / length;

                 Ac_VGM(indexOfNode, indexOfNode) = Ac_VGM(indexOfNode, indexOfNode) - dt * v_interface + dt * D_v / length;

                 b_c[indexOfNode] = length * C_v_old[indexOfNode];

              }

           } 
           else{

              for(int i = 0; i < numberOfNeighbors; i++) {

                  double radius = pointer->radii[i];

                  int indexNeighbor = pointer->neighbors[i]->index;

                  Ac_VGM(indexOfNode, indexNeighbor) = 0.0;

                  std::vector<double> coord_neighbor = pointer->neighbors[i]->coord;

                  double length = 0.0;

                  for(int j = 0; j < 3; j++) {

                      length += (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);

                  }

                  length = std::sqrt(length);

                  double p_neighbor = pointer->neighbors[i]->p_v;

                  double v_interface = -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v_k);

                  Ac_VGM(indexOfNode, indexOfNode) += length;

                  if(v_interface > 0.0){

                     Ac_VGM(indexOfNode, indexOfNode) += dt * v_interface;

                  } 
                  else{

                     Ac_VGM(indexOfNode, indexNeighbor) += dt * v_interface;

                  }

                  Ac_VGM(indexOfNode, indexOfNode) += dt * D_v / length;

                  Ac_VGM(indexOfNode, indexNeighbor) -= dt * D_v / length;

                  b_c[indexOfNode] += length * C_v_old[indexOfNode];

                  // Assemble coupling terms (nutrients)
   
                  int N_s = input.d_num_points_length;

                  int N_theta = input.d_num_points_angle;

                  std::vector<double> weights;

                  std::vector<int> id_3D_elements;

                  double length_edge = 0.0;

                  determineWeightsAndIds( N_s, N_theta, N_3D, coord, coord_neighbor, radius, h_3D, length_edge, weights, id_3D_elements);
                  
                  // Surface area of cylinder
                  double surface_area = 2.0*M_PI*0.5*length_edge*radius;

                  // Permeabilty vessel wall for nutrients
                  double L_s = pointer->L_s[i];

		  // 1D part of the coupling
		  Ac_VGM(indexOfNode, indexOfNode) += dt * L_s * surface_area;

                  // Add coupling entry to 3D3D as well as 3D1D and 1D3D matrix
                  int numberOfElements = id_3D_elements.size();

                  for(int j=0;j<numberOfElements;j++){

                      int e_id = id_3D_elements[ j ];

                      double weight = weights[ j ];

                      // get 3d pressure
		      const auto *elem = mesh.elem_ptr(e_id);
		      std::vector<unsigned int> dof_indices_pres;
		      pres_map.dof_indices(elem, dof_indices_pres, v_pres);
		      double p_t_k = pres.current_solution(dof_indices_pres[0]);

		      // get 3d nutrient
		      std::vector<unsigned int> dof_indices_nut;
		      nut_map.dof_indices(elem, dof_indices_nut, v_nut);
		      double c_t_k = nut.current_solution(dof_indices_nut[0]);

                      std::cout << "c_t_k: " << c_t_k << std::endl;

                      b_c[ indexOfNode ] += dt * L_s * surface_area * weight * c_t_k;

		      // term due to pressure difference
		      double c_transport = 0.0;

		      if(p_v_k - p_t_k >= 0.0){

		         c_transport = c_v_k;

                      }
		      else{

		         c_transport = c_t_k;

                      }

		      b_c[ indexOfNode ] -= dt * (1. - osmotic_sigma) * pointer->L_p[i] * surface_area * weight * (p_v_k - p_t_k) * c_transport;

                  }

              }

    }

    pointer = pointer->global_successor;

  }

}

void netfcfvfe::Network::solve3D1DNutrientProblem( int timeStep, double time ) {

     const auto &input = d_model_p->get_input_deck();

     int numberOfNodes = VGM.getNumberOfNodes();

     // Solver
     gmm::iteration iter(4.0e-11);

     std::cout << " " << std::endl;
     std::cout << "Assemble 3D1D nutrient matrix and right hand side" << std::endl;
     assemble3D1DSystemForNutrients();

     phi_sigma = phi_sigma_old;

     gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_nut_3D1D, 60, 1e-8);

     std::cout << " " << std::endl;
     std::cout << "Solve linear system of equations (nutrients)" << std::endl;
     gmm::bicgstab(A_nut_3D1D, phi_sigma, b_nut_3D1D, P, iter);

     auto pointer = VGM.getHead();

     std::cout << " " << std::endl;

     while( pointer ) {

            int indexOfNode = pointer->index;

            if( phi_sigma[  N_tot_3D + indexOfNode ]>1.0 ){

                phi_sigma[  N_tot_3D + indexOfNode ] = 1.0;
               
            }

            pointer->c_v = phi_sigma[  N_tot_3D + indexOfNode ];  
           // std::cout << "index: " << pointer->index << " c_v: " << pointer->c_v << " p_v: " << pointer->p_v << " coord: " << pointer->coord << std::endl;         
            pointer = pointer->global_successor;

     }

     phi_sigma_old = phi_sigma;

     // Extract nutrient concentrations
     std::cout<< " " << std::endl;
     std::cout<< "Extract nutrient concentrations" << std::endl;

     for(int i=0;i<N_tot_3D;i++){

         phi_sigma_3D[ i ] =  phi_sigma[ i ];

     }

     if( timeStep%2 == 0 ){
     
         std::cout << " " << std::endl;
         std::cout<< "Plot solutions" << std::endl;
         writeDataToVTK3D_Nutrients( phi_sigma_3D, N_3D, h_3D, timeStep );
         writeDataToVTKTimeStep_VGM( timeStep );

     }

}

void netfcfvfe::Network::solveVGMforNutrient( int timeStep, double time ) {
    
     const auto &input = d_model_p->get_input_deck();

     int numberOfNodes = VGM.getNumberOfNodes();

     // Solver
     gmm::iteration iter(5.0e-11);
     gmm::identity_matrix PR;

     std::cout << " " << std::endl;
     std::cout << "Assemble nutrient matrix and right hand side" << std::endl;
     assembleVGMSystemForNutrients();

     C_v = C_v_old;

     gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(Ac_VGM, 50, 1e-6);

     gmm::bicgstab(Ac_VGM, C_v, b_c, P, iter);

     auto pointer = VGM.getHead();

     while( pointer ) {

            int indexOfNode = pointer->index;

            pointer->c_v = C_v[ indexOfNode ];

            pointer = pointer->global_successor;

     }

     C_v_old = C_v;

     double dt_output = input.d_dt_output_interval;

     if( time/dt_output-std::floor( time/dt_output ) < 1.0e-8 ){
     
         std::cout << " " << std::endl;
         std::cout<< "Plot solutions" << std::endl;
         writeDataToVTKTimeStep_VGM( timeStep );

     }

}

double netfcfvfe::Network::getDirichletValue( std::vector< double > center_face, double L_p, double radius ){

       double dirichlet_value = 0.0;

       double dist = (center_face[ 0 ]-0.5)*(center_face[ 0 ]-0.5) + (center_face[ 1 ]-0.5)*(center_face[ 1 ]-0.5);

       dist = std::sqrt( dist );

       if(dist<radius){

          dirichlet_value = ( 1.0+center_face[ 2 ] )*L_p/( 1.0+L_p );

       }
       else{

          dirichlet_value = ( 1.0+center_face[ 2 ] )*L_p/( 1.0+L_p )*( 1.0 - radius*std::log( dist/radius ) );

       }

       return dirichlet_value;

}

double netfcfvfe::Network::getK1D( double s, double L_p, double radius ){

       double K_1D = 0.0;

       if( scenario == "test_single_vessel" ){

           K_1D = L_p/(1.0+L_p)*(1.0 + s + 0.5*s*s);

       }
       else{

           K_1D = (radius*radius*radius*radius*M_PI)/(8.0*mu);

       }

       return K_1D;

}

void netfcfvfe::Network::solve3D1DFlowProblem( int timeStep, double time ){

     const auto &input = d_model_p->get_input_deck();

     double L_x = input.d_domain_params[1];

     assemble3D1DSystemForPressure();

     // Solve linear system of equations
     std::cout<< " " << std::endl;
     std::cout<< "Solve linear system of equations (pressure)" << std::endl;

     gmm::iteration iter(10E-16);

     gmm::ilut_precond< gmm::row_matrix< gmm::wsvector<double> > > Pr(A_flow_3D1D, 50, 1e-5);

     P_3D1D = b_flow_3D1D;

     gmm::bicgstab(A_flow_3D1D, P_3D1D, b_flow_3D1D, Pr, iter);

     auto pointer = VGM.getHead();

     while( pointer ) {

            int indexOfNode = pointer->index;

            pointer->p_v = P_3D1D[  N_tot_3D + indexOfNode ];           
            pointer = pointer->global_successor;

     }

     // Extract pressures
     std::cout<< " " << std::endl;
     std::cout<< "Extract pressures" << std::endl;

     for(int i=0;i<N_tot_3D;i++){

         P_3D[ i ] = P_3D1D[ i ];

     }

     // Compute velocity field
     std::cout << "Compute velocity field" << std::endl;
     std::vector< std::vector<double> > V_3D;

     // Space directions
     std::vector<  std::vector< double > > directions;
     directions = defineDirections();

     double K_3D = input.d_tissue_flow_K;

     for(int i=0;i<N_3D;i++){ // x-loop

         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop

                 int index = i + j*N_3D + k*N_3D*N_3D;

                 std::vector<double> vel_element(3);

                 std::vector< double > center = getElementCenter( i, j, k, h_3D );

                 // Iterate over the interfaces
                 for(int face=0;face<6;face++){

                     std::vector< double > center_neighbor = getCenterNeighbor( center, directions[ face ], h_3D );

                     bool isInnerFace = isCenterInDomain( center_neighbor, L_x );

                     if( isInnerFace ){

                         int index_neighbor = getElementIndex( center_neighbor, h_3D, N_3D );

                         for(int l=0;l<3;l++){

                             vel_element[ l ] += -K_3D*( P_3D[ index_neighbor ] - P_3D[ index ] )/h_3D*directions[ face ][l];

                         }

                     }

                 }

                 V_3D.push_back( vel_element );

            }

         }

     }

     if( timeStep%2 == 0 ){
     
         std::cout << " " << std::endl;
         std::cout<< "Plot solutions" << std::endl;
         writeDataToVTK3D_Pressure( P_3D, V_3D, N_3D, h_3D, timeStep );

     }

}

void netfcfvfe::Network::writeDataToVTK3D_Pressure(std::vector<double> P_3D, std::vector< std::vector<double> > V_3D, int N_3D, double h_3D, int timeStep){

     int numberOfCells = P_3D.size();

     int numberOfPolygonData = 9*numberOfCells;

     int numberOfPoints = (N_3D+1)*(N_3D+1)*(N_3D+1);

     std::string path = scenario;
     path += "_Pressure_3D_";
     path += std::to_string(timeStep);
     path.append(".vtk");

     std::fstream filevtk;
     filevtk.open(path, std::ios::out);
     filevtk << "# vtk DataFile Version 2.0" << std::endl;
     filevtk << "Pressure 3D1D coupled problem" << std::endl;
     filevtk << "ASCII" << std::endl;
     filevtk << "DATASET UNSTRUCTURED_GRID" << std::endl;

     filevtk << "POINTS " << numberOfPoints << " float" << std::endl;

     std::vector< std::vector<int> > indices;

     for(int i=0;i<N_3D+1;i++){ // x-loop
         
         for(int j=0;j<N_3D+1;j++){ // y-loop

             for(int k=0;k<N_3D+1;k++){ // z-loop   

                 filevtk << (double) i*h_3D << " " << (double) j*h_3D << " " << (double) k*h_3D << std::endl;               

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop
         
         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop 

                 std::vector< double > center = getElementCenter( i, j, k, h_3D );
               
                 filevtk << "8" << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D ) 
                                << std::endl;

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_TYPES " << numberOfCells << std::endl;

     for(int i=0;i<numberOfCells;i++){

         filevtk << 11 << std::endl;

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_DATA " << numberOfCells << std::endl;
     filevtk << "SCALARS Pressure_(3D) float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop
         
         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop   

                 int index = i + j*N_3D + k*N_3D*N_3D;  

                 filevtk << P_3D[ index ] << std::endl;             

             }

         }

     }

     filevtk << "VECTORS velocity float" << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop
         
         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop   
 
                 int index = i + j*N_3D + k*N_3D*N_3D;  

                 filevtk << V_3D[ index ][ 0 ] << " " << V_3D[ index ][ 1 ] << " " << V_3D[ index ][ 2 ] << std::endl;             

             }

         }

     }

}

void netfcfvfe::Network::writeDataToVTK3D_Nutrients(std::vector<double> phi_sigma_3D, int N_3D, double h_3D, int timeStep){

     int numberOfCells = phi_sigma_3D.size();

     int numberOfPolygonData = 9*numberOfCells;

     int numberOfPoints = (N_3D+1)*(N_3D+1)*(N_3D+1);

     std::string path = scenario;
     path += "_Nutrients_3D_";
     path += std::to_string(timeStep);
     path.append(".vtk");

     std::fstream filevtk;
     filevtk.open(path, std::ios::out);
     filevtk << "# vtk DataFile Version 2.0" << std::endl;
     filevtk << "Pressure 3D1D coupled problem" << std::endl;
     filevtk << "ASCII" << std::endl;
     filevtk << "DATASET UNSTRUCTURED_GRID" << std::endl;

     filevtk << "POINTS " << numberOfPoints << " float" << std::endl;

     std::vector< std::vector<int> > indices;

     for(int i=0;i<N_3D+1;i++){ // x-loop
         
         for(int j=0;j<N_3D+1;j++){ // y-loop

             for(int k=0;k<N_3D+1;k++){ // z-loop   

                 filevtk << (double) i*h_3D << " " << (double) j*h_3D << " " << (double) k*h_3D << std::endl;               

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop
         
         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop 

                 std::vector< double > center = getElementCenter( i, j, k, h_3D );
               
                 filevtk << "8" << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]-0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]-0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D )
                                << " " << getIndex( center[ 0 ]-0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D ) 
                                << " " << getIndex( center[ 0 ]+0.5*h_3D, h_3D )*(N_3D+1)*(N_3D+1) + getIndex( center[ 1 ]+0.5*h_3D, h_3D )*(N_3D+1) + getIndex( center[ 2 ]+0.5*h_3D, h_3D ) 
                                << std::endl;

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_TYPES " << numberOfCells << std::endl;

     for(int i=0;i<numberOfCells;i++){

         filevtk << 11 << std::endl;

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_DATA " << numberOfCells << std::endl;
     filevtk << "SCALARS Phi_sigma_(3D) float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     for(int i=0;i<N_3D;i++){ // x-loop
         
         for(int j=0;j<N_3D;j++){ // y-loop

             for(int k=0;k<N_3D;k++){ // z-loop   

                 int index = i + j*N_3D + k*N_3D*N_3D;  

                 filevtk << phi_sigma_3D[ index ] << std::endl;             

             }

         }

     }

}

