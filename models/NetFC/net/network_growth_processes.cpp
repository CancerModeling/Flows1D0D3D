////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "netUtil.hpp"
#include "../model.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

void netfc::Network::updateNetwork(){

     std::cout << " " << std::endl;
     std::cout << "Update the network" << std::endl;

     int numberOfNodes = VGM.getNumberOfNodes();

     std::cout << " " << std::endl;
     std::cout << "Number of nodes: " << numberOfNodes << std::endl;

     std::cout << " " << std::endl;
     std::cout << "Mark nodes for growth " << std::endl;
     markApicalGrowth();

     std::cout << " " << std::endl;
     std::cout << "Process apical growth" << std::endl;
     int numberOfNewNodes = 0;
     numberOfNewNodes = processApicalGrowth();

     numberOfNodes = VGM.getNumberOfNodes();

     std::cout << " " << std::endl;
     std::cout << "Number of nodes after growing the network: " << numberOfNodes << std::endl;

     std::cout << " " << std::endl;
     std::cout << "Link terminal vessels" << std::endl;
     linkTerminalVessels();

     std::cout << " " << std::endl;
     std::cout << "Number of nodes after linking terminal vessels to the network: " << numberOfNodes << std::endl;

     std::cout << " " << std::endl;
     std::cout << "Rescale the 1D matrices and vectors" << std::endl;

     Ac_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
     b_c = std::vector<double>(numberOfNodes, 0.0);

     b_flow_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     A_nut_3D1D = gmm::row_matrix<gmm::wsvector<double>>(N_tot_3D+numberOfNodes, N_tot_3D+numberOfNodes);
     b_nut_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     phi_sigma = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);
     phi_sigma_old = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     P_3D1D = std::vector<double>(N_tot_3D+numberOfNodes, 0.0);

     for(int i=0;i<N_tot_3D;i++){

         phi_sigma[ i ]     = phi_sigma_3D[ i ];
         phi_sigma_old[ i ] = phi_sigma_3D[ i ];
         P_3D1D[ i ]        = P_3D[ i ];

     }

     auto pointer = VGM.getHead();

     while( pointer ) {

            int indexOfNode = pointer->index;

            phi_sigma[ N_tot_3D+indexOfNode ]     = pointer->c_v;
            phi_sigma_old[ N_tot_3D+indexOfNode ] = pointer->c_v;
            P_3D1D[ N_tot_3D+indexOfNode ]        = pointer->p_v;
        
            pointer = pointer->global_successor;

     }

}

void netfc::Network::linkTerminalVessels(){

     const auto &input = d_model_p->get_input_deck();

     double L_x = input.d_domain_params[1];

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     while( pointer ) {

            if( pointer->typeOfVGNode == TypeOfNode::DirichletNode ) {

                std::vector<double> coord = pointer->coord;

                std::shared_ptr<VGNode> pointer_1 = VGM.getHead();  

                double radius = pointer->radii[ 0 ];    

                int index = pointer->index;   

                std::cout << "index: " << index << "\n";

                std::vector<double> dir_term_vessel = std::vector<double>(3,0.0);  
   
                for(int i=0;i<3;i++){

                    dir_term_vessel[ i ] = coord[ i ] - pointer->neighbors[ 0 ]->coord[ i ];
          
                } 

                double length_dir = gmm::vect_norm2( dir_term_vessel );

                std::vector<double> normal_plane = std::vector<double>(3,0.0); 

                for(int i=0;i<3;i++){

                    normal_plane[ i ] = dir_term_vessel[ i ]/length_dir;
                   
                }

                while( pointer_1 ) {
/*
                       std::vector<double> coord_1 = pointer_1->coord;

                       int index_1 = pointer->index; 

                       std::cout << "index_1: " << index_1 << "\n";

                       double dist = util::dist_between_points( coord, coord_1 );

                       std::vector<double> diff, normal_plane = std::vector<double>(3,0.0);

                       for(int i=0;i<3;i++){

                           diff[ i ] = coord_1[ i ] - coord[ i ];
          
                       } 

                       double length_diff = gmm::vect_norm2( diff );
                       double dist_plane  = 0.0;

                       for(int i=0;i<3;i++){

                           dist_plane += normal_plane[ i ]*diff[ i ];

                       }

                       std::cout << "dist_plane: " << dist_plane << "\n";

                       if( dist_plane>0.02 && index != index_1 && dist < 3.0*h_3D ){

                           std::cout << " " << std::endl;
                           std::cout << "dist: " << dist << "\n";
                           std::cout << "index: " << index << "\n";
                           std::cout << "index_1: " << index_1 << "\n";
                           std::cout << "Update pointer" << "\n";
/*
                           pointer->p_boundary = 0.0;
                           pointer->c_boundary = 0.0;
                           pointer->typeOfVGNode = TypeOfNode::InnerNode;
                           pointer->apicalGrowth = false;
                           pointer->radii.push_back( radius );
                           pointer->L_p.push_back(input.d_tissue_flow_L_p);
                           pointer->L_s.push_back(input.d_tissue_nut_L_s);
                           pointer->edge_touched.push_back( true );
                           pointer->sprouting_edge.push_back( false );
                           pointer->neighbors.push_back( pointer_1 );

                           std::cout << "Update pointer 1" << "\n";

                           pointer_1->p_boundary = 0.0;
                           pointer_1->c_boundary = 0.0;
                           pointer_1->typeOfVGNode = TypeOfNode::InnerNode;
                           pointer_1->apicalGrowth = false;
                           pointer_1->radii.push_back( radius );
                           pointer_1->L_p.push_back(input.d_tissue_flow_L_p);
                           pointer_1->L_s.push_back(input.d_tissue_nut_L_s);
                           pointer_1->edge_touched.push_back( true );
                           pointer_1->sprouting_edge.push_back( false );
                           pointer_1->neighbors.push_back( pointer );                                                     

                       }
*/
                       pointer_1 = pointer_1->global_successor;

                }

            }

            pointer = pointer->global_successor;

     }

}

void netfc::Network::markApicalGrowth(){

     const auto &input = d_model_p->get_input_deck();

     double L_x = input.d_domain_params[1];

     std::shared_ptr<VGNode> pointer = VGM.getHead();
 
     std::cout << " " << std::endl;

     while( pointer ) {

            if( pointer->typeOfVGNode == TypeOfNode::DirichletNode ) {

                const auto &coord = pointer->coord;

                // if node is near the boundary, we do not process Omega = (0,L)^3 
                if( 0.01<coord[0] && coord[0]<L_x-0.01 && 0.01<coord[1] && coord[1]<L_x-0.01 && 0.01<coord[2] && coord[2]<L_x-0.01 ){

                    int index = getElementIndex( coord, h_3D, N_3D );

                    double taf_node = phi_TAF[ index ];

                    std::cout << "taf_node: " << taf_node << std::endl;

                    if( taf_node > input.d_network_update_taf_threshold ){

                        std::cout << "index: " << pointer->index << std::endl;
                          
                        pointer->apicalGrowth = true;
                     
                    }

                }
 
            }

            pointer = pointer->global_successor;

     }

}

int netfc::Network::processApicalGrowth(){

    const auto &input = d_model_p->get_input_deck();

    double L_x = input.d_domain_params[1];

    std::cout << "L_x: " << L_x << "\n";

    // Initialize random objects
    std::lognormal_distribution<> log_normal_distribution( input.d_log_normal_mean, input.d_log_normal_std_dev );
    std::random_device rd;
    std::mt19937 generator(rd());

    int numberOfNodes_old = VGM.getNumberOfNodes();

    // mark node for growth based on a certain criterion
    std::shared_ptr<VGNode> pointer = VGM.getHead();

    std::cout << "Number of nodes before: " << VGM.getNumberOfNodes() << "\n";

    while( pointer ){

           if( pointer->apicalGrowth ){

               std::cout << " " << "\n";
               std::cout << "Processing node: " << pointer->index << "\n";

               std::vector<double> coord = pointer->coord;
               std::cout << "coord: " << coord << "\n";
               std::cout << "Compute direction based on TAF" << "\n";

               int element_index = getElementIndex( coord, h_3D, N_3D );
               std::cout << "element_index: " << element_index << "\n";

               std::vector<double> node_center = getCenterFromIndex( element_index, N_3D, h_3D );
               std::cout << "node_center: " << node_center << "\n";

               std::vector<int> indicesNeighbors = getNeighboringElementIndices( coord, h_3D, L_x, N_3D );

               double TAF_node = phi_TAF[ element_index ];

               std::vector<double> TAF_neighbors;

               for(int i=0;i<indicesNeighbors.size();i++){

                   TAF_neighbors.push_back( phi_TAF[ indicesNeighbors[ i ] ] );

               }

               int max_index, max_index_2 = 0;

               double TAF_max, TAF_max_2 = 0.0;

               for(int i=0;i<indicesNeighbors.size();i++){

                   double TAF = std::abs( TAF_neighbors[ i ] );

                   //std::cout << "TAF: " << TAF << " TAF_max: " << TAF_max << " " << std::endl;

                   if( TAF > TAF_max-1.0e-8 ){

                       max_index = indicesNeighbors[ i ];

                       TAF_max = TAF;

                   }
                   else if( TAF_max-1.0e-8 > TAF && TAF > TAF_max_2-1.0e-8 ){

                       max_index_2 = indicesNeighbors[ i ];

                       TAF_max_2 = TAF;

                   }

               }

               std::cout << "max_index: " << max_index << "\n";
               std::cout << "max_index_2: " << max_index_2 << "\n";
               std::cout << "TAF_neighbors: " << TAF_neighbors << "\n";
               std::cout << "TAF_node: " << TAF_node << "\n";

               std::vector<double> new_point_1  = std::vector<double>(3,0.0);
               std::vector<double> max_center   = getCenterFromIndex( max_index, N_3D, h_3D );
               std::vector<double> max_center_2 = getCenterFromIndex( max_index_2, N_3D, h_3D );
               std::vector<double> diff         = std::vector<double>(3,0.0);

               std::vector<double> dir_term_vessel = std::vector<double>(3,0.0);
               std::vector<double> normal_plane    = std::vector<double>(3,0.0);

               for(int i=0;i<3;i++){

                   diff[ i ] = 0.5*( max_center[ i ] + max_center_2[ i ] ) - coord[ i ];
                   dir_term_vessel[ i ] = coord[ i ] - pointer->neighbors[ 0 ]->coord[ i ];
          
               } 

               double length_diff = gmm::vect_norm2( diff );
               double length_dir  = gmm::vect_norm2( dir_term_vessel );

               double prod_coord, dist_new_point = 0.0;

               for(int i=0;i<3;i++){

                   normal_plane[ i ] = dir_term_vessel[ i ]/length_dir;
                   
               }

               // lognormal distribution
               double log_dist = log_normal_distribution(generator);
               double radius_p = pointer->radii[ 0 ];

               // get length
               double length = log_dist * radius_p;

               if( length>3.0*h_3D ){

                   length = 3.0*h_3D;

               }

               for(int i=0;i<3;i++){

                   new_point_1[ i ] = coord[ i ] + ( length*diff[ i ]/length_diff );
                   dist_new_point += normal_plane[ i ]*(new_point_1[ i ]-coord[ i ]);

               }
     
               std::cout << "normal_plane: " << normal_plane << "\n";
               std::cout << "length: " << length << "\n";
               std::cout << "dir_term_vessel: " << dir_term_vessel << "\n";
               std::cout << "diff: " << diff << "\n";
               std::cout << "length_diff: " << length_diff << "\n";
               std::cout << "new_point: " << new_point_1 << "\n";
               std::cout << "dist_new_point: " << dist_new_point << "\n";

               // check if we bifurcate at this node
               bool bifurcate = false;

               double prob =  0.5 + 0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) / std::sqrt(2. * input.d_log_normal_std_dev * input.d_log_normal_std_dev));

               if( prob > 0.85 ){

                   bifurcate = true;

               }

               std::cout << "prob: " << prob << "\n";

               if( !bifurcate ){

                   if( dist_new_point>0.0 ){

                       createASingleNode( new_point_1, radius_p, pointer );

                   }

               }
               else{

		   std::cout << "Create bifuraction" << "\n";

                   double gamma = input.d_net_radius_exponent_gamma;
                   double R_c = std::pow(2.0, -1.0/gamma) * radius_p;

                   if( R_c > radius_p ){

                       R_c = radius_p;

                   }

                   // create normal distribution function
                   std::normal_distribution<> normal_distribution( R_c, R_c/35.0 );

                   double radius_b1 = normal_distribution(generator);
                   double radius_b2 = normal_distribution(generator);

                   if( radius_b1 < 5.0e-3 ){

                       radius_b1 = 5.0e-3;

                   }

                   if( radius_b2 < 5.0e-3 ){

                       radius_b2 = 5.0e-3;

                   }

                   std::vector<double> new_point_2  = std::vector<double>(3,0.0);

                   double branch_angle_1 = 0.0;
                   double branch_angle_2 = 0.0;

                   double angle_arg_1 = ( ( radius_p*radius_p*radius_p*radius_p ) + ( radius_b1*radius_b1*radius_b1*radius_b1 ) - ( radius_b2*radius_b2*radius_b2*radius_b2 ) )/ 
                                    ( 2.0*radius_p*radius_p*radius_b1*radius_b1 );

                   double angle_arg_2 = ( ( radius_p*radius_p*radius_p*radius_p ) + ( radius_b2*radius_b2*radius_b2*radius_b2 ) - ( radius_b1*radius_b1*radius_b1*radius_b1 ) )/            
                                    ( 2.0*radius_p*radius_p*radius_b2*radius_b2 );

                   if( std::abs( angle_arg_1 )<1.0 && std::abs( angle_arg_2 )<1.0 ){

                       branch_angle_1 = std::acos( angle_arg_1 );
                       branch_angle_2 = std::acos( angle_arg_2 );

                       std::cout << "radius_p: "  << radius_p << "\n";
                       std::cout << "radius_b1: " << radius_b1 << "\n";
                       std::cout << "radius_b2: " << radius_b2 << "\n";
                       std::cout << "branch_angle_1: " << branch_angle_1*180.0/M_PI << "\n";
                       std::cout << "branch_angle_2: " << branch_angle_2*180.0/M_PI << "\n";

                       double branch_angle = branch_angle_1 + branch_angle_2;

                       std::cout << "branch_angle: " << branch_angle*180.0/M_PI << "\n";

                       if( branch_angle*180.0/M_PI<100.0 && branch_angle*180.0/M_PI>20.0 ){

                           std::vector<double> rotation_axis = util::cross_prod( diff, dir_term_vessel );
                           std::vector<double> diff_2  = util::rotate( diff, branch_angle, rotation_axis );

                           double length_diff_1 = gmm::vect_norm2( diff   );
                           double length_diff_2 = gmm::vect_norm2( diff_2 );

                           // lognormal distribution
                           double log_dist = log_normal_distribution(generator);

                           // get length
                           double length_1 = log_dist * radius_b1;
                           double length_2 = log_dist * radius_b2;

                           if( length_1>3.0*h_3D ){

                               length_1 = 3.0*h_3D;

                           }

                           if( length_2>3.0*h_3D ){

                               length_2 = 3.0*h_3D;

                           }

                           for(int i=0;i<3;i++){

                               new_point_1[ i ] = coord[ i ] + ( length_1*diff[ i ]/length_diff_1 );
                               new_point_2[ i ] = coord[ i ] + ( length_2*diff_2[ i ]/length_diff_2 );
          
                           }

                           std::cout << "branch_angle: " << branch_angle << "\n";
                           std::cout << "length_1: " << length_1 << "\n";
                           std::cout << "length_2: " << length_2 << "\n";
                           std::cout << "new_point_1: " << new_point_1 << "\n";
                           std::cout << "new_point_2: " << new_point_2 << "\n";

                           createASingleNode( new_point_1, radius_b1, pointer );
                           createASingleNode( new_point_2, radius_b2, pointer );

                       }

                   }

               }

           }

           pointer = pointer->global_successor;

         }

         // reset the boolean values
         pointer = VGM.getHead();

         while( pointer ){

                int numberOfNeighbors = pointer->neighbors.size();

                pointer->apicalGrowth = false;

                for( int i=0; i<numberOfNeighbors; i++){

                     pointer->sprouting_edge[i] = false;
                     pointer->edge_touched[i] = false;

                }

                pointer = pointer->global_successor;
 
         }

         int numberOfNodes_new = VGM.getNumberOfNodes();

         std::cout << "number of nodes after growth: " << numberOfNodes_new << "\n";

         int numberOfNodes_added = numberOfNodes_new - numberOfNodes_old;

         return numberOfNodes_added;

}

void netfc::Network::createASingleNode( std::vector<double> new_point, double radius, std::shared_ptr<VGNode>& pointer ){

     const auto &input = d_model_p->get_input_deck();

     std::cout << "Create new node" << "\n";

     VGNode new_node;

     new_node.index = VGM.getNumberOfNodes();
     new_node.coord = new_point;
     new_node.p_boundary = 0.0;
     new_node.p_v = pointer->p_v;
     new_node.c_boundary = input.d_in_nutrient;
     new_node.c_v = 0.0;
     new_node.typeOfVGNode = TypeOfNode::DirichletNode;
     new_node.apicalGrowth = false;
     new_node.radii.push_back( radius );
     new_node.L_p.push_back( input.d_tissue_flow_L_p );
     new_node.L_s.push_back( input.d_tissue_nut_L_s );
     new_node.edge_touched.push_back( true );
     new_node.sprouting_edge.push_back( false );
     new_node.neighbors.push_back( pointer );

     auto sp_newNode = std::make_shared<VGNode>( new_node );
     std::cout << "New index: " << new_node.index << "\n";
     std::cout << "Neighbor index: " << pointer->index << "\n";

     std::cout << "Update old node" << "\n";

     pointer->p_boundary = 0.0;
     pointer->c_boundary = 0.0;
     pointer->typeOfVGNode = TypeOfNode::InnerNode;
     pointer->apicalGrowth = false;
     pointer->radii.push_back( radius );
     pointer->L_p.push_back(input.d_tissue_flow_L_p);
     pointer->L_s.push_back(input.d_tissue_nut_L_s);
     pointer->edge_touched.push_back( true );
     pointer->sprouting_edge.push_back( false );
     pointer->neighbors.push_back( sp_newNode );
     std::cout << "Attach new node as pointer" << "\n";

     VGM.attachPointerToNode( sp_newNode );

}


