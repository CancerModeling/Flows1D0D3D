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

void netfc::Network::readData( std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
                                std::vector<double> &radii, std::vector<std::vector<unsigned int>> &elements) {

     const auto &input = d_model_p->get_input_deck();

     std::string dgf_filename = input.d_network_init_file;

     std::cout << " " << std::endl;
     std::cout << "Read data from file: " << dgf_filename << std::endl;

     std::ostringstream out_dgf_filename;
     out_dgf_filename << dgf_filename.c_str();

     std::string file_dgf;
     file_dgf = out_dgf_filename.str().c_str();
     std::ifstream in_dgf_filename;
     in_dgf_filename.open(file_dgf.c_str());

     std::string header;
     int fileSeparator = 0, i = 0;

     while(std::getline(in_dgf_filename, header, in_dgf_filename.widen('\n'))) {

            if(header == "# " || header == "#") {

               fileSeparator += 1;
               continue;

            } 
            else if(i >= 3 && fileSeparator == 0) {

              std::vector<double> vertexInfo;

              // Split header
              std::istringstream input(header);
              std::vector<std::array<char, 50>> line;

              for (std::array<char, 50> a; input.getline(&a[0], 50, ' ');)
                  line.push_back(a);

              for (auto &a : line) {
                   vertexInfo.push_back(std::atof(&a[0]));
              }

              std::vector<double> vertex;

              for (int j = 0; j < 3; j++) {
                   vertex.push_back(vertexInfo[j]);
              }

              vertices.push_back(vertex);
              pressures.push_back( vertexInfo[3] );

            }
            else if (fileSeparator == 1) // Avoid line: SIMPLEX in the dgf file
            {
              fileSeparator += 1;
            } 
            else if (fileSeparator == 2) // Avoid line: parameters ... in the dgf file
            {
              fileSeparator += 1;
            } 
            else if (i >= 3 && fileSeparator == 3) {

              std::vector<double> segmentInfo;

              // Split header
              std::istringstream input(header);
              std::vector<std::array<char, 50>> line;

              for (std::array<char, 50> a; input.getline(&a[0], 50, ' ');)
                  line.push_back(a);

              for (auto &a : line) {
                  segmentInfo.push_back(std::atof(&a[0]));
              }

              // Insert segment
              std::vector<unsigned int> cornerIDs(2);
              cornerIDs[0] = (int)segmentInfo[0];
              cornerIDs[1] = (int)segmentInfo[1];
              double radius = segmentInfo[2];

              radii.push_back(radius);
              elements.push_back(cornerIDs);

            }

            i = i + 1;

     }

}


void netfc::Network::transferDataToVGM( std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
                                         std::vector<double> &radii, std::vector<std::vector<unsigned int>> &elements ) {

  int numberOfVertices = vertices.size();

  const auto &input = d_model_p->get_input_deck();

  double L_x = input.d_domain_params[1];

  if( scenario == "network_secomb" ){
  
      rescaleSecombData( vertices, pressures, radii, 0.08 );        

  }

  for(int i = 0; i < numberOfVertices; i++) {

      netfc::VGNode new_node;

      new_node.index = i;

      std::vector<double> coord = vertices[i];

      new_node.coord = coord;

      new_node.p_boundary = pressures[i];

      new_node.p_v = 0.0;

      new_node.c_boundary = input.d_in_nutrient;

      new_node.c_v = 0.0;

    if (new_node.p_boundary > 0.0) {

      new_node.typeOfVGNode = TypeOfNode::DirichletNode;

    } 
    else {

      new_node.typeOfVGNode = TypeOfNode::InnerNode;
    }

    VGM.attachNode(new_node);

  }

  std::shared_ptr<VGNode> pointer_1, pointer_2;

  int numberOfElements = elements.size();

  for(int i = 0; i < numberOfElements; i++){

      int index_1 = elements[i][0];

      int index_2 = elements[i][1];

      std::vector<double> coord_1 = vertices[index_1];

      std::vector<double> coord_2 = vertices[index_2];

      double length = 0.0;

      double radius = radii[i];

      pointer_1 = VGM.findNode(index_1);

      pointer_2 = VGM.findNode(index_2);

      pointer_1->radii.push_back(radius);

      pointer_1->L_p.push_back(input.d_tissue_flow_L_p);

      pointer_1->L_s.push_back(input.d_tissue_nut_L_s);

      pointer_1->neighbors.push_back(pointer_2);

      pointer_1->edge_touched.push_back(false);

      pointer_1->sprouting_edge.push_back(false);

      pointer_2->radii.push_back(radius);

      pointer_2->L_p.push_back(input.d_tissue_flow_L_p);

      pointer_2->L_s.push_back(input.d_tissue_nut_L_s);

      pointer_2->neighbors.push_back(pointer_1);

      pointer_2->edge_touched.push_back(false);

      pointer_2->sprouting_edge.push_back(false);

  }

}

void netfc::Network::printDataVGM() {

     std::cout << " " << std::endl;
     std::cout << "PrintData of network: " << std::endl;
     std::cout << " " << std::endl;

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     while (pointer) {

            std::cout << "Index of node [-]: " << pointer->index << std::endl;
            std::cout << "Type of node [-]: " << pointer->typeOfVGNode << std::endl;
            std::cout << "Boundary pressure of node [-]: " << pointer->p_boundary << std::endl;
            std::cout << "Pressure of node [-]: " << pointer->p_v << std::endl;
            std::cout << "Boundary concentration of node: " << pointer->c_boundary << std::endl;
            std::cout << "Coord: " << pointer->coord[0] << " " << pointer->coord[1] << " " << pointer->coord[2] << std::endl;

            int numberOfNeighbors = pointer->neighbors.size();

            for (int i = 0; i < numberOfNeighbors; i++) {

                 std::cout << "L_p [m/(Pa s)]: " << pointer->L_p[i] << std::endl;
                 std::cout << "L_s [m/s]: " << pointer->L_s[i] << std::endl;
                 std::cout << "Radii [m]: " << pointer->radii[i] << std::endl;
                 std::cout << "Neighbor [-]: " << pointer->neighbors[i]->index << std::endl;
            }

            std::cout << " " << std::endl;

            pointer = pointer->global_successor;

     }

}

void netfc::Network::rescaleSecombData( std::vector< std::vector<double> >& vertices, std::vector<double>& pressures, std::vector<double>& radii, double epsilon ){

     // Get min and max coordinates
     std::cout << " " << std::endl;
     std::cout << "Get min and max coordinates" << std::endl;

     int numberOfVertices = vertices.size();
     int numberOfEdges    = radii.size();

     std::vector<double> vec_min(3);
     std::vector<double> vec_max(3);

     for(int j=0;j<3;j++){

         vec_min[j] = 1.0;
         vec_max[j] = 0.0;

     }

     for(int i=0;i<numberOfVertices;i++){

         for(int j=0;j<3;j++){

             if( vertices[i][j]<vec_min[j] ){

                 vec_min[j] = vertices[i][j];

             }
 
             if( vertices[i][j]>vec_max[j] ){

                 vec_max[j] = vertices[i][j];

             }

         }

         if( pressures[ i ] < 21.0 && pressures[ i ] > 0.0 ){

             pressures[ i ] = 0.1;

         }
         else if( pressures[ i ] > 21.0 ){

             pressures[ i ] = 9.5;
   
         }

     }

     for(int j=0;j<3;j++){

         std::cout << "min: " << vec_min[j] << " max: " << vec_max[j] << std::endl;

     }

     // Get min and max coordinates
     std::cout << " " << std::endl;
     std::cout << "Get min and max coordinates" << std::endl;

     for(int i=0;i<numberOfVertices;i++){

         for(int j=0;j<3;j++){

             double new_value = 1.9*(vertices[i][j] - vec_min[j])/( vec_max[j]- vec_min[j] ) + 0.025;

             std::cout << "Value: " << new_value << std::endl;

             vertices[i][j] = new_value;

         }

     }

     std::cout << " " << std::endl;

     for(int i=0;i<numberOfEdges;i++){

         std::cout << "Radius old: " << radii[ i ] << " radius new: " << radii[ i ]/1.2e-4 << std::endl;

         radii[ i ] = radii[ i ]/1.2e-4;

     }

}

void netfc::Network::refine1DMesh() {

     const auto &input = d_model_p->get_input_deck();

     std::shared_ptr<VGNode> pointer = VGM.getHead();

     std::vector<std::shared_ptr<VGNode>> newNodes;

     int numberOfNodes = VGM.getNumberOfNodes();

     int counter = 0;

     while (pointer) {

            int numberOfNeighbors = pointer->neighbors.size();

            for (int i = 0; i < numberOfNeighbors; i++) {

                 if (!pointer->edge_touched[i]) {

                      netfc::VGNode new_node;

                      new_node.index = numberOfNodes + counter;

                      new_node.p_boundary = 0.0;

                      new_node.p_v = 0.0;

                      new_node.c_boundary = input.d_in_nutrient;

                      new_node.c_v = 0.0;

                      new_node.L_p.push_back(input.d_tissue_flow_L_p);

                      new_node.L_p.push_back(input.d_tissue_flow_L_p);

                      new_node.L_s.push_back(input.d_tissue_nut_L_s);

                      new_node.L_s.push_back(input.d_tissue_nut_L_s);

                      new_node.radii.push_back(pointer->radii[i]);

                      new_node.radii.push_back(pointer->radii[i]);

                      new_node.typeOfVGNode = TypeOfNode::InnerNode;

                      std::shared_ptr<VGNode> pointer_new_node = std::make_shared<VGNode>(new_node);

                      std::shared_ptr<VGNode> pointer_former_neighbor = pointer->neighbors[i];

                      std::vector<double> coord_pointer = pointer->coord;

                      std::vector<double> coord_former_neighbor = pointer_former_neighbor->coord;

                      std::vector<double> new_coord;

                      for (int j = 0; j < 3; j++) {

                           new_coord.push_back( 0.5 * (coord_pointer[j] + coord_former_neighbor[j]) );
        
                      }

                      pointer_new_node->coord = new_coord;

                      pointer->neighbors[i] = pointer_new_node;

                      pointer_new_node->neighbors.push_back(pointer);

                      pointer_new_node->neighbors.push_back(pointer_former_neighbor);

                      pointer_former_neighbor->replacePointerWithIndex( pointer->index, pointer_new_node );

                      pointer->edge_touched[i] = true;

                      pointer_new_node->edge_touched.push_back(true);

                      pointer_new_node->edge_touched.push_back(true);

                      pointer_new_node->sprouting_edge.push_back(false);

                      pointer_new_node->sprouting_edge.push_back(false);

                      newNodes.push_back(pointer_new_node);

                      counter = counter + 1;

                 }

            }

            pointer = pointer->global_successor;

     }

     for (int i = 0; i < newNodes.size(); i++) {

          VGM.attachPointerToNode(newNodes[i]);

     }

     pointer = VGM.getHead();

     counter = 0;

     while (pointer) {

            int numberOfNeighbors = pointer->neighbors.size();

             for (int i = 0; i < numberOfNeighbors; i++) {

                  pointer->edge_touched[i] = false;

             }

             pointer->index = counter;

             counter = counter + 1;

             pointer = pointer->global_successor;

     }

}


void netfc::Network::writeDataToVTK_3D(std::vector<double> P_3D, int N_3D, double h_3D) {

     int numberOfCells = P_3D.size();

     int numberOfPolygonData = 9 * numberOfCells;

     int numberOfPoints = (N_3D + 1) * (N_3D + 1) * (N_3D + 1);

     std::fstream filevtk;
     filevtk.open("Pressure_3D.vtk", std::ios::out);
     filevtk << "# vtk DataFile Version 2.0" << std::endl;
     filevtk << "Pressure 3D1D coupled problem" << std::endl;
     filevtk << "ASCII" << std::endl;
     filevtk << "DATASET UNSTRUCTURED_GRID" << std::endl;

     filevtk << "POINTS " << numberOfPoints << " float" << std::endl;

     std::vector<std::vector<int>> indices;

     for(int i = 0; i < N_3D + 1; i++) { // x-loop

         for(int j = 0; j < N_3D + 1; j++) { // y-loop

             for(int k = 0; k < N_3D + 1; k++) { // z-loop

                 filevtk << (double)i * h_3D << " " << (double)j * h_3D << " " << (double)k * h_3D << std::endl;

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData << std::endl;

     for(int i = 0; i < N_3D; i++) { // x-loop

         for(int j = 0; j < N_3D; j++) { // y-loop

             for(int k = 0; k < N_3D; k++) { // z-loop

                 std::vector<double> center = getElementCenter(i, j, k, h_3D);

                 filevtk << "8"  << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) + getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] - 0.5 * h_3D, h_3D) << " "  << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] - 0.5 * h_3D, h_3D) << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] - 0.5 * h_3D, h_3D) << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] - 0.5 * h_3D, h_3D) << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] + 0.5 * h_3D, h_3D) << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] + 0.5 * h_3D, h_3D) << " " << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                       getIndex(center[2] + 0.5 * h_3D, h_3D) << " " << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) * (N_3D + 1) +
                       getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) + getIndex(center[2] + 0.5 * h_3D, h_3D)  << std::endl;

             }

         }

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_TYPES " << numberOfCells << std::endl;

     for(int i = 0; i < numberOfCells; i++) {

         filevtk << 11 << std::endl;

     }

     filevtk << " " << std::endl;
     filevtk << "CELL_DATA " << numberOfCells << std::endl;
     filevtk << "SCALARS pressures(3D) float 1" << std::endl;
     filevtk << "LOOKUP_TABLE default" << std::endl;

     for(int i = 0; i < N_3D; i++) { // x-loop

         for(int j = 0; j < N_3D; j++) { // y-loop

             for(int k = 0; k < N_3D; k++) { // z-loop

                 int index = i + j * N_3D + k * N_3D * N_3D;

                 filevtk << P_3D[index] << std::endl;

             }

         }

     }

}


