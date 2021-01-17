////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "netUtil.hpp"
#include "network.hpp"

void util::unet::Network::readData(
  std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
  std::vector<double> &radii,
  std::vector<std::vector<unsigned int>> &elements) {

  const auto &input = d_model_p->get_input_deck();

  std::string dgf_filename = input.d_network_init_file;

  std::ostringstream out_dgf_filename;
  out_dgf_filename << dgf_filename.c_str();

  std::string file_dgf;
  file_dgf = out_dgf_filename.str().c_str();
  std::ifstream in_dgf_filename;
  in_dgf_filename.open(file_dgf.c_str());

  std::string header;
  int fileSeparator = 0, i = 0;

  while (std::getline(in_dgf_filename, header, in_dgf_filename.widen('\n'))) {

    if (header == "# " || header == "#") {

      fileSeparator += 1;
      continue;

    } else if (i >= 3 && fileSeparator == 0) {

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
      pressures.push_back(vertexInfo[3]);

    } else if (fileSeparator == 1) // Avoid line: SIMPLEX in the dgf file
    {
      fileSeparator += 1;
    } else if (fileSeparator == 2) // Avoid line: parameters ... in the dgf file
    {
      fileSeparator += 1;
    } else if (i >= 3 && fileSeparator == 3) {

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
      cornerIDs[0] = (int) segmentInfo[0];
      cornerIDs[1] = (int) segmentInfo[1];
      double radius = segmentInfo[2];

      radii.push_back(radius);
      elements.push_back(cornerIDs);
    }

    i = i + 1;
  }
}

void util::unet::Network::transferDataToVGM(std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
                                            std::vector<double> &radii, std::vector<std::vector<unsigned int>> &elements) {

  int numberOfVertices = vertices.size();
  const auto &input = d_model_p->get_input_deck();

  if (scenario == "network_secomb") {

    rescaleSecombData(vertices, pressures, radii, 0.08);
  }

  for (int i = 0; i < numberOfVertices; i++) {

    VGNode new_node;

    new_node.index = i;

    std::vector<double> coord = vertices[i];

    new_node.coord = coord;
    new_node.p_boundary = pressures[i];
    new_node.p_v = 0.0;
    new_node.c_boundary = input.d_in_nutrient;
    new_node.c_v = 0.0;
    new_node.is_initial_node = true;

    if (new_node.p_boundary > 0.0) {

      new_node.typeOfVGNode = TypeOfNode::DirichletNode;

    } else {

      new_node.typeOfVGNode = TypeOfNode::InnerNode;
    }

    VGM.attachNode(new_node);
  }

  std::shared_ptr<VGNode> pointer_1, pointer_2;

  int numberOfElements = elements.size();

  for (int i = 0; i < numberOfElements; i++) {

    int index_1 = elements[i][0];
    int index_2 = elements[i][1];

    std::vector<double> coord_1 = vertices[index_1];
    std::vector<double> coord_2 = vertices[index_2];

    double radius = radii[i];

    pointer_1 = VGM.findNode(index_1);
    pointer_2 = VGM.findNode(index_2);

    pointer_1->radii.push_back(radius);
    pointer_1->radii_initial.push_back(radius);
    pointer_1->L_p.push_back(input.d_tissue_flow_L_p);
    pointer_1->L_s.push_back(input.d_tissue_nut_L_s);
    pointer_1->neighbors.push_back(pointer_2);
    pointer_1->edge_touched.push_back(false);
    pointer_1->sprouting_edge.push_back(false);
    pointer_1->tau_w_initial.push_back(0.);
    pointer_1->is_given = true;

    pointer_2->radii.push_back(radius);
    pointer_2->radii_initial.push_back(radius);
    pointer_2->L_p.push_back(input.d_tissue_flow_L_p);
    pointer_2->L_s.push_back(input.d_tissue_nut_L_s);
    pointer_2->neighbors.push_back(pointer_1);
    pointer_2->edge_touched.push_back(false);
    pointer_2->sprouting_edge.push_back(false);
    pointer_2->tau_w_initial.push_back(0.);
    pointer_2->is_given = true;
  }
}

void util::unet::Network::printDataVGM() {

  std::cout << " " << std::endl;
  std::cout << "PrintData of network: " << std::endl;
  std::cout << " " << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    std::cout << "Index of node [-]: " << pointer->index << std::endl;
    std::cout << "Type of node [-]: " << pointer->typeOfVGNode << std::endl;
    std::cout << "Boundary pressure of node [-]: " << pointer->p_boundary
              << std::endl;
    std::cout << "Pressure of node [-]: " << pointer->p_v << std::endl;
    std::cout << "Boundary concentration of node: " << pointer->c_boundary
              << std::endl;
    std::cout << "Coord: " << pointer->coord[0] << " " << pointer->coord[1]
              << " " << pointer->coord[2] << std::endl;

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      std::cout << "L_p [m/(Pa s)]: " << pointer->L_p[i] << std::endl;
      std::cout << "L_s [m/s]: " << pointer->L_s[i] << std::endl;
      std::cout << "Radii [m]: " << pointer->radii[i] << std::endl;
      std::cout << "Neighbor [-]: " << pointer->neighbors[i]->index
                << std::endl;
    }

    std::cout << " " << std::endl;

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::rescaleSecombData(
  std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
  std::vector<double> &radii, double epsilon) {

  // Get min and max coordinates
  //std::cout << " " << std::endl;
  //std::cout << "Get min and max coordinates" << std::endl;

  int numberOfVertices = vertices.size();
  int numberOfEdges = radii.size();

  std::vector<double> vec_min(3);
  std::vector<double> vec_max(3);

  for (int j = 0; j < 3; j++) {

    vec_min[j] = 1.0;
    vec_max[j] = 0.0;
  }

  for (int i = 0; i < numberOfVertices; i++) {

    for (int j = 0; j < 3; j++) {

      if (vertices[i][j] < vec_min[j]) {

        vec_min[j] = vertices[i][j];
      }

      if (vertices[i][j] > vec_max[j]) {

        vec_max[j] = vertices[i][j];
      }
    }

    // TODO need to fix this
    if (pressures[i] < 21.0 && pressures[i] > 0.0) {

      pressures[i] = 1000.0;

    } else if (pressures[i] > 21.0) {

      pressures[i] = 8000.0;
    }
  }

  for (int j = 0; j < 3; j++) {

    oss << "min: " << vec_min[j] << " max: " << vec_max[j] << std::endl;
  }

  // Get min and max coordinates
  //std::cout << " " << std::endl;
  //std::cout << "Get min and max coordinates" << std::endl;

  for (int i = 0; i < numberOfVertices; i++) {

    for (int j = 0; j < 3; j++) {

      double new_value =
        1.81 * (vertices[i][j] - vec_min[j]) / (vec_max[j] - vec_min[j]) +
        0.066;

      //std::cout << "Value: " << new_value << std::endl;

      vertices[i][j] = new_value;
    }
  }

  //std::cout << " " << std::endl;

  for (int i = 0; i < numberOfEdges; i++) {

    //    std::cout << "Radius old: " << radii[i]
    //              << " radius new: " << radii[i] * 0.28 << std::endl;

    radii[i] = radii[i] * 0.28;
  }
}

void util::unet::Network::refine1DMesh() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::vector<std::shared_ptr<VGNode>> newNodes;

  int numberOfNodes = VGM.getNumberOfNodes();

  int counter = 0;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        VGNode new_node;

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

        new_node.radii_initial.push_back(pointer->radii[i]);

        new_node.radii_initial.push_back(pointer->radii[i]);

        // missed tau_w_initial?
        new_node.tau_w_initial.push_back(pointer->tau_w_initial[i]);
        new_node.tau_w_initial.push_back(pointer->tau_w_initial[i]);

        new_node.typeOfVGNode = TypeOfNode::InnerNode;

        new_node.is_given = true;

        std::shared_ptr<VGNode> pointer_new_node =
          std::make_shared<VGNode>(new_node);

        std::shared_ptr<VGNode> pointer_former_neighbor = pointer->neighbors[i];

        std::vector<double> coord_pointer = pointer->coord;

        std::vector<double> coord_former_neighbor =
          pointer_former_neighbor->coord;

        std::vector<double> new_coord;

        for (int j = 0; j < 3; j++) {

          new_coord.push_back(0.5 *
                              (coord_pointer[j] + coord_former_neighbor[j]));
        }

        pointer_new_node->coord = new_coord;

        pointer->neighbors[i] = pointer_new_node;

        pointer_new_node->neighbors.push_back(pointer);

        pointer_new_node->neighbors.push_back(pointer_former_neighbor);

        pointer_former_neighbor->replacePointerWithGlobalIndex(pointer->index, pointer_new_node);

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

void util::unet::Network::writeDataToVTK_3D(std::vector<double> P_3D,
                                            int N_3D, double h_3D) {

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

  for (int i = 0; i < N_3D + 1; i++) { // x-loop

    for (int j = 0; j < N_3D + 1; j++) { // y-loop

      for (int k = 0; k < N_3D + 1; k++) { // z-loop

        filevtk << (double) i * h_3D << " " << (double) j * h_3D << " "
                << (double) k * h_3D << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData
          << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        filevtk << "8"
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_TYPES " << numberOfCells << std::endl;

  for (int i = 0; i < numberOfCells; i++) {

    filevtk << 11 << std::endl;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfCells << std::endl;
  filevtk << "SCALARS pressures(3D) float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        filevtk << P_3D[index] << std::endl;
      }
    }
  }
}

void util::unet::Network::writeDataToVTK3D_Pressure(
  std::vector<double> P_3D, std::vector<std::vector<double>> V_3D, int N_3D,
  double h_3D, int timeStep) {

  int numberOfCells = P_3D.size();

  int numberOfPolygonData = 9 * numberOfCells;

  int numberOfPoints = (N_3D + 1) * (N_3D + 1) * (N_3D + 1);

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

  std::vector<std::vector<int>> indices;

  for (int i = 0; i < N_3D + 1; i++) { // x-loop

    for (int j = 0; j < N_3D + 1; j++) { // y-loop

      for (int k = 0; k < N_3D + 1; k++) { // z-loop

        filevtk << (double) i * h_3D << " " << (double) j * h_3D << " "
                << (double) k * h_3D << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData
          << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        filevtk << "8"
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_TYPES " << numberOfCells << std::endl;

  for (int i = 0; i < numberOfCells; i++) {

    filevtk << 11 << std::endl;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfCells << std::endl;
  filevtk << "SCALARS Pressure_(3D) float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        filevtk << P_3D[index] << std::endl;
      }
    }
  }

  filevtk << "VECTORS velocity float" << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        filevtk << V_3D[index][0] << " " << V_3D[index][1] << " "
                << V_3D[index][2] << std::endl;
      }
    }
  }
}

void util::unet::Network::writeDataToVTK3D_Nutrients(
  std::vector<double> phi_sigma_3D, int N_3D, double h_3D, int timeStep) {

  int numberOfCells = phi_sigma_3D.size();

  int numberOfPolygonData = 9 * numberOfCells;

  int numberOfPoints = (N_3D + 1) * (N_3D + 1) * (N_3D + 1);

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

  std::vector<std::vector<int>> indices;

  for (int i = 0; i < N_3D + 1; i++) { // x-loop

    for (int j = 0; j < N_3D + 1; j++) { // y-loop

      for (int k = 0; k < N_3D + 1; k++) { // z-loop

        filevtk << (double) i * h_3D << " " << (double) j * h_3D << " "
                << (double) k * h_3D << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELLS " << numberOfCells << " " << numberOfPolygonData
          << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        std::vector<double> center = getElementCenter(i, j, k, h_3D);

        filevtk << "8"
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] - 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] - 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] - 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << " "
                << getIndex(center[0] + 0.5 * h_3D, h_3D) * (N_3D + 1) *
                       (N_3D + 1) +
                     getIndex(center[1] + 0.5 * h_3D, h_3D) * (N_3D + 1) +
                     getIndex(center[2] + 0.5 * h_3D, h_3D)
                << std::endl;
      }
    }
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_TYPES " << numberOfCells << std::endl;

  for (int i = 0; i < numberOfCells; i++) {

    filevtk << 11 << std::endl;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfCells << std::endl;
  filevtk << "SCALARS Phi_sigma_(3D) float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < N_3D; i++) { // x-loop

    for (int j = 0; j < N_3D; j++) { // y-loop

      for (int k = 0; k < N_3D; k++) { // z-loop

        int index = i + j * N_3D + k * N_3D * N_3D;

        filevtk << phi_sigma_3D[index] << std::endl;
      }
    }
  }
}


void util::unet::Network::compute_elem_weights() {

  //oss << "  Computing element-weight data for network\n";

  auto pointer = VGM.getHead();

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  bool direct_find = true;
  // const auto &mesh_locator = mesh.point_locator();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    pointer->J_b_points.resize(numberOfNeighbors);

    const Point i_coords = util::to_point(pointer->coord);

    for (int j = 0; j < numberOfNeighbors; j++) {

      // get reference to J_b point data
      auto &J_b_ij = pointer->J_b_points[j];
      J_b_ij.clear();
      J_b_ij.id_seg = pointer->neighbors[j]->index;

      // const auto &j_coords = pointer->neighbors[j]->coord;
      const Point j_coords = util::to_point(pointer->neighbors[j]->coord);

      const auto &ij_rad = pointer->radii[j];

      auto dist_half_ij = 0.5 * (j_coords - i_coords).norm();

      const Point direction = (j_coords - i_coords) / (2. * dist_half_ij);

      J_b_ij.half_cyl_surf = 2. * M_PI * ij_rad * dist_half_ij;

      // discretize [0, dist_half]
      unsigned int N_length = input.d_num_points_length;
      double sum_weights = 0.;

      for (unsigned int ci = 0; ci < N_length; ci++) {

        // parametric coordinate
        double s = (double(ci) + 0.5) * dist_half_ij / (double(N_length));

        // coordinate at center of cross-section
        Point xs = i_coords + s * direction;

        // get basis vector for plane parallel to cross-section and passing
        // through origin

        // Point e1 = xs - (xs * direction) * direction;
        // Point e1 = xs;

        // TODO
        //  Compute correct vector in the plane of disk of cylinder
        //        e1(0) = 1.0;
        //        e1(1) = 1.0;
        //        e1(2) = 0.0;
        Point e1 = util::determineRotator(direction);
        e1 = e1 / e1.norm();

        // out << "e1: " << e1 << std::endl;

        // discretize [0,2pi]
        unsigned int N_theta = input.d_num_points_angle;

        for (unsigned int cj = 0; cj < N_theta; cj++) {

          double theta = (double(cj) + 0.5) * 2. * M_PI / N_theta;

          auto x = xs + util::rotate(ij_rad * e1, theta, direction);

          // out << "x: " << x << std::endl;

          // compute normalize weight of the point x
          double w = 1. / (double(N_length) * N_theta);
          sum_weights += w;

          // out << "sum_weights: " << sum_weights << std::endl;

          // get element of point x
          unsigned int elem_id =
            util::get_elem_id(x, input.d_mesh_size, input.d_num_elems, input.d_dim);
          //          else {
          //            const auto *elem_x = mesh_locator(x);
          //            elem_id = elem_x->id();
          //          }

          // add data to J_b_ij
          J_b_ij.add_unique(elem_id, w);

          //          out << "x: " << util::io::printStr(x)
          //              << ", mesh size: " << input.d_mesh_size
          //              << ", num elems: " << input.d_num_elems
          //              << ", dim: " << input.d_dim
          //              << ", elem_id: " << elem_id << "\n";

        } // loop over theta

      } // loop over length

      //       out << "Node: " << pointer->index  << ", sum weights: " << sum_weights
      //       << "\n";
    }

    pointer = pointer->global_successor;
  }
}

std::vector<util::unet::ElemWeights>
util::unet::Network::compute_elem_weights_at_node(
  std::shared_ptr<util::unet::VGNode> &pointer) const {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  int numberOfNeighbors = pointer->neighbors.size();

  std::vector<ElemWeights> J_b_points(numberOfNeighbors);

  const Point i_coords = util::to_point(pointer->coord);

  for (int j = 0; j < numberOfNeighbors; j++) {

    ElemWeights J_b_ij;

    // const auto &j_coords = pointer->neighbors[j]->coord;
    const Point j_coords = util::to_point(pointer->neighbors[j]->coord);

    const auto &ij_rad = pointer->radii[j];

    auto dist_half_ij = 0.5 * (j_coords - i_coords).norm();

    const Point direction = (j_coords - i_coords) / (2. * dist_half_ij);

    J_b_ij.half_cyl_surf = 2. * M_PI * ij_rad * dist_half_ij;

    double sum_weights = 0.;

    // discretize [0, dist_half]
    unsigned int N_length = input.d_num_points_length;

    for (unsigned int ci = 0; ci < N_length; ci++) {

      // parametric coordinate
      double s = (double(ci) + 0.5) * dist_half_ij / (double(N_length));

      // coordinate at center of cross-section
      Point xs = i_coords + s * direction;

      // get basis vector for plane parallel to cross-section and passing
      // through origin

      // Point e1 = xs - (xs * direction) * direction;
      Point e1 = util::determineRotator(direction);
      e1 = e1 / e1.norm();

      // out << "e1: " << e1 << std::endl;

      // discretize [0,2pi]
      unsigned int N_theta = input.d_num_points_angle;

      for (unsigned int cj = 0; cj < N_theta; cj++) {

        double theta = (double(cj) + 0.5) * 2. * M_PI / N_theta;

        auto x = xs + util::rotate(ij_rad * e1, theta, direction);

        // out << "x: " << x << std::endl;

        // compute normalize weight of the point x
        double w = 1. / (double(N_length) * N_theta);
        sum_weights += w;

        // out << "sum_weights: " << sum_weights << std::endl;

        // get element of point x
        unsigned int elem_id =
          util::get_elem_id(x, input.d_mesh_size, input.d_num_elems, input.d_dim);

        // add data to J_b_ij
        J_b_ij.add_unique(elem_id, w);

        // out << "elem_x->id(): " << elem_x->id() << std::endl;

      } // loop over theta

    } // loop over length

    //    out << "Node: " << pointer->index << ", sum weights: " << sum_weights
    //        << "\n";

    // add to data
    J_b_points[j] = J_b_ij;
  }

  return J_b_points;
}
