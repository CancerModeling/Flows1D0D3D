////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include "../model.hpp"
#include "nodes.hpp"
#include "utilIO.hpp"
#include "utils.hpp"
#include <random>

void netfvfe::Network::create_initial_network() {

  const auto &input = d_model_p->get_input_deck();
  d_is_network_changed = true;

  // equation system
  std::vector<std::vector<double>> vertices;
  std::vector<double> pressures;
  std::vector<double> radii;
  std::vector<std::vector<unsigned int>> elements;

  readData(vertices, pressures, radii, elements);

  transferDataToVGM(vertices, pressures, radii, elements);
  int numberOfNodes = VGM.getNumberOfNodes();
  std::cout << "Number of nodes: " << numberOfNodes << std::endl;

  // refine mesh
  //  std::cout << " " << std::endl;
  std::cout << "Refine mesh" << std::endl;

  int refinementLevel = input.d_network_init_refinement;

  for(int i = 0; i < refinementLevel; i++) {

      refine1DMesh();

  }

  // std::cout << " " << std::endl;
  numberOfNodes = VGM.getNumberOfNodes();

  std::cout << "Number of nodes: " << numberOfNodes << std::endl;

  // compute element and weights
  compute_elem_weights();

  // Print data
  //std::cout << " " << std::endl;
  //std::cout << "Print data" << std::endl;
  //printDataVGM();

  // initialize matrix and vector
  // 1D pressure: matrix, rhs, and solution
  A_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  b = std::vector<double>(numberOfNodes, 0.);
  P_v = std::vector<double>(numberOfNodes, 0.);

  // 1D nutrient: matrix, rhs, and solution
  Ac_VGM = gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  b_c = std::vector<double>(numberOfNodes, 0.);

  C_v = std::vector<double>(numberOfNodes, 0.);
  C_v_old = std::vector<double>(numberOfNodes, 0.);

  mu = input.d_init_vessel_mu;

  D_v = input.d_D_sigma_v;

  osmotic_sigma = input.d_osmotic_sigma;
}

void netfvfe::Network::readData(
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
      pressures.push_back(vertexInfo[3]
                         * 133.322); // - 100000.0 ); //*133.322 );

//      std::cout << vertexInfo[0] << " " << vertexInfo[1] << " " << vertexInfo[2]
//                << " " << vertexInfo[3] << std::endl;
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
      cornerIDs[0] = (int)segmentInfo[0];
      cornerIDs[1] = (int)segmentInfo[1];
      double radius = segmentInfo[2];

      radii.push_back(radius);
      elements.push_back(cornerIDs);

//      std::cout << cornerIDs[0] << " " << cornerIDs[1] << " " << radius << "\n";
    }

    i = i + 1;
  }
}

void netfvfe::Network::transferDataToVGM(
    std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
    std::vector<double> &radii,
    std::vector<std::vector<unsigned int>> &elements) {

  int numberOfVertices = vertices.size();
  const auto &input = d_model_p->get_input_deck();

  for (int i = 0; i < numberOfVertices; i++) {

    netfvfe::VGNode new_node;

    new_node.index = i;

    std::vector<double> coord = vertices[i];

    new_node.coord = coord;

    new_node.p_boundary = pressures[i];

    new_node.p_v = 0.0;

    new_node.c_boundary = input.d_in_nutrient;

    new_node.c_v = 0.0;

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

    double length = 0.0;

    double radius = radii[i];

    pointer_1 = VGM.findNode(index_1);

    pointer_2 = VGM.findNode(index_2);

    pointer_1->radii.push_back(radius);

    pointer_1->L_p.push_back(input.d_tissue_flow_L_p);

    pointer_1->L_s.push_back(input.d_tissue_nut_L_s);

    pointer_1->neighbors.push_back(pointer_2);

    pointer_1->edge_touched.push_back(false);

    pointer_2->radii.push_back(radius);

    pointer_2->L_p.push_back(input.d_tissue_flow_L_p);

    pointer_2->L_s.push_back(input.d_tissue_nut_L_s);

    pointer_2->neighbors.push_back(pointer_1);

    pointer_2->edge_touched.push_back(false);
  }
}

void netfvfe::Network::printDataVGM() {

  std::cout << " " << std::endl;
  std::cout << "PrintData of network: " << std::endl;
  std::cout << " " << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    std::cout << "Index of node [-]: " << pointer->index << std::endl;
    std::cout << "Type of node [-]: " << pointer->typeOfVGNode << std::endl;
    std::cout << "Boundary pressure of node [Pa]: " << pointer->p_boundary
              << std::endl;
    std::cout << "Pressure of node [Pa]: " << pointer->p_v << std::endl;
    std::cout << "Boundary concentration of node [mol/m^3]: "
              << pointer->c_boundary << std::endl;
    std::cout << "Coord [m]: " << pointer->coord[0] << " " << pointer->coord[1]
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

void netfvfe::Network::writeDataToVTK_VGM() {

  std::fstream filevtk;
  filevtk.open("ratbrain_secomb_vgm.vtk", std::ios::out);
  filevtk << "# vtk DataFile Version 2.0" << std::endl;
  filevtk << "Network ETH" << std::endl;
  filevtk << "ASCII" << std::endl;
  filevtk << "DATASET POLYDATA" << std::endl;

  filevtk << "POINTS " << VGM.getNumberOfNodes() << " float" << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    std::vector<double> coord;

    coord = pointer->coord;

    filevtk << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;

    pointer = pointer->global_successor;
  }

  int numberOfNodes = 0;

  pointer = VGM.getHead();

  while (pointer) {

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

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  int polygons = 3 * numberOfNodes;

  filevtk << " " << std::endl;
  filevtk << "LINES " << numberOfNodes << " " << polygons << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        filevtk << "2 " << pointer->index << " " << pointer->neighbors[i]->index
                << std::endl;

        pointer->edge_touched[i] = true;

        pointer->neighbors[i]->markEdge(pointer->index);
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfNodes << std::endl;
  filevtk << "SCALARS radii float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        filevtk << pointer->radii[i] << std::endl;

        pointer->edge_touched[i] = true;

        pointer->neighbors[i]->markEdge(pointer->index);
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  filevtk << " " << std::endl;
  filevtk << "POINT_DATA " << VGM.getNumberOfNodes() << std::endl;
  filevtk << "SCALARS pressure_[mmHg] float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    filevtk << pointer->p_v / 133.322 << std::endl;

    pointer = pointer->global_successor;
  }
}

void netfvfe::Network::writeDataToVTKTimeStep_VGM(int timeStep) {

  std::string path = "ratbrain_secomb_vgm_";
  path += std::to_string(timeStep);
  path.append(".vtk");

  std::fstream filevtk;
  filevtk.open(path, std::ios::out);
  filevtk << "# vtk DataFile Version 2.0" << std::endl;
  filevtk << "Network ETH" << std::endl;
  filevtk << "ASCII" << std::endl;
  filevtk << "DATASET POLYDATA" << std::endl;

  filevtk << "POINTS " << VGM.getNumberOfNodes() << " float" << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    std::vector<double> coord;

    coord = pointer->coord;

    filevtk << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;

    pointer = pointer->global_successor;
  }

  int numberOfNodes = 0;

  pointer = VGM.getHead();

  while (pointer) {

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

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  int polygons = 3 * numberOfNodes;

  filevtk << " " << std::endl;
  filevtk << "LINES " << numberOfNodes << " " << polygons << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        filevtk << "2 " << pointer->index << " " << pointer->neighbors[i]->index
                << std::endl;

        pointer->edge_touched[i] = true;

        pointer->neighbors[i]->markEdge(pointer->index);
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  filevtk << " " << std::endl;
  filevtk << "CELL_DATA " << numberOfNodes << std::endl;
  filevtk << "SCALARS radii float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        filevtk << pointer->radii[i] << std::endl;

        pointer->edge_touched[i] = true;

        pointer->neighbors[i]->markEdge(pointer->index);
      }
    }

    pointer = pointer->global_successor;
  }

  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  filevtk << " " << std::endl;
  filevtk << "POINT_DATA " << VGM.getNumberOfNodes() << std::endl;
  filevtk << "SCALARS pressure_[mmHg] float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    filevtk << pointer->p_v / 133.322 << std::endl;

    pointer = pointer->global_successor;
  }

  filevtk << "SCALARS concentration float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    filevtk << pointer->c_v << std::endl;

    pointer = pointer->global_successor;
  }
}

void netfvfe::Network::solve_system() {
  solveVGMforPressure();
  // solveVGMforNutrient();
}

void netfvfe::Network::assembleVGMSystemForPressure() {

  const auto &tum_sys = d_model_p->get_system();
  const auto &mesh = d_model_p->get_mesh();

  // 3d pressure
  const auto &pres = tum_sys.get_system<TransientLinearImplicitSystem>
      ("Pressure");
  const DofMap &pres_map = pres.get_dof_map();
  const unsigned int v_pres = pres.variable_number("pressure");

  const auto &input = d_model_p->get_input_deck();

  int numberOfNodes = VGM.getNumberOfNodes();

  // reinitialize data
  //  A_VGM =
  //      gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  if (A_VGM.nrows() != numberOfNodes) {

    A_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  } else {

    // A_VGM.clear_mat();

    for (unsigned int i=0; i<A_VGM.nrows(); i++)
      A_VGM[i].clear();
  }

  if (b.size() != numberOfNodes) {

    b.resize(numberOfNodes);
  }

  for (unsigned int i=0; i<b.size(); i++) {

    b[i] = 0.;
  }

  if (P_v.size() != numberOfNodes) {

    P_v.resize(numberOfNodes);
  }

  for (unsigned int i=0; i<P_v.size(); i++) {

    P_v[i] = 0.;
  }

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v_k = pointer->p_v;

    //std::cout << "p_v_k: " << p_v_k << std::endl;

    if (numberOfNeighbors == 1) {

      A_VGM(indexOfNode, indexOfNode) = 1.0;

      b[indexOfNode] = pointer->p_boundary;

    } else {

      // factor to enhance coupling
      double factor = 1.;
      if (pointer->L_p[0] > 1.e-18)
        factor = 1. / pointer->L_p[0];

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

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
            factor * M_PI * std::pow(radius, 4) / (8.0 * length * mu);

        A_VGM(indexOfNode, indexOfNode) -= t_seg;
        A_VGM(indexOfNode, indexNeighbor) += t_seg;

        // implicit for p_v in source
        A_VGM(indexOfNode, indexOfNode) += factor * input.d_coupling_method_theta *
                                           pointer->L_p[i] *
                                           J_b_data.half_cyl_surf;

        // loop over 3d elements
        for(unsigned int e=0; e<J_b_data.elem_id.size(); e++) {

            auto e_id = J_b_data.elem_id[e];
            auto e_w = J_b_data.elem_weight[e];

            // get 3d pressure
            const auto *elem = mesh.elem_ptr(e_id);
            std::vector<unsigned int> dof_indices_pres;
            pres_map.dof_indices(elem, dof_indices_pres, v_pres);
            Real p_t_k = pres.current_solution(dof_indices_pres[0]);

            // explicit for p_t in source
            b[indexOfNode] -=
                factor * pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                ((1. - input.d_coupling_method_theta) * p_v_k - p_t_k);
        }

      } // loop over neighbor segments

       // b[indexOfNode] = 0.0;
    }

    //    std::cout << "p_v_k: " << pointer->p_v << " s: " << pointer->coord[2]
    //              << ", b: " << b[indexOfNode]
    //              << ", L_p: " << util::io::printStr(pointer->L_p) << std::endl;

    pointer = pointer->global_successor;

  }

  //std::cout << "A_VGM: " << A_VGM << std::endl;
}

void netfvfe::Network::solveVGMforPressure() {

  std::cout << "Assemble pressure matrix and right hand side" << std::endl;
  assembleVGMSystemForPressure();

  gmm::iteration iter(10E-20);

//  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(A_VGM, 50, 1e-5);

  gmm::identity_matrix P;

  std::vector<double> P_v = b;

  gmm::bicgstab(A_VGM, P_v, b, P, iter);

  std::cout << P_v << std::endl;

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->p_v = P_v[indexOfNode];

    pointer = pointer->global_successor;
  }

  std::cout << "Write vtk" << std::endl;
  writeDataToVTK_VGM();
}

void netfvfe::Network::assembleVGMSystemForNutrient() {

  const auto &tum_sys = d_model_p->get_system();
  const auto &mesh = d_model_p->get_mesh();

  // 3d pressure
  const auto &pres = tum_sys.get_system<TransientLinearImplicitSystem>
      ("Pressure");
  const DofMap &pres_map = pres.get_dof_map();
  const unsigned int v_pres = pres.variable_number("pressure");

  // 3d nutrient
  const auto &nut = tum_sys.get_system<TransientLinearImplicitSystem>
      ("Nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  const unsigned int v_nut = nut.variable_number("nutrient");

  const auto &input = d_model_p->get_input_deck();

  int numberOfNodes = VGM.getNumberOfNodes();

  // reinitialize data
  //  Ac_VGM =
  //      gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);
  if (Ac_VGM.nrows() != numberOfNodes) {

    Ac_VGM =
        gmm::row_matrix<gmm::wsvector<double>>(numberOfNodes, numberOfNodes);

  } else {

    // Ac_VGM.clear_mat();

    for (unsigned int i=0; i<Ac_VGM.nrows(); i++)
      Ac_VGM[i].clear();
  }

  if (b_c.size() != numberOfNodes) {

    b_c.resize(numberOfNodes);
  }

  for (unsigned int i=0; i<b_c.size(); i++) {

    b_c[i] = 0.;
  }

  if (C_v.size() != numberOfNodes) {

    C_v.resize(numberOfNodes);
  }

//  for (unsigned int i=0; i<C_v.size(); i++) {
//
//    C_v[i] = 0.;
//  }

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  double dt = d_model_p->d_dt;

  while (pointer) {

    int indexOfNode = pointer->index;

    int numberOfNeighbors = pointer->neighbors.size();

    std::vector<double> coord = pointer->coord;

    double p_v_k = pointer->p_v;

    double c_v_k = pointer->c_v;

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

      double v_interface =
          -(radius * radius * M_PI) / (8.0 * length * mu) * (p_neighbor - p_v_k);

      double volume = length / 2.0 * radius * radius * M_PI;

      if (v_interface > 0.0) {

        Ac_VGM(indexOfNode, indexOfNode) = 1.0;

        // Debug
        // Remove when done debugging
        if (p_v_k < 133.322 * input.d_identify_vein_pres)
          b_c[indexOfNode] = input.d_in_nutrient_vein;
        else
          b_c[indexOfNode] = input.d_in_nutrient;

      } else {

        Ac_VGM(indexOfNode, indexOfNode) = length;

        Ac_VGM(indexOfNode, indexNeighbor) =
            dt * v_interface - dt * D_v / length;

        Ac_VGM(indexOfNode, indexOfNode) = Ac_VGM(indexOfNode, indexOfNode) -
                                           dt * v_interface + dt * D_v / length;

        b_c[indexOfNode] = length * C_v_old[indexOfNode];
      }

    } else {

      // factor to enhance coupling
      double factor = 1.;
//      if (pointer->L_s[0] > 1.e-18)
//        factor = 1. / pointer->L_s[0];

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

        double radius = pointer->radii[i];

        int indexNeighbor = pointer->neighbors[i]->index;

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

        Ac_VGM(indexOfNode, indexOfNode) += factor * length;

        if (v_interface > 0.0) {

          Ac_VGM(indexOfNode, indexOfNode) += factor * dt * v_interface;

        } else {

          Ac_VGM(indexOfNode, indexNeighbor) += factor * dt * v_interface;
        }

        Ac_VGM(indexOfNode, indexOfNode) += factor * dt * D_v / length;

        Ac_VGM(indexOfNode, indexNeighbor) -= factor * dt * D_v / length;

        b_c[indexOfNode] += factor * length * C_v_old[indexOfNode];


        // coupling between 3d and 1d nutrient
        {
          // implicit part of the coupling
          Ac_VGM(indexOfNode, indexOfNode) +=
              factor * dt * input.d_coupling_method_theta * pointer->L_s[i] *
              J_b_data.half_cyl_surf;

          // compute explicit part of the coupling
          for(unsigned int e=0; e<J_b_data.elem_id.size(); e++) {

            auto e_id = J_b_data.elem_id[e];
            auto e_w = J_b_data.elem_weight[e];

            // get 3d pressure
            const auto *elem = mesh.elem_ptr(e_id);
            std::vector<unsigned int> dof_indices_pres;
            pres_map.dof_indices(elem, dof_indices_pres, v_pres);
            Real p_t_k = pres.current_solution(dof_indices_pres[0]);

            // get 3d nutrient
            std::vector<unsigned int> dof_indices_nut;
            nut_map.dof_indices(elem, dof_indices_nut, v_nut);
            Real c_t_k = nut.current_solution(dof_indices_nut[0]);

            // explicit
            double source = dt * pointer->L_s[i] * J_b_data.half_cyl_surf * e_w *
                            ((1. - input.d_coupling_method_theta) * c_v_k -
                            c_t_k);
            b_c[indexOfNode] -= factor * source;

//            if (pointer->p_v <
//                    input.d_mmhgFactor * input.d_identify_vein_pres &&
//                i == 0 && e < 2)
//              out << "index: " << indexOfNode << ", neighbor: " << indexNeighbor
//                  << ", source: " << source << ", pressure: " << p_v_k
//                  << ", radius: " << radius << ", c_t: " << c_t_k
//                  << ", L_s: " << pointer->L_s[i] << "\n";

            // term due to pressure difference
            double c_transport = 0.;
            if (p_v_k - p_t_k >= 0.)
              c_transport = c_v_k;
            else
              c_transport = c_t_k;
            b_c[indexOfNode] -=
                factor * dt * (1. - osmotic_sigma) * pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                (p_v_k - p_t_k) * c_transport;
          }
        } // coupling
      }
    }

    pointer = pointer->global_successor;
  }
}

void netfvfe::Network::solveVGMforNutrient() {

  int numberOfNodes = VGM.getNumberOfNodes();

  // Solver
  gmm::iteration iter(5.0e-11);
  gmm::identity_matrix PR;

  std::cout << " " << std::endl;
  std::cout << "Assemble nutrient matrix and right hand side" << std::endl;
  assembleVGMSystemForNutrient();

  C_v = C_v_old;

  gmm::ilut_precond<gmm::row_matrix<gmm::wsvector<double>>> P(Ac_VGM, 50, 1e-6);

  gmm::bicgstab(Ac_VGM, C_v, b_c, P, iter);

  // std::cout << C_v << std::endl;

  auto pointer = VGM.getHead();

  while (pointer) {

    int indexOfNode = pointer->index;

    pointer->c_v = C_v[indexOfNode];

    pointer = pointer->global_successor;
  }

  C_v_old = C_v;
}

void netfvfe::Network::refine1DMesh() {

  const auto &input = d_model_p->get_input_deck();

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  std::vector<std::shared_ptr<VGNode>> newNodes;

  int numberOfNodes = VGM.getNumberOfNodes();

  int counter = 0;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (!pointer->edge_touched[i]) {

        netfvfe::VGNode new_node;

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

        pointer_former_neighbor->replacePointerWithIndex(pointer->index,
                                                         pointer_new_node);

        pointer->edge_touched[i] = true;

        pointer_new_node->edge_touched.push_back(true);

        pointer_new_node->edge_touched.push_back(true);

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

    // std::cout << "index: " << pointer->index << std::endl;

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->edge_touched[i] = false;
    }

    pointer->index = counter;

    counter = counter + 1;

    // std::cout << " " << std::endl;

    pointer = pointer->global_successor;
  }
}

void netfvfe::Network::update_network() {}

void netfvfe::Network::compute_elem_weights() {

  out << "Computing element-weight data for network\n";

  auto pointer = VGM.getHead();

  const auto &mesh_locator = d_model_p->get_mesh().point_locator();
  const auto&input = d_model_p->get_input_deck();
  const double mesh_size = input.d_mesh_size_vec[0];
  const unsigned int num_elems = input.d_num_elems;

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    pointer->J_b_points.resize(numberOfNeighbors);

    const Point i_coords =
        Point(pointer->coord[0], pointer->coord[1], pointer->coord[2]);

    for (int j = 0; j < numberOfNeighbors; j++) {

      // get reference to J_b point data
      auto &J_b_ij = pointer->J_b_points[j];

      // const auto &j_coords = pointer->neighbors[j]->coord;
      const Point j_coords = Point(pointer->neighbors[j]->coord[0],
                                   pointer->neighbors[j]->coord[1],
                                   pointer->neighbors[j]->coord[2]);

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
        Point e1 = xs;
       
        e1(0) = 1.0;
        e1(1) = 1.0;
        e1(2) = 0.0;
        e1 = e1 / e1.norm();

        //out << "e1: " << e1 << std::endl;

        // discretize [0,2pi]
        unsigned int N_theta = input.d_num_points_angle;

        for (unsigned int cj = 0; cj < N_theta; cj++) {

          double theta = (double(cj) + 0.5) * 2. * M_PI / N_theta;

          auto x = xs + util::rotate(ij_rad * e1, theta, direction);

          //out << "x: " << x << std::endl;

          // compute normalize weight of the point x
          double w = 1. / (double(N_length) * N_theta);
          sum_weights += w;

          // out << "sum_weights: " << sum_weights << std::endl;

          // get element of point x
          unsigned int elem_id = 0;
          bool direct_find = true;
          if (!direct_find) {
            const auto *elem_x = mesh_locator(x);
            elem_id = elem_x->id();
          } else {
            unsigned int i = x(0) / mesh_size;
            unsigned int j = x(1) / mesh_size;
            unsigned int k = x(2) / mesh_size;
            elem_id = k * num_elems * num_elems + j * num_elems + i;
          }


          // add data to J_b_ij
           J_b_ij.add_unique(elem_id, w);
    
          //out << "elem_x->id(): " << elem_x->id() << std::endl;

         } // loop over theta

      } // loop over length

      //std::cout << " "<<std::endl;
      // out << "Node: " << pointer->index  << ", sum weights: " << sum_weights
      // << "\n";
    }

    pointer = pointer->global_successor;
  }

}
