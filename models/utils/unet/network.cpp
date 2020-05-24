////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "network.hpp"
#include <random>

namespace {

std::ostringstream oss;

void angle_correction(const Point &parent_d, Point &child_d,
                      const double &max_angle) {

  return;

  auto child_angle = util::angle(child_d, parent_d);
  if (std::abs(child_angle) > max_angle) {

    out << "Child direction: " << child_d << ", child angle: " << child_angle
        << "\n";

    // axis for rotation
    Point axis = parent_d.cross(child_d);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    //    child_d = util::rotate(parent_d, max_angle, axis);
    out << "New child direction: " << util::rotate(parent_d, max_angle, axis)
        << "\n";
  }
}

double compute_bifurcate_child_direction(const Point &parent_d,
                                       const Point &child_d, Point &child_d_2,
                                       const double &branch_angle) {

  // get angle between child direction and parent direction and offset
  // it by branching angle
  double angle_1 = util::angle(parent_d, child_d) + branch_angle;

  // axis for rotation
  auto axis = Point();
  if (util::angle(parent_d, child_d) > 1.0e-10)
    axis = parent_d.cross(child_d);
  else
    axis = util::determineRotator(parent_d);
  axis = axis / axis.norm();

  // rotate parent_direction by -ve angle
  child_d_2 = util::rotate(parent_d, -angle_1, axis);

  return angle_1;
}

void angle_correction_bifurcation(const Point &parent_d, Point &child_d,
                                  Point &child_d_2, const double &max_angle,
                                  const double &branch_angle) {

  return;

  // check if angle of direction from parent direction is within
  // permissible range
  auto child_angle = util::angle(child_d_2, parent_d);
  if (std::abs(child_angle) > max_angle) {

    // axis for rotation
    Point axis = parent_d.cross(child_d_2);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    child_d_2 = util::rotate(parent_d, max_angle, axis);

  } else if (std::abs(child_angle) < branch_angle) {

    // axis for rotation
    Point axis = parent_d.cross(child_d_2);
    axis = axis / axis.norm();

    // rotate parent direction by allowed angle
    child_d_2 = util::rotate(parent_d, branch_angle, axis);
  }
}
} // namespace

void util::unet::Network::create_initial_network() {

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
  oss << "  Number of nodes: " << numberOfNodes << std::endl;

  // refine mesh
  //  std::cout << " " << std::endl;
  oss << "  Refine mesh" << std::endl;

  int refinementLevel = input.d_network_init_refinement;

  for (int i = 0; i < refinementLevel; i++) {

    refine1DMesh();
  }

  // std::cout << " " << std::endl;
  numberOfNodes = VGM.getNumberOfNodes();

  oss << "  Number of nodes: " << numberOfNodes << std::endl;

  // compute element and weights
  if (input.d_compute_elem_weights)
    compute_elem_weights();

  d_model_p->d_delayed_msg = oss.str();
  util::clear_oss(oss);

  // Print data
  // std::cout << " " << std::endl;
  // std::cout << "Print data" << std::endl;
  // printDataVGM();

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

  // Debug
  auto locate_point = Point(1.0683156654888104, 1.0478309232480534, 1.1339173967459324);
  std::shared_ptr<VGNode> pointer = VGM.getHead();
}

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

      //      std::cout << vertexInfo[0] << " " << vertexInfo[1] << " " <<
      //      vertexInfo[2]
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

      //      std::cout << cornerIDs[0] << " " << cornerIDs[1] << " " << radius
      //      << "\n";
    }

    i = i + 1;
  }
}

void util::unet::Network::transferDataToVGM(
    std::vector<std::vector<double>> &vertices, std::vector<double> &pressures,
    std::vector<double> &radii,
    std::vector<std::vector<unsigned int>> &elements) {

  int numberOfVertices = vertices.size();
  const auto &input = d_model_p->get_input_deck();

  for (int i = 0; i < numberOfVertices; i++) {

    VGNode new_node;

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

    pointer_1->sprouting_edge.push_back(false);

    pointer_2->radii.push_back(radius);

    pointer_2->L_p.push_back(input.d_tissue_flow_L_p);

    pointer_2->L_s.push_back(input.d_tissue_nut_L_s);

    pointer_2->neighbors.push_back(pointer_1);

    pointer_2->edge_touched.push_back(false);

    pointer_2->sprouting_edge.push_back(false);
  }
}

void util::unet::Network::printDataVGM() {

  oss << " " << std::endl;
  oss << "PrintData of network: " << std::endl;
  oss << " " << std::endl;
  d_model_p->d_log(oss);

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    oss << "Index of node [-]: " << pointer->index << std::endl;
    oss << "Type of node [-]: " << pointer->typeOfVGNode << std::endl;
    oss << "Boundary pressure of node [Pa]: " << pointer->p_boundary
              << std::endl;
    oss << "Pressure of node [Pa]: " << pointer->p_v << std::endl;
    oss << "Boundary concentration of node [mol/m^3]: "
              << pointer->c_boundary << std::endl;
    oss << "Coord [m]: " << pointer->coord[0] << " " << pointer->coord[1]
              << " " << pointer->coord[2] << std::endl;

    int numberOfNeighbors = pointer->neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      oss << "L_p [m/(Pa s)]: " << pointer->L_p[i] << std::endl;
      oss << "L_s [m/s]: " << pointer->L_s[i] << std::endl;
      oss << "Radii [m]: " << pointer->radii[i] << std::endl;
      oss << "Neighbor [-]: " << pointer->neighbors[i]->index
                << std::endl;
    }

    oss << " " << std::endl;
    d_model_p->d_log(oss);

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::writeDataToVTKTimeStep_VGM(int timeStep) {

  const auto &input = d_model_p->get_input_deck();

  std::string path = input.d_outfilename_net + "_";
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
  filevtk << "SCALARS pressure_1d float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    filevtk << pointer->p_v << std::endl;

    pointer = pointer->global_successor;
  }

  filevtk << "SCALARS nutrient_1d float 1" << std::endl;
  filevtk << "LOOKUP_TABLE default" << std::endl;

  pointer = VGM.getHead();

  while (pointer) {

    filevtk << pointer->c_v << std::endl;

    pointer = pointer->global_successor;
  }
}

void util::unet::Network::assembleVGMSystemForPressure(BaseAssembly &pres_sys) {

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

void util::unet::Network::assembleVGMSystemForNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys) {

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

void util::unet::Network::solveVGMforNutrient(BaseAssembly &pres_sys, BaseAssembly &nut_sys) {

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

void util::unet::Network::update_network(BaseAssembly &taf_sys, BaseAssembly &grad_taf_sys) {

  if (taf_sys.d_sys_name != "TAF" or grad_taf_sys.d_sys_name != "TAF_Gradient")
    libmesh_error_msg("Must pass TAF and TAF_Gradient system to update "
                      "network");

  oss << "[Mark]";
  d_model_p->d_log(oss);
  auto num_nodes_marked = markApicalGrowth("apical_taf_based", taf_sys);

  unsigned int num_new_nodes = 0;
  if (num_nodes_marked > 0) {
    oss << " -> [Process]";
    d_model_p->d_log(oss);
    num_new_nodes = processApicalGrowthTAF(taf_sys, grad_taf_sys);
  }
  oss << "\n";

  oss << "\n    ____________________\n";
  oss << "    New nodes added: " << num_new_nodes << "\n";
  d_model_p->d_log(oss);

  if (num_new_nodes > 0)
    d_update_number++;
}

void util::unet::Network::compute_elem_weights() {

  oss << "  Computing element-weight data for network\n";
  d_model_p->d_log(oss);

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

unsigned int util::unet::Network::markApicalGrowth(std::string growth_type, BaseAssembly &taf_sys) {

  const auto &tum_sys = d_model_p->get_system();
  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();
  const double domain_size = input.d_domain_params[1];

  // mark node for growth based on a certain criterion
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  unsigned int num_nodes_marked = 0;

  if (growth_type == "test") {

    while (pointer) {

      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {
        pointer->apicalGrowth = true;
        num_nodes_marked++;
      }

      pointer = pointer->global_successor;
    }
  } else if (growth_type == "apical_taf_based") {

    while (pointer) {

      //      std::cout << " " << std::endl;
      //      std::cout << "Type of growth: apical_taf_based" << std::endl;

      if ((pointer->typeOfVGNode == TypeOfNode::DirichletNode and
           d_update_number > 0) or
          (pointer->typeOfVGNode == TypeOfNode::InnerNode and
           d_update_number == 0)) {

        const auto &coord = pointer->coord;

        // if node is near the boundary, we do not process
        if (util::is_inside_box(Point(coord[0], coord[1], coord[2]),
                                domain_size, 0.01 * domain_size)) {

          // get taf
          const auto *elem_i = mesh.elem_ptr(util::get_elem_id(
              coord, input.d_mesh_size, input.d_num_elems, input.d_dim));

          // Note: TAF may be finite-element field in which case a more
          // accurate way would be to find the nearest quadrature point to
          // the coordinate of this node pointer->coord
          // But for marking, it is enough to just take taf at first vertex
          // of element
          //
          // If TAF is finite-volume than below is correct
          taf_sys.init_dof(elem_i);
          Real taf_i = taf_sys.get_current_sol(0);

          if (taf_i > input.d_network_update_taf_threshold) {
            pointer->apicalGrowth = true;
            num_nodes_marked++;
          }

        } // if inside domain

      } // check if boundary node

      pointer = pointer->global_successor;
    }
  } else {

    std::cout << " " << std::endl;
    std::cout << "Criterion not implemented!!! Program is terminated here"
              << std::endl;
    exit(1);
  }

  return num_nodes_marked;
}

unsigned int util::unet::Network::processApicalGrowthTAF(BaseAssembly &taf_sys,
                                                   BaseAssembly &grad_taf_sys) {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();
  const double domain_size = input.d_domain_params[1];
  const unsigned int dim = input.d_dim;

  // mark node for growth based on a certain criterion
  std::shared_ptr<VGNode> pointer = VGM.getHead();

  unsigned int num_nodes_old = VGM.getNumberOfNodes();

  // Initialize random objects
  std::lognormal_distribution<> log_normal_distribution(
      input.d_log_normal_mean, input.d_log_normal_std_dev);
  std::random_device rd;
  std::mt19937 generator(rd());

  // store how many new node
  unsigned int num_new_nodes_added = 0;

  while (pointer) {

    if (pointer->apicalGrowth) {

      oss << "\n      Processing node: " << pointer->index << "\n";
      d_model_p->d_log(oss);

      // compute radius, direction, length
      unsigned int num_segments = pointer->neighbors.size();
      double parent_r = 0.;
      double parent_l = 0.;
      auto parent_d = Point();
      {
        // get avg radius at the node
        auto parent_d_temp = Point();

        for (unsigned int i = 0; i < num_segments; i++) {

          parent_r += pointer->radii[i];

          parent_l += util::dist_between_points(pointer->neighbors[i]->coord,
                                                pointer->coord);

          auto neigh_d = Point();
          for (unsigned int j = 0; j < dim; j++)
            neigh_d(j) = pointer->coord[j] - pointer->neighbors[i]->coord[j];
          neigh_d = neigh_d / neigh_d.norm();

          parent_d_temp += neigh_d;
        }

        // get avg
        parent_r = parent_r / num_segments;
        parent_l = parent_l / num_segments;
        parent_d_temp = parent_d_temp / num_segments;

        // it is possible that this node has two segments in opposite direction
        // in which case parent_d vector will be very small
        if (parent_d_temp.norm_sq() < 0.0001 * domain_size) {

          // take direction along first segment
          for (unsigned int j = 0; j < dim; j++)
            parent_d_temp(j) =
                pointer->coord[j] - pointer->neighbors[0]->coord[j];
          parent_d_temp = parent_d_temp / parent_d_temp.norm();
        }

        if (pointer->typeOfVGNode == TypeOfNode::InnerNode)
          parent_d = Point();
        else
          parent_d = parent_d_temp;
      }

      // compute taf and grad taf at the node
      double parent_taf = 0.;
      Point parent_grad_taf;
      get_taf_and_gradient(taf_sys, grad_taf_sys, parent_taf, parent_grad_taf,
                           pointer->coord);

      oss << "      Parent info R: " << parent_r << ", L: " << parent_l
          << ", d: " << util::io::printStr(parent_d)
          << ", x: (" << util::io::printStr(pointer->coord) << ")\n";
      oss << "        taf: " << parent_taf
          << ", Grad(TAF): " <<  util::io::printStr(parent_grad_taf)
          << "\n";
      d_model_p->d_log(oss);

      // compute child direction
      Point child_d = input.d_net_direction_lambda_g * parent_d;
      if (parent_grad_taf.norm_sq() > 1.0e-8)
        child_d += parent_grad_taf / parent_grad_taf.norm();
      child_d = child_d / child_d.norm();

      // see if we need to correct the direction so that the angle it makes
      // is within reasonable range (Correct only if this is DirichletNode)
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode)
        angle_correction(parent_d, child_d, input.d_new_vessel_max_angle);

      // lognormal distribution
      double log_dist = log_normal_distribution(generator);

      // check if we bifurcate at this node
      bool bifurcate = false;
      if (pointer->typeOfVGNode == TypeOfNode::DirichletNode) {
        double prob =
            0.5 +
            0.5 * std::erf((std::log(log_dist) - input.d_log_normal_mean) /
                           std::sqrt(2. * input.d_log_normal_std_dev *
                                     input.d_log_normal_std_dev));

        if (prob > input.d_network_bifurcate_prob)
          bifurcate = true;

        // if this node already has 3 neighbors then we can not branch
        if (pointer->neighbors.size() >= 3)
          bifurcate = false;

        //bifurcate = false;
      }



      if (!bifurcate) {

        // get radius
        double child_r = parent_r;

        // get length
        double child_l = log_dist * child_r;

        // get end point of new vessel
        auto child_end_point = std::vector<double>(3, 0.);
        for (unsigned int i = 0; i < dim; i++)
          child_end_point[i] = pointer->coord[i] + child_l * child_d(i);

        unsigned int check_code = 0;
        auto near_node =
            check_new_node(pointer->index, pointer->coord, child_end_point,
                           0.1 * child_l, domain_size, check_code);

        oss << "      Child info R: " << child_r << ", L: " << child_l
            << ", d: " << util::io::printStr(child_d)
            << ", check point: " << check_code
            << ", x end: (" << util::io::printStr(child_end_point) << ")\n";
        d_model_p->d_log(oss);

        // check if the new point is too close to existing node
        if (check_code == 0)
          add_new_node(pointer, child_r, child_end_point);
        else if (check_code == 2) {

          // create connection between parent node and near_node_index
          add_new_node_at_existing_node(pointer, near_node);
        }

      } else {

        // get second child's direction
        Point child_d_2 = Point();
        auto angle = compute_bifurcate_child_direction(parent_d, child_d,
                                                    child_d_2,
                                          input.d_branch_angle);

        oss << "      Bifurcation angle: " << angle << "\n";
        d_model_p->d_log(oss);


        if (pointer->typeOfVGNode == TypeOfNode::DirichletNode)
          angle_correction_bifurcation(parent_d, child_d, child_d_2,
                                       input.d_new_vessel_max_angle,
                                       input.d_branch_angle);

        // compute radius to two childs
        double R_c =
            std::pow(2., -1. / input.d_net_radius_exponent_gamma) * parent_r;

        if (R_c > parent_r)
          R_c = parent_r;

        // create normal distribution function
        std::normal_distribution<> normal_distribution(R_c, R_c / 32.);

        double child_r = normal_distribution(generator);
        // method 1
        //      double R_2 = std::pow(std::pow(parent_R,
        //      input.d_net_radius_exponent_gamma) -
        //                   std::pow(R_1, input.d_net_radius_exponent_gamma),
        //                   1. / input.d_net_radius_exponent_gamma);
        // method 2
        double child_r_2 = normal_distribution(generator);

        if (child_r < 3e-6)
          child_r = 3e-6;
        if (child_r_2 < 3e-6)
          child_r_2 = 3e-6;

        // compute length
        double child_l = log_dist * child_r;
        double child_l_2 = child_l;

        // add first child
        {
          // get end point of new vessel
          auto child_end_point = std::vector<double>(3, 0.);
          for (unsigned int i = 0; i < dim; i++)
            child_end_point[i] = pointer->coord[i] + child_l * child_d(i);

          unsigned int check_code = 0;
          auto near_node =
              check_new_node(pointer->index, pointer->coord, child_end_point,
                             0.1 * child_l, domain_size, check_code);

          oss << "      Child 1 info R: " << child_r << ", L: " << child_l
              << ", d: " << util::io::printStr(child_d)
              << ", check point: " << check_code
              << ", x end: (" << util::io::printStr(child_end_point) << ")\n";
          d_model_p->d_log(oss);

          // check if the new point is too close to existing node
          if (check_code == 0)
            add_new_node(pointer, child_r, child_end_point);
          else if (check_code == 2) {

            // create connection between parent node and near_node_index
            add_new_node_at_existing_node(pointer, near_node);
          }
        }

        // add second child
        {
          // get end point of new vessel
          auto child_end_point = std::vector<double>(3, 0.);
          for (unsigned int i = 0; i < dim; i++)
            child_end_point[i] = pointer->coord[i] + child_l_2 * child_d_2(i);

          unsigned int check_code = 0;
          auto near_node =
              check_new_node(pointer->index, pointer->coord, child_end_point,
                             0.1 * child_l_2, domain_size, check_code);

          oss << "      Child 2 info R: " << child_r_2 << ", L: " << child_l_2
              << ", d: " << util::io::printStr(child_d_2)
              << ", check point: " << check_code
              << ", x end: (" << util::io::printStr(child_end_point) << ")\n";
          d_model_p->d_log(oss);

          // check if the new point is too close to existing node
          if (check_code == 0)
            add_new_node(pointer, child_r_2, child_end_point);
          else if (check_code == 2) {

            // create connection between parent node and near_node_index
            add_new_node_at_existing_node(pointer, near_node);
          }
        }
      }
    } // if apical growth

    pointer = pointer->global_successor;
  }

  // reset the boolean values
  pointer = VGM.getHead();

  while (pointer) {

    int numberOfNeighbors = pointer->neighbors.size();

    pointer->apicalGrowth = false;

    for (int i = 0; i < numberOfNeighbors; i++) {

      pointer->sprouting_edge[i] = false;

      pointer->edge_touched[i] = false;
    }

    pointer = pointer->global_successor;
  }

  return VGM.getNumberOfNodes() - num_nodes_old;
}

std::shared_ptr<util::unet::VGNode> util::unet::Network::check_new_node(
    const int &parent_index, const std::vector<double> &parent_coord,
    const std::vector<double> &child_coord, const double &dist_tol,
    const double &domain_size, unsigned int &check_code) {

  // check_code
  // 0 - no issue
  // 1 - new point is too close to boundary
  // 2 - new point is too close to one of the node
  // 3 - new segment is too close to existing segments

  // TODO
  //  In case of 3, we may want to move the child_coord to the existing
  //  nearby node. This way we can form the close loop.

  if (!util::is_inside_box(
          Point(child_coord[0], child_coord[1], child_coord[2]), domain_size,
          dist_tol)) {
    check_code = 1;
    return nullptr;
  }

  auto p_parent = util::to_point(parent_coord);
  auto p_child = util::to_point(child_coord);

  std::shared_ptr<VGNode> pointer = VGM.getHead();

  while (pointer) {

    const auto &coord = pointer->coord;

    if (util::dist_between_points(coord, child_coord) < dist_tol) {
      check_code = 2;
      return pointer;
    }

    if (pointer->index != parent_index) {

      // check the new segment might intersect or get too close to the segments
      // eminating from this node
      auto pi = util::to_point(coord);
      for (unsigned int j = 0; j < pointer->neighbors.size(); j++) {

        if (pointer->neighbors[j]->index == parent_index)
          continue;

        auto pj = util::to_point(pointer->neighbors[j]->coord);

        if (util::distance_between_segments({p_parent, p_child}, {pi, pj}) <
            dist_tol) {
          check_code = 3;
          return nullptr;
        }
      }
    }

    pointer = pointer->global_successor;
  }

  check_code = 0;
  return nullptr;
}

void util::unet::Network::add_new_node(std::shared_ptr<util::unet::VGNode> &pointer,
                                   const double &child_r,
                                   const std::vector<double> &child_end_point) {

  const auto &input = d_model_p->get_input_deck();

  VGNode new_node;

  new_node.index = VGM.getNumberOfNodes();

  new_node.typeOfVGNode = TypeOfNode::DirichletNode;

  new_node.L_p.push_back(input.d_tissue_flow_L_p);
  new_node.L_s.push_back(input.d_tissue_nut_L_s);

  new_node.radii.push_back(child_r);

  new_node.edge_touched = std::vector<bool>(1, true);

  new_node.sprouting_edge = std::vector<bool>(1, false);

  new_node.apicalGrowth = false;

  new_node.p_v = pointer->p_v;

  new_node.p_boundary = pointer->p_boundary;

  new_node.c_v = 0.0;

  new_node.c_boundary = 0.0;

  new_node.neighbors.push_back(pointer);

  new_node.coord = child_end_point;

  std::shared_ptr<VGNode> pointer_new_node = std::make_shared<VGNode>(new_node);

  pointer->neighbors.push_back(pointer_new_node);

  pointer->radii.push_back(child_r);

  pointer->L_p.push_back(input.d_tissue_flow_L_p);
  pointer->L_s.push_back(input.d_tissue_nut_L_s);

  pointer->typeOfVGNode = TypeOfNode::InnerNode;

  pointer->sprouting_edge.push_back(false);

  pointer->edge_touched.push_back(true);

  VGM.attachPointerToNode(pointer_new_node);
}

void util::unet::Network::add_new_node_at_existing_node(
    std::shared_ptr<util::unet::VGNode> &pointer, std::shared_ptr<util::unet::VGNode> &near_node) {

  const auto &input = d_model_p->get_input_deck();

  // recompute radius of child as mean radius
  double child_r = 0;
  for (unsigned int i = 0; i < pointer->neighbors.size(); i++)
    child_r += pointer->radii[i];
  child_r = child_r / pointer->neighbors.size();

  double child_r_2 = 0;
  for (unsigned int i = 0; i < near_node->neighbors.size(); i++)
    child_r_2 += near_node->radii[i];
  child_r_2 = child_r_2 / near_node->neighbors.size();

  child_r = 0.5 * child_r + 0.5 * child_r_2;

  // update parent pointer
  pointer->neighbors.push_back(near_node);

  pointer->radii.push_back(child_r);

  pointer->L_p.push_back(input.d_tissue_flow_L_p);
  pointer->L_s.push_back(input.d_tissue_nut_L_s);

  pointer->typeOfVGNode = TypeOfNode::InnerNode;

  pointer->sprouting_edge.push_back(false);

  pointer->edge_touched.push_back(true);

  // update near node
  near_node->neighbors.push_back(pointer);

  near_node->radii.push_back(child_r);

  near_node->L_p.push_back(input.d_tissue_flow_L_p);
  near_node->L_s.push_back(input.d_tissue_nut_L_s);

  near_node->typeOfVGNode = TypeOfNode::InnerNode;

  near_node->sprouting_edge.push_back(false);

  near_node->edge_touched.push_back(true);
}

void util::unet::Network::get_taf_and_gradient(BaseAssembly &taf_sys,
                                         BaseAssembly &grad_taf_sys,
                                         double &taf_val, Point &grad_taf_val,
                                         const std::vector<double> &coord) const {

  const auto &mesh = d_model_p->get_mesh();
  const auto &input = d_model_p->get_input_deck();

  // reset
  taf_val = 0.;
  grad_taf_val = Point(0., 0., 0.);

  // compute analytically
  bool compute_analytically = true;
  if (compute_analytically and (input.d_test_name == "test_taf" or
       input.d_test_name == "test_taf_2") and !input.d_taf_source_type.empty()) {

    auto x = util::to_point(coord);

    for (int i=0; i<input.d_taf_source_type.size(); i++) {

      // get center
      auto xc = util::to_point(input.d_taf_source_center[i]);

      if ((xc - x).norm() < input.d_taf_source_radius[i]) {
        taf_val += 0.;
        grad_taf_val += Point();
      } else {

        taf_val += (xc - x).norm();
        grad_taf_val += (xc - x) / taf_val;
      }
    }

    return;
  }

  // get element, update dof maps, update quadrature data
  const auto *elem_i = mesh.elem_ptr(
      util::get_elem_id(coord, input.d_mesh_size, input.d_num_elems,
                        input.d_dim));

  // initialize dofs at the element
  taf_sys.init_dof(elem_i);
  grad_taf_sys.init_dof(elem_i);

  // compute values depending on the model type
  if (input.d_model_name == "NetFV") {

    taf_val = taf_sys.get_current_sol(0);
    for (unsigned int var = 0; var < input.d_dim; var++)
      grad_taf_val(var) = grad_taf_sys.get_current_sol_var(0, var);

  } else if (input.d_model_name == "NetFVFE") {

    // taf and grad_taf are fe field so find the quadrature point closest to
    // coord and evaluate fields at that point
    unsigned int qp_near = 0;
    double dist = input.d_mesh_size;
    for (unsigned int i = 0; i < taf_sys.d_qpoints.size(); i++) {
      Point dx = taf_sys.d_qpoints[i] - util::to_point(coord);

      if (dist > dx.norm()) {
        dist = dx.norm();
        qp_near = i;
      }
    }

    // evaluate fields
    for (unsigned int l = 0; l < taf_sys.d_phi.size(); l++) {

      taf_val += taf_sys.d_phi[l][qp_near] * taf_sys.get_current_sol(l);

      // use interpolation function of tar even for grad_taf
      // Since grad_taf is system with multiple variables (equal to dimension
      // of the problem), we use get_current_sol_var
      for (unsigned int i = 0; i < input.d_dim; i++)
        grad_taf_val(i) +=
            taf_sys.d_phi[l][qp_near] * grad_taf_sys.get_current_sol_var(l, i);
    }
  }
}