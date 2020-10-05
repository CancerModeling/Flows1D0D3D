////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "utilIO.hpp"
#include "utils.hpp"
#include <cmath>

bool util::io::read_network_file(const std::string &filename,
                                 std::vector<Point> &nodes,
                                 std::vector<std::vector<unsigned int>> &elems,
                                 std::vector<unsigned int> &node_boundary_flag,
                                 std::vector<bool> &node_branch_flag,
                                 std::vector<double> &elem_radius,
                                 std::vector<double> &elem_viscosity) {

  // read node and element from file
  std::ifstream filein(filename);

  if (!filein) {
    libmesh_error_msg("Error: Can not open network init file " +
                      filename);
    exit(1);
  }

  std::string line;

  nodes.clear();
  elems.clear();
  node_boundary_flag.clear();
  node_branch_flag.clear();
  elem_radius.clear();
  elem_viscosity.clear();

  //
  // Structure of input file should be as follows
  // 1. Coordinate of nodes should be first, then element-node connectivity
  //, and then element data in order: radius, viscosity, etc
  //
  while (std::getline(filein, line)) {
    std::string tag = line;

    // read node
    if (tag == "$Nodes") {

      // get next line
      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_nodes = 0;
      iss >> num_nodes;

      nodes.resize(num_nodes);

      // start reading nodes
      for (unsigned int i = 0; i < num_nodes; i++) {
        double x, y, z;
        unsigned int id;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> x >> y >> z;

        nodes[id] = Point(x, y, z);
      }

      // read end of node # Nodes end
      std::getline(filein, line);
    } else if (tag == "$Elements") {

      // get next line
      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_elems = 0;
      iss >> num_elems;

      elems.resize(num_elems);

      // start reading nodes
      for (unsigned int i = 0; i < num_elems; i++) {

        unsigned int id, type;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> type;

        // we currently support 2-node line element, however, later we can
        // implement other line elements
        if (type == EDGE2) {

          unsigned int n1, n2;
          iss >> n1 >> n2;

          elems[id] = {n1, n2};
        }
      }

      // read end of node # Elements end
      std::getline(filein, line);
    } else if (tag == "$Nodes Boundary Flags") {

      // get next line
      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_nodes = 0;
      iss >> num_nodes;

      node_boundary_flag.resize(num_nodes);

      // start reading nodes
      for (unsigned int i = 0; i < num_nodes; i++) {
        unsigned int id, flag;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> flag;

        node_boundary_flag[id] = flag;
      }

      // read end of node # Nodes Boundary Flags end
      std::getline(filein, line);
    } else if (tag == "$Nodes Branch Flags") {

      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_nodes = 0;
      iss >> num_nodes;

      node_branch_flag.resize(num_nodes, false);

      // start reading nodes
      for (unsigned int i = 0; i < num_nodes; i++) {
        unsigned int id, flag;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> flag;

        node_branch_flag[id] = flag != 0;
      }

      // read end of node # Nodes Branch Flags end
      std::getline(filein, line);
    } else if (tag == "$Elements Radius") {

      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_elems = 0;
      iss >> num_elems;

      elem_radius.resize(num_elems);

      // start reading nodes
      for (unsigned int i = 0; i < num_elems; i++) {

        unsigned int id;
        double r;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> r;

        elem_radius[id] = r;
      }

      // read end of node # Elements Radius end
      std::getline(filein, line);
    } else if (tag == "$Elements Viscosity") {

      std::getline(filein, line);
      std::istringstream iss(line);

      // get number of nodes first
      unsigned int num_elems = 0;
      iss >> num_elems;

      elem_viscosity.resize(num_elems);

      // start reading nodes
      for (unsigned int i = 0; i < num_elems; i++) {

        unsigned int id;
        double mu;

        std::getline(filein, line);
        std::istringstream iss(line);

        iss >> id >> mu;

        elem_viscosity[id] = mu;
      }

      // read end of node # Elements Viscosity end
      std::getline(filein, line);
    }

    if (filein.eof()) break;
  }

  filein.close();

  return true;
}