////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNET_SEG_FV_H
#define UTIL_UNET_SEG_FV_H

#include "list_structure.hpp"
#include "utils.hpp"

// definitions
#define UNET_FREE_MASK 0x0000
#define UNET_PRES_BDRY_DIRIC 0x0001
#define UNET_PRES_BDRY_INNER 0x0002
#define UNET_PRES_INNER 0x0004
#define UNET_NUT_BDRY_ARTERY_INLET 0x0008
#define UNET_NUT_BDRY_VEIN_INLET 0x0010
#define UNET_NUT_BDRY_INNER_INLET 0x0020
#define UNET_NUT_BDRY_ARTERY_OUTLET 0x0040
#define UNET_NUT_BDRY_VEIN_OUTLET 0x0080
#define UNET_NUT_BDRY_INNER_OUTLET 0x0100
#define UNET_NUT_INNER 0x0200

namespace util {

namespace unet {

struct ElemWeights {

  unsigned id_seg;
  double half_cyl_surf;
  std::vector<unsigned int> elem_id;
  std::vector<double> elem_weight;

  ElemWeights(unsigned int id = 0) : id_seg(id), half_cyl_surf(0.){};

  void add_unique(const unsigned int &elem, const double &weight) {

    for (unsigned int i = 0; i < elem_id.size(); i++) {
      if (elem_id[i] == elem) {
        elem_weight[i] += weight;
        return;
      }
    }

    elem_id.push_back(elem);
    elem_weight.push_back(weight);
  }

  void clear() {
    id_seg = 0;
    half_cyl_surf = 0.;
    elem_id.clear();
    elem_weight.clear();
  }
};

enum TypeOfSegment { DirBoundary,
                     Inner };

class SegFV {

public:
  /*!
   * @brief Constructor
   */

  SegFV() : index(0), typeOfSegment(Inner), length(0.0),
            radius(0.0), L_p(0.0), p_boundary_1(0.0), p_boundary_2(0.0),
            p_v(0.0), global_successor(NULL),
            index_1(0), index_2(0), coord_1(0.0), coord_2(0.0) {}

  std::vector<std::shared_ptr<SegFV>> neighbors_1;

  std::vector<std::shared_ptr<SegFV>> neighbors_2;

  int index, index_1, index_2;

  std::vector<double> coord_1, coord_2;

  TypeOfSegment typeOfSegment;

  double length, radius, L_p, p_boundary_1, p_boundary_2, t_seg, p_v, mu;

  std::shared_ptr<SegFV> global_successor;

  double getTransmissibilty() {

    double t = 0.0;

    t = (2.0 * M_PI * radius * radius * radius * radius) / (8.0 * length * mu);

    return t;
  }
};

enum TypeOfNode { DirichletNode,
                  InnerNode,
                  NeumannNode };

class VGNode {

public:
  /*!
   * @brief Constructor
   */
  VGNode() : index(0), p_v(0.0), c_v(0.0), p_boundary(0.0),
             c_boundary(0.0), edge_touched(false), sprouting_edge(false), apicalGrowth(false),
             coord(0.0), radii(0.0), radii_initial(0.0), tau_w_initial(0.0), L_p(0.0), notUpdated(0), is_initial_node(false), is_sprouting_node(false),
             node_marked(false), is_given(false) {}

  int index, notUpdated;

  bool apicalGrowth, is_sprouting_node, is_given;

  double p_v, c_v, p_boundary, c_boundary;

  std::vector<double> coord, radii, L_p, L_s, radii_initial, tau_w_initial;

  std::vector<bool> edge_touched;

  std::vector<bool> sprouting_edge;

  std::vector<std::shared_ptr<VGNode>> neighbors;

  TypeOfNode typeOfVGNode;

  std::shared_ptr<VGNode> global_successor, global_predecessor;

  std::vector<ElemWeights> J_b_points;

  /*! @brief Marks nodes which have been there from the beginning. */
  bool is_initial_node;

  /*! @brief Generic attribute to mark the node to simplify algorithms on the graph. */
  bool node_marked;

  /*! @brief Marks the edge belonging to the neighbor node with the given global index. */
  void markEdge(int globalIndex) {
    edge_touched[getLocalNeighborIndex(globalIndex)] = true;
  }

  /*! @brief Returns the local index of the neighbor with the given global index. */
  int getLocalNeighborIndex(int globalIndex) {
    for (int i = 0; i < neighbors.size(); i++) {
      if (globalIndex == neighbors[i]->index) {
        return i;
      }
    }
    throw std::runtime_error("the node with index " + std::to_string(index) + " has no neighbor with index " + std::to_string(globalIndex));
  }

  void markEdgeLocalIndex(int localIndex) {

    edge_touched[localIndex] = true;
  }

  /*! @brief Returns the local index of the given neighbor, and throws an exception if it is not present. */
  int getLocalIndexOfNeighbor(std::shared_ptr<VGNode> neighbor) {
    for (int i = 0; i < neighbors.size(); i++) {
      if (util::dist_between_points(neighbor->coord, neighbors[i]->coord) < 1.0e-12)
        return i;
    }
    throw std::runtime_error("could not find neighbor");
  }

  void replacePointerWithIndex(int index_new, std::shared_ptr<VGNode> new_pointer) {

    int numberOfNeighbors = neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (index_new == neighbors[i]->index) {

        neighbors[i] = new_pointer;

        edge_touched[i] = true;
      }
    }
  }


  void attachNeighbor(std::shared_ptr<VGNode> new_pointer) {

    neighbors.push_back(new_pointer);

    typeOfVGNode = InnerNode;
  }

  double getTotalVolume() {

    double totalVolume = 0.0;

    int numberOfNeighbors = neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      std::vector<double> coord_neighbor = neighbors[i]->coord;

      double length = 0.0;

      for (int j = 0; j < 3; j++) {

        length += (coord[j] - coord_neighbor[j]) * (coord[j] - coord_neighbor[j]);
      }

      length = std::sqrt(length);

      totalVolume = totalVolume + (M_PI * radii[i] * radii[i] * length / 2.0);
    }

    return totalVolume;
  }

  void markEdgeForSprouting(int edgeNumber) {

    sprouting_edge[edgeNumber] = true;
  }

  void markNodeForApicalGrowth() {

    apicalGrowth = true;
  }

  void printInformationOfNode() {

    std::cout << " " << std::endl;
    std::cout << "Data of Node: " << index << std::endl;
    std::cout << "notUpdated: " << notUpdated << std::endl;
    std::cout << "apicalGrowth: " << apicalGrowth << std::endl;
    std::cout << "p_v: " << p_v << std::endl;
    std::cout << "c_v: " << c_v << std::endl;
    std::cout << "p_boundary: " << p_boundary << std::endl;
    std::cout << "c_boundary: " << c_boundary << std::endl;
    std::cout << "coord: " << coord << std::endl;
    std::cout << "radii: " << radii << std::endl;
    std::cout << "L_p: " << L_p << std::endl;
    std::cout << "L_s: " << L_s << std::endl;
    std::cout << "edge_touched: " << edge_touched << std::endl;
    std::cout << "sprouting_edge: " << sprouting_edge << std::endl;
    std::cout << "typeOfVGNode: " << typeOfVGNode << std::endl;
    std::cout << "radii_initial: " << radii_initial << std::endl;
    std::cout << "tau_w_initial: " << tau_w_initial << std::endl;
    std::cout << "is_sprouting_node: " << is_sprouting_node << std::endl;

    if (global_successor) {

      std::cout << "global_successor: " << global_successor->index << std::endl;
    }

    if (global_predecessor) {

      std::cout << "global_predecessor: " << global_predecessor->index << std::endl;
    }

    int numberOfNeighbors = neighbors.size();
    std::cout << "numberOfNeighbors: " << numberOfNeighbors << std::endl;

    std::cout << " " << std::endl;

    for (int i = 0; i < numberOfNeighbors; i++) {

      std::cout << "index_neighbor: " << neighbors[i]->index << std::endl;
      std::cout << "coord_neighbor: " << neighbors[i]->coord << std::endl;
      std::cout << "radius: " << radii[i] << std::endl;
    }
  }

  void remove(const std::function<bool(const VGNode &)> &predicate) {
    std::vector<double> radii_new, L_p_new, L_s_new, radii_initial_new, tau_w_initial_new;
    std::vector<bool> edge_touched_new, sprouting_edge_new;
    std::vector<std::shared_ptr<VGNode>> neighbors_new;

    for (int i = 0; i < neighbors.size(); i++) {
      if (!predicate(*neighbors[i])) {
        radii_new.push_back(radii[i]);
        L_p_new.push_back(L_p[i]);
        L_s_new.push_back(L_s[i]);
        radii_initial_new.push_back(radii_initial[i]);
        tau_w_initial_new.push_back(tau_w_initial[i]);
        edge_touched_new.push_back(edge_touched[i]);
        sprouting_edge_new.push_back(sprouting_edge[i]);
        neighbors_new.push_back(neighbors[i]);
      }
    }

    radii = radii_new;
    L_p = L_p_new;
    L_s = L_s_new;
    radii_initial = radii_initial_new;
    tau_w_initial = tau_w_initial_new;
    edge_touched = edge_touched_new;
    sprouting_edge = sprouting_edge_new;
    neighbors = neighbors_new;
  }

  void removeComponent(int comp) {

    int numberOfComponents = neighbors.size();

    std::vector<double> radii_new, L_p_new, L_s_new, radii_initial_new, tau_w_initial_new;

    std::vector<bool> edge_touched_new, sprouting_edge_new;

    std::vector<std::shared_ptr<VGNode>> neighbors_new;

    for (int i = 0; i < numberOfComponents; i++) {

      if (i != comp) {

        radii_new.push_back(radii[i]);
        L_p_new.push_back(L_p[i]);
        L_s_new.push_back(L_s[i]);
        radii_initial_new.push_back(radii_initial[i]);
        tau_w_initial_new.push_back(tau_w_initial[i]);
        edge_touched_new.push_back(edge_touched[i]);
        sprouting_edge_new.push_back(sprouting_edge[i]);
        neighbors_new.push_back(neighbors[i]);
      }
    }

    radii = radii_new;
    L_p = L_p_new;
    L_s = L_s_new;
    radii_initial = radii_initial_new;
    tau_w_initial = tau_w_initial_new;
    edge_touched = edge_touched_new;
    sprouting_edge = sprouting_edge_new;
    neighbors = neighbors_new;
  }

  void removePointer(int comp) {

    int numberOfComponents = neighbors.size();

    std::vector<std::shared_ptr<VGNode>> neighbors_new;

    for (int i = 0; i < numberOfComponents; i++) {

      if (i != comp) {

        neighbors_new.push_back(neighbors[i]);
      }
    }

    neighbors = neighbors_new;
  }

  void removeComponents(std::vector<int> components) {

    int numberOfComponents = neighbors.size();

    int numberOfCompToRemove = components.size();

    std::vector<double> radii_new, L_p_new, L_s_new, radii_initial_new, tau_w_initial_new;

    std::vector<bool> edge_touched_new, sprouting_edge_new;

    std::vector<std::shared_ptr<VGNode>> neighbors_new;

    for (int i = 0; i < numberOfComponents; i++) {

      bool removeComponent = false;

      for (int j = 0; j < numberOfCompToRemove; j++) {

        if (i == components[j]) {

          removeComponent = true;
        }
      }

      if (!removeComponent) {

        radii_new.push_back(radii[i]);
        L_p_new.push_back(L_p[i]);
        L_s_new.push_back(L_s[i]);
        radii_initial_new.push_back(radii_initial[i]);
        tau_w_initial_new.push_back(tau_w_initial[i]);
        edge_touched_new.push_back(edge_touched[i]);
        sprouting_edge_new.push_back(sprouting_edge[i]);
        neighbors_new.push_back(neighbors[i]);
      }
    }

    radii = radii_new;
    L_p = L_p_new;
    L_s = L_s_new;
    radii_initial = radii_initial_new;
    tau_w_initial = tau_w_initial_new;
    edge_touched = edge_touched_new;
    sprouting_edge = sprouting_edge_new;
    neighbors = neighbors_new;
  }

  void replacePointerWithLocalIndex(int index_new, std::shared_ptr<VGNode> new_pointer) {

    neighbors[index_new] = new_pointer;
  }

  void replacePointerWithGlobalIndex(int index_new, std::shared_ptr<VGNode> new_pointer) {

    int numberOfNeighbors = neighbors.size();

    for (int i = 0; i < numberOfNeighbors; i++) {

      if (index_new == neighbors[i]->index) {

        neighbors[i] = new_pointer;

        edge_touched[i] = true;
      }
    }
  }
};

} // namespace unet

} // namespace util

#endif // UTIL_UNET_SEG_FV_H
