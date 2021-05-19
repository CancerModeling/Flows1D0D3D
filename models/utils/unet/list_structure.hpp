////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
/////////////////////////////////////////////////////////////////////////////

#ifndef UTIL_UNET_LISTSTRUCTURE_H
#define UTIL_UNET_LISTSTRUCTURE_H

#include <boost/shared_ptr.hpp>

namespace util {

namespace unet {

template<class Node>
class ListStructure {

  std::shared_ptr<Node> head, tail;

  int numberOfNodes;

public:
  ListStructure() : numberOfNodes(0) {

    head = tail = NULL;
  }

  ~ListStructure() = default;

  bool isEmpty() {

    return (head == NULL);
  }

  /*! Removes the given node from list by splicing it out. */
  void remove(std::shared_ptr<Node> pointer) {
    auto predecessor = pointer->global_predecessor;
    auto successor = pointer->global_successor;

    // sanity check for predecessor:
    if (predecessor){
      if( predecessor->global_successor.get() != pointer.get() )
        throw std::runtime_error("ListStructure corrupted! Did you try to remove the same node twice?");
    } else {
      if ( head.get() != pointer.get() )
        throw std::runtime_error("ListStructure corrupted! Did you try to remove the same node twice?");
    }
    // sanity check for successor:
    if (successor){
      if( successor->global_predecessor.get() != pointer.get() )
        throw std::runtime_error("ListStructure corrupted! Did you try to remove the same node twice?");
    } else {
      if ( tail.get() != pointer.get() )
        throw std::runtime_error("ListStructure corrupted! Did you try to remove the same node twice?");
    }

    // decrement the total number of nodes in the list
    numberOfNodes -= 1;

    // case: we were the last node in the graph
    if (predecessor == nullptr && successor == nullptr) {
      setHead(nullptr);
      setTail(nullptr);
    }

    if (predecessor != nullptr) {
      predecessor->global_successor = successor;
    }
    // case: we were the first node
    else {
      setHead(successor);
    }

    if (successor != nullptr) {
      successor->global_predecessor = predecessor;
    }
    // case: we were the last node
    else {
      setTail(predecessor);
    }

    pointer->global_predecessor = nullptr;
    pointer->global_successor = nullptr;
  }

  void attachNode(Node newNode) {

    auto sp_newNode = std::make_shared<Node>(newNode);

    if (isEmpty()) {

      tail = head = sp_newNode;

    } else {

      tail->global_successor = sp_newNode;
      sp_newNode->global_predecessor = tail;
      tail = sp_newNode;
    }

    numberOfNodes = numberOfNodes + 1;
  }

  void attachPointerToNode(std::shared_ptr<Node> pointer) {

    if (isEmpty()) {

      tail = head = pointer;

    } else {

      tail->global_successor = pointer;
      pointer->global_predecessor = tail;
      tail = pointer;
    }

    numberOfNodes = numberOfNodes + 1;
  }

  std::shared_ptr<Node> getHead() {
    return head;
  }

  std::shared_ptr<const Node> getHead() const {
    return head;
  }

  std::shared_ptr<Node> getTail() {
    return tail;
  }

  int getNumberOfNodes() {
    return numberOfNodes;
  }

  void determineNumberOfNodes() {
    std::shared_ptr<Node> pointer = head;

    auto oldNumberOfNodes = numberOfNodes;

    numberOfNodes = 0;

    while (pointer) {
      numberOfNodes = numberOfNodes + 1;
      pointer = pointer->global_successor;
    }

    if (oldNumberOfNodes != numberOfNodes)
      throw std::runtime_error("ListStructure corrupted!");
  }

  void findVertex(int numberOfVertex, std::vector<double> &coord, bool &allVerticesFound) {

    std::shared_ptr<Node> pointer = head;

    bool vertexFound = false;

    while (pointer) {

      int index_1 = pointer->index_1;
      int index_2 = pointer->index_2;

      if (index_1 == numberOfVertex) {

        coord = pointer->coord_1;
        vertexFound = true;
        break;

      } else if (index_2 == numberOfVertex) {

        coord = pointer->coord_2;
        vertexFound = true;
        break;
      }

      if (vertexFound == true) {

        break;
      }

      pointer = pointer->global_successor;
    }

    if (vertexFound == false) {

      allVerticesFound = true;
    }
  }

  std::shared_ptr<Node> findNode(int indexOfNode) {

    std::shared_ptr<Node> pointer = head;

    while (pointer) {

      int index = pointer->index;

      if (index == indexOfNode) {

        return pointer;
      }

      pointer = pointer->global_successor;
    }

    throw std::runtime_error("could not find node with index " + std::to_string(indexOfNode));
  }

private:
  void setHead(const std::shared_ptr<Node> &_head) {
    head = _head;
  }

  void setTail(const std::shared_ptr<Node> &_tail) {
    tail = _tail;
  }
};

} // namespace unet

} // namespace util

#endif // UTIL_UNET_LISTSTRUCTURE_H
