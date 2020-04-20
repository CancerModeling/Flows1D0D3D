////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
/////////////////////////////////////////////////////////////////////////////

#ifndef NETFC_LISTSTRUCTURE_H
#define NETFC_LISTSTRUCTURE_H

namespace netfc {

template<class Node>
class ListStructure{

std::shared_ptr<Node> head, tail;

int numberOfNodes;

public:

ListStructure(): numberOfNodes(0){

   head = tail = NULL;

}

~ListStructure(){}

bool isEmpty(){

     return (head == NULL) ? true : false;

}

void attachNode(Node newNode){

     auto sp_newNode = std::make_shared<Node>( newNode );

     if( isEmpty() ){

         tail = head = sp_newNode;

     }
     else{

         tail->global_successor = sp_newNode;
         tail = sp_newNode;

    }

    numberOfNodes = numberOfNodes+1;

}

void attachPointerToNode( std::shared_ptr<Node> pointer ){

    if( isEmpty() ){

        tail = head = pointer;

    }
    else{

        tail->global_successor = pointer;
        tail = pointer;

    }

    numberOfNodes = numberOfNodes+1;

}

std::shared_ptr<Node> getHead(){

     return head;

}

const std::shared_ptr<Node> getHead() const{

      return head;

}


std::shared_ptr<Node> getTail(){

    return tail;

}

int getNumberOfNodes(){

    return numberOfNodes;

}

void resetNumberOfNodes(){

    numberOfNodes = 0;

}

void determineNumberOfNodes(){

     std::shared_ptr<Node> pointer = head;

     numberOfNodes = 0;

     while( pointer ){

            numberOfNodes = numberOfNodes+1;

            pointer = pointer->global_successor;

     }

}

void findVertex(int numberOfVertex, std::vector<double>& coord, bool& allVerticesFound){

    std::shared_ptr<Node> pointer = head;

    bool vertexFound = false;

    while( pointer ){

           int index_1 = pointer->index_1;
           int index_2 = pointer->index_2;

           if( index_1 == numberOfVertex ){

               coord = pointer->coord_1;
               vertexFound = true;
               break;

           }
           else if( index_2 == numberOfVertex ){

               coord = pointer->coord_2;
               vertexFound = true;
               break;

           }

           if( vertexFound == true ){

               break;

           }

           pointer = pointer->global_successor;

    }

    if( vertexFound == false ){

        allVerticesFound = true;

    }
    
}

std::shared_ptr<Node> findNode(int indexOfNode){

    std::shared_ptr<Node> pointer = head;

    while( pointer ){

           int index = pointer->index;

           if( index == indexOfNode ){

               return pointer;

           }

           pointer = pointer->global_successor;

    }

}


};

}

#endif // NETFC_LISTSTRUCTURE_H
