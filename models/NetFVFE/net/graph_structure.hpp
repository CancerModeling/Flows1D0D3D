////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
/////////////////////////////////////////////////////////////////////////////


namespace netfvfe {

template<class Node>
class GraphStructure{

std::shared_ptr<Node> head, tail;

public:

GraphStructure(){

   head = tail = NULL;

}

~GraphStructure(){}

bool isEmpty(){

     return (head == NULL) ? true : false;

}

void attachNode(Node newNode){

    auto sp_newNode = std::make_shared<Node>( newNode );

    if( isEmpty() ){
        tail = head = sp_newNode;
    }
    else{
        tail->successor = sp_newNode;
        tail = sp_newNode;
    }

}

};

}
