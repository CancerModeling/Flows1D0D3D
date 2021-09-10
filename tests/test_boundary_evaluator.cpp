////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "catch2/catch.hpp"
#include <utility>

#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/petsc/petsc_vec.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/edge_boundary_evaluator.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "mpi.h"

namespace mc = macrocirculation;

TEST_CASE("BoundaryEvaluator", "[BoundaryEvaluator]") {
  auto comm = MPI_COMM_WORLD;

  if (mc::mpi::size(comm) != 2)
    throw std::runtime_error("test only supports 2 ranks");

  const size_t num_micro_edges = 4;

  auto physical_data = mc::PhysicalData::set_from_data(0, 0, 0, 0, 0, 1.);

  auto graph_A = std::make_shared<mc::GraphStorage>();
  std::vector< std::shared_ptr< mc::Vertex > > vertices_A;
  for (size_t k = 0; k < 9; k+=1)
    vertices_A.push_back(graph_A->create_vertex());
  std::vector< std::shared_ptr< mc::Edge > > edges_A;
  for (size_t k = 0; k < vertices_A.size()-1; k+=1)
    edges_A.push_back(graph_A->connect(*vertices_A.at(k), *vertices_A.at(k+1), num_micro_edges));
  for (const auto& edge: edges_A)
    edge->add_physical_data(physical_data);

  auto graph_B = std::make_shared<mc::GraphStorage>();
  std::vector< std::shared_ptr< mc::Vertex > > vertices_B;
  for (size_t k = 0; k < 9; k+=1)
    vertices_B.push_back(graph_B->create_vertex());
  std::vector< std::shared_ptr< mc::Edge > > edges_B;
  for (size_t k = 0; k < vertices_B.size()-1; k+=1)
    edges_B.push_back(graph_B->connect(*vertices_B.at(k), *vertices_B.at(k+1), num_micro_edges));
  for (const auto& edge: edges_B)
    edge->add_physical_data(physical_data);

  // connect the graphs
  mc::Vertex::connect(graph_A, *vertices_A.back(), graph_B, *vertices_B.front());

  graph_A->finalize_bcs();
  graph_B->finalize_bcs();

  for (size_t k=0; k<edges_A.size()/2; k += 1)
    graph_A->assign_edge_to_rank(*edges_A[k], 0);
  for (size_t k=edges_A.size()/2; k<edges_A.size(); k += 1)
    graph_A->assign_edge_to_rank(*edges_A[k], 1);

  for (size_t k=0; k<edges_B.size()/2; k += 1)
    graph_B->assign_edge_to_rank(*edges_B[k], 0);
  for (size_t k=edges_B.size()/2; k<edges_B.size(); k += 1)
    graph_B->assign_edge_to_rank(*edges_B[k], 1);

  auto dof_map_A = std::make_shared< mc::DofMap > (*graph_A);
  dof_map_A->create(comm, *graph_A, 1, 0, true);
  auto dof_map_B = std::make_shared< mc::DofMap > (*graph_B);
  dof_map_B->create(comm, *graph_B, 1, 0, true);

  mc::PetscVec vec_A ("vecA", dof_map_A->num_owned_dofs());
  for (size_t k = 0; k < dof_map_A->num_owned_dofs(); k+= 1)
    vec_A.set(k + dof_map_A->first_owned_global_dof(), 1. + k + dof_map_A->first_owned_global_dof());

  mc::PetscVec vec_B ("vecB", dof_map_B->num_owned_dofs());
  for (size_t k = 0; k < dof_map_B->num_owned_dofs(); k+= 1)
    vec_B.set(k + dof_map_B->first_owned_global_dof(), 200. + k + dof_map_B->first_owned_global_dof());

  mc::EdgeBoundaryEvaluator evaluator_A(comm, graph_A, dof_map_A, 0);
  mc::EdgeBoundaryEvaluator evaluator_B(comm, graph_B, dof_map_B, 0);

  evaluator_A.init(vec_A);
  evaluator_B.init(vec_B);

  if (mc::mpi::rank(comm) == 0)
  {
    std::vector< double > values (1, 0);
    evaluator_A(*vertices_A.back(), values);
    REQUIRE( values[0] == Approx(edges_A.size() * num_micro_edges));
  }

  if (mc::mpi::rank(comm) == 1)
  {
    std::vector< double > values (1, 0);
    evaluator_B(*vertices_B.front(), values);
    REQUIRE( values[0] == Approx(200.) );
  }
}
