////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "implicit_transport_solver.hpp"

#include <utility>

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "explicit_nonlinear_flow_solver.hpp"
#include "fe_type.hpp"
#include "graph_storage.hpp"
#include "implicit_linear_flow_solver.hpp"
#include "linearized_flow_upwind_evaluator.hpp"
#include "nonlinear_flow_upwind_evaluator.hpp"
#include "petsc/petsc_ksp.hpp"
#include "petsc/petsc_vec.hpp"
#include "petsc_assembly_blocks.hpp"

namespace macrocirculation {

ConstantUpwindProvider::ConstantUpwindProvider(double speed)
    : d_speed(speed) {}

ConstantUpwindProvider::~ConstantUpwindProvider() = default;

void ConstantUpwindProvider::get_values_at_qp(double t,
                                              const Edge &edge,
                                              size_t micro_edge,
                                              const QuadratureFormula &qf,
                                              std::vector<double> &v_qp) const {
  assert(qf.size() == v_qp.size());

  std::uniform_real_distribution<double> distribution(0.5, 1.);

  for (size_t k = 0; k < qf.size(); k += 1)
    v_qp[k] = d_speed;
}

/*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
void ConstantUpwindProvider::get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const {
  assert(v_qp.size() == edge.num_micro_vertices());

  std::uniform_real_distribution<double> distribution(0.5, 1.);

  for (size_t k = 0; k < edge.num_micro_vertices(); k += 1)
    v_qp[k] = d_speed;
}

void ConstantUpwindProvider::get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const {
  assert(v.get_edge_neighbors().size() == A.size());
  assert(v.get_edge_neighbors().size() == Q.size());

  for (size_t k = 0; k < A.size(); k += 1) {
    A[k] = 1.;
    Q[k] = d_speed;
  }
}

UpwindProviderNonlinearFlow::UpwindProviderNonlinearFlow(std::shared_ptr<NonlinearFlowUpwindEvaluator> evaluator, std::shared_ptr<ExplicitNonlinearFlowSolver> solver)
    : d_evaluator(std::move(evaluator)),
      d_solver(std::move(solver)) {}

void UpwindProviderNonlinearFlow::init(double t, const std::vector<double> &u) {
  d_evaluator->init(t, u);
}

void UpwindProviderNonlinearFlow::get_values_at_qp(double t,
                                                   const Edge &edge,
                                                   size_t micro_edge,
                                                   const QuadratureFormula &qf,
                                                   std::vector<double> &v_qp) const {
  assert(v_qp.size() == qf.size());

  FETypeNetwork fe(qf, d_solver->get_degree());
  auto &ldof_map = d_solver->get_dof_map().get_local_dof_map(edge);
  std::vector<size_t> dof_indices(ldof_map.num_basis_functions());
  std::vector<double> dof_values(ldof_map.num_basis_functions());

  std::vector<double> values_A(qf.size());
  std::vector<double> values_Q(qf.size());

  ldof_map.dof_indices(micro_edge, d_solver->A_component, dof_indices);
  extract_dof(dof_indices, d_solver->get_solution(), dof_values);
  fe.evaluate_dof_at_quadrature_points(dof_values, values_A);

  ldof_map.dof_indices(micro_edge, d_solver->Q_component, dof_indices);
  extract_dof(dof_indices, d_solver->get_solution(), dof_values);
  fe.evaluate_dof_at_quadrature_points(dof_values, values_Q);

  for (size_t k = 0; k < qf.size(); k += 1)
    v_qp[k] = values_Q[k] / values_A[k];
}

/*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
void UpwindProviderNonlinearFlow::get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const {
  std::vector<double> Q_up(v_qp.size());
  std::vector<double> A_up(v_qp.size());
  d_evaluator->get_fluxes_on_macro_edge(t, edge, d_solver->get_solution(), Q_up, A_up);
  for (size_t k = 0; k < v_qp.size(); k += 1)
    v_qp[k] = Q_up[k] / A_up[k];
}

void UpwindProviderNonlinearFlow::get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const {
  d_evaluator->get_fluxes_on_nfurcation(t, v, Q, A);
}

UpwindProviderLinearizedFlow::UpwindProviderLinearizedFlow(std::shared_ptr<GraphStorage> graph, std::shared_ptr<LinearizedFlowUpwindEvaluator> evaluator, std::shared_ptr<ImplicitLinearFlowSolver> solver)
    : d_graph(std::move(graph)),
      d_evaluator(std::move(evaluator)),
      d_solver(std::move(solver)) {}

void UpwindProviderLinearizedFlow::init(double t, const std::vector<double> &u) {
  d_evaluator->init(t, u);
}

void UpwindProviderLinearizedFlow::init(double t, const PetscVec &u) {
  d_evaluator->init(t, u);
}

void UpwindProviderLinearizedFlow::get_values_at_qp(double t,
                                                    const Edge &edge,
                                                    size_t micro_edge,
                                                    const QuadratureFormula &qf,
                                                    std::vector<double> &v_qp) const {
  assert(v_qp.size() == qf.size());

  FETypeNetwork fe(qf, d_solver->get_degree());
  auto &ldof_map = d_solver->get_dof_map().get_local_dof_map(edge);
  std::vector<size_t> dof_indices(ldof_map.num_basis_functions());
  std::vector<double> dof_values(ldof_map.num_basis_functions());

  std::vector<double> values_q(qf.size());

  ldof_map.dof_indices(micro_edge, d_solver->q_component, dof_indices);
  extract_dof(dof_indices, d_solver->get_solution(), dof_values);
  fe.evaluate_dof_at_quadrature_points(dof_values, values_q);

  auto &param = edge.get_physical_data();

  for (size_t k = 0; k < qf.size(); k += 1)
    v_qp[k] = values_q[k] / param.A0;
}

/*! @brief Returns the upwinded values for Q and A for a whole macro-edge at the micro-edge boundaries. */
void UpwindProviderLinearizedFlow::get_upwinded_values(double t, const Edge &edge, std::vector<double> &v_qp) const {
  assert(v_qp.size() == edge.num_micro_vertices());
  std::vector<double> p_up(edge.num_micro_vertices());
  std::vector<double> q_up(edge.num_micro_vertices());
  d_evaluator->get_fluxes_on_macro_edge(t, edge, d_solver->get_solution(), p_up, q_up);
  assert(edge.has_physical_data());
  auto A0 = edge.get_physical_data().A0;
  for (size_t k = 0; k < v_qp.size(); k += 1)
    v_qp[k] = q_up[k] / A0;
}

void UpwindProviderLinearizedFlow::get_upwinded_values(double t, const Vertex &v, std::vector<double> &A, std::vector<double> &Q) const {
  std::vector<double> p_up(v.get_edge_neighbors().size());
  d_evaluator->get_fluxes_on_nfurcation(t, v, p_up, Q);

  std::vector<double> A0;
  for (size_t k = 0; k < v.get_edge_neighbors().size(); k += 1) {
    auto &edge = *d_graph->get_edge(v.get_edge_neighbors()[k]);
    assert(edge.has_physical_data());
    A[k] = edge.get_physical_data().A0;
  }
}

ImplicitTransportSolver::ImplicitTransportSolver(MPI_Comm comm,
                                                 std::shared_ptr<GraphStorage> graph,
                                                 std::shared_ptr<DofMap> dof_map,
                                                 std::shared_ptr<UpwindProvider> upwind_provider,
                                                 size_t degree)
    : ImplicitTransportSolver(comm,
                              std::vector<std::shared_ptr<GraphStorage>>({std::move(graph)}),
                              std::vector<std::shared_ptr<DofMap>>({std::move(dof_map)}),
                              std::vector<std::shared_ptr<UpwindProvider>>({std::move(upwind_provider)}), degree) {}


ImplicitTransportSolver::ImplicitTransportSolver(MPI_Comm comm,
                                                 std::vector<std::shared_ptr<GraphStorage>> graph,
                                                 std::vector<std::shared_ptr<DofMap>> dof_map,
                                                 std::vector<std::shared_ptr<UpwindProvider>> upwind_provider,
                                                 size_t degree)
    : d_comm(comm),
      d_graph(std::move(graph)),
      d_dof_map(std::move(dof_map)),
      d_upwind_provider(std::move(upwind_provider)),
      d_degree(degree),
      u(std::make_shared<PetscVec>("u", d_dof_map)),
      rhs(std::make_shared<PetscVec>("rhs", d_dof_map)),
      A(std::make_shared<PetscMat>("A", d_dof_map)),
      mass(std::make_shared<PetscVec>("mass", d_dof_map)),
      linear_solver(PetscKsp::create_with_pc_ilu(*A)) {

  if (d_graph.size() != d_dof_map.size())
    throw std::runtime_error("ImplicitTransportSolver::ImplicitTransportSolver: number of graphs and dof maps must coincide");
  if (d_graph.size() != d_upwind_provider.size())
    throw std::runtime_error("ImplicitTransportSolver::ImplicitTransportSolver: number of graphs and upwind providers must coincide");

  // TODO: preallocate the nonzeros properly!
  MatSetOption(A->get_mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  assemble_mass(d_comm, d_graph, d_dof_map, *mass);
  u->zero();
  rhs->zero();
  // initialize zero pattern of matrix with velocity 1 vector:
  // otherwise in some vessels the velocity might be zero and thus the sparsity pattern would change
  sparsity_pattern();
}

void ImplicitTransportSolver::solve(double tau, double t) {
  assemble(tau, t);
  linear_solver->solve(*rhs, *u);
}

void ImplicitTransportSolver::assemble(double tau, double t) {
  assemble_matrix(tau, t);
  assemble_rhs(tau, t);
}

void ImplicitTransportSolver::sparsity_pattern() {
  ConstantUpwindProvider upwind_provider_plus(1.);
  ConstantUpwindProvider upwind_provider_minus(-42.);
  const double tau = 1e-3;
  const double t = 1.;
  A->zero();
  for (size_t k = 0; k < d_graph.size(); k += 1) {
    implicit_transport::additively_assemble_matrix(d_comm, tau, t, upwind_provider_plus, *d_dof_map[k], *d_graph[k], *A);
    implicit_transport::additively_assemble_matrix(d_comm, tau, t, upwind_provider_minus, *d_dof_map[k], *d_graph[k], *A);
  }
  A->assemble();
}

void ImplicitTransportSolver::assemble_matrix(double tau, double t) {
  A->zero();
  for (size_t k = 0; k < d_graph.size(); k += 1) {
    implicit_transport::additively_assemble_matrix(d_comm, tau, t, *d_upwind_provider[k], *d_dof_map[k], *d_graph[k], *A);
  }
  assemble_characteristics(tau, t);
  A->assemble();
}

void ImplicitTransportSolver::assemble_characteristics(double tau, double t) {
  for (size_t graph_index = 0; graph_index < d_graph.size(); graph_index += 1) {
    auto graph = d_graph[graph_index];

    for (const auto &v_id : graph->get_active_vertex_ids(mpi::rank(d_comm))) {
      auto &vertex = *graph->get_vertex(v_id);

      // only vertices which are coupled
      if (!(vertex.is_nonlinear_characteristic_inflow() || vertex.is_linear_characteristic_inflow()))
        continue;

      auto edge = graph->get_edge(vertex.get_edge_neighbors()[0]);

      auto conn = vertex.get_inter_graph_connections();
      auto &adjacent_graph = conn[0].get_graph();
      auto &adjacent_vertex = conn[0].get_vertex();
      auto adjacent_edge = adjacent_graph.get_edge(adjacent_vertex.get_edge_neighbors()[0]);

      auto adjacent_graph_index = -1;
      for (size_t j = 0; j < d_graph.size(); j += 1)
        if (d_graph[j].get() == &adjacent_graph)
          adjacent_graph_index = j;

      // we only assemble once with the smaller index!
      if (adjacent_graph_index < graph_index)
        continue;

      std::vector<Edge const *> edges = {edge.get(), adjacent_edge.get()};

      std::vector<double> sigma = {
        edge->is_pointing_to(vertex.get_id()) ? +1. : -1.,
        adjacent_edge->is_pointing_to(adjacent_vertex.get_id()) ? +1. : -1.};

      // get upwinded values
      std::vector<double> A_up;
      std::vector<double> Q_up;
      std::vector<double> A_up_value(1), Q_up_value(1);
      d_upwind_provider[graph_index]->get_upwinded_values(t, vertex, A_up_value, Q_up_value);
      A_up.push_back(A_up_value[0]);
      Q_up.push_back(Q_up_value[0]);
      d_upwind_provider[adjacent_graph_index]->get_upwinded_values(t, adjacent_vertex, A_up_value, Q_up_value);
      A_up.push_back(A_up_value[0]);
      Q_up.push_back(Q_up_value[0]);

      // get boundary point type
      using BPT = BoundaryPointType;
      std::vector<BPT> boundary_type;
      boundary_type.push_back(edge->is_pointing_to(vertex.get_id()) ? BPT::Right : BPT::Left);
      boundary_type.push_back(adjacent_edge->is_pointing_to(adjacent_vertex.get_id()) ? BPT::Right : BPT::Left);

      std::vector<std::vector<size_t>> dof_indices_list;
      {
        auto &local_dof_map = d_dof_map[graph_index]->get_local_dof_map(*edge);
        auto micro_edge_id = edge->is_pointing_to(vertex.get_id()) ? edge->num_micro_edges() - 1 : 0;
        std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
        local_dof_map.dof_indices(micro_edge_id, 0, dof_indices);
        dof_indices_list.push_back(dof_indices);
      }
      {
        auto &local_dof_map = d_dof_map[adjacent_graph_index]->get_local_dof_map(*adjacent_edge);
        auto micro_edge_id = adjacent_edge->is_pointing_to(adjacent_vertex.get_id()) ? adjacent_edge->num_micro_edges() - 1 : 0;
        std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
        local_dof_map.dof_indices(micro_edge_id, 0, dof_indices);
        dof_indices_list.push_back(dof_indices);
      }

      implicit_transport::assembly_kernel_nfurcation(d_comm, tau, dof_indices_list, edges, sigma, Q_up, A_up, boundary_type, *A);
    }
  }
}

void ImplicitTransportSolver::assemble_rhs(double tau, double t) {
  rhs->zero();
  implicit_transport::additively_assemble_rhs_cells(*mass, *u, *rhs);
  for (size_t k = 0; k < d_graph.size(); k += 1) {
    implicit_transport::additively_assemble_rhs_inflow(d_comm, tau, t, *d_upwind_provider[k], *d_dof_map[k], *d_graph[k], inflow_function, *rhs);
  }
  rhs->assemble();
}

Eigen::MatrixXd create_QA_phi_grad_psi(const FETypeNetwork &fe,
                                       const LocalEdgeDofMap &local_dof_map,
                                       const std::vector<double> &v_qp) {
  const auto &phi = fe.get_phi();
  const auto &dphi = fe.get_dphi();
  const auto &JxW = fe.get_JxW();

  Eigen::MatrixXd k_loc(local_dof_map.num_basis_functions(), local_dof_map.num_basis_functions());
  for (int j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
    for (int i = 0; i < local_dof_map.num_basis_functions(); i += 1) {
      k_loc(j, i) = 0;
      for (int qp = 0; qp < phi[i].size(); qp += 1)
        k_loc(j, i) += v_qp[qp] * phi[i][qp] * dphi[j][qp] * JxW[qp];
    }
  }
  return k_loc;
}


void ImplicitTransportSolver::applySlopeLimiter(){
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);
    const auto &param = macro_edge->get_physical_data();
    
    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions());    
    std::vector<std::size_t> dof_indices_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_right(local_dof_map.num_basis_functions());
      
    for (size_t micro_vertex_id = 0; micro_vertex_id < macro_edge->num_micro_vertices(); micro_vertex_id += 1) {
    
      auto edge_id = micro_vertex_id;
      
      local_dof_map.dof_indices(edge_id, 0, dof_indices); 
      std::vector<double> dof_edge(local_dof_map.num_basis_functions());
      extract_dof(dof_indices, *u, dof_edge);
	     	
      if( micro_vertex_id == 0 ){                  
          auto &vertex = *d_graph->get_vertex(micro_vertex_id);          
          auto right_edge_id = micro_vertex_id + 1;        
          local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);        
          std::vector<double> dof_right_edge(local_dof_map.num_basis_functions());
          extract_dof(dof_indices_right, *u, dof_right_edge); 
          
	  for( int j=dof_indices.size()-1;j>0;j--){
	    double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	    double diff_left  = gamma * ( dof_edge[ j-1 ] - 0.0 ); 
	    double diff_right = gamma * ( dof_right_edge[ j-1 ] - dof_edge[ j-1 ] ); 
	    double dof_central = dof_edge[ j ];
	    double new_dof = 0.0;

	    if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	    }
		  
	    if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	    }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);  
          } 
             
      }
      else if( micro_vertex_id == macro_edge->num_micro_vertices() - 1 ){    
      
          auto left_edge_id  = micro_vertex_id - 1;        
          local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);        
          std::vector<double> dof_left_edge(local_dof_map.num_basis_functions());
          extract_dof(dof_indices_left, *u, dof_left_edge);                     

          for(int j=dof_indices.size()-1;j>0;j--){	      
	      double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	      double diff_left  = gamma * ( dof_edge[ j-1 ] - dof_left_edge[ j-1 ] ); 
	      double diff_right = gamma * ( 0.0 - dof_edge[ j-1 ] ); 
	      double dof_central = dof_edge[ j ];
	      double new_dof = 0.0;

	      if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		  new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	      }
		  
	      if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	      }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);  		         
          }        
      
      }
      else{ 

        auto left_edge_id  = micro_vertex_id - 1;
        auto right_edge_id = micro_vertex_id + 1;
      
	local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);
	local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);
	     
	std::vector<double> dof_left_edge(local_dof_map.num_basis_functions());
	std::vector<double> dof_right_edge(local_dof_map.num_basis_functions());  
	
	extract_dof(dof_indices_left, *u, dof_left_edge);
	extract_dof(dof_indices_right, *u, dof_right_edge);      
	      
	for(int j=dof_indices.size()-1;j>0;j--){	      
	    double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	    double diff_left  = gamma * ( dof_edge[ j-1 ] - dof_left_edge[ j-1 ] ); 
	    double diff_right = gamma * ( dof_right_edge[ j-1 ] - dof_edge[ j-1 ] ); 
	    double dof_central = dof_edge[ j ];
	    double new_dof = 0.0;

	    if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	    }
		  
	    if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	    }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);		         
        } 	      
      }  
    }      
  }
}


namespace implicit_transport {

void additively_assemble_rhs_cells(PetscVec &mass, PetscVec &u, PetscVec &rhs) {
  CHKERRABORT(PETSC_COMM_WORLD, VecPointwiseMult(rhs.get_vec(), u.get_vec(), mass.get_vec()));
}

void additively_assemble_matrix(MPI_Comm comm,
                                double tau,
                                double t,
                                const UpwindProvider &upwind_provider,
                                const DofMap &dof_map,
                                const GraphStorage &graph,
                                PetscMat &A) {
  additively_assemble_matrix_cells(comm, tau, t, upwind_provider, dof_map, graph, A);
  additively_assemble_matrix_inner_boundaries(comm, tau, t, upwind_provider, dof_map, graph, A);
  additively_assemble_matrix_nfurcations(comm, tau, t, upwind_provider, dof_map, graph, A);
  additively_assemble_matrix_outflow(comm, tau, t, upwind_provider, dof_map, graph, A);
}

void additively_assemble_matrix_cells(MPI_Comm comm,
                                      double tau,
                                      double t,
                                      const UpwindProvider &upwind_provider,
                                      const DofMap &dof_map,
                                      const GraphStorage &graph,
                                      PetscMat &A) {

  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto macro_edge = graph.get_edge(e_id);

    const auto &local_dof_map = dof_map.get_local_dof_map(*macro_edge);
    const auto &param = macro_edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    std::vector<std::size_t> dof_indices_gamma(local_dof_map.num_basis_functions());

    const auto qf = create_gauss4();
    FETypeNetwork fe(qf, local_dof_map.num_basis_functions() - 1);
    fe.reinit(h);

    std::vector<double> v_qp(qf.size());

    auto m_loc{create_mass(fe, local_dof_map)};

    for (const auto &edge : macro_edge->micro_edges()) {
      local_dof_map.dof_indices(edge, 0, dof_indices_gamma);

      upwind_provider.get_values_at_qp(t, *macro_edge, edge.get_local_id(), fe.get_quadrature_formula(), v_qp);

      auto k_loc{create_QA_phi_grad_psi(fe, local_dof_map, v_qp)};

      Eigen::MatrixXd mat = m_loc - tau * k_loc;
      A.add(dof_indices_gamma, dof_indices_gamma, mat);
    }
  }
}

void applySlopeLimiter(){
  for (const auto &e_id : d_graph->get_active_edge_ids(mpi::rank(d_comm))) {
    const auto macro_edge = d_graph->get_edge(e_id);

    const auto &local_dof_map = d_dof_map->get_local_dof_map(*macro_edge);
    const auto &param = macro_edge->get_physical_data();
    
    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions());    
    std::vector<std::size_t> dof_indices_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_right(local_dof_map.num_basis_functions());
      
    for (size_t micro_vertex_id = 0; micro_vertex_id < macro_edge->num_micro_vertices(); micro_vertex_id += 1) {
    
      auto edge_id = micro_vertex_id;
      
      local_dof_map.dof_indices(edge_id, 0, dof_indices); 
      std::vector<double> dof_edge(local_dof_map.num_basis_functions());
      extract_dof(dof_indices, *u, dof_edge);
	     	
      if( micro_vertex_id == 0 ){                  
          auto &vertex = *d_graph->get_vertex(micro_vertex_id);          
          auto right_edge_id = micro_vertex_id + 1;        
          local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);        
          std::vector<double> dof_right_edge(local_dof_map.num_basis_functions());
          extract_dof(dof_indices_right, *u, dof_right_edge); 
          
	  for( int j=dof_indices.size()-1;j>0;j--){
	    double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	    double diff_left  = gamma * ( dof_edge[ j-1 ] - 0.0 ); 
	    double diff_right = gamma * ( dof_right_edge[ j-1 ] - dof_edge[ j-1 ] ); 
	    double dof_central = dof_edge[ j ];
	    double new_dof = 0.0;

	    if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	    }
		  
	    if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	    }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);  
          } 
             
      }
      else if( micro_vertex_id == macro_edge->num_micro_vertices() - 1 ){    
      
          auto left_edge_id  = micro_vertex_id - 1;        
          local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);        
          std::vector<double> dof_left_edge(local_dof_map.num_basis_functions());
          extract_dof(dof_indices_left, *u, dof_left_edge);                     

          for(int j=dof_indices.size()-1;j>0;j--){	      
	      double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	      double diff_left  = gamma * ( dof_edge[ j-1 ] - dof_left_edge[ j-1 ] ); 
	      double diff_right = gamma * ( 0.0 - dof_edge[ j-1 ] ); 
	      double dof_central = dof_edge[ j ];
	      double new_dof = 0.0;

	      if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		  new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	      }
		  
	      if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	      }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);  		         
          }        
      
      }
      else{ 

        auto left_edge_id  = micro_vertex_id - 1;
        auto right_edge_id = micro_vertex_id + 1;
      
	local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);
	local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);
	     
	std::vector<double> dof_left_edge(local_dof_map.num_basis_functions());
	std::vector<double> dof_right_edge(local_dof_map.num_basis_functions());  
	
	extract_dof(dof_indices_left, *u, dof_left_edge);
	extract_dof(dof_indices_right, *u, dof_right_edge);      
	      
	for(int j=dof_indices.size()-1;j>0;j--){	      
	    double gamma = 1.0/(2.0*( 2.0*(double) j -1.0 ));
	    double diff_left  = gamma * ( dof_edge[ j-1 ] - dof_left_edge[ j-1 ] ); 
	    double diff_right = gamma * ( dof_right_edge[ j-1 ] - dof_edge[ j-1 ] ); 
	    double dof_central = dof_edge[ j ];
	    double new_dof = 0.0;

	    if( dof_central > 0.0 && diff_left > 0.0  && diff_right > 0.0 ){          
		new_dof = std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) ); 
  	    }
		  
	    if( dof_central < 0.0 && diff_left < 0.0  && diff_right < 0.0 ){          
		new_dof = -std::min( std::min( std::abs( dof_central), std::abs( diff_left ) ), std::abs( diff_right ) );        
	    }
	    
	    //if( std::abs( dof_central-new_dof )>1.0e-10 ) 
	    	 u->set(dof_indices[j],new_dof);		         
        } 	      
      }  
    }      
  }
}

void additively_assemble_matrix_inner_boundaries(MPI_Comm comm,
                                                 double tau,
                                                 double t,
                                                 const UpwindProvider &upwind_provider,
                                                 const DofMap &dof_map,
                                                 const GraphStorage &graph,
                                                 PetscMat &A) {
  std::vector<double> v_up;

  for (const auto &e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    const auto macro_edge = graph.get_edge(e_id);

    v_up.resize(macro_edge->num_micro_vertices());

    const auto &local_dof_map = dof_map.get_local_dof_map(*macro_edge);

    std::vector<std::size_t> dof_indices_left(local_dof_map.num_basis_functions());
    std::vector<std::size_t> dof_indices_right(local_dof_map.num_basis_functions());

    using BPT = BoundaryPointType;
    auto pattern_ll = create_boundary(local_dof_map, BPT::Left, BPT::Left);
    auto pattern_lr = create_boundary(local_dof_map, BPT::Left, BPT::Right);
    auto pattern_rl = create_boundary(local_dof_map, BPT::Right, BPT::Left);
    auto pattern_rr = create_boundary(local_dof_map, BPT::Right, BPT::Right);

    upwind_provider.get_upwinded_values(t, *macro_edge, v_up);

    for (size_t micro_vertex_id = 1; micro_vertex_id < macro_edge->num_micro_vertices() - 1; micro_vertex_id += 1) {
      auto left_edge_id = micro_vertex_id - 1;
      auto right_edge_id = micro_vertex_id;

      local_dof_map.dof_indices(left_edge_id, 0, dof_indices_left);
      local_dof_map.dof_indices(right_edge_id, 0, dof_indices_right);

      const double v = v_up[micro_vertex_id];

      if (v > 0) {
        Eigen::MatrixXd mat_ll = +tau * v * pattern_rr;
        Eigen::MatrixXd mat_rl = -tau * v * pattern_lr;

        A.add(dof_indices_left, dof_indices_left, mat_ll);
        A.add(dof_indices_right, dof_indices_left, mat_rl);
      } else {
        Eigen::MatrixXd mat_lr = +tau * v * pattern_rl;
        Eigen::MatrixXd mat_rr = -tau * v * pattern_ll;

        A.add(dof_indices_left, dof_indices_right, mat_lr);
        A.add(dof_indices_right, dof_indices_right, mat_rr);
      }
    }
  }
}

void additively_assemble_rhs_inflow(MPI_Comm comm,
                                    double tau,
                                    double t,
                                    const UpwindProvider &upwind_provider,
                                    const DofMap &dof_map,
                                    const GraphStorage &graph,
                                    const std::function<double(double)> &inflow_function,
                                    PetscVec &rhs) {
  std::vector<double> Q_up(1);
  std::vector<double> A_up(1);

  for (const auto &v_id : graph.get_active_vertex_ids(mpi::rank(comm))) {
    auto &vertex = *graph.get_vertex(v_id);
    if (!vertex.is_inflow())
      continue;

    auto &neighbor_edge = *graph.get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = dof_map.get_local_dof_map(neighbor_edge);
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_id) ? neighbor_edge.num_micro_edges() - 1 : 0;
    const double sigma = neighbor_edge.is_pointing_to(v_id) ? +1 : -1;

    upwind_provider.get_upwinded_values(t, vertex, A_up, Q_up);

    const double v_up = Q_up[0] / A_up[0];

    // TODO: generalize this to allow different functions
    const auto c_in = A_up[0] * inflow_function(t);

    std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, 0, dof_indices);
    std::vector<double> rhs_values(local_dof_map.num_basis_functions());
    for (size_t j = 0; j < local_dof_map.num_basis_functions(); j += 1) {
      // L^{-1} * tau * q_in(t) * phi_j(-1) = L^{-1} * tau * q_in(t) * (-1)^{j}
      rhs_values[j] = tau * (-sigma * v_up * c_in) * std::pow(sigma, j);
    }

    rhs.add(dof_indices, rhs_values);
  }
}


void additively_assemble_matrix_outflow(MPI_Comm comm,
                                        double tau,
                                        double t,
                                        const UpwindProvider &upwind_provider,
                                        const DofMap &dof_map,
                                        const GraphStorage &graph,
                                        PetscMat &A) {
  std::vector<double> Q_up(1);
  std::vector<double> A_up(1);

  for (const auto &v_id : graph.get_active_vertex_ids(mpi::rank(comm))) {
    auto &vertex = *graph.get_vertex(v_id);

    if (!vertex.is_leaf())
      continue;
    if (vertex.is_inflow())
      continue;

    if (vertex.is_linear_characteristic_inflow())
      continue;
    if (vertex.is_nonlinear_characteristic_inflow())
      continue;

    auto &neighbor_edge = *graph.get_edge(vertex.get_edge_neighbors()[0]);
    auto &local_dof_map = dof_map.get_local_dof_map(neighbor_edge);
    auto micro_edge_idx = neighbor_edge.is_pointing_to(v_id) ? neighbor_edge.num_micro_edges() - 1 : 0;

    std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
    local_dof_map.dof_indices(micro_edge_idx, 0, dof_indices);

    const double sigma = neighbor_edge.is_pointing_to(v_id) ? +1 : -1;

    auto pattern = neighbor_edge.is_pointing_to(v_id)
                     ? create_boundary(local_dof_map, BoundaryPointType::Right, BoundaryPointType::Right)
                     : create_boundary(local_dof_map, BoundaryPointType::Left, BoundaryPointType::Left);

    upwind_provider.get_upwinded_values(t, vertex, A_up, Q_up);

    const double v_up = Q_up[0] / A_up[0];

    Eigen::MatrixXd mat = sigma * tau * v_up * pattern;

    A.add(dof_indices, dof_indices, mat);
  }
}

void additively_assemble_matrix_nfurcations(MPI_Comm comm,
                                            double tau,
                                            double t,
                                            const UpwindProvider &upwind_provider,
                                            const DofMap &dof_map,
                                            const GraphStorage &graph,
                                            PetscMat &A) {
  for (auto v_idx : graph.get_active_vertex_ids(mpi::rank(comm))) {
    auto &vertex = *graph.get_vertex(v_idx);

    if (!vertex.is_bifurcation())
      continue;

    const auto num_edges = vertex.get_edge_neighbors().size();

    // collect edges:
    std::vector<Edge const *> edges;
    for (auto e_id : vertex.get_edge_neighbors())
      edges.push_back(graph.get_edge(e_id).get());

    // collect normals
    std::vector<double> sigma;
    for (auto edge : edges)
      sigma.push_back(edge->is_pointing_to(v_idx) ? +1. : -1.);

    // get upwinded values
    std::vector<double> A_up(num_edges);
    std::vector<double> Q_up(num_edges);
    upwind_provider.get_upwinded_values(t, vertex, A_up, Q_up);

    // get boundary point type
    using BPT = BoundaryPointType;
    std::vector<BPT> boundary_type;
    for (auto edge : edges)
      boundary_type.push_back(edge->is_pointing_to(v_idx) ? BPT::Right : BPT::Left);

    std::vector<std::vector<size_t>> dof_indices_list;
    for (int i = 0; i < num_edges; i += 1) {
      auto edge = edges[i];
      auto &local_dof_map = dof_map.get_local_dof_map(*edge);
      auto micro_edge_id = edge->is_pointing_to(v_idx) ? edge->num_micro_edges() - 1 : 0;
      std::vector<size_t> dof_indices(local_dof_map.num_basis_functions());
      local_dof_map.dof_indices(micro_edge_id, 0, dof_indices);
      dof_indices_list.push_back(dof_indices);
    }

    assembly_kernel_nfurcation(comm, tau, dof_indices_list, edges, sigma, Q_up, A_up, boundary_type, A);
  }
}

void assembly_kernel_nfurcation(MPI_Comm comm,
                                double tau,
                                const std::vector<std::vector<size_t>> &dof_indices_list,
                                const std::vector<Edge const *> &edges,
                                const std::vector<double> &sigma,
                                const std::vector<double> &Q_up,
                                const std::vector<double> &A_up,
                                const std::vector<BoundaryPointType> &boundary_type,
                                PetscMat &A) {
  const size_t num_edges = dof_indices_list.size();

  if (num_edges != dof_indices_list.size() || num_edges != sigma.size() || num_edges != Q_up.size() || num_edges != A_up.size() || num_edges != boundary_type.size())
    throw std::runtime_error("assembly_kernel_nfurcation: inconsistent sizes");

  // divide into inflow and outflow edges
  std::vector<size_t> inflows;
  std::vector<size_t> outflows;
  for (size_t k = 0; k < num_edges; k += 1) {
    // if the flux is very small (mainly at the beginning), we have a bias to make it an inflow
    if (Q_up[k] / A_up[k] * sigma[k] >= 0 || std::abs(Q_up[k]) < 1e-14)
      inflows.push_back(k);
    else
      outflows.push_back(k);
  }

  // get total inflow flux
  double Q_in = 0;
  for (auto k : inflows)
    Q_in += Q_up[k] * sigma[k];

  // assemble the inflow edges
  for (auto i : inflows) {
    auto &edge = *edges[i];

    // we only assemble rows which belong to us:
    if (edge.rank() != mpi::rank(comm))
      continue;

    const auto &dof_indices_i = dof_indices_list[i];
    auto pattern = create_boundary(dof_indices_i.size(), boundary_type[i], boundary_type[i]);
    Eigen::MatrixXd mat = tau * Q_up[i] / A_up[i] * sigma[i] * pattern;
    A.add(dof_indices_i, dof_indices_i, mat);
  }

  // assemble the outflow edges
  for (auto i : outflows) {
    auto &edge_i = *edges[i];

    // we only assemble rows which belong to us:
    if (edge_i.rank() != mpi::rank(comm))
      continue;

    const auto &dof_indices_i = dof_indices_list[i];

    for (auto j : inflows) {
      const auto &dof_indices_j = dof_indices_list[j];
      auto pattern = create_boundary(dof_indices_j.size(), boundary_type[i], boundary_type[j]);
      const double q_factor = std::abs(Q_in) < 1e-16 ? 0.0 : ((Q_up[i] * Q_up[j]) / Q_in);
      Eigen::MatrixXd mat = tau * sigma[i] * q_factor * sigma[j] / A_up[j] * pattern;
      A.add(dof_indices_i, dof_indices_j, mat);
    }
  }
}

} // namespace implicit_transport

} // namespace macrocirculation
