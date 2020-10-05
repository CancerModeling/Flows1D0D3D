////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"
#include <utilIO.hpp>

// local anonymous namespace
namespace {

double diri_bc_val = 0.;
double force_c = 1.;
double coeff_a = 1.;

void random_init() { srand(time(nullptr)); }

// random number between 0 and 1
double get_random_number() { return double(std::rand()) / double(RAND_MAX); }

// random ppint with each element between 0 and 1
Point get_random_point() {
  return {double(std::rand()) / double(RAND_MAX),
          double(std::rand()) / double(RAND_MAX),
          double(std::rand()) / double(RAND_MAX)};
}

/*!
 * @brief Segment of blood network. It consists of two ids of end nodes, and
 * ids of other segments it is associated to
 */
struct Segment {
  unsigned int d_n1;
  unsigned int d_n2;

  std::vector<unsigned int> d_conn;

  Segment() : d_n1(0), d_n2(0){};
  Segment(const unsigned int &n1, const unsigned int &n2)
      : d_n1(n1), d_n2(n2){};
};

/*!
 * @brief Square matrix of general size
 */
struct Matrix {

  std::vector<std::vector<double>> d_data;
  unsigned int d_size;

  Matrix(const unsigned int &size, bool identity = false) : d_size(size) {
    for (unsigned int i = 0; i < d_size; i++)
      d_data.push_back(std::vector<double>(d_size, 0.));

    if (identity)
      for (unsigned int i = 0; i < d_size; i++)
        d_data[i][i] = 1.;
  }

  void add(Matrix &m, double factor = 1.) {
    if (m.d_size != d_size) {
      libmesh_error_msg("Matrix must be of same size\n");
      exit(1);
    }
    for (unsigned int i = 0; i < d_size; i++)
      for (unsigned int j = 0; j < d_size; j++)
        (*this)(i, j) += factor * m(i, j);
  }

  Matrix mult(double a) {

    Matrix n(d_size);
    for (unsigned int i = 0; i < d_size; i++)
      for (unsigned int j = 0; j < d_size; j++)
        n(i, j) = a * (*this)(i, j);

    return n;
  }

  Matrix mult(Matrix &m) {
    if (m.d_size != d_size) {
      libmesh_error_msg("Matrix must be of same size\n");
      exit(1);
    }

    Matrix n(d_size);
    for (unsigned int i = 0; i < d_size; i++)
      for (unsigned int j = 0; j < d_size; j++)
        for (unsigned int k = 0; k < d_size; k++)
          n(i, j) += (*this)(i, k) * m(k, j);

    return n;
  }

  std::vector<double> mult(std::vector<double> &b) {
    if (b.size() != d_size) {
      libmesh_error_msg("Matrix and vector must be of same size\n");
      exit(1);
    }

    std::vector<double> r(b.size(), 0.);
    for (unsigned int i = 0; i < d_size; i++)
      for (unsigned int j = 0; j < d_size; j++)
        r[i] += (*this)(i, j) * b[j];

    return r;
  }

  double det() const {
    if (d_size == 1)
      return (*this)(0, 0);
    else if (d_size == 2)
      return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
    else if (d_size == 3)
      return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) -
                              (*this)(2, 1) * (*this)(1, 2)) -
             (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) -
                              (*this)(2, 0) * (*this)(1, 2)) +
             (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) -
                              (*this)(2, 0) * (*this)(1, 1));

    libmesh_error_msg("Determinant of matrix above size 3 is not implemented\n");
    exit(1);
  }

  Matrix inv() {
    if (d_size == 1) {
      Matrix n(d_size);
      n(0, 0) = 1. / (*this)(0, 0);

      return n;
    } else if (d_size == 2) {

      auto det_inv = 1. / this->det();

      Matrix n(d_size);
      n(0, 0) = det_inv * (*this)(1, 1);
      n(1, 1) = det_inv * (*this)(0, 0);

      n(0, 1) = -det_inv * (*this)(0, 1);
      n(1, 0) = -det_inv * (*this)(1, 0);

      return n;
    } else if (d_size == 3) {

      auto det_inv = 1. / this->det();

      Matrix n(d_size);

      n(0, 0) = det_inv *
                ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2));
      n(0, 1) = -det_inv *
                ((*this)(0, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(0, 2));
      n(0, 2) = det_inv *
                ((*this)(0, 1) * (*this)(1, 2) - (*this)(1, 1) * (*this)(0, 2));

      n(1, 0) = -det_inv *
                ((*this)(1, 0) * (*this)(2, 2) - (*this)(2, 0) * (*this)(1, 2));
      n(1, 1) = det_inv *
                ((*this)(0, 0) * (*this)(2, 2) - (*this)(2, 0) * (*this)(0, 2));
      n(1, 2) = -det_inv *
                ((*this)(0, 0) * (*this)(1, 2) - (*this)(1, 0) * (*this)(0, 2));

      n(2, 0) = det_inv *
                ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1));
      n(2, 1) = -det_inv *
                ((*this)(0, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(0, 1));
      n(2, 2) = det_inv *
                ((*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1));

      return n;
    }

    libmesh_error_msg("Inverse of matrix above size 3 is not implemented\n");
    exit(1);
  }

  Matrix pow(int a) {

    if (a < 0) {
      libmesh_error_msg("Can not compute negative power of matrix\n");
      exit(1);
    }

    if (a == 0)
      return Matrix(d_size, true);
    else if (a == 1) {
      Matrix n(d_size);
      this->copy(n);

      return n;
    } else {
      Matrix n(d_size);
      this->copy(n);
      for (int p = 2; p <= a; p++)
        n = this->mult(n);

      return n;
    }
  }

  void copy(Matrix &n) {

    n.d_size = d_size;
    n.d_data = d_data;
  }

  double &operator()(size_t i, size_t j) { return d_data[i][j]; }
  const double &operator()(size_t i, size_t j) const { return d_data[i][j]; }

  /*!
   * @brief Prints the information
   *
   * @param nt Number of tabs to append before printing
   * @param lvl Information level (higher means more information)
   */
  std::string printStr(int nt = 0, int lvl = 0) const {

    auto tabS = util::io::getTabS(nt);
    std::ostringstream oss;

    oss << tabS << "------- Matrix --------" << std::endl
        << std::endl;
    oss << tabS << "Size = " << d_size << std::endl;
    oss << tabS << "Elements:" << std::endl;
    oss << tabS << "[";
    for (unsigned int i = 0; i < d_size; i++) {
      if (i > 0)
        oss << tabS;

      oss << util::io::printStr(d_data[i], 1) << std::endl;
    }
    oss << tabS << "]" << std::endl;

    if (lvl > 0) {
      oss << tabS << "Determinant = " << this->det() << std::endl;
    }

    oss << std::endl;

    return oss.str();
  }

  void print(int nt = 0, int lvl = 0) const { std::cout << printStr(nt, lvl); }
};

unsigned int create_node(Point p, std::vector<unsigned int> &nodes,
                         ReplicatedMesh &mesh) {

  auto node = mesh.add_point(p);
  nodes.push_back(node->id());

  return node->id();
}

unsigned int create_elem(unsigned int node_1, unsigned int node_2,
                         std::vector<unsigned int> &elems,
                         ReplicatedMesh &mesh) {

  auto elem = Elem::build(EDGE2).release();
  elem->set_node(0) = mesh.node_ptr(node_1);
  elem->set_node(1) = mesh.node_ptr(node_2);
  auto add_elem = mesh.add_elem(elem);
  elems.push_back(add_elem->id());

  return add_elem->id();
}

void assemble_net(EquationSystems &es, const std::string &system_name) {

  // flow system
  auto &flow = es.get_system<TransientLinearImplicitSystem>("Flow");
  const unsigned int v_flow = flow.variable_number("pressure");
  const DofMap &flow_map = flow.get_dof_map();
  std::vector<unsigned int> dof_indices_flow;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = flow.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Parameters

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    flow_map.dof_indices(elem, dof_indices_flow, v_flow);

    const unsigned int n_dofs = dof_indices_flow.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number flow_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++)
        flow_cur += phi[l][qp] * flow.current_solution(dof_indices_flow[l]);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += JxW[qp] * force_c * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += JxW[qp] * coeff_a * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    flow_map.heterogenously_constrain_element_matrix_and_vector(
      Ke, Fi, dof_indices_flow);
    flow.matrix->add_matrix(Ke, dof_indices_flow);
    flow.rhs->add_vector(Fi, dof_indices_flow);
  }
}

void assemble_net_2(EquationSystems &es, const std::string &system_name) {

  // flow system
  auto &flow = es.get_system<TransientLinearImplicitSystem>("Flow");
  const unsigned int v_flow = flow.variable_number("pressure");
  const DofMap &flow_map = flow.get_dof_map();
  std::vector<unsigned int> dof_indices_flow;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = flow.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Parameters

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    flow_map.dof_indices(elem, dof_indices_flow, v_flow);

    const unsigned int n_dofs = dof_indices_flow.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number flow_cur = 0.;

      for (unsigned int l = 0; l < phi.size(); l++)
        flow_cur += phi[l][qp] * flow.current_solution(dof_indices_flow[l]);

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += JxW[qp] * force_c * phi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          // now we have derivative of shape function

          Ke(i, j) += JxW[qp] * coeff_a * dphi[j][qp] * dphi[i][qp];
        }
      }
    } // loop over quadrature points

    flow_map.heterogenously_constrain_element_matrix_and_vector(
      Ke, Fi, dof_indices_flow);
    flow.matrix->add_matrix(Ke, dof_indices_flow);
    flow.rhs->add_vector(Fi, dof_indices_flow);
  }
}

void assemble_net_transient(EquationSystems &es,
                            const std::string &system_name) {

  // flow system
  auto &flow = es.get_system<TransientLinearImplicitSystem>("Flow");
  const unsigned int v_flow = flow.variable_number("pressure");
  const DofMap &flow_map = flow.get_dof_map();
  std::vector<unsigned int> dof_indices_flow;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = flow.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Parameters
  const Real dt = es.parameters.get<Real>("time_step");

  // Arranging matrix
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fi;

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    flow_map.dof_indices(elem, dof_indices_flow, v_flow);

    const unsigned int n_dofs = dof_indices_flow.size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fi.resize(n_dofs);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

      // Computing solution
      Number flow_cur = 0.;
      Number flow_old = 0.;
      Gradient flow_grad;

      for (unsigned int l = 0; l < phi.size(); l++) {
        flow_cur += phi[l][qp] * flow.current_solution(dof_indices_flow[l]);
        flow_old += phi[l][qp] * flow.old_solution(dof_indices_flow[l]);

        flow_grad.add_scaled(dphi[l][qp],
                             flow.old_solution(dof_indices_flow[l]));
      }

      // Assembling matrix
      for (unsigned int i = 0; i < phi.size(); i++) {

        Fi(i) += JxW[qp] * (flow_old + dt * force_c) * phi[i][qp];

        Fi(i) -= JxW[qp] * dt * coeff_a * flow_grad * dphi[i][qp];

        for (unsigned int j = 0; j < phi.size(); j++) {

          Ke(i, j) += JxW[qp] * phi[j][qp] * phi[i][qp];
        }
      }
    } // loop over quadrature points

    flow_map.heterogenously_constrain_element_matrix_and_vector(
      Ke, Fi, dof_indices_flow);
    flow.matrix->add_matrix(Ke, dof_indices_flow);
    flow.rhs->add_vector(Fi, dof_indices_flow);
  }
}

Number ic_net(const Point &p, const Parameters &es,
              const std::string &system_name, const std::string &var_name) {

  if (var_name == "pressure") {

    return 0.;
  }

  return 0.;
}

void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  sys.project_solution(ic_net, nullptr, es.parameters);
}

void boundary_condition(EquationSystems &eq_sys) {

  auto &sys = eq_sys.get_system<TransientLinearImplicitSystem>("Flow");

  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(0);

  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("pressure"));

  // constant function
  ConstFunction<Number> twof(diri_bc_val);

  DirichletBoundary diri_bc(boundary_ids, vars, &twof);

  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

void solve_eq_sys(EquationSystems &eq_sys, std::vector<unsigned int> &nodes,
                  std::vector<unsigned int> &elems, std::string mesh_file,
                  std::string sim_file) {

  // get flow system
  auto &flow = eq_sys.get_system<TransientLinearImplicitSystem>("Flow");

  // get mesh
  MeshBase &mesh = eq_sys.get_mesh();

  // get parameters
  const Real dt = eq_sys.parameters.get<Real>("time_step");
  const Real init_time = eq_sys.parameters.get<Real>("init_time");
  const unsigned int num_steps = eq_sys.parameters.get<unsigned int>("num_steps");

  eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = 10;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln(flow.solution->clone());

  size_t t_step = 0;
  double current_time = init_time;
  do {

    out << "Solving for time step = " << t_step << std::endl;

    // if dt is modified (say due to adaptive time discretization) then
    // update it in the parameter list
    eq_sys.parameters.set<Real>("time_step") = dt;

    // Prepare time step
    t_step++;
    current_time += dt;
    flow.time = current_time;
    eq_sys.parameters.set<Real>("time") = flow.time;

    *flow.old_local_solution = *flow.current_local_solution;

    // solve
    flow.solve();

    // write
    {
      // write solution to file
      char sim_file_i[400];
      sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
      ExodusII_IO exo(mesh);
      exo.write_timestep(sim_file_i, eq_sys, t_step, flow.time);
    }

  } while (t_step < num_steps);
}

void verify_solution(EquationSystems &eq_sys, std::vector<unsigned int> &nodes,
                     std::vector<unsigned int> &elems, std::string test_case,
                     double p_in, double seg_length) {

  out << "Verifying solution " << std::endl;

  // get flow system
  auto &flow = eq_sys.get_system<TransientLinearImplicitSystem>("Flow");

  // get mesh
  MeshBase &mesh = eq_sys.get_mesh();

  // verify solution
  if (test_case == "one_segment") {

    double p2_exact = seg_length * seg_length * 0.5 + p_in;

    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));

    out << "Exact sol node 2 = (" << p2_exact << "), "
        << "Approx sol node 2 = (" << p2_approx << ")" << std::endl;

    if (std::abs(p2_approx - p2_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;
      libmesh_error();
    }
  } else if (test_case == "two_segments") {

    double p2_exact = 3. * seg_length * seg_length * 0.5 + p_in;
    double p3_exact = 2. * seg_length * seg_length + p_in;

    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));
    double p3_approx =
      flow.current_solution(mesh.node_ptr(nodes[2])->dof_number(0, 0, 0));

    out << "Exact sols node 2 and 3 = (" << p2_exact << ", " << p3_exact
        << "), "
        << "Approx sols node 2 and 3 = (" << p2_approx << ", " << p3_approx
        << ")" << std::endl;

    if (std::abs(p2_approx - p2_exact) > 1.0E-10 ||
        std::abs(p3_approx - p3_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;
      libmesh_error();
    }
  } else if (test_case == "three_segments" ||
             test_case == "three_segments_branch") {

    double p2_exact = 5. * seg_length * seg_length * 0.5 + p_in;
    double p3_exact = 4. * seg_length * seg_length + p_in;
    double p4_exact = 9. * seg_length * seg_length * 0.5 + p_in;
    if (test_case == "three_segments_branch") {
      p3_exact = 3. * seg_length * seg_length + p_in;
      p4_exact = 3. * seg_length * seg_length + p_in;
    }

    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));
    double p3_approx =
      flow.current_solution(mesh.node_ptr(nodes[2])->dof_number(0, 0, 0));
    double p4_approx =
      flow.current_solution(mesh.node_ptr(nodes[3])->dof_number(0, 0, 0));

    out << "Exact sols node 2, 3, 4 = (" << p2_exact << ", " << p3_exact << ", "
        << p4_exact << "), "
        << "Approx sols node 2, 3, 4 = (" << p2_approx << ", " << p3_approx
        << ", " << p4_approx << ")" << std::endl;

    if (std::abs(p2_approx - p2_exact) > 1.0E-10 ||
        std::abs(p3_approx - p3_exact) > 1.0E-10 ||
        std::abs(p4_approx - p4_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;
      libmesh_error();
    }
  }
}

void verify_transient_solution(EquationSystems &eq_sys,
                               std::vector<unsigned int> &nodes,
                               std::vector<unsigned int> &elems,
                               std::string test_case,
                               double p_in,
                               double seg_length) {

  // get flow system
  auto &flow = eq_sys.get_system<TransientLinearImplicitSystem>("Flow");

  // get mesh
  MeshBase &mesh = eq_sys.get_mesh();

  // get parameters
  const Real dt = eq_sys.parameters.get<Real>("time_step");
  const Real init_time = eq_sys.parameters.get<Real>("init_time");
  const unsigned int num_steps =
    eq_sys.parameters.get<unsigned int>("num_steps");

  // verify solution
  if (test_case == "one_segment") {

    double p2_exact = 0.;
    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));

    // compute p2_exact at time step = num_steps
    double alpha = 3. * dt / (seg_length * seg_length);
    Matrix A(1);
    A(0, 0) = 1. - alpha;
    std::vector<double> b(1, 0.);
    b[0] = 3. * dt / 2. + alpha * p_in;

    // get matrix
    Matrix M(1);
    for (int i = 0; i < num_steps; i++) {
      auto Mi = A.pow(i);
      M.add(Mi);
    }

    auto Mb = M.mult(b);
    p2_exact = Mb[0];

    out << "Exact sol node 2 = (" << p2_exact << "), "
        << "Approx sol node 2 = (" << p2_approx << ")" << std::endl;

    if (std::abs(p2_approx - p2_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;
      libmesh_error();
    }
  } else if (test_case == "two_segments") {

    double p2_exact = 0.;
    double p3_exact = 0.;
    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));
    double p3_approx =
      flow.current_solution(mesh.node_ptr(nodes[2])->dof_number(0, 0, 0));

    // compute p2_exact and p3_exact at time step = num_steps
    double alpha = 3. * dt / (seg_length * seg_length);

    // mass matrix and its inverse
    Matrix M(2);
    M(0, 0) = 2.;
    M(0, 1) = 0.5;
    M(1, 0) = 0.5;
    M(1, 1) = 1.;

    Matrix Mbar = M.inv();

    // stiffness matrix
    Matrix K(2);
    K(0, 0) = 2. * alpha;
    K(0, 1) = -alpha;
    K(1, 0) = -alpha;
    K(1, 1) = alpha;

    // A matrix
    Matrix Id(2, true);
    Matrix A_temp = Mbar.mult(K);
    auto A = A_temp.mult(-1.);
    A.add(Id);

    // get force vector
    std::vector<double> b_temp(2, 0.);
    b_temp[0] = 3. * dt + alpha * p_in;
    b_temp[1] = 3. * dt / 2.;

    auto b = Mbar.mult(b_temp);

    // get factor matrix
    Matrix D(2);
    for (int i = 0; i < num_steps; i++) {
      auto Di = A.pow(i);
      D.add(Di);
    }

    auto Db = D.mult(b);
    p2_exact = Db[0];
    p3_exact = Db[1];

    out << "Exact sols node 2 and 3 = (" << p2_exact << ", " << p3_exact
        << "), "
        << "Approx sols node 2 and 3 = (" << p2_approx << ", " << p3_approx
        << ")" << std::endl;

    if (std::abs(p2_approx - p2_exact) > 1.0E-10 ||
        std::abs(p3_approx - p3_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;

      libmesh_error();
    }
  } else if (test_case == "three_segments_branch") {

    double p2_exact = 0.;
    double p3_exact = 0.;
    double p4_exact = 0.;
    double p2_approx =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));
    double p3_approx =
      flow.current_solution(mesh.node_ptr(nodes[2])->dof_number(0, 0, 0));
    double p4_approx =
      flow.current_solution(mesh.node_ptr(nodes[3])->dof_number(0, 0, 0));

    // compute p2_exact and p3_exact at time step = num_steps
    double alpha = 3. * dt / (seg_length * seg_length);

    // mass matrix and its inverse
    Matrix M(3);
    M(0, 0) = 3.;
    M(0, 1) = 0.5;
    M(0, 2) = 0.5;
    M(1, 0) = 0.5;
    M(1, 1) = 1.;
    M(1, 2) = 0.;
    M(2, 0) = 0.5;
    M(2, 1) = 0.;
    M(2, 2) = 1.;

    Matrix Mbar = M.inv();

    // stiffness matrix
    Matrix K(3);
    K(0, 0) = 3. * alpha;
    K(0, 1) = -alpha;
    K(0, 2) = -alpha;
    K(1, 0) = -alpha;
    K(1, 1) = alpha;
    K(1, 2) = 0.;
    K(2, 0) = -alpha;
    K(2, 1) = 0.;
    K(2, 2) = alpha;

    // A matrix
    Matrix Id(3, true);
    Matrix A_temp = Mbar.mult(K);
    auto A = A_temp.mult(-1.);
    A.add(Id);

    // get force vector
    std::vector<double> b_temp(3, 0.);
    b_temp[0] = 9. * dt / 2. + alpha * p_in;
    b_temp[1] = 3. * dt / 2.;
    b_temp[2] = 3. * dt / 2.;

    auto b = Mbar.mult(b_temp);

    // get factor matrix
    Matrix D(3);
    for (int i = 0; i < num_steps; i++) {
      auto Di = A.pow(i);
      D.add(Di);
    }

    auto Db = D.mult(b);
    p2_exact = Db[0];
    p3_exact = Db[1];
    p4_exact = Db[2];

    out << "Exact sols node 2, 3, 4 = (" << p2_exact << ", " << p3_exact << ", "
        << p4_exact << "), "
        << "Approx sols node 2, 3, 4 = (" << p2_approx << ", " << p3_approx
        << ", " << p4_approx << ")" << std::endl;

    // for debug
    //    out << "M" << std::endl;
    //    out << M.printStr(1, 1);
    //    out << "Mbar" << std::endl;
    //    out << Mbar.printStr(1);
    //    out << "K" << std::endl;
    //    out << K.printStr(1);
    //    out << "A" << std::endl;
    //    out << A.printStr(1);
    //    out << "D" << std::endl;
    //    out << D.printStr(1);

    if (std::abs(p2_approx - p2_exact) > 1.0E-10 ||
        std::abs(p3_approx - p3_exact) > 1.0E-10 ||
        std::abs(p4_approx - p4_exact) > 1.0E-10) {
      out << "Approx solution does not match exact solution for test = "
          << test_case << std::endl;


      libmesh_error();
    }
  }
}

} // namespace

//
// Test: Same as add_node_elem_eq_sys_test() but now we add
// aseembly function for simple 1-d flow equation on network and solve the
// 1-d equation on changing mesh
//
void test::mesh::eq_sys_assemble(int argc, char **argv,
                                 Parallel::Communicator *comm) {

  out << "********** TestMesh 3 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t3/mesh";
  std::string sim_file = "t3/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<unsigned int> nodes;
  std::vector<unsigned int> elems;
  std::vector<Segment> segments;

  // create two nodes at point (0,0,0) and (1,1,1) and element between these
  // nodes
  {
    auto add_node = mesh.add_point(Point(0., 0., 0.));
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(Point(1., 1., 1.));
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[0]);
    elem->set_node(1) = mesh.node_ptr(nodes[1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[0], nodes[1]);

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // assign boundary ids
  {
    auto &bd_info = mesh.get_boundary_info();
    // auto elem_0 = mesh.elem(elems[0]);
    // flow.get_dof_map().constrain_p_dofs(0, elem_0, 0, elem_0->level());

    // add 1st node to boundary 0 and 2nd to 1
    bd_info.add_node(nodes[0], 0);
    bd_info.add_node(nodes[1], 1);

    // add element's left side to boundary
    bd_info.add_side(elems[0], 0, 0);

    // add element's right side to boundary
    bd_info.add_side(elems[0], 1, 1);
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;

  // set some parameters
  force_c = 1.0;
  coeff_a = 1.0;

  // attach initial condition
  flow.attach_init_function(initial_condition);

  // attach assembly function
  flow.attach_assemble_function(assemble_net);

  // Add boundary condition
  size_t bc_case = 4;
  diri_bc_val = bc_case;
  boundary_condition(eq_sys);

  // init
  eq_sys.init();

  // print extra information
  if (false) {
    out << "old sol at node 1 = "
        << flow.old_solution(mesh.node_ptr(nodes[0])->dof_number(0, 0, 0))
        << ", current sol at node 1 = "
        << flow.current_solution(mesh.node_ptr(nodes[0])->dof_number(0, 0, 0))
        << std::endl;

    // auto sol = flow.current_local_solution.release();
    // sol->print(out);

    // print dof info
    flow.get_dof_map().print_info(out);
    if (flow.get_dof_map().is_constrained_dof(0))
      out << "0 dof is constrained" << std::endl;
    else
      out << "0 dof is not constrained" << std::endl;

    // print bd info
    auto &bd_info = mesh.get_boundary_info();
    bd_info.print_info(out);
  }

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // solve above system with just 2 nodes and 1 edge elements and compare the
  // results with exact solution
  {
    eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") =
      10;

    // to compute the nonlinear convergence
    UniquePtr<NumericVector<Number>> last_nonlinear_soln(
      flow.solution->clone());

    size_t t_step = 1;
    double current_time = 0.;
    do {

      out << "Solving for time step = " << t_step << std::endl;

      // if dt is modified (say due to adaptive time discretization) then
      // update it in the parameter list
      eq_sys.parameters.set<Real>("time_step") = 1.;

      // Prepare time step
      t_step++;
      current_time += t_step;
      flow.time = current_time;
      eq_sys.parameters.set<Real>("time") = flow.time;

      *flow.old_local_solution = *flow.current_local_solution;

      // solve
      flow.solve();

      // write
      {
        // write solution to file
        char sim_file_i[400];
        sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
        ExodusII_IO exo(mesh);
        exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
      }

    } while (t_step < 2);

    // verify that solution at node 2 is equal to exact numerical solution
    double exact_sol = 1.5 - diri_bc_val / 2.;

    double approx_sol =
      flow.current_solution(mesh.node_ptr(nodes[1])->dof_number(0, 0, 0));

    out << "Exact sol at node 2 = " << exact_sol
        << ", Approx sol at node 2 = " << approx_sol << std::endl;

    if (std::abs(approx_sol - exact_sol) > 1.0E-10) {
      out << "Approx solution does not match exact solution" << std::endl;
      libmesh_error();
    }
  }

  // end
}

//
// Test: Same as add_node_elem_eq_sys_test() but now we add
// aseembly function for simple 1-d flow equation on network and solve the
// 1-d equation on changing mesh
//
void test::mesh::eq_sys_assemble_2(int argc, char **argv, double p_in,
                                   Parallel::Communicator *comm) {

  out << "********** TestMesh 4 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t4/mesh";
  std::string sim_file = "t4/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<unsigned int> nodes;
  std::vector<unsigned int> elems;
  std::vector<Segment> segments;
  double seg_length = std::sqrt(3.);

  // create two nodes at point (0,0,0) and (1,1,1) and element between these
  // nodes
  {
    auto add_node = mesh.add_point(Point(0., 0., 0.));
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(Point(1., 1., 1.));
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[0]);
    elem->set_node(1) = mesh.node_ptr(nodes[1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[0], nodes[1]);

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // assign boundary ids
  {
    auto &bd_info = mesh.get_boundary_info();
    // auto elem_0 = mesh.elem(elems[0]);
    // flow.get_dof_map().constrain_p_dofs(0, elem_0, 0, elem_0->level());

    // add 1st node to boundary 0 and 2nd to 1
    bd_info.add_node(nodes[0], 0);
    bd_info.add_node(nodes[1], 1);

    // add element's left side to boundary
    bd_info.add_side(elems[0], 0, 0);

    // add element's right side to boundary
    bd_info.add_side(elems[0], 1, 1);
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;

  // set some parameters
  force_c = 1.0;
  coeff_a = 1.0;

  // attach initial condition
  flow.attach_init_function(initial_condition);

  // attach assembly function
  flow.attach_assemble_function(assemble_net_2);

  // Add boundary condition
  diri_bc_val = p_in;
  boundary_condition(eq_sys);

  // add parameters
  eq_sys.parameters.set<Real>("init_time") = 0.;
  eq_sys.parameters.set<Real>("time_step") = 1. / 1;
  eq_sys.parameters.set<unsigned int>("num_steps") = 1;
  flow.time = 0.;

  // init
  eq_sys.init();

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // Step 1:
  // Solve the final system
  {
    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "one_segment", p_in, seg_length);
  }

  // Step 2:
  // add third vertex to the system and also add second edge element between
  // second and third vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(1.5, 1. - std::sqrt(5. / 2.), 0.5), nodes, mesh);

    auto elem_i = create_elem(nodes[nodes.size() - 2], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "two_segments", p_in, seg_length);
  }

  // Step 3:
  // add fourth vertex to the system and also add third edge element between
  // third and fourth vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(2., 1. - 2. * std::sqrt(5. / 2.), 1.), nodes, mesh);

    auto elem_i = create_elem(nodes[nodes.size() - 2], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "three_segments", p_in, seg_length);
  }

  // end
}

//
// Test: Same as eq_sys_assemble_2() but now the third step has branching
// network
//
void test::mesh::eq_sys_assemble_3(int argc, char **argv, double p_in,
                                   Parallel::Communicator *comm) {

  out << "********** TestMesh 5 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t5/mesh";
  std::string sim_file = "t5/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<unsigned int> nodes;
  std::vector<unsigned int> elems;
  std::vector<Segment> segments;
  double seg_length = std::sqrt(3.);

  // create two nodes at point (0,0,0) and (1,1,1) and element between these
  // nodes
  {
    auto add_node = mesh.add_point(Point(0., 0., 0.));
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(Point(1., 1., 1.));
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[0]);
    elem->set_node(1) = mesh.node_ptr(nodes[1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[0], nodes[1]);

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // assign boundary ids
  {
    auto &bd_info = mesh.get_boundary_info();
    // auto elem_0 = mesh.elem(elems[0]);
    // flow.get_dof_map().constrain_p_dofs(0, elem_0, 0, elem_0->level());

    // add 1st node to boundary 0 and 2nd to 1
    bd_info.add_node(nodes[0], 0);
    bd_info.add_node(nodes[1], 1);

    // add element's left side to boundary
    bd_info.add_side(elems[0], 0, 0);

    // add element's right side to boundary
    bd_info.add_side(elems[0], 1, 1);
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;

  // set some parameters
  force_c = 1.0;
  coeff_a = 1.0;

  // attach initial condition
  flow.attach_init_function(initial_condition);

  // attach assembly function
  flow.attach_assemble_function(assemble_net_2);

  // Add boundary condition
  diri_bc_val = p_in;
  boundary_condition(eq_sys);

  // add parameters
  eq_sys.parameters.set<Real>("init_time") = 0.;
  eq_sys.parameters.set<Real>("time_step") = 1. / 1;
  eq_sys.parameters.set<unsigned int>("num_steps") = 1;
  flow.time = 0.;

  // init
  eq_sys.init();

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // Step 1:
  // Solve the final system
  {
    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "one_segment", p_in, seg_length);
  }

  // Step 2:
  // add third vertex to the system and also add second edge element between
  // second and third vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(1.5, 1. - std::sqrt(5. / 2.), 0.5), nodes, mesh);

    auto elem_i = create_elem(nodes[nodes.size() - 2], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "two_segments", p_in, seg_length);
  }

  // Step 3:
  // add fourth vertex to the system and also add third edge element between
  // third and fourth vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(0.5, 1. + std::sqrt(5. / 2.), 1.5), nodes, mesh);

    // branching at second vertex (note the id of node_1)
    auto elem_i = create_elem(nodes[nodes.size() - 3], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_solution(eq_sys, nodes, elems, "three_segments_branch", p_in,
                    seg_length);
  }

  // end
}

//
// Test: Transient system
//
void test::mesh::eq_sys_assemble_transient(int argc, char **argv, std::string test_case,
                                           unsigned int num_steps, double p_in,
                                           Parallel::Communicator *comm) {

  out << "********** TestMesh 6 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t6/mesh";
  std::string sim_file = "t6/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<unsigned int> nodes;
  std::vector<unsigned int> elems;
  std::vector<Segment> segments;
  double seg_length = std::sqrt(3.);

  // create two nodes at point (0,0,0) and (1,1,1) and element between these
  // nodes
  {
    auto add_node = mesh.add_point(Point(0., 0., 0.));
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(Point(1., 1., 1.));
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[0]);
    elem->set_node(1) = mesh.node_ptr(nodes[1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[0], nodes[1]);

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // assign boundary ids
  {
    auto &bd_info = mesh.get_boundary_info();
    // auto elem_0 = mesh.elem(elems[0]);
    // flow.get_dof_map().constrain_p_dofs(0, elem_0, 0, elem_0->level());

    // add 1st node to boundary 0 and 2nd to 1
    bd_info.add_node(nodes[0], 0);
    bd_info.add_node(nodes[1], 1);

    // add element's left side to boundary
    bd_info.add_side(elems[0], 0, 0);

    // add element's right side to boundary
    bd_info.add_side(elems[0], 1, 1);
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;

  // set some parameters
  force_c = 1.0;
  coeff_a = 1.0;

  // attach initial condition
  flow.attach_init_function(initial_condition);

  // attach assembly function
  flow.attach_assemble_function(assemble_net_transient);

  // Add boundary condition
  diri_bc_val = p_in;
  boundary_condition(eq_sys);

  // add parameters
  eq_sys.parameters.set<Real>("init_time") = 0.;
  eq_sys.parameters.set<Real>("time_step") = 0.01 / num_steps;
  eq_sys.parameters.set<unsigned int>("num_steps") = num_steps;
  flow.time = 0.;

  // init
  eq_sys.init();

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // Step 1:
  // Solve the final system
  if (test_case == "one_segment") {
    // solve the system
    solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

    // verify solution
    verify_transient_solution(eq_sys, nodes, elems, "one_segment", p_in,
                              seg_length);

    return;
  }


  // Step 2:
  // add third vertex to the system and also add second edge element between
  // second and third vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(1.5, 1. - std::sqrt(5. / 2.), 0.5), nodes, mesh);

    auto elem_i = create_elem(nodes[nodes.size() - 2], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    if (test_case == "two_segments") {
      // solve the system
      solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

      // verify solution
      verify_transient_solution(eq_sys, nodes, elems, test_case, p_in,
                                seg_length);

      return;
    }
  }

  // Step 3:
  // add fourth vertex to the system and also add third edge element between
  // third and fourth vertices.
  //
  // Solve the final system
  {
    auto node_i =
      create_node(Point(0.5, 1. + std::sqrt(5. / 2.), 1.5), nodes, mesh);

    // branching at second vertex (note the id of node_1)
    auto elem_i = create_elem(nodes[nodes.size() - 3], nodes[nodes.size() - 1],
                              elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // reinit the system
    eq_sys.reinit();

    if (test_case == "three_segments_branch") {

      // solve the system
      solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

      // verify solution
      verify_transient_solution(eq_sys, nodes, elems, test_case, p_in,
                                seg_length);

      return;
    }
  }

  // end
}