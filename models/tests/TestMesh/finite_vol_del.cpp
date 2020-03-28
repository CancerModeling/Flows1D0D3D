////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"
#include <utilIO.hpp>
#include <utils.hpp>

// local anonymous namespace
namespace {

double diri_bc_val = 0.;
double force_c = 1.;
double coeff_a = 1.;
double p_in = 1.;

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
  bool d_state;
  unsigned int d_n1;
  unsigned int d_n2;
  std::vector<unsigned int> d_elems;

  Segment() : d_state(false), d_n1(0), d_n2(0){};
  Segment(const unsigned int &n1, const unsigned int &n2)
      : d_state(true), d_n1(n1), d_n2(n2){};
};

struct NetNode {

  bool d_state;
  std::vector<unsigned int> d_elems;

  NetNode() : d_state(false) {};
  explicit NetNode(bool state) : d_state(state) {};
};

unsigned int create_node(const Point& p, std::vector<NetNode> &nodes,
                         ReplicatedMesh &mesh) {

  auto node = mesh.add_point(p);
  unsigned int id = node->id();
  if (nodes.size() <= id)
    nodes.resize(id+1);
  nodes[id] = NetNode(true);

  return id;
}

unsigned int create_elem(unsigned int node_1, unsigned int node_2,
                         std::vector<Segment> &elems,
                         ReplicatedMesh &mesh) {

  auto elem = Elem::build(EDGE2).release();
  elem->set_node(0) = mesh.node_ptr(node_1);
  elem->set_node(1) = mesh.node_ptr(node_2);
  auto add_elem = mesh.add_elem(elem);

  unsigned int id = add_elem->id();
  if (elems.size() <= id)
    elems.resize(id+1);

  elems[id] = Segment(node_1, node_2);

  return id;
}

void add_unique(unsigned int i, std::vector<unsigned int> &list) {

  for (const auto &l : list)
    if (l == i)
      return;

  list.push_back(i);
}

void add_unique(unsigned int dof, Real val,
                std::vector<unsigned int> &list, std::vector<Real>
                &list_val) {

  for (unsigned int i=0; i<list.size(); i++) {

    if (list[i] == dof) {
      list_val[i] += val;

      return;
    }
  }

  list.push_back(dof);
  list_val.push_back(val);
}

void update_connectivity(std::vector<NetNode> &nodes, std::vector<Segment>
    &elems, ReplicatedMesh &mesh) {

  // clear old connectivity
  for (auto &node: nodes)
    node.d_elems.clear();

  for (size_t e=0; e<elems.size(); e++) {

    auto & elem = elems[e];

    if (!elem.d_state)
      continue;

    add_unique(e, nodes[elem.d_n1].d_elems);
    add_unique(e, nodes[elem.d_n2].d_elems);
  }
}

void assemble_net(EquationSystems &es, const std::string &system_name) {

  // flow system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();

  // Parameters
  const auto *net_nodes =
      es.parameters.get<std::vector<NetNode> *>("net_nodes");
  const auto *net_elems =
      es.parameters.get<std::vector<Segment> *>("net_elems");

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Number> bc_vals;

  // Looping through elements
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {

    const Elem *elem = *el;
    pres_map.dof_indices(elem, dof_indices_pres, v_pres);

    // Arranging matrix
    std::vector<Real> Ke_val_col;
    std::vector<unsigned int> Ke_dof_col;
    std::vector<unsigned int> Ke_dof_row;
    Ke_dof_row.push_back(dof_indices_pres[0]);

    // loop over vertex of element
    const auto &net_e = (*net_elems)[elem->id()];
    for (unsigned int I=0; I<2; I++) {

      unsigned int node_I = net_e.d_n1;
      if (I == 1)
        node_I = net_e.d_n2;

      // add contribution 1 to Ke
      add_unique(dof_indices_pres[0], 1., Ke_dof_col, Ke_val_col);

      // loop over neighboring elements
      // compute sum of coefficients. by design, coefficient for each element
      // is 1 so sum will be 1 times number of neighbors
      const auto &neighbors = (*net_nodes)[node_I].d_elems;
      Real bar_t = neighbors.size() * 1.;

      out << "  elem = " << elem->id()
          << ", node = " << node_I
          << ", num elems = " << neighbors.size() << "\n";

      // loop over elements
      for (const auto neigh_e : neighbors) {

        // get dof id
        std::vector<unsigned int> dof_indices_pres_neigh;
        const Elem *neigh_elem = mesh.elem_ptr(neigh_e);
        pres_map.dof_indices(neigh_elem, dof_indices_pres_neigh, v_pres);

        // compute quantity to be added
        Real a = -1. * 1. / (bar_t);

        add_unique(dof_indices_pres_neigh[0], a, Ke_dof_col, Ke_val_col);
      }
    }

    // add to matrix
    DenseMatrix<Number> Ke;
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i=0; i<Ke_dof_col.size(); i++)
      Ke(0,i) = Ke_val_col[i];

    pres.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // check for bc
    if (elem->id() == 0) {
      bc_rows.push_back(dof_indices_pres[0]);
      bc_vals.push_back(p_in);
    }
  }

  out << "Adding BC constraint\n";
  pres.matrix->close();
  pres.rhs->close();

  out << "\n ++++ \n";
  pres.matrix->print(out);

  // apply bc constraint
  pres.matrix->zero_rows(bc_rows, 1.);
  for (unsigned int i=0; i<bc_rows.size(); i++)
    pres.rhs->set(bc_rows[i], bc_vals[i]);

  pres.matrix->close();
  pres.rhs->close();

  out << "\n ++++ \n";
  pres.matrix->print(out);
}

void boundary_condition(EquationSystems &es) {

  return;

  auto &sys = es.get_system<TransientLinearImplicitSystem>("Pressure");

  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(0);

  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("pressure"));

  // constant function
  ConstFunction<Number> twof(diri_bc_val);

  DirichletBoundary diri_bc(boundary_ids, vars, &twof);

  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

void solve_eq_sys(EquationSystems &eq_sys, std::vector<NetNode> &nodes,
                  std::vector<Segment> &elems, std::string mesh_file,
                  std::string sim_file) {

  // get flow system
  auto &pres = eq_sys.get_system<TransientLinearImplicitSystem>("Pressure");

  // get mesh
  MeshBase &mesh = eq_sys.get_mesh();

  // get parameters
  const Real dt = eq_sys.parameters.get<Real>("time_step");
  const Real init_time = eq_sys.parameters.get<Real>("init_time");
  const unsigned int num_steps = eq_sys.parameters.get<unsigned int>("num_steps");

  eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = 10;

  // to compute the nonlinear convergence
  UniquePtr<NumericVector<Number>> last_nonlinear_soln(pres.solution->clone());

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
    pres.time = current_time;
    eq_sys.parameters.set<Real>("time") = pres.time;

    *pres.old_local_solution = *pres.current_local_solution;

    // solve
    pres.solve();

    // write
    {
      // write solution to file
      char sim_file_i[400];
      sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
      ExodusII_IO exo(mesh);
      exo.write_timestep(sim_file_i, eq_sys, t_step+1, pres.time);
    }

  } while (t_step < num_steps);
}

} // namespace

//
// Test
//
void test::fv::add_and_delete(int argc, char **argv, Parallel::Communicator *comm,
    double pressure_in, bool branch) {

  out << "********** Finite Volume: Add and Delete **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "fv2/mesh";
  std::string sim_file = "fv2/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<NetNode> nodes;
  std::vector<Segment> elems;

  // create two nodes
  {
    auto n1 = create_node(Point(0., 0., 0.), nodes, mesh);
    auto n2 = create_node(Point(2., 0., 0.), nodes, mesh);
    create_elem(n1, n2, elems, mesh);

    n1 = n2;
    n2 = create_node(Point(4., 0., 0.), nodes, mesh);
    create_elem(n1, n2, elems, mesh);

    n1 = n2;
    n2 = create_node(Point(6., 0., 0.), nodes, mesh);
    create_elem(n1, n2, elems, mesh);

    if (!branch) {
      n1 = n2;
      n2 = create_node(Point(8., 0., 0.), nodes, mesh);
      create_elem(n1, n2, elems, mesh);
    } else {

      unsigned int n_branch = n1;
      n2 = create_node(Point(4., 2., 0.), nodes, mesh);
      create_elem(n_branch, n2, elems, mesh);
    }

    // update connectivity
    update_connectivity(nodes, elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // assign boundary ids
  {
    auto &bd_info = mesh.get_boundary_info();

    // add 1st node to boundary 0 and 2nd to 1
    bd_info.add_node(mesh.node_ptr(0), 0);

    // add element's left side to boundary
    bd_info.add_side(mesh.elem_ptr(0), 0, 0);
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &pres = eq_sys.add_system<TransientLinearImplicitSystem>("Pressure");
  pres.add_variable("pressure", CONSTANT, MONOMIAL);
  pres.time = 0.0;

  auto &flux = eq_sys.add_system<TransientLinearImplicitSystem>("Flux");
  flux.add_variable("flux", FIRST);
  flux.time = 0.0;

  // set some parameters
  force_c = 1.0;
  coeff_a = 1.0;

  // attach assembly function
  pres.attach_assemble_function(assemble_net);

  // Add boundary condition
  p_in = pressure_in;
  boundary_condition(eq_sys);

  // add parameters
  eq_sys.parameters.set<Real>("init_time") = 0.;
  int num_steps = 4;
  eq_sys.parameters.set<Real>("time_step") = 0.01 / num_steps;
  eq_sys.parameters.set<unsigned int>("num_steps") = num_steps;
  eq_sys.parameters.set<std::vector<NetNode> *>("net_nodes") = &nodes;
  eq_sys.parameters.set<std::vector<Segment> *>("net_elems") = &elems;

  // init
  eq_sys.init();

  // set Petsc matrix option to suppress the error
  PetscMatrix<Number> *pet_mat =
      dynamic_cast<PetscMatrix<Number>*>(pres.matrix);
  MatSetOption( pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, pres.time);
  }

  // Solve the final system
  solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);


  //
  // Add new element at mid point of element 3 (id 2 for zero based counting)
  //       n1        2         n2
  //       o-------------------o    <-- parent element 3
  //
  //
  //      n1   2   n2     3      n3
  //       o--------x----------o
  //                |
  //                |   4
  //                |
  //                |
  //                x n4
  //
  {
    // get element 3
    auto &parent_elem = elems[2];

    // get node n1 and n2
    auto n1 = parent_elem.d_n1;
    auto n2 = parent_elem.d_n2;

    // create node n3 and n4
    // get mid point of the element
    auto p_n1 = mesh.point(parent_elem.d_n1);
    auto p_n2_old = mesh.point(parent_elem.d_n2);
    auto p_n2_new = util::point_on_line(p_n1, p_n2_old, 0.5);

    // move node n2 from old position to new position
    mesh.node(n2) += p_n2_new - p_n2_old;

    auto n3 = create_node(p_n2_old, nodes, mesh);
    auto n4 = create_node(p_n2_new + Point(1., 0., 1.), nodes, mesh);

    // add new segment between n2 - n3 and n2 - n4
    create_elem(n2, n3, elems, mesh);
    create_elem(n2, n4, elems, mesh);

    // update connectivity
    update_connectivity(nodes, elems, mesh);

    // clear old boundary info in mesh
    // assign boundary ids
    if (false) {
      auto &bd_info = mesh.get_boundary_info();
      bd_info.clear();

      // add 1st node to boundary 0 and 2nd to 1
      bd_info.add_node(mesh.node_ptr(0), 0);

      // add element's left side to boundary
      bd_info.add_side(mesh.elem_ptr(0), 0, 0);
    }

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // new files for output
    mesh_file = "fv2/mesh_2";
    sim_file = "fv2/sim_2";

    // write this mesh to file
    mesh.write(mesh_file + ".e");

    // update equation system
    eq_sys.reinit();

    // set Petsc matrix option to suppress the error
    PetscMatrix<Number> *pet_mat =
        dynamic_cast<PetscMatrix<Number>*>(pres.matrix);
    MatSetOption( pet_mat->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }

  // Solve the final system
  std::cout << "\nSolving on new network mesh\n";
  //  exit(1);

  solve_eq_sys(eq_sys, nodes, elems, mesh_file, sim_file);

  // end
}