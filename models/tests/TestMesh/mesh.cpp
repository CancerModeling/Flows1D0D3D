////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"
#include "utils.hpp"

static double diri_bc_val = 0.;

// local anonymous namespace
namespace {

void random_init() { srand(time(nullptr)); }

// random number between 0 and 1
double get_random_number() { return double(std::rand()) / double(RAND_MAX); }

// random ppint with each element between 0 and 1
Point get_random_point(double scale = 1.) {
  return {scale * double(std::rand()) / double(RAND_MAX),
          scale * double(std::rand()) / double(RAND_MAX),
          scale * double(std::rand()) / double(RAND_MAX)};
}

std::string print_pt(const Point &p) {

  std::ostringstream oss;
  oss << "(" << p(0) << ", " << p(1) << ", " << p(2) << ")";

  return oss.str();
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

Number ic_net(const Point &p, const Parameters &es,
              const std::string &system_name, const std::string &var_name) {

  if (var_name == "pressure") {

    if (p.norm() < 1.0E-8)
      return 2.;
    else
      return 0.;
  }

  return 0.;
}

unsigned int create_node(Point p, std::vector<unsigned int> &nodes,
    ReplicatedMesh &mesh) {

  auto node = mesh.add_point(p);
  nodes.push_back(node->id());

  return node->id();
}

unsigned int create_elem(unsigned int node_1, unsigned int node_2,
    std::vector<unsigned int> &elems, ReplicatedMesh &mesh) {

  auto elem = Elem::build(EDGE2).release();
  elem->set_node(0) = mesh.node_ptr(node_1);
  elem->set_node(1) = mesh.node_ptr(node_2);
  auto add_elem = mesh.add_elem(elem);
  elems.push_back(add_elem->id());

  return add_elem->id();
}

} // namespace

//
// Test: Start with empty mesh and successively add new nodes and elements
//
void test::mesh::add_node_elem_test(int argc, char **argv,
                                    Parallel::Communicator *comm) {

  out << "********** TestMesh 1 **************\n";

  // test 1 (start with empty mesh and successively add new nodes and
  // elements)

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t1/mesh";
  std::string sim_file = "t1/sim";

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

  // add new nodes and segments
  for (int i = 1; i <= 4; i++) {
    auto add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[1], nodes[2]);

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);
  }

  // attach a branch to the ending node
  size_t el_size = elems.size();
  for (int i = el_size; i < el_size + 1; i++) {

    auto add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 2]);
    add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[1], nodes[2]);

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);
  }
  // end
}

//
// Test: Same as add_node_elem_test() but now we add equation
// system for pressure
//
void test::mesh::add_node_elem_eq_sys_test(int argc, char **argv,
                                           Parallel::Communicator *comm) {

  out << "********** TestMesh 2 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t2/mesh";
  std::string sim_file = "t2/sim";

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

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;
  eq_sys.init();

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // add new nodes and segments
  for (int i = 1; i <= 4; i++) {
    auto add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[1], nodes[2]);

    // reinit equation system
    eq_sys.reinit();

    // increment time
    flow.time = i;
    // set the values of pressure at new dofs
    for (const auto &node : mesh.local_node_ptr_range()) {

      if (i != 1 && node->id() != nodes[i + 1])
        continue;

      auto j = node->dof_number(0, 0, 0);
      flow.solution->set(j, (j + 1) * 10.);
    }
    flow.solution->close();
    flow.update();

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);

    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), i + 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // attach a branch to the ending node
  size_t el_size = elems.size();
  for (int i = el_size; i < el_size + 1; i++) {

    auto add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    add_node = mesh.add_point(get_random_point());
    nodes.push_back(add_node->id());

    auto elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 1]);
    auto add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    elem = Elem::build(EDGE2).release();
    elem->set_node(0) = mesh.node_ptr(nodes[i]);
    elem->set_node(1) = mesh.node_ptr(nodes[i + 2]);
    add_elem = mesh.add_elem(elem);
    elems.push_back(add_elem->id());

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[1], nodes[2]);

    // reinit equation system
    eq_sys.reinit();

    // increment time
    flow.time = i;
    // set the values of pressure at new dofs
    for (const auto &node : mesh.local_node_ptr_range()) {

      if (node->id() != nodes[i + 1] || node->id() != nodes[i + 2])
        continue;

      auto j = node->dof_number(0, 0, 0);
      flow.solution->set(j, (j + 1) * 10.);
    }
    flow.solution->close();
    flow.update();

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);

    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), i + 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // end
}

//
// Test: Same as add_node_elem_test() but now we add equation
// system for pressure
//
void test::mesh::add_node_elem_eq_sys_test_2(int argc, char **argv,
                                           Parallel::Communicator *comm) {

  out << "********** TestMesh 2 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t2/mesh";
  std::string sim_file = "t2/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);

  // create list of nodes and segments
  std::vector<unsigned int> nodes;
  std::vector<unsigned int> elems;
  std::vector<Segment> segments;

  // create two nodes at point (0,0,0) and (1,1,1) and element between these
  // nodes
  {
    auto node_1 = create_node(Point(0., 0., 0.), nodes, mesh);
    auto node_2 = create_node(Point(1., 1., 1.), nodes, mesh);

    auto elem_1 = create_elem(node_1, node_2, elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(node_1, node_2);

    // write this mesh to file
    mesh.write(mesh_file + "_1.e");
  }

  // add equation system for pressure
  EquationSystems eq_sys(mesh);
  auto &flow = eq_sys.add_system<TransientLinearImplicitSystem>("Flow");
  flow.add_variable("pressure", FIRST);
  flow.time = 0.0;
  eq_sys.init();

  {
    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // add new nodes and segments
  for (int i = 1; i <= 4; i++) {

    auto node_i = create_node(get_random_point(), nodes, mesh);

    auto elem_i = create_elem(nodes[i], nodes[i+1], elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[i], nodes[i+1]);

    // reinit equation system
    eq_sys.reinit();

    // increment time
    flow.time = i;
    // set the values of pressure at new dofs
    for (const auto &node : mesh.local_node_ptr_range()) {

      if (i != 1 && node->id() != nodes[i + 1])
        continue;

      auto j = node->dof_number(0, 0, 0);
      flow.solution->set(j, (j + 1) * 10.);
    }
    flow.solution->close();
    flow.update();

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);

    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), i + 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // attach a branch to the ending node
  size_t el_size = elems.size();
  for (int i = el_size; i < el_size + 1; i++) {

    auto node_i = create_node(get_random_point(), nodes, mesh);
    auto node_j = create_node(get_random_point(), nodes, mesh);

    auto elem_i = create_elem(nodes[i], nodes[i+1], elems, mesh);
    auto elem_j = create_elem(nodes[i], nodes[i+2], elems, mesh);

    // finish adding nodes and elements
    mesh.prepare_for_use();

    // create a segment
    segments.emplace_back(nodes[1], nodes[2]);

    // reinit equation system
    eq_sys.reinit();

    // increment time
    flow.time = i;
    // set the values of pressure at new dofs
    for (const auto &node : mesh.local_node_ptr_range()) {

      if (node->id() != nodes[i + 1] || node->id() != nodes[i + 2])
        continue;

      auto j = node->dof_number(0, 0, 0);
      flow.solution->set(j, (j + 1) * 10.);
    }
    flow.solution->close();
    flow.update();

    // write this mesh to file
    std::string mesh_file_i = mesh_file + "_" + std::to_string(i + 1) + ".e";
    mesh.write(mesh_file_i);

    // write solution to file
    char sim_file_i[400];
    sprintf(sim_file_i, "%s_%d.e", sim_file.c_str(), i + 1);
    ExodusII_IO exo(mesh);
    exo.write_timestep(sim_file_i, eq_sys, 1, flow.time);
  }

  // end
}

//
// Test: Same as add_node_elem_test() but now we add equation
// system for pressure
//
void test::mesh::elem_id_numbering(int argc, char **argv,
                                           Parallel::Communicator *comm) {

  out << "********** TestMesh 10 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t10/mesh";
  std::string sim_file = "t10/sim";

  // dummy mesh
  ReplicatedMesh mesh(*comm);
  unsigned int N = 4;
  double h = 1./ N;
  MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                        1., 0., 1., HEX8);

  // Looping through elements
  out << "Element information:\n";
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    const Point xc = elem->centroid();
    unsigned int id_guess = util::get_elem_id(xc, h, N, 3);

    out << "id: " << elem->id()
        << ", id_guess: " << id_guess
        << ", center: " << print_pt(xc);

    // try 4 perturbations of point
    out << ", rdm pertb. guess: ";
    for (unsigned int i=0; i<4; i++) {
      Point x = xc + get_random_point(0.8 * h);
      unsigned int id_guess = util::get_elem_id(xc, h, N, 3);

      out << id_guess << " ";
    }
    out << "\n\n";
  }

  // Looping over nodes
  out << "Node information:\n";
  for (const auto &node : mesh.node_ptr_range()) {

    const Point& x = mesh.point(node->id());
    unsigned int id_guess = 0;
    {
      unsigned int i = x(0) / h;
      unsigned int j = x(1) / h;
      unsigned int k = x(2) / h;
      id_guess = k * N * N + j * N + i;
    }

    out << "id: " << node->id()
        << ", id_guess: " << id_guess
        << ", coord: " << print_pt(x)
        << "\n";
  }

}


