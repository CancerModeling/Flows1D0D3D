////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"
#include "utilLibs.hpp"
#include <utils.hpp>
// using namespace libMesh;

static double diri_bc_val = 0.;

// local anonymous namespace
namespace {

void assemble_pres(EquationSystems &es, const std::string &system_name) {

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();

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
    const unsigned int n_dofs = dof_indices_pres.size();

    // Arranging matrix
    std::vector<Real> Ke_val_col;
    std::vector<unsigned int> Ke_dof_col;
    std::vector<unsigned int> Ke_dof_row;
    Ke_dof_row.push_back(dof_indices_pres[0]);
    DenseVector<Number> Fe;
    Fe.resize(n_dofs);

    // loop over sides of the element
    const auto elem_center = elem->centroid();
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);
        const auto neighbor_center = neighbor->centroid();
        const auto dx = (elem_center - neighbor_center).norm();

        // get dof id
        std::vector<unsigned int> dof_indices_pres_neigh;
        pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);

        // get coefficient
        const Real a = 1. / dx;

        // add contribution
        // +a to (e,e) where e is id of this element
        // -a to (e, n_e) where n_e is id of neighbor of element
        util::add_unique(dof_indices_pres[0], a, Ke_dof_col, Ke_val_col);
        util::add_unique(dof_indices_pres_neigh[0], -a, Ke_dof_col, Ke_val_col);
      }
    } // loop over faces

    // add to matrix
    DenseMatrix<Number> Ke;
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    pres.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    pres.rhs->add_vector(Fe, dof_indices_pres);
  } // element loop

  // finish
  pres.matrix->close();
  pres.rhs->close();
}

class Model;

class PressureAssembly : public System::Assembly {
public:
  PressureAssembly(Model *model, const std::string system_name)
      : d_model_p(model), d_sys_name(system_name) {}
  void assemble() override;

private:
  std::string d_sys_name;
  Model *d_model_p;
};

class Model {

public:
  Model(int argc, char **argv, Parallel::Communicator *comm);

  Parallel::Communicator *d_comm_p;
  ReplicatedMesh d_mesh;
  EquationSystems d_eq_sys;
  PressureAssembly d_pres_assembly;
};

Model::Model(int argc, char **argv, Parallel::Communicator *comm)
    : d_comm_p(comm), d_mesh(ReplicatedMesh(*d_comm_p)), d_eq_sys(d_mesh),
      d_pres_assembly(this, "Pressure") {

  out << "********** TestMesh 11 **************\n";

  // for debug output
  std::ostringstream oss;

  std::string mesh_file = "t11/mesh";
  std::string sim_file = "t11/sim";

  // dummy mesh
  unsigned int num_elems = 30;
  MeshTools::Generation::build_cube(d_mesh, num_elems, num_elems, num_elems, 0.,
                                    1., 0., 1., 0., 1., HEX8);

  // create equation system
  auto &pres = d_eq_sys.add_system<TransientLinearImplicitSystem>("Pressure");
  pres.add_variable("pressure", CONSTANT, MONOMIAL);

  // method 1
  pres.attach_assemble_function(assemble_pres);

  // method 2
  //pres.attach_assemble_object(d_pres_assembly);

  d_eq_sys.init();

  // solve
  {
    *pres.old_local_solution = *pres.current_local_solution;
    d_eq_sys.parameters.set<Real>("linear solver tolerance") = 1.0e-6;

    for (unsigned int l = 0; l < 1000; ++l) {

      out << "step: " << l << "\n";
      pres.solve();
    }
  }
}

void PressureAssembly::assemble() {

  // return;
  assemble_pres(d_model_p->d_eq_sys, "Pressure");
}
} // namespace

void test::mesh::memory_leak(int argc, char **argv,
                             Parallel::Communicator *comm) {

  // auto model = Model(argc, argv, comm);
  auto model = Model(argc, argv, comm);

  return;
}
