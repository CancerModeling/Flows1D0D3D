////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"

namespace {

double get_exact_source(const Point &p) {

  return std::sin(2. * M_PI * p(0)) * std::sin(2. * M_PI * p(1)) *
         std::sin(2. * M_PI * p(2));
}
} // namespace

Number netfc::initial_condition_pres(const Point &p, const Parameters &es,
                                     const std::string &system_name,
                                     const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Pressure");

  if (var_name == "pressure") {

    const auto *deck = es.get<netfc::InputDeck *>("input_deck");

    return deck->d_pressure_ic_val;
  }

  return 0.;
}

void netfc::boundary_condition_pres(EquationSystems &es) {

  const auto *deck = es.parameters.get<netfc::InputDeck *>("input_deck");

  std::set<boundary_id_type> ids;
  if (deck->d_nutrient_bc_north)
    ids.insert(2);
  if (deck->d_nutrient_bc_north)
    ids.insert(0);
  if (deck->d_nutrient_bc_north)
    ids.insert(1);
  if (deck->d_nutrient_bc_north)
    ids.insert(3);

  // get variable ids
  auto &sys = es.get_system<TransientLinearImplicitSystem>("Pressure");
  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("pressure"));

  // create constant function for bc and apply
  ConstFunction<Number> cf(deck->d_pressure_bc_val);
  DirichletBoundary diri_bc(ids, vars, &cf);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

void netfc::PressureAssembly::assemble() {

  const auto &deck = d_model_p->get_input_deck();

  assemble_1();
}

void netfc::PressureAssembly::assemble_1() {

  /*  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // FEM parameters
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Finite element type for tumor variable
  FEType fe_type_tum = tum.variable_type(v_tum[0]);
  std::unique_ptr<FEBase> fe_elem_face_tum(FEBase::build(dim, fe_type_tum));
  QGauss qface(dim - 1, fe_type_tum.default_quadrature_order());
  fe_elem_face_tum->attach_quadrature_rule(&qface);

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> &phi_face = fe_elem_face_tum->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi_face =
      fe_elem_face_tum->get_dphi();
  const std::vector<Real> &JxW_face = fe_elem_face_tum->get_JxW();
  const std::vector<Point> &qface_normals = fe_elem_face_tum->get_normals();
  const std::vector<Point> &qface_points = fe_elem_face_tum->get_xyz();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();
  Real mat_coeff = deck.d_tissue_flow_rho * deck.d_tissue_flow_K / deck.d_tissue_flow_mu;
  const double factor_p_t = deck.d_coupling_factor_p_t;

  const Real mesh_size = deck.d_mesh_size_vec[0];
  Real face_area = mesh_size * mesh_size;
  if (mesh.mesh_dimension() == 2)
    face_area = mesh_size;

  Real elem_vol = face_area * mesh_size;

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Number> bc_vals;

  // Looping through elements
  for (const auto & elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

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
        const Real a = factor_p_t * mat_coeff * face_area / dx;

        // add contribution
        // +a to (e,e) where e is id of this element
        // -a to (e, n_e) where n_e is id of neighbor of element
        util::add_unique(dof_indices_pres[0], a, Ke_dof_col, Ke_val_col);
        util::add_unique(dof_indices_pres_neigh[0], -a, Ke_dof_col, Ke_val_col);

        //
        // compute force term due to div(chem_tum * grad(tum)) term
        //
        // we need to compute
        // int_{face} chem_tum * grad(tum) . normal(face) dS
        // Get quadrature points on the face

        // Reinitialize shape functions on the element side
        fe_elem_face_tum->reinit(elem, side);

        // loop over quadrature point
        for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

          Number chem_tum_cur = 0.;
          Gradient tum_grad;

          for (unsigned int l = 0; l < phi_face.size(); l++) {

            chem_tum_cur += phi_face[l][qp] *
                            tum.current_solution(dof_indices_tum_var[1][l]);
            tum_grad.add_scaled(
                dphi_face[l][qp],
                tum.current_solution(dof_indices_tum_var[0][l]));
          }

          Fe(0) += JxW_face[qp] * factor_p_t * mat_coeff * chem_tum_cur *
                   tum_grad * qface_normals[qp];
        }

      } // elem neighbor is not null
      if (false) {
        // element side is in the boundary

        // check if this is -x side

        // build elem side
        std::unique_ptr<const Elem> elem_side (elem->build_side_ptr(side));

        // get center of the element
        const auto elem_side_center = elem_side->centroid();

        auto dx = (elem_side_center - elem_center).norm();

        // get coefficient
        const Real a = factor_p_t * mat_coeff * face_area / dx;

        util::add_unique(dof_indices_pres[0], a, Ke_dof_col, Ke_val_col);

        if (false) {
          double bc = 0.0;

          // std::cout << "elem_side_center: " << elem_side_center << std::endl;

          double dist = std::sqrt(
              ((elem_side_center(0) - 0.5) * (elem_side_center(0) - 0.5)) +
              ((elem_side_center(1) - 0.5) * (elem_side_center(1) - 0.5)));

          // std::cout << "dist: " << dist << std::endl;

          if (dist < 0.1) {

            bc = 0.1 / 1.1 * (1.0 + elem_side_center(2));

          } else {

            bc = 0.1 / 1.1 * (1.0 + elem_side_center(2)) *
                 (1.0 - 0.1 * std::log(dist / 0.1));
          }

          // std::cout << "bc: " << bc << std::endl;

          Fe(0) += bc * a;
          //}
        }
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

  //
  // Handle coupling between tissue pressure and vessel pressure
  //
  {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    while (pointer) {
      
    std::cout << " " <<std::endl;
    std::cout << " " <<std::endl;
    std::cout << " " <<std::endl;
    std::cout << "index: " << pointer->index <<std::endl;

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;

      // get element data at points on cylinder surface
      const auto J_b_points; // = network.compute_elem_weights_at_node(pointer);

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = J_b_points[i];

       // std::cout << "Neighbor: " << i <<std::endl;

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = mesh.elem_ptr(e_id);
          pres_map.dof_indices(elem, dof_indices_pres, v_pres);
          double p_t_k = pres.current_solution(dof_indices_pres[0]);

          std::vector<unsigned int> Ke_dof_row;
          Ke_dof_row.push_back(dof_indices_pres[0]);
          DenseMatrix<Number> Ke(1,1);
          DenseVector<Number> Fe(1);

          
          std::cout << "e_id: " << e_id <<std::endl;
          std::cout << "L_p: " << pointer->L_p[i] <<std::endl;
          std::cout << "J_b_data.half_cyl_surf: " << J_b_data.half_cyl_surf <<std::endl;
          std::cout << "e_w: " << e_w <<std::endl;

          // implicit for p_t in source
          Ke(0,0) += factor_p_t * deck.d_coupling_method_theta * pointer->L_p[i] *
              J_b_data.half_cyl_surf *
              e_w;

          // explicit for p_v in source
          Fe(0) += factor_p_t * pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                   (p_v_k - (1. - deck.d_coupling_method_theta) * p_t_k);

          // update matrix
          pres.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);
          pres.rhs->add_vector(Fe, dof_indices_pres);
        }
      } // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }

  // finish
  pres.matrix->close();
  pres.rhs->close();

  if (mesh.n_elem() < 5) {
    out << "Diri rows\n";
    out << util::io::printStr(bc_rows);

    out << "number of elements = " << mesh.n_elem() << "\n";
    pres.rhs->print(out);
    out << "matrix:\n";
    pres.matrix->print(out);
  } */
}
