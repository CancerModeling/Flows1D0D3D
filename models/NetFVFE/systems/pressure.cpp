////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {

double get_exact_source(const Point &p) {

  return std::sin(2. * M_PI * p(0)) * std::sin(2. * M_PI * p(1)) *
         std::sin(2. * M_PI * p(2));
}
} // namespace

Number netfvfe::initial_condition_pres(const Point &p, const Parameters &es,
                                       const std::string &system_name,
                                       const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Pressure");

  if (var_name == "pressure") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_pressure_ic_val;
  }

  return 0.;
}

void netfvfe::boundary_condition_pres(EquationSystems &es) {

  const auto *deck = es.parameters.get<InpDeck *>("input_deck");

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

// Assembly class

void netfvfe::PressureAssembly::assemble() {
  assemble_1d_coupling();
  assemble_face();
  assemble_1();
}

void netfvfe::PressureAssembly::assemble_1d_coupling() {

  // Get required system alias
  // auto &pres = d_model_p->get_pres_assembly();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();
  const double factor_p = deck.d_assembly_factor_p_t;

  // coupling between tissue pressure and vessel pressure
  {
    const auto &network = d_model_p->get_network();

    DenseMatrix<Number> Ke(1, 1);
    DenseVector<Number> Fe(1);

    double h_3D = network.h_3D;
    int N_3D = network.N_3D;
    double L_p = network.L_p;

    int node_proc = 0;
    int node_neigh = 1;
    std::vector<unsigned int> nodes(2, 0);
    std::vector<std::vector<double>> coords;
    coords.emplace_back(3, 0.);
    coords.emplace_back(3, 0.);
    double radius = 0.;
    double length = 0.;
    double surface_area = 0.;
    std::vector<double> weights;
    std::vector<int> id_3D_elements;
    int numberOfElements = 0;
    double p_v = 0.;
    unsigned int assembly_cases;

    for (unsigned int i = 0; i < network.d_numSegments; i++) {

      nodes[0] = network.d_segments[2 * i + 0];
      nodes[1] = network.d_segments[2 * i + 1];
      radius = network.d_segmentData[i];
      for (unsigned int j = 0; j < 3; j++) {
        coords[0][j] = network.d_vertices[3 * nodes[0] + j];
        coords[1][j] = network.d_vertices[3 * nodes[1] + j];
      }
      length = util::dist_between_points(coords[0], coords[1]);

      for (unsigned int j = 0; j < 2; j++) {

        node_proc = 0;
        node_neigh = 1;
        if (j == 1) {
          node_proc = 1;
          node_neigh = 0;
        }

        p_v = network.P_v[nodes[node_proc]];
        assembly_cases = network.d_vertexBdFlag[nodes[node_proc]];

        if (!(assembly_cases & UNET_PRES_BDRY_DIRIC)) {

          // Surface area of cylinder
          surface_area = 2.0 * M_PI * (0.5 * length) * radius;
          util::unet::determineWeightsAndIds(
            deck.d_num_points_length, deck.d_num_points_angle, N_3D, coords[node_proc],
            coords[node_neigh], radius, h_3D, 0.5 * length, weights,
            id_3D_elements,
            deck.d_coupling_3d1d_integration_method, d_mesh, true);

          // Add coupling entry
          numberOfElements = id_3D_elements.size();

          for (int k = 0; k < numberOfElements; k++) {

            if (id_3D_elements[k] > -1) {

              const auto *elem = d_mesh.elem_ptr(id_3D_elements[k]);
              if (elem->processor_id() == d_model_p->get_comm()->rank()) {
                init_dof(elem);

                // implicit for p_t in source
                Ke(0, 0) = factor_p * L_p * surface_area * weights[k];

                // explicit for p_v in source
                Fe(0) = factor_p * L_p * surface_area * weights[k] * p_v;

                // update matrix
                d_sys.matrix->add_matrix(Ke, d_dof_indices_sys,
                                         d_dof_indices_sys);
                d_sys.rhs->add_vector(Fe, d_dof_indices_sys);
              }
            }
          } // loop over 3D elements
        }   // if not dirichlet
      }     // segment's node loop

    } // loop over segments
  }
}

void netfvfe::PressureAssembly::assemble_face() {

  // Get required system alias
  auto &pro = d_model_p->get_pro_assembly();
  auto &hyp = d_model_p->get_hyp_assembly();
  auto &nut = d_model_p->get_nut_assembly();
  auto &ecm = d_model_p->get_ecm_assembly();

  // Parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real factor_p = deck.d_assembly_factor_p_t;

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Real> bc_vals;

  // to store pair of column dof and row-column matrix value
  std::vector<unsigned int> Ke_dof_row(1, 0);
  std::vector<Real> Ke_val_col;
  std::vector<unsigned int> Ke_dof_col;

  // local matrix and vector
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe(1);

  // store neighboring element's dof information
  std::vector<unsigned int> dof_indices_pres_neigh;

  // Store current and old solution
  Real nut_cur = 0.;
  Real ecm_cur = 0.;
  Real chem_pro_cur = 0.;
  Real chem_hyp_cur = 0.;

  Real nut_proj = 0.;
  Real ecm_proj = 0.;

  Gradient pro_grad = 0.;
  Gradient hyp_grad = 0.;

  Gradient Sp = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    pro.init_dof(elem);
    hyp.init_dof(elem);
    nut.init_dof(elem);
    ecm.init_dof(elem);

    // reset matrix anf force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get finite-volume quantities
    nut_cur = nut.get_current_sol(0);
    nut_proj = util::project_concentration(nut_cur);

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);

        // get dof id
        // pres
        init_dof(neighbor, dof_indices_pres_neigh);

        // get coefficient
        const Real a_diff = factor_p * deck.d_tissue_flow_coeff * deck.d_face_by_h;

        // diffusion
        util::add_unique(get_global_dof_id(0), a_diff, Ke_dof_col, Ke_val_col);
        util::add_unique(dof_indices_pres_neigh[0], -a_diff, Ke_dof_col,
                         Ke_val_col);

        // div(chem_tum * grad(tum)) term
        // requires integration over face of an element
        pro.d_fe_face->reinit(elem, side);

        // loop over quadrature points
        for (unsigned int qp = 0; qp < pro.d_qrule_face.n_points(); qp++) {

          chem_pro_cur = 0.;
          chem_hyp_cur = 0.;
          pro_grad = 0.;
          hyp_grad = 0.;
          ecm_cur = 0.;
          ecm_proj = 0.;
          for (unsigned int l = 0; l < pro.d_phi_face.size(); l++) {

            chem_pro_cur +=
              pro.d_phi_face[l][qp] * pro.get_current_sol_var(l, 1);

            pro_grad.add_scaled(pro.d_dphi_face[l][qp],
                                pro.get_current_sol_var(l, 0));

            chem_hyp_cur +=
              pro.d_phi_face[l][qp] * hyp.get_current_sol_var(l, 1);

            hyp_grad.add_scaled(pro.d_dphi_face[l][qp],
                                hyp.get_current_sol_var(l, 0));


            ecm_cur +=
              pro.d_phi_face[l][qp] * ecm.get_current_sol(l);
          }

          ecm_proj = util::project_concentration(ecm_cur);

          Sp = (chem_pro_cur + deck.d_chi_c * nut_cur +
                deck.d_chi_h * ecm_cur) *
                 pro_grad +
               (chem_hyp_cur + deck.d_chi_c * nut_cur +
                deck.d_chi_h * ecm_cur) *
                 hyp_grad;

          // add to force
          Fe(0) += -factor_p * pro.d_JxW_face[qp] * deck.d_tissue_flow_coeff * Sp * pro.d_qface_normals[qp];
        } // loop over quadrature points on face

      } // elem neighbor is not null
    }   // loop over faces

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  } // element loop
}

void netfvfe::PressureAssembly::assemble_1() {

  // Get required system alias
  // auto &pres = d_model_p->get_pres_assembly();

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}