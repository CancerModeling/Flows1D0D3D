////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"

namespace {
double get_nut_source(const std::string &test_name, const Point &x,
                      const std::vector<double> &x0, const double &r) {

  if (test_name != "test_tum_2")
    return 0.;

  // For cylinder, take z-coord as 0
  const Point xc = util::to_point({x0[0], x0[1], 0.});

  // spherical source
  Point dx = Point(x(0), x(1), 0.) - xc;
  if (dx.norm() < r)
    return 1.;

  return 0.;
}
} // namespace

Number netpresnut::initial_condition_nut(const Point &p, const Parameters &es,
                                      const std::string &system_name,
                                      const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<InpDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void netpresnut::boundary_condition_nut(EquationSystems &es) {

  const auto *deck = es.parameters.get<InpDeck *>("input_deck");

  std::set<boundary_id_type> ids;
  if (deck->d_nutrient_bc_north)
    ids.insert(2);
  if (deck->d_nutrient_bc_south)
    ids.insert(0);
  if (deck->d_nutrient_bc_east)
    ids.insert(1);
  if (deck->d_nutrient_bc_west)
    ids.insert(3);

  // get variable ids
  auto &sys = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  std::vector<unsigned int> vars;
  vars.push_back(sys.variable_number("nutrient"));

  // create constant function for bc and apply
  ConstFunction<Number> one(1);
  DirichletBoundary diri_bc(ids, vars, &one);
  sys.get_dof_map().add_dirichlet_boundary(diri_bc);
}

// Assembly class
void netpresnut::NutAssembly::assemble() {
  assemble_1d_coupling();
  assemble_face();
  assemble_1();
}

void netpresnut::NutAssembly::assemble_1d_coupling() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

  // coupling between tissue pressure and vessel pressure
  {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    DenseMatrix<Number> Ke(1, 1);
    DenseVector<Number> Fe(1);

    double h_3D = network.h_3D;
    int N_3D = network.N_3D;
    double L_p = network.L_p;
    double L_s = network.L_s;
    double osmotic_sigma = network.osmotic_sigma;

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
    double c_v = 0.;
    double p_t = 0.;
    unsigned int assembly_cases;
    bool dirichlet_fixed = false;

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
        c_v = network.C_v[nodes[node_proc]];
        assembly_cases = network.d_vertexBdFlag[nodes[node_proc]];

        dirichlet_fixed = false;
        if (assembly_cases & UNET_NUT_BDRY_ARTERY_INLET or
            assembly_cases & UNET_NUT_BDRY_VEIN_INLET)
          dirichlet_fixed = true;

        if (dirichlet_fixed == false) {

          // Surface area of cylinder
          surface_area = 2.0 * M_PI * (0.5 * length) * radius;
          util::unet::determineWeightsAndIds(
              deck.d_num_points_length, deck.d_num_points_angle, N_3D,
              coords[node_proc], coords[node_neigh], radius, h_3D, 0.5 * length,
              weights, id_3D_elements,
              deck.d_coupling_3d1d_integration_method, d_mesh, true);

          // Add coupling entry
          numberOfElements = id_3D_elements.size();

          for (int j = 0; j < numberOfElements; j++) {

            if (id_3D_elements[j] > -1) {

              const auto *elem = d_mesh.elem_ptr(id_3D_elements[j]);
              if (elem->processor_id() == d_model_p->get_comm()->rank()) {

                // get 3d pressure
                pres.init_dof(elem);
                p_t = pres.get_current_sol(0);

                // init nutrient dof
                init_dof(elem);

                // implicit for c_t in source
                Ke(0, 0) = dt * factor_nut * L_s * surface_area * weights[j];

                // explicit for c_v in source
                Fe(0) = dt * factor_nut * L_s * surface_area * weights[j] * c_v;

                // osmotic reflection term
                if (p_v - p_t > 0.0) {

                  // 3D equation
                  // 2pi R (p_v - p_t) phi_v term in right hand side of 3D
                  // equation
                  Fe(0) += dt * factor_nut * (1. - osmotic_sigma) * L_p *
                           surface_area * weights[j] * (p_v - p_t) * c_v;

                } else {

                  // 3D equation
                  // 2pi R (p_v - p_t) phi_sigma term in right hand side of 3D
                  // equation
                  Ke(0, 0) += -dt * factor_nut * (1. - osmotic_sigma) * L_p *
                              surface_area * weights[j] * (p_v - p_t);
                }

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

void netpresnut::NutAssembly::assemble_face() {

  // Get required system alias
  // auto &nut = d_model_p->get_nut_assembly();
  auto &pres = d_model_p->get_pres_assembly();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = d_model_p->d_dt;
  const Real factor_nut = deck.d_assembly_factor_c_t;

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
  std::vector<unsigned int> dof_indices_nut_neigh;
  std::vector<unsigned int> dof_indices_pres_neigh;

  // Store current and old solution
  Real pres_cur = 0.;

  // Store current and old solution of neighboring element
  Real pres_neigh_cur = 0.;

  // Looping through elements
  for (const auto &elem : d_mesh.active_local_element_ptr_range()) {

    init_dof(elem);
    pres.init_dof(elem);

    // reset matrix and force
    Ke_dof_col.clear();
    Ke_val_col.clear();

    Ke_dof_row[0] = get_global_dof_id(0);
    Fe(0) = 0.;

    // get solution in this element
    pres_cur = pres.get_current_sol(0);

    // loop over sides of the element
    for (auto side : elem->side_index_range()) {

      if (elem->neighbor_ptr(side) != nullptr) {

        const Elem *neighbor = elem->neighbor_ptr(side);

        // get dof id
        // nut
        init_dof(neighbor, dof_indices_nut_neigh);

        // pres
        pres.init_dof(neighbor, dof_indices_pres_neigh);
        pres_neigh_cur = pres.get_current_sol(0, dof_indices_pres_neigh);

        // diffusion
        const Real a_diff = factor_nut * dt * deck.d_D_sigma * deck.d_face_by_h;
        util::add_unique(get_global_dof_id(0), a_diff, Ke_dof_col, Ke_val_col);
        util::add_unique(dof_indices_nut_neigh[0], -a_diff, Ke_dof_col,
                         Ke_val_col);

        // advection (pressure gradient term)
        Real v = deck.d_tissue_flow_coeff * deck.d_face_by_h *
                 (pres_cur - pres_neigh_cur);

        // upwinding
        if (v >= 0.)
          util::add_unique(get_global_dof_id(0), factor_nut * dt * v,
                           Ke_dof_col, Ke_val_col);
        else
          util::add_unique(dof_indices_nut_neigh[0], factor_nut * dt * v,
                           Ke_dof_col, Ke_val_col);

      } // elem neighbor is not null
    }   // loop over faces

    // add to matrix
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    d_sys.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    d_sys.rhs->add_vector(Fe, Ke_dof_row);
  }
}

void netpresnut::NutAssembly::assemble_1() {

  // finish
  d_sys.matrix->close();
  d_sys.rhs->close();
}