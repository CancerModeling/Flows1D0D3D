////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "../model.hpp"
#include "systems.hpp"

Number netfvfe::initial_condition_nut(const Point &p, const Parameters &es,
                                     const std::string &system_name,
                                     const std::string &var_name) {

  libmesh_assert_equal_to(system_name, "Nutrient");

  if (var_name == "nutrient") {

    const auto *deck = es.get<netfvfe::InputDeck *>("input_deck");

    return deck->d_nut_ic_value;
  }

  return 0.;
}

void netfvfe::boundary_condition_nut(EquationSystems &es) {

  const auto *deck = es.parameters.get<netfvfe::InputDeck *>("input_deck");

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

void netfvfe::NutAssembly::assemble() {

    const auto &deck = d_model_p->get_input_deck();

    if (deck.d_assembly_method == 1)
        assemble_1();
    else if (deck.d_assembly_method == 2)
        assemble_2();
    else if (deck.d_assembly_method == 3)
      assemble_3();
}

void netfvfe::NutAssembly::assemble_1() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Nutrient system
  auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // TAF system
  auto &taf = es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // ECM system
  auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
  const unsigned int v_ecm = ecm.variable_number("ecm");
  const DofMap &ecm_map = ecm.get_dof_map();
  std::vector<unsigned int> dof_indices_ecm;

  // MDE system
  auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
  const unsigned int v_mde = mde.variable_number("mde");
  const DofMap &mde_map = mde.get_dof_map();
  std::vector<unsigned int> dof_indices_mde;

  // Velocity system
  auto &vel = es.get_system<TransientLinearImplicitSystem>("Velocity");
  std::vector<unsigned int> v_vel(2);
  v_vel[0] = vel.variable_number("velocity_x");
  v_vel[1] = vel.variable_number("velocity_y");
  if (dim > 2)
    v_vel.push_back(vel.variable_number("velocity_z"));

  const DofMap &vel_map = vel.get_dof_map();
  std::vector<unsigned int> dof_indices_vel;
  std::vector<std::vector<dof_id_type>> dof_indices_vel_var(2);
  if (dim > 2)
    dof_indices_vel_var.resize(3);

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // FEM parameters
  FEType fe_tum_type = tum.variable_type(0);
  UniquePtr<FEBase> fe_tum(FEBase::build(dim, fe_tum_type));
  QGauss qrule(dim, fe_tum_type.default_quadrature_order());
  fe_tum->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe_tum->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe_tum->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_tum->get_dphi();

  // Finite element type for tumor variable
  std::unique_ptr<FEBase> fe_elem_face_tum(FEBase::build(dim, fe_tum_type));
  QGauss qface(dim - 1, fe_tum_type.default_quadrature_order());
  fe_elem_face_tum->attach_quadrature_rule(&qface);

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> &phi_face = fe_elem_face_tum->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi_face =
      fe_elem_face_tum->get_dphi();
  const std::vector<Real> &JxW_face = fe_elem_face_tum->get_JxW();
  const std::vector<Point> &qface_normals = fe_elem_face_tum->get_normals();
  const std::vector<Point> &qface_points = fe_elem_face_tum->get_xyz();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
  const Real mesh_size = deck.d_mesh_size_vec[0];
  Real face_area = mesh_size * mesh_size;
  if (mesh.mesh_dimension() == 2)
    face_area = mesh_size;

  Real elem_vol = face_area * mesh_size;

  double factor_nut = 1.;
  if (deck.d_tissue_nut_L_s > 1.e-18)
    factor_nut = 1. / deck.d_tissue_nut_L_s;

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Number> bc_vals;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);
    taf_map.dof_indices(elem, dof_indices_taf, v_taf);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    vel_map.dof_indices(elem, dof_indices_vel);
    for (unsigned int var = 0; var < dim; var++)
      vel_map.dof_indices(elem, dof_indices_vel_var[var], v_vel[var]);

    pres_map.dof_indices(elem, dof_indices_pres, v_pres);

    const unsigned int n_dofs = dof_indices_nut.size();

    fe_tum->reinit(elem);

    std::vector<Real> Ke_val_col;
    std::vector<unsigned int> Ke_dof_col;
    std::vector<unsigned int> Ke_dof_row;
    Ke_dof_row.push_back(dof_indices_nut[0]);
    DenseVector<Number> Fe;
    Fe.resize(n_dofs);

    // Finite volume contribution
    {
      const auto elem_center = elem->centroid();

      double elem_p_t_k = pres.current_solution(dof_indices_pres[0]);
      double elem_nut_old = nut.current_solution(dof_indices_nut[0]);

      // loop over sides of the element

      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);
          const auto neighbor_center = neighbor->centroid();
          const auto dx = (elem_center - neighbor_center).norm();

          // get dof id
          std::vector<unsigned int> dof_indices_nut_neigh;
          nut_map.dof_indices(neighbor, dof_indices_nut_neigh, v_nut);

          std::vector<unsigned int> dof_indices_pres_neigh;
          pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);
          double neighbor_p_t_k =
              pres.current_solution(dof_indices_pres_neigh[0]);

          // get coefficient
          const Real a = factor_nut * dt * deck.d_D_sigma * face_area / dx;
          //          if (elem->id() > 100 and elem->id() < 120)
          //            out << "Elem: " << elem->id() << ", neighbor: " << neighbor->id()
          //                << ", diffusion: " << a
          //                << ", dt: " << dt << ", D_s: " << deck.d_D_sigma << "\n";

          // add contribution
          // +a to (e,e) where e is id of this element
          // -a to (e, n_e) where n_e is id of neighbor of element
          util::add_unique(dof_indices_nut[0], a, Ke_dof_col, Ke_val_col);
          util::add_unique(dof_indices_nut_neigh[0], -a, Ke_dof_col,
                           Ke_val_col);

          // advection
          double v = -deck.d_tissue_flow_K * (neighbor_p_t_k - elem_p_t_k) /
                     (deck.d_tissue_flow_mu * mesh_size);
          if (v >= 0.)
            util::add_unique(dof_indices_nut[0],
                             factor_nut * dt * v * face_area, Ke_dof_col,
                             Ke_val_col);
          else
            util::add_unique(dof_indices_nut_neigh[0],
                             factor_nut * dt * v * face_area, Ke_dof_col,
                             Ke_val_col);

          // mass matrix
          util::add_unique(dof_indices_nut[0], factor_nut * elem_vol,
                           Ke_dof_col, Ke_val_col);

          Fe(0) += factor_nut * elem_nut_old * elem_vol;

          // chemotaxi term
          // Reinitialize shape functions on the element side
          fe_elem_face_tum->reinit(elem, side);

          // loop over quadrature point
          for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

            Gradient tum_grad;

            for (unsigned int l = 0; l < phi_face.size(); l++) {

              tum_grad.add_scaled(
                  dphi_face[l][qp],
                  tum.current_solution(dof_indices_tum_var[0][l]));
            }

            Fe(0) += factor_nut * JxW_face[qp] * deck.d_chi_c * tum_grad * qface_normals[qp];
          }
        } // elem neighbor is not null
      }   // loop over faces
    }

    // contribution from other species which are linear fields
    Number nut_cur = nut.current_solution(dof_indices_nut[0]);
    {
      for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

        // Computing solution
        Number tum_cur = 0.;
        Number hyp_cur = 0.;
        Number nec_cur = 0.;
        Number taf_cur = 0.;
        Number ecm_cur = 0.;
        Number mde_cur = 0.;

        for (unsigned int l = 0; l < phi.size(); l++) {

          tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
          hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
          nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
          taf_cur += phi[l][qp] * taf.current_solution(dof_indices_taf[l]);
          ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
          mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
        }

        // add source
        Fe(0) +=
            factor_nut * JxW[qp] * (dt * deck.d_lambda_A * (tum_cur - nec_cur) +
                       dt * deck.d_lambda_ECM_D * ecm_cur * mde_cur);

        // handle all coupling with nutrient as source term
        double compute_mat =
            JxW[qp] * dt *
            (deck.d_lambda_P * (tum_cur - hyp_cur - nec_cur) +
             deck.d_lambda_Ph * hyp_cur +
             deck.d_lambda_ECM_P * (1. - ecm_cur) *
                 util::heaviside(ecm_cur - deck.d_bar_phi_ECM_P));

        util::add_unique(dof_indices_nut[0], factor_nut * compute_mat,
                           Ke_dof_col, Ke_val_col);
      } // loop over quadrature points
    }

    // add to matrix
    DenseMatrix<Number> Ke;
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    nut.rhs->add_vector(Fe, dof_indices_nut);
  }

  //
  // Handle coupling between tissue nutrient and vessel nutrient
  //
  if (true) {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;
      double c_v_k = pointer->c_v;

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = mesh.elem_ptr(e_id);
          pres_map.dof_indices(elem, dof_indices_pres, v_pres);
          double p_t_k = pres.current_solution(dof_indices_pres[0]);

          // get 3d nutrient
          nut_map.dof_indices(elem, dof_indices_nut, v_nut);
          double c_t_k = nut.current_solution(dof_indices_nut[0]);

          std::vector<unsigned int> Ke_dof_row;
          Ke_dof_row.push_back(dof_indices_nut[0]);
          DenseMatrix<Number> Ke(1, 1);
          DenseVector<Number> Fe(1);

          // implicit for c_t in source
          Ke(0, 0) += dt * factor_nut * deck.d_coupling_method_theta *
                      pointer->L_s[i] * J_b_data.half_cyl_surf * e_w;

          // explicit for c_v in source
          Fe(0) += dt * factor_nut * pointer->L_s[i] * J_b_data.half_cyl_surf *
                   e_w * (c_v_k - (1. - deck.d_coupling_method_theta) * c_t_k);

          // term due to pressure difference
          double c_transport = 0.;
          if (p_v_k - p_t_k >= 0.)
            c_transport = c_v_k;
          else
            c_transport = c_t_k;
          Fe(0) += dt * factor_nut * (1. - deck.d_osmotic_sigma) *
                   pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                   (p_v_k - p_t_k) * c_transport;

          // update matrix
          nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);
          nut.rhs->add_vector(Fe, dof_indices_nut);
        }
      } // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }

  // finish
  nut.matrix->close();
  nut.rhs->close();
}

void netfvfe::NutAssembly::assemble_2() {

    // get tumor equation system
    EquationSystems &es = d_model_p->get_system();

    // Mesh
    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Tumor system
    auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
    std::vector<unsigned int> v_tum(2);
    v_tum[0] = tum.variable_number("tumor");
    v_tum[1] = tum.variable_number("chemical_tumor");

    const DofMap &tum_map = tum.get_dof_map();
    std::vector<unsigned int> dof_indices_tum;
    std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

    // Nutrient system
    auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
    const unsigned int v_nut = nut.variable_number("nutrient");
    const DofMap &nut_map = nut.get_dof_map();
    std::vector<unsigned int> dof_indices_nut;

    // Hypoxic system
    auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
    const unsigned int v_hyp = hyp.variable_number("hypoxic");
    const DofMap &hyp_map = hyp.get_dof_map();
    std::vector<unsigned int> dof_indices_hyp;

    // Necrotic system
    auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
    const unsigned int v_nec = nec.variable_number("necrotic");
    const DofMap &nec_map = nec.get_dof_map();
    std::vector<unsigned int> dof_indices_nec;

    // TAF system
    auto &taf = es.get_system<TransientLinearImplicitSystem>("TAF");
    const unsigned int v_taf = taf.variable_number("taf");
    const DofMap &taf_map = taf.get_dof_map();
    std::vector<unsigned int> dof_indices_taf;

    // ECM system
    auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
    const unsigned int v_ecm = ecm.variable_number("ecm");
    const DofMap &ecm_map = ecm.get_dof_map();
    std::vector<unsigned int> dof_indices_ecm;

    // MDE system
    auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
    const unsigned int v_mde = mde.variable_number("mde");
    const DofMap &mde_map = mde.get_dof_map();
    std::vector<unsigned int> dof_indices_mde;

    // Velocity system
    auto &vel = es.get_system<TransientLinearImplicitSystem>("Velocity");
    std::vector<unsigned int> v_vel(2);
    v_vel[0] = vel.variable_number("velocity_x");
    v_vel[1] = vel.variable_number("velocity_y");
    if (dim > 2)
        v_vel.push_back(vel.variable_number("velocity_z"));

    const DofMap &vel_map = vel.get_dof_map();
    std::vector<unsigned int> dof_indices_vel;
    std::vector<std::vector<dof_id_type>> dof_indices_vel_var(2);
    if (dim > 2)
        dof_indices_vel_var.resize(3);

    // Pressure system
    auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
    const unsigned int v_pres = pres.variable_number("pressure");
    const DofMap &pres_map = pres.get_dof_map();
    std::vector<unsigned int> dof_indices_pres;

    // FEM parameters
    FEType fe_tum_type = tum.variable_type(0);
    UniquePtr<FEBase> fe_tum(FEBase::build(dim, fe_tum_type));
    QGauss qrule(dim, fe_tum_type.default_quadrature_order());
    fe_tum->attach_quadrature_rule(&qrule);
    const std::vector<Real> &JxW = fe_tum->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe_tum->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe_tum->get_dphi();

    // Finite element type for tumor variable
    std::unique_ptr<FEBase> fe_elem_face_tum(FEBase::build(dim, fe_tum_type));
    QGauss qface(dim - 1, fe_tum_type.default_quadrature_order());
    fe_elem_face_tum->attach_quadrature_rule(&qface);

    // Data for surface integrals on the element boundary
    const std::vector<std::vector<Real>> &phi_face = fe_elem_face_tum->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi_face =
            fe_elem_face_tum->get_dphi();
    const std::vector<Real> &JxW_face = fe_elem_face_tum->get_JxW();
    const std::vector<Point> &qface_normals = fe_elem_face_tum->get_normals();
    const std::vector<Point> &qface_points = fe_elem_face_tum->get_xyz();

    // Model parameters
    const auto &deck = d_model_p->get_input_deck();
    const Real dt = es.parameters.get<Real>("time_step");
    const Real mesh_size = deck.d_mesh_size_vec[0];
    Real face_area = mesh_size * mesh_size;
    if (mesh.mesh_dimension() == 2)
        face_area = mesh_size;

    Real elem_vol = face_area * mesh_size;

    double factor_nut = 1.;
    if (deck.d_tissue_nut_L_s > 1.e-18)
        factor_nut = 1. / deck.d_tissue_nut_L_s;

    // store boundary condition constraints
    std::vector<unsigned int> bc_rows;
    std::vector<Number> bc_vals;

    // Looping through elements
    for (const auto &elem : mesh.active_local_element_ptr_range()) {

        tum_map.dof_indices(elem, dof_indices_tum);
        for (unsigned int var = 0; var < 2; var++)
            tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

        nut_map.dof_indices(elem, dof_indices_nut, v_nut);
        hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
        nec_map.dof_indices(elem, dof_indices_nec, v_nec);
        taf_map.dof_indices(elem, dof_indices_taf, v_taf);
        ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
        mde_map.dof_indices(elem, dof_indices_mde, v_mde);

        vel_map.dof_indices(elem, dof_indices_vel);
        for (unsigned int var = 0; var < dim; var++)
            vel_map.dof_indices(elem, dof_indices_vel_var[var], v_vel[var]);

        pres_map.dof_indices(elem, dof_indices_pres, v_pres);

        const unsigned int n_dofs = dof_indices_nut.size();

        fe_tum->reinit(elem);

        std::vector<Real> Ke_val_col;
        std::vector<unsigned int> Ke_dof_col;
        std::vector<unsigned int> Ke_dof_row;
        Ke_dof_row.push_back(dof_indices_nut[0]);
        DenseVector<Number> Fe;
        Fe.resize(n_dofs);

        // Finite volume contribution
        {
            const auto elem_center = elem->centroid();

            double elem_p_t_k = pres.current_solution(dof_indices_pres[0]);
            double elem_nut_old = nut.current_solution(dof_indices_nut[0]);

            // loop over sides of the element

            for (auto side : elem->side_index_range()) {

                if (elem->neighbor_ptr(side) != nullptr) {

                    const Elem *neighbor = elem->neighbor_ptr(side);
                    const auto neighbor_center = neighbor->centroid();
                    const auto dx = (elem_center - neighbor_center).norm();

                    // get dof id
                    std::vector<unsigned int> dof_indices_nut_neigh;
                    nut_map.dof_indices(neighbor, dof_indices_nut_neigh, v_nut);

                    std::vector<unsigned int> dof_indices_pres_neigh;
                    pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);
                    double neighbor_p_t_k =
                            pres.current_solution(dof_indices_pres_neigh[0]);

                    // get coefficient
                    const Real a = factor_nut * dt * deck.d_D_sigma * face_area / dx;
                    //          if (elem->id() > 100 and elem->id() < 120)
                    //            out << "Elem: " << elem->id() << ", neighbor: " << neighbor->id()
                    //                << ", diffusion: " << a
                    //                << ", dt: " << dt << ", D_s: " << deck.d_D_sigma << "\n";

                    // add contribution
                    // +a to (e,e) where e is id of this element
                    // -a to (e, n_e) where n_e is id of neighbor of element
                    util::add_unique(dof_indices_nut[0], a, Ke_dof_col, Ke_val_col);
                    util::add_unique(dof_indices_nut_neigh[0], -a, Ke_dof_col,
                                     Ke_val_col);

                    // advection
                    double v = -deck.d_tissue_flow_K * (neighbor_p_t_k - elem_p_t_k) /
                               (deck.d_tissue_flow_mu * mesh_size);
                    if (v >= 0.)
                        util::add_unique(dof_indices_nut[0],
                                         factor_nut * dt * v * face_area, Ke_dof_col,
                                         Ke_val_col);
                    else
                        util::add_unique(dof_indices_nut_neigh[0],
                                         factor_nut * dt * v * face_area, Ke_dof_col,
                                         Ke_val_col);

                    // mass matrix
                    util::add_unique(dof_indices_nut[0], factor_nut * elem_vol,
                                     Ke_dof_col, Ke_val_col);

                    Fe(0) += factor_nut * elem_nut_old * elem_vol;

                    // chemotaxi term
                    // Reinitialize shape functions on the element side
                    fe_elem_face_tum->reinit(elem, side);

                    // loop over quadrature point
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

                        Gradient tum_grad;

                        for (unsigned int l = 0; l < phi_face.size(); l++) {

                            tum_grad.add_scaled(
                                    dphi_face[l][qp],
                                    tum.current_solution(dof_indices_tum_var[0][l]));
                        }

                        Fe(0) += factor_nut * JxW_face[qp] * deck.d_chi_c * tum_grad * qface_normals[qp];
                    }
                } // elem neighbor is not null
            }   // loop over faces
        }

        // contribution from other species which are linear fields
        {
            for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

                // Computing solution
                Number tum_cur = 0.;
                Number hyp_cur = 0.;
                Number nec_cur = 0.;
                Number ecm_cur = 0.;
                Number mde_cur = 0.;

                for (unsigned int l = 0; l < phi.size(); l++) {

                    tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
                    hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
                    nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
                    ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
                    mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
                }

                // get projected values of species
                Number tum_proj = util::project_concentration(tum_cur);
                Number hyp_proj = util::project_concentration(hyp_cur);
                Number nec_proj = util::project_concentration(nec_cur);
                Number ecm_proj = util::project_concentration(ecm_cur);
                Number mde_proj = util::project_concentration(mde_cur);

                // add source
                Fe(0) +=
                        factor_nut * JxW[qp] * dt * (deck.d_lambda_A * (tum_proj - nec_proj) +
                                                deck.d_lambda_ECM_D * ecm_proj * mde_proj);

                // handle all coupling with nutrient as source term
                double compute_mat =
                        JxW[qp] * dt *
                        (deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) +
                         deck.d_lambda_Ph * hyp_proj +
                         deck.d_lambda_ECM_P * (1. - ecm_proj) *
                         util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P));

                util::add_unique(dof_indices_nut[0], factor_nut * compute_mat,
                                 Ke_dof_col, Ke_val_col);
            } // loop over quadrature points
        }

        // add to matrix
        DenseMatrix<Number> Ke;
        Ke.resize(1, Ke_dof_col.size());
        for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
            Ke(0, i) = Ke_val_col[i];

        nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

        // add to vector
        nut.rhs->add_vector(Fe, dof_indices_nut);
    }

    //
    // Handle coupling between tissue nutrient and vessel nutrient
    //
    if (true) {
        const auto &network = d_model_p->get_network();
        auto pointer = network.get_mesh().getHead();

        while (pointer) {

            int numberOfNeighbors = pointer->neighbors.size();

            double p_v_k = pointer->p_v;
            double c_v_k = pointer->c_v;

            for (int i = 0; i < numberOfNeighbors; i++) {

                const auto &J_b_data = pointer->J_b_points[i];

                // loop over 3d elements
                for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

                    auto e_id = J_b_data.elem_id[e];
                    auto e_w = J_b_data.elem_weight[e];

                    // get 3d pressure
                    const auto *elem = mesh.elem_ptr(e_id);
                    pres_map.dof_indices(elem, dof_indices_pres, v_pres);
                    double p_t_k = pres.current_solution(dof_indices_pres[0]);

                    // get 3d nutrient
                    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
                    double c_t_k = nut.current_solution(dof_indices_nut[0]);
                    Number c_t_k_proj = util::project_concentration(c_t_k);

                    std::vector<unsigned int> Ke_dof_row;
                    Ke_dof_row.push_back(dof_indices_nut[0]);
                    DenseMatrix<Number> Ke(1, 1);
                    DenseVector<Number> Fe(1);

                    // implicit for c_t in source
                    Ke(0, 0) += dt * factor_nut * deck.d_coupling_method_theta *
                                pointer->L_s[i] * J_b_data.half_cyl_surf * e_w;

                    // explicit for c_v in source
                    Fe(0) += dt * factor_nut * pointer->L_s[i] * J_b_data.half_cyl_surf *
                             e_w * (c_v_k - (1. - deck.d_coupling_method_theta) * c_t_k_proj);

                    // term due to pressure difference
                    double c_transport = 0.;
                    if (p_v_k - p_t_k >= 0.)
                        c_transport = c_v_k;
                    else
                        c_transport = c_t_k_proj;
                    Fe(0) += dt * factor_nut * (1. - deck.d_osmotic_sigma) *
                             pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                             (p_v_k - p_t_k) * c_transport;

                    // update matrix
                    nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);
                    nut.rhs->add_vector(Fe, dof_indices_nut);
                }
            } // loop over neighbor segments

            pointer = pointer->global_successor;
        } // loop over vertex in 1-d
    }

    // finish
    nut.matrix->close();
    nut.rhs->close();
}

void netfvfe::NutAssembly::assemble_3() {

  // get tumor equation system
  EquationSystems &es = d_model_p->get_system();

  // Mesh
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Tumor system
  auto &tum = es.get_system<TransientLinearImplicitSystem>("Tumor");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("tumor");
  v_tum[1] = tum.variable_number("chemical_tumor");

  const DofMap &tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices_tum;
  std::vector<std::vector<dof_id_type>> dof_indices_tum_var(2);

  // Nutrient system
  auto &nut = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = nut.variable_number("nutrient");
  const DofMap &nut_map = nut.get_dof_map();
  std::vector<unsigned int> dof_indices_nut;

  // Hypoxic system
  auto &hyp = es.get_system<TransientLinearImplicitSystem>("Hypoxic");
  const unsigned int v_hyp = hyp.variable_number("hypoxic");
  const DofMap &hyp_map = hyp.get_dof_map();
  std::vector<unsigned int> dof_indices_hyp;

  // Necrotic system
  auto &nec = es.get_system<TransientLinearImplicitSystem>("Necrotic");
  const unsigned int v_nec = nec.variable_number("necrotic");
  const DofMap &nec_map = nec.get_dof_map();
  std::vector<unsigned int> dof_indices_nec;

  // TAF system
  auto &taf = es.get_system<TransientLinearImplicitSystem>("TAF");
  const unsigned int v_taf = taf.variable_number("taf");
  const DofMap &taf_map = taf.get_dof_map();
  std::vector<unsigned int> dof_indices_taf;

  // ECM system
  auto &ecm = es.get_system<TransientLinearImplicitSystem>("ECM");
  const unsigned int v_ecm = ecm.variable_number("ecm");
  const DofMap &ecm_map = ecm.get_dof_map();
  std::vector<unsigned int> dof_indices_ecm;

  // MDE system
  auto &mde = es.get_system<TransientLinearImplicitSystem>("MDE");
  const unsigned int v_mde = mde.variable_number("mde");
  const DofMap &mde_map = mde.get_dof_map();
  std::vector<unsigned int> dof_indices_mde;

  // Velocity system
  auto &vel = es.get_system<TransientLinearImplicitSystem>("Velocity");
  std::vector<unsigned int> v_vel(2);
  v_vel[0] = vel.variable_number("velocity_x");
  v_vel[1] = vel.variable_number("velocity_y");
  if (dim > 2)
    v_vel.push_back(vel.variable_number("velocity_z"));

  const DofMap &vel_map = vel.get_dof_map();
  std::vector<unsigned int> dof_indices_vel;
  std::vector<std::vector<dof_id_type>> dof_indices_vel_var(2);
  if (dim > 2)
    dof_indices_vel_var.resize(3);

  // Pressure system
  auto &pres = es.get_system<TransientLinearImplicitSystem>("Pressure");
  const unsigned int v_pres = pres.variable_number("pressure");
  const DofMap &pres_map = pres.get_dof_map();
  std::vector<unsigned int> dof_indices_pres;

  // FEM parameters
  FEType fe_tum_type = tum.variable_type(0);
  UniquePtr<FEBase> fe_tum(FEBase::build(dim, fe_tum_type));
  QGauss qrule(dim, fe_tum_type.default_quadrature_order());
  fe_tum->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe_tum->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe_tum->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_tum->get_dphi();

  // Finite element type for tumor variable
  std::unique_ptr<FEBase> fe_elem_face_tum(FEBase::build(dim, fe_tum_type));
  QGauss qface(dim - 1, fe_tum_type.default_quadrature_order());
  fe_elem_face_tum->attach_quadrature_rule(&qface);

  // Data for surface integrals on the element boundary
  const std::vector<std::vector<Real>> &phi_face = fe_elem_face_tum->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi_face =
      fe_elem_face_tum->get_dphi();
  const std::vector<Real> &JxW_face = fe_elem_face_tum->get_JxW();
  const std::vector<Point> &qface_normals = fe_elem_face_tum->get_normals();
  const std::vector<Point> &qface_points = fe_elem_face_tum->get_xyz();

  // Model parameters
  const auto &deck = d_model_p->get_input_deck();
  const Real dt = es.parameters.get<Real>("time_step");
  const Real mesh_size = deck.d_mesh_size_vec[0];
  Real face_area = mesh_size * mesh_size;
  if (mesh.mesh_dimension() == 2)
    face_area = mesh_size;

  Real elem_vol = face_area * mesh_size;

  double factor_nut = 1.;
  if (deck.d_tissue_nut_L_s > 1.e-18)
    factor_nut = 1. / deck.d_tissue_nut_L_s;

  // store boundary condition constraints
  std::vector<unsigned int> bc_rows;
  std::vector<Number> bc_vals;

  // Looping through elements
  for (const auto &elem : mesh.active_local_element_ptr_range()) {

    tum_map.dof_indices(elem, dof_indices_tum);
    for (unsigned int var = 0; var < 2; var++)
      tum_map.dof_indices(elem, dof_indices_tum_var[var], v_tum[var]);

    nut_map.dof_indices(elem, dof_indices_nut, v_nut);
    hyp_map.dof_indices(elem, dof_indices_hyp, v_hyp);
    nec_map.dof_indices(elem, dof_indices_nec, v_nec);
    taf_map.dof_indices(elem, dof_indices_taf, v_taf);
    ecm_map.dof_indices(elem, dof_indices_ecm, v_ecm);
    mde_map.dof_indices(elem, dof_indices_mde, v_mde);

    vel_map.dof_indices(elem, dof_indices_vel);
    for (unsigned int var = 0; var < dim; var++)
      vel_map.dof_indices(elem, dof_indices_vel_var[var], v_vel[var]);

    pres_map.dof_indices(elem, dof_indices_pres, v_pres);

    const unsigned int n_dofs = dof_indices_nut.size();

    fe_tum->reinit(elem);

    std::vector<Real> Ke_val_col;
    std::vector<unsigned int> Ke_dof_col;
    std::vector<unsigned int> Ke_dof_row;
    Ke_dof_row.push_back(dof_indices_nut[0]);
    DenseVector<Number> Fe;
    Fe.resize(n_dofs);

    // Finite volume contribution
    {
      const auto elem_center = elem->centroid();

      double elem_p_t_k = pres.current_solution(dof_indices_pres[0]);
      double elem_nut_old = nut.current_solution(dof_indices_nut[0]);

      // loop over sides of the element

      for (auto side : elem->side_index_range()) {

        if (elem->neighbor_ptr(side) != nullptr) {

          const Elem *neighbor = elem->neighbor_ptr(side);
          const auto neighbor_center = neighbor->centroid();
          const auto dx = (elem_center - neighbor_center).norm();

          // get dof id
          std::vector<unsigned int> dof_indices_nut_neigh;
          nut_map.dof_indices(neighbor, dof_indices_nut_neigh, v_nut);

          std::vector<unsigned int> dof_indices_pres_neigh;
          pres_map.dof_indices(neighbor, dof_indices_pres_neigh, v_pres);
          double neighbor_p_t_k =
              pres.current_solution(dof_indices_pres_neigh[0]);

          // get coefficient
          const Real a = factor_nut * dt * deck.d_D_sigma * face_area / dx;
          //          if (elem->id() > 100 and elem->id() < 120)
          //            out << "Elem: " << elem->id() << ", neighbor: " << neighbor->id()
          //                << ", diffusion: " << a
          //                << ", dt: " << dt << ", D_s: " << deck.d_D_sigma << "\n";

          // add contribution
          // +a to (e,e) where e is id of this element
          // -a to (e, n_e) where n_e is id of neighbor of element
          util::add_unique(dof_indices_nut[0], a, Ke_dof_col, Ke_val_col);
          util::add_unique(dof_indices_nut_neigh[0], -a, Ke_dof_col,
                           Ke_val_col);

          // advection
          double v = -deck.d_tissue_flow_K * (neighbor_p_t_k - elem_p_t_k) /
                     (deck.d_tissue_flow_mu * mesh_size);
          if (v >= 0.)
            util::add_unique(dof_indices_nut[0],
                             factor_nut * dt * v * face_area, Ke_dof_col,
                             Ke_val_col);
          else
            util::add_unique(dof_indices_nut_neigh[0],
                             factor_nut * dt * v * face_area, Ke_dof_col,
                             Ke_val_col);

          // mass matrix
          util::add_unique(dof_indices_nut[0], factor_nut * elem_vol,
                           Ke_dof_col, Ke_val_col);

          Fe(0) += factor_nut * elem_nut_old * elem_vol;

          // chemotaxi term
          // Reinitialize shape functions on the element side
          fe_elem_face_tum->reinit(elem, side);

          // loop over quadrature point
          for (unsigned int qp = 0; qp < qface.n_points(); qp++) {

            Gradient tum_grad;

            for (unsigned int l = 0; l < phi_face.size(); l++) {

              tum_grad.add_scaled(
                  dphi_face[l][qp],
                  tum.current_solution(dof_indices_tum_var[0][l]));
            }

            Fe(0) += factor_nut * JxW_face[qp] * deck.d_chi_c * tum_grad * qface_normals[qp];
          }
        } // elem neighbor is not null
      }   // loop over faces
    }

    // contribution from other species which are linear fields
    {
      Number nut_cur = nut.current_solution(dof_indices_nut[0]);
      Number nut_proj = util::project_concentration(nut_cur);
      for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {

        // Computing solution
        Number tum_cur = 0.;
        Number hyp_cur = 0.;
        Number nec_cur = 0.;
        Number ecm_cur = 0.;
        Number mde_cur = 0.;

        for (unsigned int l = 0; l < phi.size(); l++) {

          tum_cur += phi[l][qp] * tum.current_solution(dof_indices_tum_var[0][l]);
          hyp_cur += phi[l][qp] * hyp.current_solution(dof_indices_hyp[l]);
          nec_cur += phi[l][qp] * nec.current_solution(dof_indices_nec[l]);
          ecm_cur += phi[l][qp] * ecm.current_solution(dof_indices_ecm[l]);
          mde_cur += phi[l][qp] * mde.current_solution(dof_indices_mde[l]);
        }

        // get projected values of species
        Number tum_proj = util::project_concentration(tum_cur);
        Number hyp_proj = util::project_concentration(hyp_cur);
        Number nec_proj = util::project_concentration(nec_cur);
        Number ecm_proj = util::project_concentration(ecm_cur);
        Number mde_proj = util::project_concentration(mde_cur);

        // add source
        Fe(0) +=
            factor_nut * JxW[qp] * dt * (deck.d_lambda_A * (tum_proj - nec_proj) +
                                         deck.d_lambda_ECM_D * ecm_proj * mde_proj);

        // handle all coupling with nutrient as source term
        Fe(0) +=
            factor_nut * JxW[qp] * dt *
            (-deck.d_lambda_P * (tum_proj - hyp_proj - nec_proj) * nut_proj -
             deck.d_lambda_Ph * hyp_proj * nut_proj -
             deck.d_lambda_ECM_P * (1. - ecm_proj) *
             util::heaviside(ecm_proj - deck.d_bar_phi_ECM_P) * nut_proj);
      } // loop over quadrature points
    }

    // add to matrix
    DenseMatrix<Number> Ke;
    Ke.resize(1, Ke_dof_col.size());
    for (unsigned int i = 0; i < Ke_dof_col.size(); i++)
      Ke(0, i) = Ke_val_col[i];

    nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_col);

    // add to vector
    nut.rhs->add_vector(Fe, dof_indices_nut);
  }

  //
  // Handle coupling between tissue nutrient and vessel nutrient
  //
  if (true) {
    const auto &network = d_model_p->get_network();
    auto pointer = network.get_mesh().getHead();

    while (pointer) {

      int numberOfNeighbors = pointer->neighbors.size();

      double p_v_k = pointer->p_v;
      double c_v_k = pointer->c_v;

      for (int i = 0; i < numberOfNeighbors; i++) {

        const auto &J_b_data = pointer->J_b_points[i];

        // loop over 3d elements
        for (unsigned int e = 0; e < J_b_data.elem_id.size(); e++) {

          auto e_id = J_b_data.elem_id[e];
          auto e_w = J_b_data.elem_weight[e];

          // get 3d pressure
          const auto *elem = mesh.elem_ptr(e_id);
          pres_map.dof_indices(elem, dof_indices_pres, v_pres);
          double p_t_k = pres.current_solution(dof_indices_pres[0]);

          // get 3d nutrient
          nut_map.dof_indices(elem, dof_indices_nut, v_nut);
          double c_t_k = nut.current_solution(dof_indices_nut[0]);
          Number c_t_k_proj = util::project_concentration(c_t_k);

          std::vector<unsigned int> Ke_dof_row;
          Ke_dof_row.push_back(dof_indices_nut[0]);
          DenseMatrix<Number> Ke(1, 1);
          DenseVector<Number> Fe(1);

          // implicit for c_t in source
          Ke(0, 0) += dt * factor_nut * deck.d_coupling_method_theta *
                      pointer->L_s[i] * J_b_data.half_cyl_surf * e_w;

          // explicit for c_v in source
          Fe(0) += dt * factor_nut * pointer->L_s[i] * J_b_data.half_cyl_surf *
                   e_w * (c_v_k - (1. - deck.d_coupling_method_theta) * c_t_k_proj);

          // term due to pressure difference
          double c_transport = 0.;
          if (p_v_k - p_t_k >= 0.)
            c_transport = c_v_k;
          else
            c_transport = c_t_k_proj;
          Fe(0) += dt * factor_nut * (1. - deck.d_osmotic_sigma) *
                   pointer->L_p[i] * J_b_data.half_cyl_surf * e_w *
                   (p_v_k - p_t_k) * c_transport;

          // update matrix
          nut.matrix->add_matrix(Ke, Ke_dof_row, Ke_dof_row);
          nut.rhs->add_vector(Fe, dof_indices_nut);
        }
      } // loop over neighbor segments

      pointer = pointer->global_successor;
    } // loop over vertex in 1-d
  }

  // finish
  nut.matrix->close();
  nut.rhs->close();
}

