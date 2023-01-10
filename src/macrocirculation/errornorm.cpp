////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "errornorm.hpp"

#include "communication/mpi.hpp"
#include "dof_map.hpp"
#include "fe_type.hpp"
#include "gmm_legacy_facade.hpp"

namespace macrocirculation {

double errornorm(const MPI_Comm comm,
                 const GraphStorage &graph,
                 const DofMap &map,
                 std::size_t component,
                 const std::vector<double> &u_h,
                 const EvaluationFunction &u) {
  double local_sum = 0;

  const auto qf = create_gauss4();

  for (auto e_id : graph.get_active_edge_ids(mpi::rank(comm))) {
    auto edge = graph.get_edge(e_id);

    auto &local_dof_map = map.get_local_dof_map(*edge);

    const auto degree = local_dof_map.num_basis_functions() - 1;

    FETypeNetwork fe(qf, degree);
    QuadraturePointMapper qpm(qf);

    std::vector<std::size_t> dof_indices(local_dof_map.num_basis_functions());
    std::vector<double> u_h_local(local_dof_map.num_basis_functions());
    std::vector<double> u_h_qp(qf.size());
    std::vector<double> u_qp(qf.size());
    std::vector<double> diff_qp(qf.size());

    const auto &param = edge->get_physical_data();
    const double h = param.length / local_dof_map.num_micro_edges();

    fe.reinit(h);

    for (std::size_t micro_edge_id = 0; micro_edge_id < local_dof_map.num_micro_edges(); micro_edge_id += 1) {
      qpm.reinit(micro_edge_id * h, static_cast<double>(micro_edge_id + 1) * h);

      local_dof_map.dof_indices(micro_edge_id, component, dof_indices);
      extract_dof(dof_indices, u_h, u_h_local);
      fe.evaluate_dof_at_quadrature_points(u_h_local, u_h_qp);
      u(qpm.get_quadrature_points(), u_qp);
      gmm::add(u_h_qp, gmm::scaled(u_qp, -1), diff_qp);

      for (std::size_t qp = 0; qp < qf.size(); qp += 1)
        local_sum += std::pow(diff_qp[qp], 2) * fe.get_JxW()[qp];
    }
  }

  double global_sum = 0;

  CHECK_MPI_SUCCESS(MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm));

  return std::sqrt(global_sum);
}

} // namespace macrocirculation
