////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "errornorm.hpp"

#include "fe_type_network.hpp"
#include "dof_map_network.hpp"
#include "gmm.h"

namespace macrocirculation {

template < std::size_t degree >
double errornorm(const GraphStorage & graph,
                 const DofMapNetwork & map,
                 std::size_t component,
                 const std::vector<double>& u_h,
                 EvaluationFunction u)
{
  const auto qf = create_gauss4();
  FETypeNetwork<degree> fe(qf);
  QuadraturePointMapper qpm(qf);

  std::vector<std::size_t> dof_indices(degree+1);
  std::vector<double> u_h_local(degree+1);
  std::vector<double> u_h_qp(fe.get_phi()[0].size());
  std::vector<double> u_qp(fe.get_phi()[0].size());
  std::vector<double> diff_qp(fe.get_phi()[0].size());

  double sum = 0;

  for (auto e_id : graph.get_edge_ids()) {
    auto edge = graph.get_edge(e_id);

    fe.reinit(*edge);
    qpm.reinit(*edge);

    map.dof_indices(*edge, dof_indices, component);
    extract_dof(dof_indices, u_h, u_h_local);
    fe.evaluate_dof_at_quadrature_points(u_h_local, u_h_qp);
    u(qpm.get_quadrature_points(), u_qp);
    gmm::add(u_h_qp, gmm::scaled(u_qp, -1), diff_qp);

    for (std::size_t qp = 0; qp<qf.size(); qp+=1)
      sum += std::pow(diff_qp[qp], 2) * fe.get_JxW()[qp];
  }

  return std::sqrt(sum);
}


template double errornorm<0>(const GraphStorage &, const DofMapNetwork &, std::size_t , const std::vector<double> &, EvaluationFunction );
template double errornorm<1>(const GraphStorage &, const DofMapNetwork &, std::size_t , const std::vector<double> &, EvaluationFunction );
template double errornorm<2>(const GraphStorage &, const DofMapNetwork &, std::size_t , const std::vector<double> &, EvaluationFunction );
template double errornorm<3>(const GraphStorage &, const DofMapNetwork &, std::size_t , const std::vector<double> &, EvaluationFunction );

}
