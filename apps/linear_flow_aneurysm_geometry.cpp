////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cxxopts.hpp>
#include <memory>
#include <petsc.h>
#include <utility>

#include "macrocirculation/0d_boundary_conditions.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/petsc/petsc_ksp.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/quantities_of_interest.hpp"

namespace mc = macrocirculation;

double inflow_q(const double t) {
  return 1.1400127037858874*sin(2*M_PI*t) + 0.10000157404719523*sin(4*M_PI*t) - 0.11286043261149026*sin(6*M_PI*t) - 0.10272017307948093*sin(8*M_PI*t) - 0.0054169100548564046*sin(10*M_PI*t) - 0.02557859579598561*sin(12*M_PI*t) - 0.078047431190483588*sin(14*M_PI*t) - 0.022033974667631729*sin(16*M_PI*t) - 0.034738483688161952*sin(18*M_PI*t) - 0.045615688241977252*sin(20*M_PI*t) + 0.018999552643855375*sin(22*M_PI*t) - 0.0084982548911405227*sin(24*M_PI*t) + 0.0066704505306768112*sin(26*M_PI*t) + 0.0205896141914854*sin(28*M_PI*t) + 0.00056547376330811749*sin(30*M_PI*t) + 0.0079161507767753492*sin(32*M_PI*t) + 0.011957921955760809*sin(34*M_PI*t) + 0.0019131544051249399*sin(36*M_PI*t) + 0.0048081271514646045*sin(38*M_PI*t) - 0.48681323281506877*cos(2*M_PI*t) - 0.53309141008559258*cos(4*M_PI*t) - 0.20633570858479031*cos(6*M_PI*t) - 0.036670810452806915*cos(8*M_PI*t) - 0.014049767327029486*cos(10*M_PI*t) - 0.066663082737686313*cos(12*M_PI*t) - 0.015523473625300599*cos(14*M_PI*t) + 0.023855278352215264*cos(16*M_PI*t) - 0.016520131034039841*cos(18*M_PI*t) + 0.04615603374764362*cos(20*M_PI*t) + 0.031258225120331377*cos(22*M_PI*t) + 0.0052060140772440802*cos(24*M_PI*t) + 0.030280763300256783*cos(26*M_PI*t) + 0.00202669624626461*cos(28*M_PI*t) - 0.00026986726295702868*cos(30*M_PI*t) + 0.0091328960066964764*cos(32*M_PI*t) - 0.0025995007712682956*cos(34*M_PI*t) - 0.003740301110564225*cos(36*M_PI*t) + 0.0020915964477192118*cos(38*M_PI*t) + 3.0675015000000001;
}


int main(int argc, char *argv[]) {
  // initialize petsc
  CHKERRQ(PetscInitialize(&argc, &argv, nullptr, "solves linear flow problem"));
  Eigen::setNbThreads(1);
  {
    const std::size_t degree = 1;

    cxxopts::Options options(argv[0], "Linearized solver for breast geometry");
    options.add_options()                                                                                                                   //
      ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output/"))                           //
      // ("mesh", "filepath to the given mesh", cxxopts::value<std::string>()->default_value("./data/1d-meshes/coarse-breast1-geometry.json")) //
      ("mesh", "filepath to the given mesh", cxxopts::value<std::string>()->default_value("./data/1d-meshes/Graph0.json")) //
      ("tau", "time step size", cxxopts::value<double>()->default_value("1e-2"))                                                            //
      ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-1"))                                         //
      ("t-end", "Endtime for simulation", cxxopts::value<double>()->default_value("10"))                                                    //
      ("h,help", "print usage");
    options.allow_unrecognised_options(); // for petsc
    auto args = options.parse(argc, argv);
    if (args.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    if (!args.unmatched().empty()) {
      std::cout << "The following arguments were unmatched: " << std::endl;
      for (auto &it : args.unmatched())
        std::cout << " " << it;
      std::cout << "\nAre they part of petsc or a different auxillary library?" << std::endl;
    }

    const double t_end = args["t-end"].as<double>();

    const auto tau = args["tau"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();

    std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);

    // vessel parameters

    // create the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();

    mc::EmbeddedGraphReader graph_reader;
    graph_reader.append(args["mesh"].as<std::string>(), *graph);


    //const std::string path_inflow_pressures = "data/1d-input-pressures/from-33-vessels-with-small-extension.json";
    //auto inflow_pressure = mc::read_input_pressures(path_inflow_pressures);
    //auto inflow_function = mc::piecewise_linear_source_function(inflow_pressure[0].t, inflow_pressure[0].p, inflow_pressure[0].periodic);
    //graph->find_vertex_by_name("Inflow")->set_to_inflow_with_fixed_pressure(inflow_function);
    graph->find_vertex_by_name("Inflow")->set_to_inflow_with_fixed_flow(inflow_q);
    // mc::set_0d_tree_boundary_conditions(graph, "Outflow");
    // graph_reader.set_boundary_data("./data/meshes/boundary-combined-geometry-linear-part.json", *graph);

    {
      double A_total = 0;
      for (auto v : graph->find_vertices_by_name_prefix("Outflow"))
      {
        // MCA
        if (v->get_id() > 24)
          continue;
        A_total += graph->get_edge( v->get_edge_neighbors()[0] )->get_physical_data().A0;
      }
      const double C_tot = 0.0116;
      const double R_tot = 59.7;
     /*
      const double C_tot = 0.0082;
      const double R_tot = 84.8;
     */
      for (auto v : graph->find_vertices_by_name_prefix("Outflow"))
      {
        if (v->get_id() > 24)
          continue;
        const double A_i = graph->get_edge( v->get_edge_neighbors()[0] )->get_physical_data().A0;
        const double ratio = A_i / A_total;
        v->set_to_windkessel_outflow(R_tot/ratio, ratio*C_tot);
      }
    }

    {
      double A_total = 0;
      for (auto v : graph->find_vertices_by_name_prefix("Outflow"))
      {
        if (v->get_id() <= 24)
          continue;
        A_total += graph->get_edge( v->get_edge_neighbors()[0] )->get_physical_data().A0;
      }

      const double C_tot = 0.0082;
      const double R_tot = 84.8;
     /*
      const double C_tot = 0.0116;
      const double R_tot = 59.7;
     */
      for (auto v : graph->find_vertices_by_name_prefix("Outflow"))
      {
        if (v->get_id() <= 24)
          continue;
        const double A_i = graph->get_edge( v->get_edge_neighbors()[0] )->get_physical_data().A0;
        const double ratio = A_i / A_total;
        v->set_to_windkessel_outflow(R_tot/ratio, ratio*C_tot);
      }
    }

    graph->finalize_bcs();

    // mc::naive_mesh_partitioner(*graph, PETSC_COMM_WORLD);
    mc::flow_mesh_partitioner(PETSC_COMM_WORLD, *graph, degree);

    auto dof_map = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map->create(PETSC_COMM_WORLD, *graph, 2, degree, true);

    if (mc::mpi::rank(PETSC_COMM_WORLD) == 0)
      std::cout << "num_1d_dof = " << dof_map->num_dof() << std::endl;
    std::cout << mc::mpi::rank(PETSC_COMM_WORLD) << " owns " << dof_map->num_owned_dofs() << " dof" << std::endl;

    mc::ImplicitLinearFlowSolver solver(PETSC_COMM_WORLD, graph, dof_map, degree);
    solver.setup(tau);
    //solver.use_named_solver("ilf_");

    std::vector<mc::Point> points;
    std::vector<double> vessel_A0;
    std::vector<double> vessel_ids;
    std::vector<double> vertex_ids;
    mc::fill_with_vessel_A0(MPI_COMM_WORLD, *graph, points, vessel_A0);
    mc::fill_with_vessel_id(MPI_COMM_WORLD, *graph, points, vessel_ids);
    mc::fill_with_vertex_id(MPI_COMM_WORLD, *graph, points, vertex_ids);

    mc::GraphPVDWriter writer(PETSC_COMM_WORLD, args["output-directory"].as<std::string>(), "breast_geometry_linearized");

    double t = 0;
    const auto t_max_idx = static_cast<size_t>(std::ceil(t_end / tau));
    const auto begin_t = std::chrono::steady_clock::now();
    for (size_t t_idx = 0; t_idx < t_max_idx; t_idx += 1) {
      t += tau;
      solver.solve(tau, t);

      if (t_idx % output_interval == 0) {
        std::cout << "iter = " << t_idx << ", t = " << t << std::endl;
        std::cout << "solver iterations = " << solver.get_solver().get_iterations() << std::endl;

        std::vector<mc::Point> points;
        std::vector<double> p_vertex_values;
        std::vector<double> q_vertex_values;
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.p_component, solver.get_solution(), points, p_vertex_values);
        interpolate_to_vertices(PETSC_COMM_WORLD, *graph, *dof_map, solver.q_component, solver.get_solution(), points, q_vertex_values);

        std::vector<double> r_vertex_values;
        calculate_linearized_r(PETSC_COMM_WORLD, *graph, p_vertex_values, points, r_vertex_values);

        writer.set_points(points);
        writer.add_vertex_data("p", p_vertex_values);
        writer.add_vertex_data("q", q_vertex_values);
        writer.add_vertex_data("a0", vessel_A0);
        writer.add_vertex_data("r", r_vertex_values);
        writer.add_vertex_data("vessel_ids", vessel_ids);
        writer.add_vertex_data("vertex_ids", vertex_ids);
        writer.write(t);

        int num = 0;
        for (auto &v_id : graph->get_active_vertex_ids(mc::mpi::rank(PETSC_COMM_WORLD))) {
          auto &vertex = *graph->get_vertex(v_id);
          if (vertex.is_vessel_tree_outflow()) {
            const auto &vertex_dof_map = dof_map->get_local_dof_map(vertex);
            const auto &vertex_dofs = vertex_dof_map.dof_indices();
            const auto &u = solver.get_solution();
            std::cout << num << std::endl;
            num += 1;

            std::cout << "rank = " << mc::mpi::rank(PETSC_COMM_WORLD) << std::endl;
            std::cout << "vertex = " << vertex.get_name() << "\n";
            std::cout << " p = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << u.get(vertex_dofs[k]) << ", ";
            std::cout << std::endl;
            std::cout << " R = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << vertex.get_vessel_tree_data().resistances[k] << ", ";
            std::cout << std::endl;
            std::cout << " C = ";
            for (size_t k = 0; k < vertex_dofs.size(); k += 1)
              std::cout << vertex.get_vessel_tree_data().capacitances[k] << ", ";
            std::cout << std::endl;
          } else {
            // std::cout << "vertex = " << vertex.get_name() << " (no tree outflow)\n";
          }
        }

        std::cout << "iter = " << t_idx << ", t = " << t << std::endl;
      }
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
  }

  CHKERRQ(PetscFinalize());
}