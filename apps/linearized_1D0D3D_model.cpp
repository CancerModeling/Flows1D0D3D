////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner, Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <fmt/format.h>
#include <utility>

#include "petsc.h"

#include "macrocirculation/0d_boundary_conditions.hpp"
#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/coupled_explicit_implicit_1d_solver.hpp"
#include "macrocirculation/csv_vessel_tip_writer.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_csv_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/heart_to_breast_1d_solver.hpp"
#include "macrocirculation/heart_to_breast_3d_solver.hpp"
#include "macrocirculation/implicit_linear_flow_solver.hpp"
#include "macrocirculation/implicit_transport_solver.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/libmesh_utils.hpp"
#include "macrocirculation/linearized_flow_upwind_evaluator.hpp"
#include "macrocirculation/nonlinear_linear_coupling.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/vessel_tree_flow_integrator.hpp"

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

namespace macrocirculation {

class LinearizedHeartToBreast1DSolver {
public:
  explicit LinearizedHeartToBreast1DSolver(MPI_Comm comm)
      : d_comm(comm),
        d_degree(2),
        d_is_setup(false),
        graph{std::make_shared<GraphStorage>()},
        d_tau_flow{NAN},
        solver{nullptr} {}

  void set_path_inflow_pressures(const std::string &path) {
    if (d_is_setup)
      throw std::runtime_error("HeartToBreast1DSolver::set_path_inflow_pressures: cannot be called after setup.");

    path_inflow_pressures = path;
  }

  void set_path_geometry(const std::string &path) {
    if (d_is_setup)
      throw std::runtime_error("HeartToBreast1DSolver::set_path_geometry: cannot be called after setup.");

    path_linear_geometry = path;
  }

  void setup(size_t degree, double tau) {
    d_is_setup = true;
    setup_graphs();
    setup_solver(degree, tau);
    setup_output();
  }

  void solve_flow(double tau, double t) {
    if (std::abs(tau - d_tau_flow) > 1e-12)
      throw std::runtime_error("changing time step width not supported yet.");

    solver->solve(tau, t);
  }

  void solve_transport(double tau, double t) {
    upwind_evaluator->init(t, solver->get_solution());
    transport_solver->solve(tau, t);
  }

  void apply_slope_limiter_transport(double t) {
    transport_solver->apply_slope_limiter(t);
  }

  void write_output(double t) {
    if (solver == nullptr)
      throw std::runtime_error("solver must be initialized before output");

    auto &dof_map = solver->get_dof_map();

    auto &dof_map_transport = *transport_solver->get_dof_maps_transport().at(0);

    csv_writer->add_data("p", solver->get_solution());
    csv_writer->add_data("q", solver->get_solution());
    csv_writer->add_data("c", transport_solver->get_solution());
    csv_writer->write(t);

    interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map, solver->p_component, solver->get_solution(), points, p_vertex_values);
    interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map, solver->q_component, solver->get_solution(), points, q_vertex_values);
    interpolate_to_vertices(MPI_COMM_WORLD, *graph, dof_map_transport, 0, transport_solver->get_solution(), points, c_vertex_values);

    graph_pvd_writer->set_points(points);
    graph_pvd_writer->add_vertex_data("p", p_vertex_values);
    graph_pvd_writer->add_vertex_data("q", q_vertex_values);
    graph_pvd_writer->add_vertex_data("c", c_vertex_values);
    graph_pvd_writer->add_vertex_data("vessel_id", vessel_ids);
    graph_pvd_writer->add_vertex_data("r", vessel_radii);
    graph_pvd_writer->add_vertex_data("A", vessel_A);
    graph_pvd_writer->write(t);

    vessel_tip_writer->write(t, {solver->get_solution(), transport_solver->get_solution(), transport_solver->get_volumes()});
  }

  void set_output_folder(std::string output_dir) {
    output_folder_name = std::move(output_dir);
  }

  std::vector<VesselTipCurrentCouplingData> get_vessel_tip_pressures() {
    auto &u_flow = solver->get_solution();
    auto pressure_values = get_vessel_tip_dof_values(d_comm, *graph, *d_dof_map, u_flow);

    auto &dof_map_transport = *transport_solver->get_dof_maps_transport().back();
    auto &u_transport = transport_solver->get_solution();
    auto concentration_values = get_vessel_tip_dof_values(d_comm, *graph, dof_map_transport, u_transport);

    std::vector<VesselTipCurrentCouplingData> results;

    for (auto v_id : graph->get_vertex_ids()) {
      auto &v = *graph->get_vertex(v_id);

      if (v.is_vessel_tree_outflow()) {
        auto &e = *graph->get_edge(v.get_edge_neighbors()[0]);

        auto &R = v.get_vessel_tree_data().resistances;

        if (!e.has_embedding_data())
          throw std::runtime_error("cannot determine coupling data for an unembedded graph");

        Point p = e.is_pointing_to(v_id) ? e.get_embedding_data().points.back() : e.get_embedding_data().points.front();

        // 1e3 since we have to convert kg -> g:
        auto p_out = pressure_values[v_id] * 1e3;

        // 1e3 since we have to convert kg -> g:
        auto R2 = R.back() * 1e3;

        auto radius = calculate_edge_tree_parameters(e).radii.back();

        auto concentration = concentration_values[v_id];

        results.push_back({p, v.get_id(), p_out, concentration, R2, radius, R.size()});
      }
    }

    return results;
  }

  void update_vessel_tip_pressures(const std::map<size_t, double> &pressures_at_outlets) {
    for (auto v_id : graph->get_active_vertex_ids(mpi::rank(d_comm))) {
      auto &v = *graph->get_vertex(v_id);
      if (v.is_vessel_tree_outflow())
        // convert [Ba] to [kg / cm / s^2]
        v.update_vessel_tip_pressures(pressures_at_outlets.at(v_id) * 1e-3);
    }
  }

  void update_vessel_tip_concentrations(const std::map<size_t, double> &concentrations_at_outlets) {
    for (auto v_id : graph->get_active_vertex_ids(mpi::rank(d_comm))) {
      auto &v = *graph->get_vertex(v_id);
      if (v.is_vessel_tree_outflow())
        transport_solver->set_inflow_value(*graph, v, concentrations_at_outlets.at(v_id));
    }
  }

private:
  MPI_Comm d_comm;

  size_t d_degree;

  bool d_is_setup;

  std::shared_ptr<GraphStorage> graph;

  std::shared_ptr<DofMap> d_dof_map;

  double d_tau_flow;

  std::string path_linear_geometry{"data/1d-meshes/coarse-breast-geometry-with-extension.json"};

  std::string path_inflow_pressures{"data/1d-input-pressures/from-33-vessels-with-small-extension.json"};

  std::string output_folder_name{"output"};
  std::string filename_csv{"linearized_heart_to_breast_1d_solution"};
  std::string filename_csv_tips{"linearized_heart_to_breast_1d_solution_tips"};
  std::string filename_pvd{"linearized_heart_to_breast_1d_solution"};

  std::shared_ptr<ImplicitLinearFlowSolver> solver;

  std::shared_ptr<GraphCSVWriter> csv_writer;
  std::shared_ptr<CSVVesselTipWriter> vessel_tip_writer;
  std::shared_ptr<GraphPVDWriter> graph_pvd_writer;

  std::vector<Point> points;
  std::vector<double> p_vertex_values;
  std::vector<double> q_vertex_values;
  std::vector<double> c_vertex_values;

  std::vector<double> vessel_ids;
  std::vector<double> vessel_radii;
  std::vector<double> vessel_A;

  std::shared_ptr<LinearizedFlowUpwindEvaluator> upwind_evaluator;
  std::shared_ptr<UpwindProviderLinearizedFlow> upwind_provider;
  std::shared_ptr<ImplicitTransportSolver> transport_solver;

  void setup_graphs() {
    EmbeddedGraphReader graph_reader;

    graph_reader.append(path_linear_geometry, *graph);

    auto inflow_pressures = read_input_pressures(path_inflow_pressures);
    for (const auto &inflow_pressure : inflow_pressures) {
      auto &v_in = *graph->find_vertex_by_name(inflow_pressure.name);
      v_in.set_to_inflow_with_fixed_pressure(piecewise_linear_source_function(inflow_pressure.t, inflow_pressure.p, inflow_pressure.periodic));
    }

    set_0d_tree_boundary_conditions(graph);

    graph->finalize_bcs();

    flow_mesh_partitioner(d_comm, *graph, d_degree);
  }

  void setup_solver(size_t degree, double tau) {
    setup_solver_flow(degree, tau);
    setup_solver_transport(degree);
  }

  void setup_solver_flow(size_t degree, double tau) {
    d_degree = degree;
    d_tau_flow = tau;
    d_dof_map = std::make_shared<DofMap>(*graph);
    d_dof_map->create(d_comm, *graph, 2, d_degree, true);
    solver = std::make_shared<ImplicitLinearFlowSolver>(d_comm, graph, d_dof_map, d_degree);
    solver->setup(tau);
  }

  void setup_solver_transport(size_t degree) {
    // upwind evaluators:
    upwind_evaluator = std::make_shared<LinearizedFlowUpwindEvaluator>(MPI_COMM_WORLD, graph, d_dof_map);
    upwind_provider = std::make_shared<UpwindProviderLinearizedFlow>(graph, upwind_evaluator, solver);

    auto dof_map_transport = std::make_shared<DofMap>(*graph);
    DofMap::create_for_transport(MPI_COMM_WORLD, {graph}, {dof_map_transport}, degree);

    transport_solver = std::make_shared<ImplicitTransportSolver>(MPI_COMM_WORLD,
                                                                 std::vector<std::shared_ptr<GraphStorage>>({graph}),
                                                                 std::vector<std::shared_ptr<DofMap>>({dof_map_transport}),
                                                                 std::vector<std::shared_ptr<UpwindProvider>>({upwind_provider}),
                                                                 degree);
  }


  void setup_output() {
    if (solver == nullptr)
      throw std::runtime_error("solver must be initialized before output");

    auto dof_map_transport = transport_solver->get_dof_maps_transport().at(0);

    csv_writer = std::make_shared<GraphCSVWriter>(d_comm, output_folder_name, filename_csv, graph);
    csv_writer->add_setup_data(d_dof_map, solver->p_component, "p");
    csv_writer->add_setup_data(d_dof_map, solver->q_component, "q");
    csv_writer->add_setup_data(dof_map_transport, 0, "c");
    csv_writer->setup();

    graph_pvd_writer = std::make_shared<GraphPVDWriter>(d_comm, output_folder_name, filename_pvd);

    vessel_tip_writer = std::make_shared<CSVVesselTipWriter>(MPI_COMM_WORLD,
                                                             output_folder_name, filename_csv_tips,
                                                             graph,
                                                             std::vector<std::shared_ptr<DofMap>>({d_dof_map, dof_map_transport, transport_solver->get_dof_maps_volume().back()}),
                                                             std::vector<std::string>({"p", "c", "V"}));

    // vessels ids and radii do not change, thus we can precalculate them
    fill_with_vessel_id(d_comm, *graph, points, vessel_ids);
    fill_with_radius(d_comm, *graph, points, vessel_radii);
    fill_with_vessel_A0(d_comm, *graph, points, vessel_A);
  }
};

} // namespace macrocirculation


int main(int argc, char *argv[]) {
  // Libmesh init
  lm::LibMeshInit init(argc, argv);
  lm::Parallel::Communicator *comm = &init.comm();

  char petsc_version[1000];
  PetscGetVersion(petsc_version, 1000);
  std::cout << petsc_version << std::endl;

  cxxopts::Options options(argv[0], "Fully coupled 1D-0D-3D solver.");
  options.add_options()                                                                                                                 //
    ("tau", "time step size for the 1D model", cxxopts::value<double>()->default_value("1.5625e-05"))                                   //
    ("tau-transport", "time step size for the 1D transport", cxxopts::value<double>()->default_value("1.5625e-04"))                     //
    ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("5e-1"))                                       //
    ("tau-coup", "time step size for updating the coupling", cxxopts::value<double>()->default_value("5e-3"))                           //
    ("t-coup-start", "The time when the 3D coupling gets activated", cxxopts::value<double>()->default_value("2"))                      //
    ("t-end", "Simulation period for simulation", cxxopts::value<double>()->default_value("50."))                                       //
    ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output_full_1d0d3d_pkj/"))         //
    ("time-step", "time step size", cxxopts::value<double>()->default_value("0.01"))                                                    //
    ("mesh-size", "mesh size", cxxopts::value<double>()->default_value("0.02"))                                                         //
    ("mesh-file", "mesh filename", cxxopts::value<std::string>()->default_value("data/3d-meshes/test_full_1d0d3d_cm.e"))                //
    ("deactivate-3d-1d-coupling", "deactivates the 3d-1d coupling", cxxopts::value<bool>()->default_value("false"))                     //
    ("no-slope-limiter", "Disables the slope limiter", cxxopts::value<bool>()->default_value("false"))                                  //
    ("input-file", "input filename for parameters", cxxopts::value<std::string>()->default_value(""))                                   //
    ("input-file-pressures", "file containing the pressures at the input", cxxopts::value<std::string>()->default_value(""))            //
    ("input-file-1d-geometry", "file containing the linearized part of our geometry", cxxopts::value<std::string>()->default_value("")) //
    ("h,help", "print usage");
  options.allow_unrecognised_options(); // for petsc

  auto args = options.parse(argc, argv);

  if (args.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const bool slope_limiter_active = !args["no-slope-limiter"].as<bool>();

  // setup 1D solver
  const double t_end = args["t-end"].as<double>();
  const double t_coup_start = args["t-coup-start"].as<double>();
  const std::size_t max_iter = 160000000;

  const auto tau = args["tau"].as<double>();
  const auto tau_out = args["tau-out"].as<double>();
  const auto tau_coup = args["tau-coup"].as<double>();
  const auto tau_transport = args["tau-transport"].as<double>();
  auto out_dir = args["output-directory"].as<std::string>();

  const auto input_file_pressures = args["input-file-pressures"].as<std::string>();
  const auto input_file_geometry = args["input-file-1d-geometry"].as<std::string>();

  if (tau > tau_coup || tau > tau_transport) {
    std::cerr << "flow time step width too large" << std::endl;
    exit(0);
  }

  const auto activate_3d_1d_coupling = !args["deactivate-3d-1d-coupling"].as<bool>();

  const auto output_interval = static_cast<std::size_t>(std::round(tau_out / tau));
  const auto coupling_interval = static_cast<std::size_t>(std::round(tau_coup / tau));
  const auto transport_update_interval = static_cast<std::size_t>(std::round(tau_transport / tau));

  std::cout << output_interval << " " << coupling_interval << " " << transport_update_interval << std::endl;

  if (std::abs(coupling_interval * tau - tau_coup) > 1e-3 * tau) {
    std::cerr << "coupling time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  if (std::abs(transport_update_interval * tau - tau_transport) > 1e-3 * tau) {
    std::cerr << "transport time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  if (std::abs(output_interval * tau - tau_out) > 1e-3 * tau) {
    std::cerr << "output time step width must be a multiple of tau" << std::endl;
    exit(0);
  }

  mc::LinearizedHeartToBreast1DSolver solver_1d(MPI_COMM_WORLD);
  // customize the 1d solver from the respective input parameters
  if (!input_file_pressures.empty()) {
    std::cout << "Using input pressures from " << input_file_pressures << std::endl;
    solver_1d.set_path_inflow_pressures(input_file_pressures);
  }
  if (!input_file_geometry.empty()) {
    std::cout << "Using custom linear geometry at " << input_file_geometry << std::endl;
    solver_1d.set_path_geometry(input_file_geometry);
  }
  solver_1d.set_output_folder(out_dir);
  solver_1d.setup(degree, tau);

  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "updating transport every " << transport_update_interval << " interval" << std::endl;

  // create logger
  mc::Logger log(out_dir + "sim", comm->rank());

  // setup 3D solver
  log("setting up 3D solver\n");
  auto filename = args["input-file"].as<std::string>();
  // read input parameters
  mc::HeartToBreast3DSolverInputDeck input(filename);
  if (filename.empty()) {
    input.d_T = t_end;
    input.d_dt = args["time-step"].as<double>();
    input.d_h = args["mesh-size"].as<double>();
    input.d_mesh_file = args["mesh-file"].as<std::string>();
    input.d_out_dir = out_dir;
    input.d_debug_lvl = 1;
    input.d_perf_regularized = false;
    input.d_perf_fn_type = "const";
    input.d_perf_regularized = true;
    input.d_perf_fn_type = "linear";
    input.d_perf_neigh_size = std::make_pair(4., 10.);
  }
  log("input data \n" + input.print_str() + "\n");

  // create mesh
  log("creating mesh\n");
  lm::ReplicatedMesh mesh(*comm);
  if (!input.d_mesh_file.empty()) {
    mesh.read(input.d_mesh_file);
    //input.d_h = mc::get_min_nodal_spacing(mesh);
    input.d_h = mc::get_mesh_size_estimate_using_element_volume(mesh);
    log(fmt::format("mesh size = {}\n", input.d_h));
  } else {
    long N = long(1. / input.d_h);
    lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                          1., 0., 1., lm::HEX8);
  }

  // create equation system
  log("creating equation system\n");
  lm::EquationSystems eq_sys(mesh);
  eq_sys.parameters.set<mc::HeartToBreast3DSolverInputDeck *>("input_deck") = &input;
  eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
  auto &p_cap = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Capillary_Pressure");
  p_cap.add_variable("p_cap", lm::FIRST);
  auto &p_tis = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tissue_Pressure");
  p_tis.add_variable("p_tis", lm::FIRST);
  auto &nut_cap = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Capillary_Nutrient");
  nut_cap.add_variable("nut_cap", lm::FIRST);
  auto &nut_tis = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tissue_Nutrient");
  nut_tis.add_variable("nut_tis", lm::FIRST);
  auto &tum = eq_sys.add_system<lm::TransientLinearImplicitSystem>("Tumor");
  tum.add_variable("tum", lm::FIRST);
  tum.add_variable("mu_tum", lm::FIRST);

  // create spatial field of hydraulic conductivity
  auto &K_tis = eq_sys.add_system<lm::ExplicitSystem>("Tissue_K");
  K_tis.add_variable("k_tis", lm::CONSTANT, lm::MONOMIAL);
  auto &Dnut_tis_field = eq_sys.add_system<lm::ExplicitSystem>("Tissue_D_Nut");
  Dnut_tis_field.add_variable("Dnut_tis", lm::CONSTANT, lm::MONOMIAL);
  auto &N_bar_cap_field = eq_sys.add_system<lm::ExplicitSystem>("Avg_Capillary_Surf_Area");
  N_bar_cap_field.add_variable("n_bar_cap", lm::CONSTANT, lm::MONOMIAL);
  auto &N_bar_surf_cap_field = eq_sys.add_system<lm::ExplicitSystem>("Avg_Capillary_Cross_Section_Area");
  N_bar_surf_cap_field.add_variable("n_bar_surf_cap", lm::CONSTANT, lm::MONOMIAL);


  // create model that holds all essential variables
  log("creating model\n");
  auto solver_3d = mc::HeartToBreast3DSolver(MPI_COMM_WORLD, comm,
                                             input, mesh, eq_sys, p_cap, p_tis,
                                             nut_cap, nut_tis, tum,
                                             K_tis, Dnut_tis_field,
                                             N_bar_cap_field, N_bar_surf_cap_field,
                                             log);
  eq_sys.init();
  solver_3d.set_conductivity_fields();

  // setup the 1D pressure data in 3D solver
  log("setting 1D-3D coupling data in 3D solver\n");
  {
    auto data_1d = solver_1d.get_vessel_tip_pressures();
    solver_3d.setup_1d3d(data_1d);
  }

  // finalize 3D solver setup
  log("finalizing setup of 3D solver\n");
  solver_3d.setup();
  solver_3d.write_output();

  // NOTE to get relevant values from 3D system to solve 1D system
  // call get_vessel_tip_data_3d()
  // data_3d contains vector of coefficients a and b and also weighted avg of 3D pressure
  auto data_3d = solver_3d.get_vessel_tip_data_3d();
  log("initial 3D data at outlet tips");
  for (const auto &a : data_3d)
    log(fmt::format("avg_p = {}, avg_nut = {}\n", a.d_p_3d_w, a.d_nut_3d_w));

  double t = 0;
  double t_transport = 0;

  // time integration
  const auto begin_t = std::chrono::steady_clock::now();
  for (std::size_t it = 0; it < max_iter; it += 1) {
    // solve flow:
    solver_1d.solve_flow(tau, t);

    t += tau;

    // solve transport:
    if (it % transport_update_interval == 0) {
      solver_1d.solve_transport(tau_transport, t_transport + tau_transport);
      if (slope_limiter_active)
        solver_1d.apply_slope_limiter_transport(t_transport + transport_update_interval * tau);
      t_transport += tau_transport;
    }

    if (it % output_interval == 0) {
      if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
        std::cout << "iter = " << it << ", t = " << t << std::endl;

      solver_1d.write_output(t);
    }

    // solve 3D system:
    if ((t >= t_coup_start - 1e-12) && (it % coupling_interval == 0)) {
      std::cout << "calculates coupling " << std::endl;
      auto data_1d = solver_1d.get_vessel_tip_pressures();

      for (auto &d : data_1d) {
        // just return the values for now:
        if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
          std::cout << "vertex with id = " << d.vertex_id << ", "
                    << "coordinates = (" << d.p.x << ", " << d.p.y << ", " << d.p.z << "), "
                    << "p = " << d.pressure << ", "
                    << "c = " << d.concentration << ", "
                    << "R = " << d.R2 << ", "
                    << "r = " << d.radius << std::endl;
      }

      // Some condition to solve the 3D system
      {
        log("update 1d data in 3d solver\n");
        solver_3d.update_1d_data(data_1d);

        log("solve 3d systems\n");
        solver_3d.solve();

        log("update 3d data for 1d systems\n");
        solver_3d.update_3d_data();

        if (it % coupling_interval == 0)
          solver_3d.write_output();

        // recompute avg 3d values at outlet tips
        data_3d = solver_3d.get_vessel_tip_data_3d();
      }

      // update the boundary conditions of the 1D system:
      {
        std::map<size_t, double> new_tip_pressures;
        std::map<size_t, double> new_tip_concentrations;

        if (activate_3d_1d_coupling) {
          std::cout << "size of 3D coupling data is " << data_3d.size() << ", size of 1D coupling data is " << data_1d.size() << std::endl;
        }

        for (size_t k = 0; k < data_1d.size(); k += 1) {
          auto &d = data_1d[k];
          if (activate_3d_1d_coupling) {
            new_tip_pressures[d.vertex_id] = data_3d.at(k).d_p_3d_w;
            new_tip_concentrations[d.vertex_id] = data_3d.at(k).d_nut_3d_w; // FIXME: Gets units right
          } else {
            // constant 30 mmHg pressures
            new_tip_pressures[d.vertex_id] = 30 * 1333.3;
            new_tip_concentrations[d.vertex_id] = 0.;
          }
        }
        solver_1d.update_vessel_tip_pressures(new_tip_pressures);
        solver_1d.update_vessel_tip_concentrations(new_tip_concentrations);
      }
    }

    // break
    if (t > t_end + 1e-12)
      break;
  }

  const auto end_t = std::chrono::steady_clock::now();
  const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
  if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
}