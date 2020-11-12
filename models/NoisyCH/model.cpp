////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "model.hpp"
#include "rw/vtk_io.hpp"

namespace noisych {

/*!
 * @brief set initial condition
 *
 * @param es Equation system
 * @param system_name Name of system
 */
void initial_condition(EquationSystems &es, const std::string &system_name) {

  auto &sys = es.get_system<TransientLinearImplicitSystem>(system_name);

  GetPot input();

  if (system_name == "CahnHilliard") {
    auto initial_solution_type = es.parameters.get<std::string>("initial_solution");
    if (initial_solution_type == "random") {
      sys.project_solution(noisych::initial_condition_cahnhilliard_random, nullptr, es.parameters);
    } else if (initial_solution_type == "circle") {
      sys.project_solution(noisych::initial_condition_cahnhilliard_circle, nullptr, es.parameters);
    } else {
      throw std::runtime_error("unknown initial solution type " + initial_solution_type);
    }
  } else {
    return;
  }
}

// Model setup and run
void model_setup_run(const std::string &filename, Parallel::Communicator *comm) {

  auto sim_begin = steady_clock::now();

  GetPot input(filename);

  const std::size_t nelx = 1u << static_cast<unsigned int>(input("size_x", 6));
  const std::size_t nely = 1u << static_cast<unsigned int>(input("size_y", 6));
  Mesh mesh(*comm);
  MeshTools::Generation::build_square(mesh, nelx, nely, 0., 1, 0., 1, QUAD4);

  EquationSystems sys(mesh);

  sys.parameters.set<std::string>("config_filename") = filename;
  sys.parameters.set<std::string>("initial_solution") = input("initial_solution", "undefined");

  // Add systems, variables and assemble
  auto &ch = sys.add_system<TransientLinearImplicitSystem>("CahnHilliard");

  ch.add_variable("concentration", FIRST);
  ch.add_variable("potential", FIRST);
  ch.attach_init_function(initial_condition);

  auto model = Model(filename, comm, mesh, sys, ch);

  // run model
  model.run();
}

// Model class
Model::Model(const std::string &filename,
             Parallel::Communicator *comm,
             Mesh &mesh,
             EquationSystems &sys,
             TransientLinearImplicitSystem &ch)
    : d_sys(sys),
      d_mesh(mesh),
      d_ch("CahnHilliard", mesh, ch, CahnHilliardConfig::from_parameter_file(filename)),
      d_steps(100),
      d_step(0),
      d_time(0),
      d_dt(1e-2) {
  sys.parameters.set<unsigned int>("rank") = comm->rank();
  ch.attach_assemble_object(d_ch);
  sys.init();

  GetPot input(filename);
  d_dt = input("dt", 1e-2);
  const auto final_time = input("final_time", 100*d_dt);
  d_steps = std::ceil(final_time / d_dt);
}

void Model::run() {
  write_system(0);

  do {
    // Prepare time step
    d_step++;
    d_time += d_dt;

    std::cout << "Time step: " << std::to_string(d_step) << ", time: " << std::to_string(d_time) << std::endl;

    // Prepare solution for the current timestep
    *(d_ch.d_sys.old_local_solution) = *(d_ch.d_sys.current_local_solution);

    // stochastic contributions from the cylindrical Wiener process
    d_ch.calculate_new_stochastic_coefficients(d_dt);

    // solve cahn-hilliard system, where the old solution is our initial guess
    (*d_ch.d_sys.current_local_solution) = *d_ch.d_sys.old_local_solution;
    d_ch.solve();

    // write solution to disk
    write_system(d_step);

  } while (d_step < d_steps);
}

void Model::write_system(const unsigned int &t_step) {
  rw::VTKIO(d_mesh).write_equation_systems("output/cahnhilliard_" + std::to_string(t_step) + ".pvtu", d_sys);
}

} // namespace noisych
