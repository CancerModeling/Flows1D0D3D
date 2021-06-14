////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cxxopts.hpp>
#include <memory>

#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "macrocirculation/assembly_system.hpp"
#include "macrocirculation/base_model.hpp"

namespace mc = macrocirculation;

// define input, systems, model
namespace darcy3d {
// read and store input parameters
    struct InputDeck {
        InputDeck(const std::string &filename = "") : d_K(1.), d_T(1.),
                                                      d_dt(0.01), d_h(0.1) {
            if (!filename.empty() and filename != "")
                read_parameters(filename);
        }

        void read_parameters(const std::string &filename) {
            GetPot input(filename);
            d_K = input("K", 1.);
            d_T = input("T", 1.);
            d_dt = input("dt", 0.01);
            d_h = input("d_h", 0.1);
        }

        double d_K;
        double d_T;
        double d_dt;
        double d_h;
    };

// forward declare model class (specific definition requires first defining systems)
    class Model;

// define pressure system

// bc
    void bc_p(lm::EquationSystems &es) {
        std::set<lm::boundary_id_type> ids;
        ids.insert(2);
        auto &sys = es.get_system<lm::TransientLinearImplicitSystem>("P");
        std::vector<unsigned int> vars(1, sys.variable_number("p"));

        lm::ConstFunction<lm::Number> cf(100.);
        lm::DirichletBoundary bc(ids, vars, &cf);
        sys.get_dof_map().add_dirichlet_boundary(bc);
    }

// ic
    lm::Number ic_p(const lm::Point &p, const lm::Parameters &es,
                    const std::string &system_name,
                    const std::string &var_name) { return 1000.; }

    void ic(lm::EquationSystems &es, const std::string &system_name) {
        auto &sys = es.get_system<lm::TransientLinearImplicitSystem>(
                system_name);
        if (system_name == "P")
            sys.project_solution(ic_p, nullptr,
                                 es.parameters);
    }

// assembly
    class Pres : public mc::BaseAssembly {
    public:
        Pres(Model *model, lm::MeshBase &mesh,
             lm::TransientLinearImplicitSystem &sys)
                : mc::BaseAssembly("P", mesh, sys, 1,
                                   {sys.variable_number("p")}),
                  d_model_p(model) {
            sys.attach_assemble_object(
                    *this); // attach this element assembly object
            sys.attach_init_function(ic);      // add ic
            bc_p(sys.get_equation_systems());  // add bc
        }

        void assemble() override {};
        Model *d_model_p;
    };

// complete definition of model
    class Model : public mc::BaseModel {
    public:
        Model(lm::Parallel::Communicator *comm,
              InputDeck &input,
              lm::ReplicatedMesh &mesh,
              lm::EquationSystems &eq_sys,
              lm::TransientLinearImplicitSystem &pres,
              mc::Logger &log)
                : mc::BaseModel(comm, mesh, eq_sys, log, "Darcy_3D"),
                  d_input(input), d_pres(this, d_mesh, pres) {
            d_log("model created");
        };

        Pres &get_pres() { return d_pres; }

        void run() override {};

        void write_system(const unsigned int &t_step) override {};

        void solve_system() override {};

        void compute_qoi() override {};

        InputDeck &d_input;
        Pres d_pres;
    };

} // namespace darcy3d

int main(int argc, char *argv[]) {
    auto sim_begin = std::chrono::steady_clock::now();

    lm::LibMeshInit init(argc, argv);
    lm::Parallel::Communicator *comm = &init.comm();

    cxxopts::Options options(argv[0], "Darcy's flow in 3D tissue domain");
    options.add_options()("input-file", "path to the input file",
                          cxxopts::value<std::string>()->default_value(""))(
            "final-time", "final simulation time", cxxopts::value<double>()
                    ->default_value("1."))(
            "time-step", "time step size", cxxopts::value<double>()
                    ->default_value("0.01"))(
            "mesh-size", "mesh size", cxxopts::value<double>()->default_value
            ("0.1"))("hyd-cond",
                                                                "hydraulic conductivity",
                                                                cxxopts::value<double>()->default_value
                                                                        ("1."))(
            "output-directory", "directory for the output",
            cxxopts::value<std::string>()->default_value("./output/"))("h,help",
                                                                       "print usage");
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
        std::cout
                << "\nAre they part of petsc or a different auxiliary library?"
                << std::endl;
    }

    auto filename = args["input-file"].as<std::string>();
    auto out_dir = args["output-directory"].as<std::string>();

    // read input parameters
    auto input = darcy3d::InputDeck(filename);
    if (filename == "") {
        input.d_T = args["final-time"].as<double>();
        input.d_dt = args["time-step"].as<double>();
        input.d_h = args["mesh-size"].as<double>();
        input.d_K = args["hyd-cond"].as<double>();
    }

    // create logger
    mc::Logger log(out_dir + "sim.log", comm->rank());

    // create mesh
    log("creating mesh\n");
    lm::ReplicatedMesh mesh(*comm);
    long N = long(1. / input.d_h);
    lm::MeshTools::Generation::build_cube(mesh, N, N, N, 0., 1., 0.,
                                          1., 0., 1., lm::HEX8);

    // create equation system
    log("creating equation system\n");
    lm::EquationSystems eq_sys(mesh);
    eq_sys.parameters.set<darcy3d::InputDeck *>("input_deck") = &input;
    eq_sys.parameters.set<lm::Real>("time_step") = input.d_dt;
    auto &pres = eq_sys.add_system<lm::TransientLinearImplicitSystem>("P");
    pres.add_variable("p", lm::FIRST);

    // create model that holds all essential variables
    log("creating model\n");
    auto model = darcy3d::Model(comm, input, mesh, eq_sys, pres, log);

    return 0;
}