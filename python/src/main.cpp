#include <pybind11/pybind11.h>
#include "macrocirculation/heart_to_breast_1d_solver.hpp"
#include "petsc.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
namespace mc = macrocirculation;

static int len = 0;
static char * program_name = "";
static char ** argv = &program_name;

void initialize (){
  CHKERRABORT(MPI_COMM_WORLD, PetscInitialize(&len, &argv, nullptr, "flows1d0d3d library"));
};

void finalize(){
  CHKERRABORT(MPI_COMM_WORLD, PetscFinalize());
}

PYBIND11_MODULE(_core, m) {
  m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: flows1d0d3d
        .. autosummary::
           :toctree: _generate
           HeartToBreast1DSolver
    )pbdoc";

  m.def("initialize", initialize, "initializes the flow1d0d3d library");
  m.def("finalize", finalize, "finalizes the flow1d0d3d library");

  py::class_<mc::HeartToBreast1DSolver>(m, "HeartToBreast1DSolver")
    .def(py::init<>())
    .def("set_output_folder", &mc::HeartToBreast1DSolver::set_output_folder)
    .def("set_path_inflow_pressures", &mc::HeartToBreast1DSolver::set_path_inflow_pressures)
    .def("set_path_nonlinear_geometry", &mc::HeartToBreast1DSolver::set_path_nonlinear_geometry)
    .def("set_path_linear_geometry", &mc::HeartToBreast1DSolver::set_path_linear_geometry)
    .def("set_path_coupling_conditions", &mc::HeartToBreast1DSolver::set_path_coupling_conditions)
    .def("setup", &mc::HeartToBreast1DSolver::setup)
    .def("solve_flow", &mc::HeartToBreast1DSolver::solve_flow)
    .def("write_output", &mc::HeartToBreast1DSolver::write_output)
    ;

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}