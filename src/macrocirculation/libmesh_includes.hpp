////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Prashant K. Jha.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef LIBMESH_INCLUDES_H
#define LIBMESH_INCLUDES_H

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/point.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/transient_system.h"
#include "libmesh/analytic_function.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/error_vector.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/getpot.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/perf_log.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system_norm.h"
#include "libmesh/utility.h"
#include "libmesh/zero_function.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/communicator.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/enum_xdr_mode.h"

namespace lm = libMesh;

#endif // LIBMESH_INCLUDES_H
