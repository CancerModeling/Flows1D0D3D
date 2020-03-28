#!/bin/bash
PWD=$(pwd)

# locate executible
execsrc="../bin/NetFV"

# libmesh options
libmesh_opts=" -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10"

# run test
mpiexec -np 1 "$execsrc" input.in "$libmesh_opts" 2>&1 | tee output.txt
