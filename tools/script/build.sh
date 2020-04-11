#!/bin/bash
MY_PWD=$pwd
local="${HOME}/Softwares/local_libmesh"

# libraries
libmesh="$local/libmesh/1.5.0-opt"
petsc="$local/petsc/3.12.1-opt"

# target directory where code will be built
target_build=$pwd

# source of code (TGMOdels directory)
source="../../."

# cmake (can use different versions of cmake)
cmake_c="cmake"

"$cmake_c" -DLIBMESH_DIR="$libmesh" \
					 -DPETSC_DIR="$petsc" \
					 -DVTK_DIR="/usr/lib/cmake/vtk-7.1" \
					 -DLIBTG_BUILD_FLAG="Build_Prashant" \
					 -DCMAKE_INSTALL_PREFIX="$target_build" \
			  	 -DCMAKE_BUILD_TYPE=Debug \
			  	 "$source"

# compile
make -j 8

