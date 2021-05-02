#!/bin/bash

## Usage
## Create directory <path to TumorModels>/build/test and cd to it
## Copy this script to <path to TumorModels>/build/test
## Set the paths (see below) where libmesh and petsc are installed 
## Run: ./build.sh 

MY_PWD=$pwd
local="${HOME}/Softwares/local_libmesh"

# libraries (set paths appropriately)
libmesh="0"
petsc="0"
if [[ $libmesh -eq "0" -or $petsc -eq "0" ]]; then
  echo "Specify the paths where libmesh and petsc are installed."
  exit
fi

# build
target_build=$(pwd)
source="../../."
cmake -DLIBMESH_DIR="$libmesh" \
			-DPETSC_LIB="$petsc/lib" \
			-DCMAKE_BUILD_TYPE=Release \
			-DCMAKE_C_COMPILER=/usr/bin/gcc \
			-DCMAKE_CXX_COMPILER=/usr/bin/g++
			"$source"

make -j 4

