#!/bin/bash

## Usage
## Create directory <path to TumorModels>/build/test and cd to it
## Copy this script to <path to TumorModels>/build/test
## Set the paths (see below) where libmesh and petsc are installed 
## Run: ./build.sh 

# libraries (set paths appropriately)
libmesh="0"
petsc="0"
if [[ $libmesh -eq "0" -or $petsc -eq "0" ]]; then
  echo "Specify the paths where libmesh and petsc are installed."
  exit
fi

# build
source="../../."
cmake -DLIBMESH_DIR="$libmesh" \
			-DPETSC_DIR="$petsc" \
			-DCMAKE_BUILD_TYPE=Release \
  	  -DEnable_Tests=ON \
			"$source"

make -j 4

