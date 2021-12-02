#!/bin/bash

## Usage
## Create directory <path to TumorModels>/build/test and cd to it
## Copy this script to <path to TumorModels>/build/test
## Set the path (see below) where petsc is installed
## Run: ./build.sh 

# libraries (set paths appropriately)
petsc="0"
if [[ $petsc -eq "0" ]]; then
  echo "Specify the path where petsc is installed."
  exit
fi

# build
source="../../."
cmake -DPETSC_DIR="$petsc" \
			-DCMAKE_BUILD_TYPE=Release \
  	  -DEnable_Tests=ON \
			"$source"

make -j 4

