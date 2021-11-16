# Flows1D0D3D

[![CircleCI](https://circleci.com/gh/CancerModeling/Flows1D0D3D.svg?style=shield)](https://circleci.com/gh/CancerModeling/Flows1D0D3D) [![GitHub release](https://img.shields.io/github/release/CancerModeling/Flows1D0D3D.svg)](https://GitHub.com/CancerModeling/Flows1D0D3D/releases/) [![GitHub license](https://img.shields.io/github/license/CancerModeling/Flows1D0D3D.svg)](https://github.com/CancerModeling/Flows1D0D3D/blob/main/LICENSE) [![GitHub issues](https://img.shields.io/github/issues/CancerModeling/Flows1D0D3D.svg)](https://github.com/CancerModeling/Flows1D0D3D/issues) [![GitHub repo size](https://img.shields.io/github/repo-size/CancerModeling/Flows1D0D3D.svg)](https://GitHub.com/CancerModeling/Flows1D0D3D/)

<p align="center"> <img src="https://github.com/CancerModeling/Flows1D0D3D/blob/main/assets/logo/logo.png" width="300"> </p>


## Table of contents

  - [Introduction](#Introduction)
  - [Directory structure](#Directory-structure)
  - [Examples](#Examples)
  - [Installation](#Installation)
    * [Dependencies](#Dependencies)
    * [Building the code](#Building-the-code)
    * [Recommendations for quick build](#Recommendations-for-quick-build)
  - [Run unittests](#Run-unittests)

## Introduction
Application which (currently) simulates bloodflows with coupled 1D, 0D and 3D models.
Its future goal is to simulate and study breast cancer models.

## Directory structure
  - `apps`: contains the demo applications (executibles are created inside `<build path>/bin/macrocirculation`)
  - `assets`: contains asset files for doxygen
  - `cmake`: contains cmake modules to locate petsc and libmesh
  - `docs`: doxygen related files
  - `external`: keep external libraries
  - `data`: contains all the data needed for the simulations
    - `1d-boundary`: Contains data for the calibrated RCR-models for different meshes.
    - `1d-coupling`: Contains data coupling two geometries together. 
    - `1d-input-pressures`: Contains input pressures for other models to decoupled the simulations.
    - `1d-meshes`: Contains the 1d networks.
    - `3d-meshes`: Contains the 3d domains.
  - `tools`: contains utility scripts
    - `mesh_creation`: scripts for creating the 1D mesh 
    - `script/install`: install scripts to facilitate setup
    - `visualization`: scripts for visualizing our output data
  - `src`: the source code of the library 
  - `tests`: a few catch2 unittests for critical components
  - `python`: the python bindings
    - `flows1d0d3d`: the parts of our library written in python
    - `examples`: contains example applications in python calling our (*already installed*) library.

## Installation

### Dependencies
Core dependencies are:
  - [cmake](https://cmake.org/) (3.12 or above) 
  - [vtk](https://vtk.org/) (7.1.1)
    * recommend to install using `apt-get`
  - [petsc](https://github.com/petsc/petsc) (3.13.3 or above)
    * required to build libmesh
    * see further below on how to build petsc
  - [libmesh](https://github.com/libMesh/libmesh) (1.5.0 or above)
    * core library
    * see further below on how to build this
  - [aixlog](https://github.com/badaix/aixlog) (1.2.4)
    * included as external library in the code
    * for logging
  - [fast-cpp-csv-parser](https://github.com/ben-strasser/fast-cpp-csv-parser) (1.2.4)
    * included as external library in the code
    * for reading csv file
  - [json](https://github.com/nlohmann/json) (3.7.3)
    * directly pulled using CMake `FetchContent` module
    * for reading json file
  - [cxxopts](https://github.com/jarro2783/cxxopts) (2.2.1)
    * directly pulled using CMake `FetchContent` module
    * for reading command line options
  - [fmt](https://github.com/fmtlib/fmt) (7.1.3)
    * directly pulled using CMake `FetchContent` module
    * for easy formatting and output
  - [nanoflann](https://github.com/jlblancoc/nanoflann) (1.3.2)
    * directly pulled using CMake `FetchContent` module
    * for tree search
  - [Catch2](https://github.com/catchorg/Catch2.git) (2.13.1)
    * directly pulled using CMake `FetchContent` module
    * for ctest
  - [gmm](http://getfem.org/project/libdesc_gmm.html)
    * included as external library in the code
    * provides framework to solve linear systems associated to the 1D networks

Dependencies for running the examples:
  - [python3](https://www.python.org/)
    * required to run the test python scripts
  - [numpy](https://numpy.org/)
    * required to run the test python scripts

### Building the code
Assuming all dependencies are installed in global path (`/usr/local/` etc), we build the code using
```sh
git clone https://github.com/CancerModeling/Flows1D0D3D.git

cd Flows1D0D3D && mkdir build && cd build

cmake   -DEnable_Documentation=ON \
        -DEnable_Tests=ON \
        -DCMAKE_BUILD_TYPE=Release \
        ..

make -j 4

ctest --verbose
```

If libmesh and petsc are installed at custom paths, we will do
```sh
git clone https://github.com/CancerModeling/Flows1D0D3D.git

cd Flows1D0D3D && mkdir build && cd build

cmake   -DLIBMESH_DIR="<libmesh install path>" \
        -DPETSC_DIR="<petsc install path>" \
        -DEnable_Documentation=ON \
        -DEnable_Tests=ON \
        -DCMAKE_BUILD_TYPE=Release \
        .. 

make -j 4

ctest --verbose
```

### Recommendations for quick build
1. Install most of the dependencies using `apt-get`:
```sh
sudo apt-get update 
  
sudo apt-get install -y build-essential ubuntu-dev-tools \
  wget curl lzip \
  cmake gfortran \
  libopenmpi-dev openmpi-bin \
  libboost-all-dev libvtk7-dev \
  liblapack-dev libblas-dev \
  doxygen doxygen-latex graphviz ghostscript \
  python3-pip 

# pyvista and pandas are not required, so they can be excluded
pip3 install numpy
```

2. Build petsc and libmesh. To learn how to build petsc and libmesh on ubuntu 18.04 and 20.04, you can follow the [dockerfile](https://github.com/prashjha/dockerimages/blob/main/angio-base-bionic/Dockerfile). [Shell script](tools/script/install/install_petsc_libmesh.sh) could also be helpful.

3. Use instructions in previous section to build `Flows1D0D3D` (provide right paths to libmesh and petsc). 

### Docker
For `circle-ci` testing, we use docker images `prashjha/angio-base-bionic` and `prashjha/angio-base-focal` of ubuntu 18.04 and 20.04 with petsc and libmesh installed. The associated dockerfiles can be found [here](https://github.com/prashjha/dockerimages). 

In [Packages](https://github.com/orgs/CancerModeling/packages?repo_name=Flows1D0D3D), docker images of `Flows1D0D3D` is provided. 

### Installing the python bindings

For installation execute
```bash
LIBMESH_DIR=<REPLACE_WITH_LIBMESH_DIRECTORY> PETSC_DIR=<REPLACE_WITH_PETSC_DIRECTORY> python3 -m pip install .
```
in the directory containing the `setup.py` file. Test if the installation was successful with
```bash
python3 -c "import flows1d0d3d"
```

## Run unittests

The unittests can be run right after building the project with ctest
```
cmake <cmake-options-here> ..
make 
ctest
```

## Developers
  - [Andreas Wagner](wagneran@ma.tum.de)
  - [Tobias Koeppl](koepplto@ma.tum.de)
  - [Prashant K. Jha](pjha.sci@gmail.com)
  - [Marvin Fritz](marvin.fritz@ma.tum.de)
  - [Chengyue Wu](cw35926@utexas.edu)
