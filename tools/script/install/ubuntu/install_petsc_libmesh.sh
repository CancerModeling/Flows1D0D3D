#!/bin/bash

## usage
## Set the paths PETSC_INSTALL_PATH and LIBMESH_INSTALL_PATH where you want to install petsc and libmesh
## Run: ./install_petsc_libmesh.sh Release

(

SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
SOURCEDIR="$SCRIPTPATH/source/"
BUILDDIR="$SCRIPTPATH/build/"
BUILDTHREADS="1"
if [ ! -d "$SOURCEDIR" ]; then
    mkdir -p "$SOURCEDIR"
fi
if [ ! -d "$BUILDDIR" ]; then
    mkdir -p "$BUILDDIR"
fi

echo "SCRIPTPATH = $SCRIPTPATH"
echo "SOURCEDIR = $SOURCEDIR"
echo "BUILDDIR = $BUILDDIR"

## process input
BUILD_TYPE="$1"
if [ "$1" = "" ]; then
    echo "Build Type Not Specified"
    exit -1
fi

if [ "$1" != "Debug" ] && [ "$1" != "Release" ] && [ "$1" != "RelWithDebInfo" ]; then
    echo "Build Type Is Not Correct"
    exit -1
fi

echo "Shell: $(which sh)"

CMAKE_EXE="cmake"

echo "<<<<<<<<<<< >>>>>>>>>>>"
echo "apt-get installs"
echo "<<<<<<<<<<< >>>>>>>>>>>"
# sudo apt-get install libboost-all-dev libvtk7-dev libblas-dev liblapack-dev

echo "<<<<<<<<<<< >>>>>>>>>>>"
echo "PETSC"
echo "<<<<<<<<<<< >>>>>>>>>>>"
PETSC_VER="3.13.3"
PETSC_TAR_FILE="petsc-lite-${PETSC_VER}.tar.gz"
PETSC_INSTALL_PATH="0" # <<--- fix install path (e.g. /usr/local/)
PETSC_SOURCE_DIR=$SOURCEDIR/petsc/$PETSC_VER

if [[ $PETSC_INSTALL_PATH -eq "0" ]]; then
  echo "Need to fix the petsc install path. You can install either at global directory /usr/local/ or local directory of your choice."
  exit
fi

petsc_build="1"
if [[ $petsc_build -eq "1" ]]; then
  # download library
  cd $SOURCEDIR
  wget "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PETSC_TAR_FILE}"
  
  mkdir -p $PETSC_SOURCE_DIR
  tar -zxf $PETSC_TAR_FILE -C $PETSC_SOURCE_DIR --strip-components=1

  # build library
  cd $PETSC_SOURCE_DIR
  ./configure --prefix=$PETSC_INSTALL_PATH \
              --with-blas-lib=libblas.a \
              --with-lapack-lib=liblapack.a \
              --download-hypre  \
              --download-scalapack \
              --download-mumps \
              --download-metis \
              --download-parmetis \
              --download-superlu \
              --download-superlu_dist \
              --with-debugging=0 \
              COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3' 

  make VERBOSE=1 PETSC_DIR=`pwd` -j $BUILDTHREADS all
  make PETSC_DIR=`pwd` install 
fi

echo "<<<<<<<<<<< >>>>>>>>>>>"
echo "LIBMESH"
echo "<<<<<<<<<<< >>>>>>>>>>>"
LIBMESH_VER="1.5.0"
LIBMESH_TAR_FILE="v$LIBMESH_VER.tar.gz"
LIBMESH_INSTALL_PATH="0" # <<--- fix install path (e.g. /usr/local/)
LIBMESH_SOURCE_DIR=$SOURCEDIR/libmesh/$LIBMESH_VER

if [[ $LIBMESH_INSTALL_PATH -eq "0" ]]; then
  echo "Need to fix the libmesh install path. You can install either at global directory /usr/local/ or local directory of your choice."
  exit
fi

libmesh_build="1"
if [[ $libmesh_build -eq "1" ]]; then
  # download library
  cd $SOURCEDIR
  wget "https://github.com/libMesh/libmesh/archive/refs/tags/$LIBMESH_TAR_FILE"

  mkdir -p $LIBMESH_SOURCE_DIR
  tar -zxf $LIBMESH_TAR_FILE -C $LIBMESH_SOURCE_DIR --strip-components=1

  # build library
  cd "$LIBMESH_SOURCE_DIR"
  export PETSC_DIR="$PETSC_INSTALL_PATH"
  ./configure --with-methods="opt" \
              --prefix="$LIBMESH_INSTALL_PATH" \
              --with-metis=PETSc
  
  make -j $BUILDTHREADS
  make install
  make check
fi

) |& tee "build.log"
