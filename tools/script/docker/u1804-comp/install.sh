#!/bin/bash

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
echo "HPX"
echo "<<<<<<<<<<< >>>>>>>>>>>"
HPX_VER="1.3.0"
HPX_TAR_FILE="${HPX_VER}.tar.gz"
HPX_INSTALL_PATH="/usr/local/"
HPX_BUILD_PATH=$BUILDDIR/hpx/$HPX_VER/$BUILD_TYPE/
HPX_SOURCE_DIR=$SOURCEDIR/hpx/$HPX_VER

hpx_build="1"
if [[ $hpx_build -eq "1" ]]; then
  # download library
  cd $SOURCEDIR
  if [ ! -f "$HPX_TAR_FILE" ]; then
    wget "https://github.com/STEllAR-GROUP/hpx/archive/${HPX_TAR_FILE}"
  fi

  if [ ! -d "$HPX_SOURCE_DIR" ]; then
    mkdir -p $HPX_SOURCE_DIR
    tar -zxf $HPX_TAR_FILE -C $HPX_SOURCE_DIR --strip-components=1
  fi

  # build library
  cd $BUILDDIR

  if [ ! -d "$HPX_BUILD_PATH" ]; then
    mkdir -p "$HPX_BUILD_PATH"

    cd "$HPX_BUILD_PATH"

    $CMAKE_EXE -DCMAKE_BUILD_TYPE=$BUILD_TYPE   \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
          -DHPX_WITH_THREAD_IDLE_RATES=ON \
          -DHPX_WITH_DISABLED_SIGNAL_EXCEPTION_HANDLERS=ON \
          -DHPX_WITH_MALLOC=jemalloc \
          -DHPX_WITH_EXAMPLES=OFF \
          -DHPX_WITH_NETWORKING=OFF \
          -DCMAKE_INSTALL_PREFIX=$HPX_INSTALL_PATH \
          $HPX_SOURCE_DIR
          
          # -DBOOST_ROOT=$BOOST_INSTALL_PATH \
          # -DHWLOC_ROOT=$HWLOC_INSTALL_PATH \
          # -DJEMALLOC_ROOT=$JEMALLOC_INSTALL_PATH \
  fi

  cd "$HPX_BUILD_PATH"
  make -j -l$BUILDTHREADS
  make install
fi

echo "<<<<<<<<<<< >>>>>>>>>>>"
echo "PETSC"
echo "<<<<<<<<<<< >>>>>>>>>>>"
PETSC_VER="3.12.1"
PETSC_TAR_FILE="petsc-${PETSC_VER}.tar.gz"
PETSC_INSTALL_PATH="/usr/local/"
PETSC_BUILD_PATH=$BUILDDIR/petsc/$PETSC_VER/$BUILD_TYPE/
PETSC_SOURCE_DIR=$SOURCEDIR/petsc/$PETSC_VER

petsc_build="1"
if [[ $petsc_build -eq "1" ]]; then
  # download library
  cd $SOURCEDIR
  if [ ! -f "$PETSC_TAR_FILE" ]; then
    wget "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PETSC_TAR_FILE}"
  fi

  if [ ! -d "$PETSC_SOURCE_DIR" ]; then
    mkdir -p $PETSC_SOURCE_DIR
    tar -zxf $PETSC_TAR_FILE -C $PETSC_SOURCE_DIR --strip-components=1
  fi

  # build library
  cd $BUILDDIR

  petsc_arch="arch-linux2-c-debug"
  if [ ! -d "$PETSC_BUILD_PATH" ]; then
    mkdir -p "$PETSC_BUILD_PATH"

    cd "$PETSC_BUILD_PATH"

    ${PETSC_SOURCE_DIR}/configure \
        --prefix="$PETSC_INSTALL_PATH" \
        --COPTFLAGS='-O3' \
        --CXXOPTFLAGS='-O3' \
        --FOPTFLAGS='-O3' \
        --download-fblaslapack=1 \
        --with-mumps=1 --download-mumps=1 \
        --with-metis=1 --download-metis=1 \
        --with-mumps=1 --download-mumps=1 \
        --with-blacs=1 --download-blacs=1 \
        --with-hypre=1 --download-hypre=1 \
        --with-parmetis=1 --download-parmetis=1 \
        --with-scalapack=1 --download-scalapack=1 \
        --with-superlu_dist=1 --download-superlu_dist=1 \
        --with-superlu=1 --download-superlu=1 \
        --download-hdf5=ifneeded
  fi

  cd "$PETSC_BUILD_PATH"

  make PETSC_DIR="$PETSC_BUILD_PATH" PETSC_ARCH="$petsc_arch" all -j 14
  make PETSC_DIR="$PETSC_BUILD_PATH" PETSC_ARCH="$petsc_arch" install
  make PETSC_DIR="$PETSC_INSTALL_PATH" PETSC_ARCH="" test
  make PETSC_DIR="$PETSC_INSTALL_PATH" PETSC_ARCH="" streams
fi

echo "<<<<<<<<<<< >>>>>>>>>>>"
echo "LIBMESH"
echo "<<<<<<<<<<< >>>>>>>>>>>"
LIBMESH_VER="1.5.0"
LIBMESH_TAR_FILE="v$LIBMESH_VER.tar.gz"
LIBMESH_INSTALL_PATH="/usr/local/"
LIBMESH_BUILD_PATH=$BUILDDIR/libmesh/$LIBMESH_VER/$BUILD_TYPE/
LIBMESH_SOURCE_DIR=$SOURCEDIR/libmesh/$LIBMESH_VER

libmesh_build="1"
if [[ $libmesh_build -eq "1" ]]; then
  # download library
  cd $SOURCEDIR
  if [ ! -f "$LIBMESH_TAR_FILE" ]; then
    wget "https://github.com/libMesh/libmesh/archive/refs/tags/$LIBMESH_TAR_FILE"
  fi

  if [ ! -d "$LIBMESH_SOURCE_DIR" ]; then
    mkdir -p $LIBMESH_SOURCE_DIR
    tar -zxf $LIBMESH_TAR_FILE -C $LIBMESH_SOURCE_DIR --strip-components=1
  fi

  # build library
  cd $BUILDDIR

  if [ ! -d "$LIBMESH_BUILD_PATH" ]; then
    mkdir -p "$LIBMESH_BUILD_PATH"

    cd "$LIBMESH_BUILD_PATH"

    export PETSC_DIR="$PETSC_INSTALL_PATH"
    unset PETSC_ARCH
    ${LIBMESH_SOURCE_DIR}/configure \
          --prefix="$LIBMESH_INSTALL_PATH" \
          --with-metis=PETSc
  fi
  cd "$LIBMESH_BUILD_PATH"
  make -j -l$BUILDTHREADS
  make install
  make check
fi

) |& tee "build.log"
