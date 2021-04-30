#////////////////////////////////////////////////////////////////////////////////
#//  Copyright (c) 2020-2021 Prashant K. Jha, Tobias Koeppl, Andreas Wagner
#//
#//  Distributed under the Boost Software License, Version 1.0. (See accompanying
#//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#////////////////////////////////////////////////////////////////////////////////

FROM ubuntu:bionic

# Set metadata
LABEL maintainer="Prashant K. Jha <pjha.sci@gmail.com>"
LABEL description="Ubuntu 18.04 with HPX, petsc, libmesh"

# Set different working directory
WORKDIR /lib_builds/

# Set noninteractive mode (otherwise the build stucks at step where ubuntu asks for geographical location)
ENV DEBIAN_FRONTEND=nonintercative

# Install essentials
RUN apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y build-essential ubuntu-dev-tools \
  git wget curl lzip tzdata ssh gcc g++ gdb clang rsync tar python \
  cmake autoconf libtool pkg-config \
  rpm gcovr ruby-coveralls \
  liblapack-dev libblas-dev libopenmpi-dev \
  doxygen doxygen-latex graphviz ghostscript \
  gfortran libmpfr-dev libgmp-dev \
  libhwloc-dev libjemalloc-dev libboost-all-dev libyaml-cpp-dev \
  libvtk7-dev gmsh libflann-dev python3-pip && \
  apt-get autoremove -y && \
  apt-get autoclean -y

RUN pip3 install numpy pyvista pandas

# Copy the install script
COPY install.sh install.sh

# Run install script
RUN /bin/bash install.sh Release

# Copy the build logs
RUN cp build.log /build.log

# Change to root
WORKDIR /

# Clean the build directories
RUN rm -rf /lib_builds



