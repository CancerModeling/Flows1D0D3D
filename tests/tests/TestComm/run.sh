#!/bin/bash
MYPWD=$(pwd)

# locate executible
execsrc="../bin/Test"

# libmesh options
libmesh_opts=""

mpiexec -np 2 "$execsrc" input.in "$libmesh_opts" 2>&1 | tee output.txt
# "$execsrc" input.in "$libmesh_opts" 2>&1 | tee output.txt
