#!/bin/bash
MYPWD=$(pwd)

# locate executible
execsrc="../bin/Test"

# libmesh options
libmesh_opts=""

# create directories if not present
declare -a dirs=("t1" "t2" "t3" "t4")

# loop over directories
for c_dir in ${dirs[@]}; do

	if [[ ! -d $c_dir ]]; then
		mkdir $c_dir
	fi
done


mpiexec -np 1 "$execsrc" input.in "$libmesh_opts" 2>&1 | tee output.txt
# "$execsrc" input.in "$libmesh_opts" 2>&1 | tee output.txt
