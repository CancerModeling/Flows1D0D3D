#!/bin/bash
MY_PWD=$(pwd)

## executible
EXEC_DIR="../../../../bin/"

## libmesh options
LIBMESH_OPTS=" -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10"

(
echo "Setting up simulation ... "
python3 -B setup.py
echo "Running simulation ... "
f_inp="input$f_suf.yaml"
"${EXEC_DIR}/NetFVFE" input.in $LIBMESH_OPTS
) |& tee output.log

# check if we have produced 'output_10.vtu' file
if [[ -f "output_10.vtu" ]]; then
    exit 0
else
    exit 1
fi
