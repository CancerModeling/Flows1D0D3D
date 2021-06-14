#!/bin/bash
MYPWD=$(pwd)

# locate executible
execsrc="/home/prashant/work/tum_works/CLRemote_Mac_TumorModels/build/rem-dbg/bin/"

# libmesh options
libmesh_opts=" -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10"

if [[ $# -lt 3 ]]; then
	echo "Require three arguments"
	echo "Arg1: model name. E.g. NetFV or NetFVFE"
	echo "Arg2: pp tag. E.g. dbg_1, t1, etc"
	echo "Arg3: 0 or 1"
	exit
fi

model_name="$1"
pp_tag="$2"
script_flag="$3"

if [[ $model_name != "NetFV" ]] && [[ $model_name != "NetFVFE" ]] && [[ $model_name != "NetFCFVFE" ]] && [[ $model_name != "NetFC" ]] && [[ $model_name != "TwoSpecies" ]]; then
	echo "Check model name or modify the script to include it"
	exit
fi


if [[ $script_flag == "0" ]]; then

	echo "script_flag = $script_flag"

	# go to output directory
	if [[ ! -d "$pp_tag" ]]; then
		echo "mkdir dir = $pp_tag"
		mkdir "$pp_tag" 
	else
		echo "re-make dir = $pp_tag"
		rm -r "$pp_tag"
		mkdir "$pp_tag"
	fi
elif [[ $script_flag == "1" ]]; then
	
	cd "$pp_tag"
	mpiexec -n 1 "$execsrc$model_name" "input.in" "$libmesh_opts"	
else
	echo "check first argument passed to run.sh script"
fi
