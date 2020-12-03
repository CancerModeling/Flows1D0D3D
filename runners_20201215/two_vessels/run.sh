#!/bin/bash
MY_PWD=$(pwd)

## Fixed variable which should be changed less often
EXEC_DIR="/home/wagneran/tumor-model-stack/TumorModels/build_new_model_dev_tobias/bin"

## libmesh options
# LIBMESH_OPTS=" -log_view -info myinfo.txt -ksp_type minres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10"
# LIBMESH_OPTS=" -log_view -ksp_view -ksp_type gmres -ksp_monitor -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 0"
#LIBMESH_OPTS="--solver-system-names -log_view -ksp_view -hypoxic_ksp_view -hypoxic_ksp_type chebyshev -prolific_ksp_type gmres -Prolific_ksp_type gmres -ksp_type minres"
#LIBMESH_OPTS="--solver-system-names -ksp_type gmres -pc_type bjacobi -prolific_ksp_view -prolific_ksp_type gmres"
#LIBMESH_OPTS="--solver-system-names -ksp_type gmres -pc_type bjacobi -prolific_ksp_type preonly -pc_type lu"

#LIBMESH_OPTS2="--solver-system-names "
#LIBMESH_OPTS2+="--ksp_monitor_true_residual "
#LIBMESH_OPTS2+="-prolific_ksp_monitor_true_residual -hypoxic_ksp_monitor_true_residual -necrotic_ksp_monitor_true_residual -taf_ksp_monitor_true_residual -mde_ksp_monitor_true_residual "
#LIBMESH_OPTS2+="-prolific_ksp_type gmres -prolific_pc_type fieldsplit -prolific_pc_fieldsplit_type schur -prolific_pc_fieldsplit_schur_fact_type lower "
#LIBMESH_OPTS2+="-prolific_ksp_typepreonly -prolific_pc_type lu -prolific_pc_factor_mat_solver_type mumps "
#LIBMESH_OPTS2+="-hypoxic_ksp_type gmres -hypoxic_pc_type fieldsplit -hypoxic_pc_fieldsplit_type schur -hypoxic_pc_fieldsplit_schur_fact_type lower "
#LIBMESH_OPTS2+="-prolific_ksp_type gmres -prolific_pc_type bjacobi -prolific_sub_pc_type ilu -prolific_sub_pc_factor_levels 0 "
#LIBMESH_OPTS2+="-hypoxic_ksp_type gmres -hypoxic_pc_type bjacobi -hypoxic_sub_pc_type ilu -hypoxic_sub_pc_factor_levels 0 "
#LIBMESH_OPTS2+="-nutrient_ksp_type gmres -nutrient_pc_type bjacobi -nutrient_sub_pc_type ilu -nutrient_sub_pc_factor_levels 0 "
#LIBMESH_OPTS2+="-necrotic_ksp_type gmres -necrotic_pc_type bjacobi -necrotic_sub_pc_type ilu -necrotic_sub_pc_factor_levels 0 "
#LIBMESH_OPTS2+="-taf_ksp_type gmres -taf_pc_type bjacobi "
#LIBMESH_OPTS2+="-ecm_ksp_type gmres -ecm_pc_type bjacobi "
#LIBMESH_OPTS2+="-mde_ksp_type gmres -mde_pc_type bjacobi "

#LIBMESH_OPTS="-ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 0 "
#LIBMESH_OPTS=$LIBMESH_OPTS2

# LIBMESH_OPTS+="-prefix_ksp_type gmres -prefix_pc_type bjacobi -prefix_sub_pc_type ilu -prefix_sub_pc_factor_levels 0 "
# LIBMESH_OPTS=" -log_view -ksp_monitor -ksp_view -info myinfo.txt -ksp_type minres -pc_type bjacobi -sub_pc_type ilu"
# LIBMESH_OPTS="-ksp_type gmres -pc_type bjacob -sub_pc_type ilu -sub_pc_factor_levels 10"

TUMOR_MODEL_LIBMESH_OPTS=${TUMOR_MODEL_LIBMESH_OPTS:-" -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 0"}

model_name="$1"
sim_dir="$2"
n_mpi="$3"
run_screen="$4"
pp_tag=$5

(

# if sim dir exist, clean it first
if [[ -d $sim_dir ]]; then
	rm $sim_dir/*vtu > /dev/null 2>&1
	rm $sim_dir/*txt > /dev/null 2>&1
	rm $sim_dir/*vtk > /dev/null 2>&1
	rm $sim_dir/*log > /dev/null 2>&1
fi

cd $sim_dir

echo "$pwd"

if [[ $n_mpi == "0" ]]; then

    if [[ $run_screen == "0" ]]; then
        "$EXEC_DIR/$model_name" input.in $LIBMESH_OPTS
    else
        screen -dmS $pp_tag-$model_name "$EXEC_DIR/$model_name" input.in $LIBMESH_OPTS
    fi
else

    if [[ $run_screen == "0" ]]; then
        mpiexec -np $n_mpi "$EXEC_DIR/$model_name" input.in $LIBMESH_OPTS
    else
        screen -dmS $pp_tag-$model_name mpiexec -np $n_mpi "$EXEC_DIR/$model_name" input.in $LIBMESH_OPTS
	echo 'starting in shell ' + ${pp_tag-$model_name}
    fi
fi

) |& tee "$sim_dir/sim.log"
