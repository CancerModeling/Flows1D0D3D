import os


def call_executibles(dps, run_screen=True):
    print('calling executible ... ')
    for i,dp in enumerate(dps):
        run_screen = 1 if run_screen else 0  # 1 - true, 0 - false
        os.system('./run.sh ' + dp.model_name + ' ' + dp.run_path + ' ' + str(dp.n_mpi) + ' ' + str(run_screen) + ' ' + dp.pp_tag)

