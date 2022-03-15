import argparse
import os
import flows1d0d3d as f
import dolfin as df
from _utils import open_input_pressures 


def cli():
    parser = argparse.ArgumentParser(description='Full 1d0d3d solver.')
    parser.add_argument("--use-fully-coupled", action="store_true", help="Should the fully coupled solver be used?")
    parser.add_argument("--tip-pressures-input-file", type=str, help="Should the pressures be initialized?", required=False)
    parser.add_argument("--tip-pressures-output-file", type=str, help="Should the pressures be initialized?", required=False)
    parser.add_argument('--t-preiter', type=float, help='How long should we preiterate?', default=6)
    parser.add_argument('--tau-log2', type=int, help='Defines the time step width 1/2^k where k is the parameter', required=False)
    parser.add_argument('--t-end', type=float, help='When should the simulation stop', default=80)
    parser.add_argument("--output-folder", type=str, help="Into which directory should we write?", default='./tmp_flow')
    parser.add_argument("--data-folder", type=str, help="Into which directory should we write?", default='../../data')

    args = parser.parse_args()
    return args


def run():
    args = cli()

    data_folder = args.data_folder 
    output_folder = args.output_folder 

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau_out = 1. / 2 ** 4 
    t_end = args.t_end
    t = 0
    gid = 2

    if args.use_fully_coupled:
        tau = 1. / 2 ** args.tau_log2 if args.tau_log2 is not None else 1. / 2 ** 16
        solver1d = f.FullyCoupledHeartToBreast1DSolver()
        solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"))
        solver1d.set_path_linear_geometry(os.path.join(data_folder, f"1d-meshes/coarse-breast{gid}-geometry-with-extension.json"))
        solver1d.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"))
    else:
        tau = 1. / 2 ** args.tau_log2 if args.tau_log2 is not None else 1. / 2 ** 4
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, f"1d-meshes/coarse-breast{gid}-geometry-with-extension.json"))
    solver1d.set_output_folder(output_folder)
    solver1d.setup(degree, tau)

    # if an input pressure file is present we read it in
    if args.tip_pressures_input_file is not None:
        tip_pressures = open_input_pressures(args.tip_pressures_input_file)
        solver1d.update_vessel_tip_pressures(tip_pressures)
    
    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    for i in range(0, int(t_end / tau_out)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('iter = {}, t = {}'.format(i, t))
        t = solver1d.solve_flow(tau, t, int(tau_out / tau))
        solver1d.write_output(t)
        vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('vessel tip pressures', vessel_tip_pressures)


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()