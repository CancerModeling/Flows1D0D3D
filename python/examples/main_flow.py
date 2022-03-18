import time
import argparse
import json
import os
import flows1d0d3d as f
import dolfin as df
import numpy as np
from _utils import read_mesh, open_input_pressures, AverageQuantityWriter
from implicit_pressure_solver import ImplicitPressureSolver


def cli():
    parser = argparse.ArgumentParser(description='Full 1d0d3d solver.')
    parser.add_argument("--use-fully-coupled", action="store_true", help="Should the fully coupled solver be used?")
    parser.add_argument("--tip-pressures-input-file", type=str, help="Should the pressures be initialized?", required=False)
    parser.add_argument("--tip-pressures-output-file", type=str, help="Should the pressures be initialized?", required=False)
    parser.add_argument('--t-3dcoup-start', type=float, help='When should the 3d coupling be active?', default=6)
    parser.add_argument('--t-preiter', type=float, help='How long should we preiterate?', default=6)
    parser.add_argument('--tau-log2', type=int, help='Defines the time step width 1/2^k where k is the parameter', required=False)
    parser.add_argument('--tau-coup-log2', type=int, help='Defines the time step width 1/2^k for the coupling where k is the parameter', required=False)
    parser.add_argument('--t-end', type=float, help='When should the simulation stop', default=81)
    parser.add_argument("--output-folder", type=str, help="Into which directory should we write?", default='./tmp_flow')
    parser.add_argument("--disable-output", action='store_true')
    parser.add_argument("--data-folder", type=str, help="From which directory should get the input files?", default='../../data')
    parser.add_argument("--geometry-id", type=int, help="Which breast geometry should we use", default=1)

    args = parser.parse_args()
    return args


def run():
    args = cli()

    data_folder = args.data_folder 
    output_folder = args.output_folder 

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    t_end = args.t_end
    t = 0
    t_coup_start = args.t_3dcoup_start 
    t_preiter = args.t_preiter 
    gid = args.geometry_id

    mesh_3d_filename = os.path.join(data_folder, f'3d-meshes/test{gid}_full_1d0d3d_cm.xdmf')

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

    tau_out = max(1. / 2 ** 6, tau)
    tau_coup = max(1. / 2 ** args.tau_coup_log2, tau)

    # if an input pressure file is present we read it in
    if args.tip_pressures_input_file is not None:
        tip_pressures = open_input_pressures(args.tip_pressures_input_file)
        solver1d.update_vessel_tip_pressures(tip_pressures)
    
    coupling_interval_0d3d = int(round(tau_coup / tau_out))

    t = solver1d.solve_flow(tau, t, int(t_preiter / tau))

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    q_writer = AverageQuantityWriter(vessel_tip_pressures)

    mesh = read_mesh(mesh_3d_filename)
    solver3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    # initialize the 3D solver: 
    solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
    solver3d.solve()
    if not args.disable_output:
        solver3d.write_solution(0.)
        solver3d.write_subdomains(output_folder)

    new_pressures = solver3d.get_pressures()

    flow_solver3d_times = []
    flow_solver1d_times = []

    start_simulation = time.time()

    for i in range(0, int(t_end / tau_out)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('iter = {}, t = {}'.format(i, t))
        steps = int(tau_out / tau)
        start = time.time()
        t = solver1d.solve_flow(tau, t, steps)
        flow_solver1d_times.append((time.time() - start) / steps)

        if not args.disable_output:
            solver1d.write_output(t)
        if (t > t_coup_start) and (i % coupling_interval_0d3d == 0):
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D pressures')
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)

            start = time.time()
            solver3d.solve()
            flow_solver3d_times.append(time.time() - start)

            if not args.disable_output:
                solver3d.write_solution(t)
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D pressures')
            new_pressures = solver3d.get_pressures()
            solver1d.update_vessel_tip_pressures(new_pressures)

            if not args.disable_output:
                max_value = solver3d.current.vector().max()
                min_value = solver3d.current.vector().min()
                if df.MPI.rank(df.MPI.comm_world) == 0:
                    print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))

        print(f'3d solver (avg {np.mean(flow_solver3d_times)}): {flow_solver3d_times}')
        print(f'1d solver (avg {np.mean(flow_solver1d_times)}, steps {steps}): {flow_solver1d_times}')

        if not args.disable_output:
            q_writer.update(t, new_pressures, solver3d.get_pressures(1))
            q_writer.write(os.path.join(output_folder, "average_quantities.json"))
        
        if args.tip_pressures_output_file is not None and df.MPI.rank(df.MPI.comm_world) == 0 and not args.disable_output:
            print('writing tip pressures to {}'.format(args.tip_pressures_output_file))
            with open(args.tip_pressures_output_file, 'w') as file:
                file.write(json.dumps(new_pressures))

    elapsed = time.time() - start_simulation
    print(f'elapsed simulation time: {elapsed}')


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()