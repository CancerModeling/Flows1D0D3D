import time
import argparse
import json
import os
import flows1d0d3d as f
import dolfin as df
import numpy as np
from _utils import read_mesh, open_input_pressures, AverageQuantityWriter, StopWatch
from implicit_pressure_solver import ImplicitPressureSolver, estimate_coeffs_ca, get_volumes
from parameters import FlowModelParameters


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
    parser.add_argument("--use-integrated-average-pressures", action='store_true')
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

    compare_pressure_multiplier_integral = False

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

    print('start 1d preiter')
    t = solver1d.solve_flow(tau, t, int(t_preiter / tau))
    print('end 1d preiter')

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    q_writer = AverageQuantityWriter(vessel_tip_pressures)

    mesh = read_mesh(mesh_3d_filename)
    volumes, volume = get_volumes(vessel_tip_pressures, mesh)

    flow_config = FlowModelParameters()
    factor = estimate_coeffs_ca(flow_config, vessel_tip_pressures, volumes) / 1.05214002e-06
    print(f'factor {factor}')
    flow_config.L_cv *= factor
    flow_config.L_ct *= factor
    flow_config.L_tl *= factor
    print(f'L_cv = {flow_config.L_cv}, L_ct = {flow_config.L_ct}, L_tl = {flow_config.L_tl}')

    solver3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures, flow_config)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    # initialize the 3D solver:
    solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
    print('start 3d preiter')
    solver3d.solve()
    print('end 3d preiter')
    if not args.disable_output:
        solver3d.write_solution(0.)
        solver3d.write_subdomains(output_folder)

    new_pressures = solver3d.get_pressures()

    watch_simulation = StopWatch()
    watch_1d = StopWatch()
    watch_3d = StopWatch()
    watch_1d_coup_out = StopWatch()
    watch_3d_coup_out = StopWatch()
    watch_1d_coup_in = StopWatch()
    watch_3d_coup_in = StopWatch()
    watch_io = StopWatch()

    watch_simulation.start()
    for i in range(0, int(t_end / tau_out)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('iter = {}, t = {}'.format(i, t))
        steps = int(tau_out / tau)
        watch_1d.start()
        t = solver1d.solve_flow(tau, t, steps)
        watch_1d.end(num_steps=steps)

        if not args.disable_output:
            watch_io.start()
            solver1d.write_output(t)
            watch_io.end()
        if (t > t_coup_start) and (i % coupling_interval_0d3d == 0):
            watch_1d_coup_out.start()
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            watch_1d_coup_out.end()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D pressures')
            watch_3d_coup_in.start()
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
            watch_3d_coup_in.end()

            watch_3d.start()
            solver3d.solve()
            watch_3d.end()

            if not args.disable_output:
                watch_io.start()
                solver3d.write_solution(t)
                watch_io.end()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D pressures')
            watch_3d_coup_out.start()
            new_pressures = solver3d.get_pressures() if args.use_integrated_average_pressures else solver3d.get_pressures_from_multiplier()
            if compare_pressure_multiplier_integral:
                new_pressures_old = solver3d.get_pressures()
                new_pressures_mul = solver3d.get_pressures_from_multiplier()
                for x in new_pressures_old.keys():
                    print(f'{x} {new_pressures_old[x]} {new_pressures_mul[x]} {np.abs(new_pressures_old[x] - new_pressures_mul[x])}')
            watch_3d_coup_out.end()
            watch_1d_coup_in.start()
            solver1d.update_vessel_tip_pressures(new_pressures)
            watch_1d_coup_in.end()

            if not args.disable_output:
                max_value = solver3d.current.vector().max()
                min_value = solver3d.current.vector().min()
                if df.MPI.rank(df.MPI.comm_world) == 0:
                    print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))

        if not args.disable_output:
            watch_io.start()
            q_writer.update(t, new_pressures, solver3d.get_pressures(1))
            q_writer.write(os.path.join(output_folder, "average_quantities.json"))
            watch_io.end()

        if args.tip_pressures_output_file is not None and df.MPI.rank(df.MPI.comm_world) == 0 and not args.disable_output:
            print('writing tip pressures to {}'.format(args.tip_pressures_output_file))
            with open(args.tip_pressures_output_file, 'w') as file:
                watch_io.start()
                file.write(json.dumps(new_pressures))
                watch_io.end()

        print(f'3d solver       {watch_3d}')
        print(f'1d solver       {watch_1d}')
        print(f'1d coupling out {watch_1d_coup_out}')
        print(f'1d coupling in  {watch_1d_coup_in}')
        print(f'3d coupling out {watch_3d_coup_out}')
        print(f'3d coupling in  {watch_3d_coup_in}')
        print(f'io              {watch_io}')

    watch_simulation.end()

    print(f'elapsed simulation time: {watch_simulation.average()}')


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()