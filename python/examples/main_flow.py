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
    args = parser.parse_args()
    return args


def run():
    args = cli()

    data_folder = '../../data'
    output_folder = './tmp_flow'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau_out = 1. / 2 ** 6
    #tau_out = 1.
    #tau_coup = 2.
    tau_coup = 1. / 2 ** 3
    t_end = 80.
    t = 0
    t_coup_start = 2.

    if args.use_fully_coupled:
        tau = 1. / 2 ** 16
        solver1d = f.FullyCoupledHeartToBreast1DSolver()
        solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"))
        solver1d.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))
        solver1d.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"))
    else:
        tau = 1. / 2 ** 4
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))
    solver1d.set_output_folder(output_folder)
    solver1d.setup(degree, tau)

    # if an input pressure file is present we read it in
    if args.tip_pressures_input_file is not None:
        tip_pressures = open_input_pressures(args.tip_pressures_input_file)
        solver1d.update_vessel_tip_pressures(tip_pressures)
    
    coupling_interval_0d3d = int(round(tau_coup / tau_out))

    t = solver1d.solve_flow(tau, t, int(6 / tau))

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    q_writer = AverageQuantityWriter(vessel_tip_pressures)

    mesh = read_mesh(mesh_3d_filename)
    solver3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    # initialize the 3D solver: 
    solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
    solver3d.solve()
    solver3d.write_solution(0.)
    solver3d.write_subdomains(output_folder)

    for i in range(0, int(t_end / tau_out)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('iter = {}, t = {}'.format(i, t))
        t = solver1d.solve_flow(tau, t, int(tau_out / tau))
        solver1d.write_output(t)
        if (t > t_coup_start) and (i % coupling_interval_0d3d == 0):
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D pressures')
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
            solver3d.solve()
            solver3d.write_solution(t)
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D pressures')
            new_pressures = solver3d.get_pressures()
            solver1d.update_vessel_tip_pressures(new_pressures)

            max_value = solver3d.current.vector().max()
            min_value = solver3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))

        q_writer.update(t, new_pressures, solver3d.get_pressures(1))
        q_writer.write(os.path.join(output_folder, "average_quantities.json"))
        
        if args.tip_pressures_output_file is not None and df.MPI.rank(df.MPI.comm_world) == 0:
            print('writing tip pressures to {}'.format(args.tip_pressures_output_file))
            with open(args.tip_pressures_output_file, 'w') as file:
                file.write(json.dumps(new_pressures))


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()