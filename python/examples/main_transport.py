import json
import argparse
from dataclasses import dataclass
import os
import flows1d0d3d as f
import dolfin as df
import numpy as np
from _utils import read_mesh, open_input_pressures
from implicit_pressure_solver import ImplicitPressureSolver
from transport_solver import TransportSolver 


def cli():
    parser = argparse.ArgumentParser(description='Full 1d0d3d transport solver.')
    parser.add_argument("--use-fully-coupled", action="store_true", help="Should the fully coupled solver be used?")
    parser.add_argument("--tip-pressures-input-file", type=str, help="Should the pressures be initialized?", required=True)
    args = parser.parse_args()
    return args


def run():
    args = cli()

    data_folder = '../../data'
    output_folder = './tmp'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau_out = 1. / 2 ** 3
    tau_coup = tau_out
    t_end = 80.
    t = 0
    t_coup_start = 2.
    t_preiter = 6 

    if args.use_fully_coupled:
        tau = 1. / 2 ** 16
        tau_transport = 1. / 2 ** 7
        solver1d = f.FullyCoupledHeartToBreast1DSolver()
        solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"))
        solver1d.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))
        solver1d.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"))
        solver1d.set_start_time_transport(t_preiter)
    else:
        tau = 1. / 2 ** 4
        tau_transport = 1. / 2 ** 4
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))

    solver1d.set_output_folder(output_folder)
    solver1d.setup(degree, tau)

    tip_pressures = open_input_pressures(args.tip_pressures_input_file)
    solver1d.update_vessel_tip_pressures(tip_pressures)
    
    coupling_interval_0d3d = int(round(tau_coup / tau_out))

    for j in range(int((t_preiter-t)/tau_transport)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('preiter = {}, t = {}'.format(j, t))
        t_old = t
        t = solver1d.solve_flow(tau, t_old, int(tau_transport/ tau))
        assert np.isclose(t_old + tau_transport, t)
        solver1d.solve_transport(tau_transport, t)

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    mesh = read_mesh(mesh_3d_filename)
    flow_solver_3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(flow_solver_3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), flow_solver_3d.V.sub(0).collapse()))
    df.assign(flow_solver_3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), flow_solver_3d.V.sub(1).collapse()))
    # initialize the 3D solver: 
    flow_solver_3d.update_vessel_tip_pressures(vessel_tip_pressures)
    flow_solver_3d.solve()
    flow_solver_3d.write_solution(0.)
    flow_solver_3d.write_subdomains(output_folder)

    transport_solver_3d = TransportSolver(mesh, output_folder, flow_solver_3d, vessel_tip_pressures)

    for i in range(0, int((t_end-t_preiter) / tau_out)):
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('iter = {}, t = {}'.format(i, t))

        for j in range(int(tau_out/tau_transport)):
            t_old = t
            t = solver1d.solve_flow(tau, t_old, int(tau_transport/ tau))
            assert np.isclose(t_old + tau_transport, t)
            solver1d.solve_transport(tau_transport, t)
        solver1d.write_output(t)

        if (t > t_coup_start) and (i % coupling_interval_0d3d == 0):
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D pressures')
            flow_solver_3d.update_vessel_tip_pressures(vessel_tip_pressures)
            flow_solver_3d.solve()
            flow_solver_3d.write_solution(t)
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D pressures')
            new_pressures = flow_solver_3d.get_pressures()
            solver1d.update_vessel_tip_pressures(new_pressures)

            # fix me with real concentrations 
            # for v in vessel_tip_pressures:
            #     v.concentration = 1
            # vessel_tip_pressures[0].concentration = 1

            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D transport')
            transport_solver_3d.update_average_pressures(new_pressures)
            transport_solver_3d.update_vessel_tip_concentrations(vessel_tip_pressures)
            transport_solver_3d.solve()
            transport_solver_3d.write_solution(t)
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D transport')
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print ('new_pressures', new_pressures)
                print ('vessel_tip_pressures = [')
                for v in vessel_tip_pressures:
                    print("VesselTipData(")
                    print('p = Point(x={}, y={}, z={}),'.format(v.p.x, v.p.y, v.p.z))
                    print('vertex_id = {},'.format(v.vertex_id))
                    print('pressure = {},'.format(v.pressure))
                    print('concentration = {},'.format(v.concentration))
                    print('R2 = {},'.format(v.R2))
                    print('radius_first = {},'.format(v.radius_first))
                    print('radius_last = {},'.format(v.radius_last))
                    print('level = {}'.format(v.level))
                    print("),")
                print (']')

            max_value = transport_solver_3d.current.vector().max()
            min_value = transport_solver_3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('transport max = {}, min = {}'.format(max_value, min_value))

            if transport_solver_3d.current.vector().min() < -10:
                os.sys.exit()

            max_value = flow_solver_3d.current.vector().max()
            min_value = flow_solver_3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('flow max = {}, min = {}'.format(max_value / 1333, min_value / 1333))


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()