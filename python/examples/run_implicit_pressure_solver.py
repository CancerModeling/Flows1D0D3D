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
    parser.add_argument("--tip-pressures-input-file", type=str, help="The 1D pressures.", required=True)
    parser.add_argument("--geometry-id", type=int, help="Which breast geometry should we use", default=1)
    parser.add_argument("--data-folder", type=str, help="From which directory should get the input files?", default='../../data')
    parser.add_argument("--output-folder", type=str, help="Into which directory should we write?", default='./tmp_flow')
    parser.add_argument("--refinements", type=int, help="Number of refinements", default=0)

    args = parser.parse_args()
    return args


def run():
    args = cli()

    data_folder = args.data_folder
    output_folder = args.output_folder

    os.makedirs(output_folder, exist_ok=True)

    gid = args.geometry_id

    mesh_3d_filename = os.path.join(data_folder, f'3d-meshes/test{gid}_full_1d0d3d_cm.xdmf')

    tip_pressures = open_input_pressures(args.tip_pressures_input_file)
    print(tip_pressures)

    # dummy 1d solver to get the information from the mesh:
    solver1d = f.LinearizedHeartToBreast1DSolver()
    solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
    solver1d.set_path_geometry(os.path.join(data_folder, f"1d-meshes/coarse-breast{gid}-geometry-with-extension.json"))
    solver1d.set_output_folder(output_folder)
    solver1d.setup(0, 0)
    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
    for v in vessel_tip_pressures:
        v.pressure = tip_pressures[v.vertex_id]

    # 3d solver:
    mesh = read_mesh(mesh_3d_filename)
    for i in range(args.refinements):
        mesh = df.refine(mesh)
    solver3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    # initialize the 3D solver:

    # prerun:
    solver3d.solve()
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))

    print('s')
    start_simulation = time.time()
    solver3d.solve()
    elapsed = time.time() - start_simulation
    print(f'time = {elapsed} s')


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
