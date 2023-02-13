import time
import os
import flows1d0d3d as f
import math


def cli():
    import argparse
    parser = argparse.ArgumentParser(description='Full 1d0d3d solver.')
    parser.add_argument("--t-max", type=float, help="Maximum simulation time", default=1)
    parser.add_argument("--dt", type=float, help="Time-step width", default=1e-5)
    parser.add_argument("--dt-out", type=float, help="Time-step width for output", default=1e-3)
    parser.add_argument("--output-directory", type=str, help="Output directory", default='output')

    args = parser.parse_args()
    return args


def path_to_mesh():
    import pathlib
    directory = pathlib.Path(__file__).parent.resolve()
    return os.path.join(directory, '../../data/1d-meshes/3-vessels.json')


def run():
    args = cli()

    # dummy 1d solver to get the information from the mesh:
    mesh_path = path_to_mesh()
    os.makedirs(args.output_directory, exist_ok=True)
    output_path = os.path.join(args.output_directory, '3-vessels' )

    solver1d = f.SimpleLinearizedSolver(mesh_path, output_path, args.dt)

    num_steps = int(math.ceil(args.t_max / args.dt))
    output_interval = int(math.ceil(args.dt_out / args.dt))

    print('s')
    start_simulation = time.time()
    for i in range(num_steps):
        solver1d.solve()
        if i % output_interval == 0:
            print(f'iter = {i}')
            result = solver1d.get_result(0)
            print(f'vertex 0: a = {result.a}, p = {result.p}, q = {result.q}')
            result = solver1d.get_result(1)
            print(f'vertex 1: a = {result.a}, p = {result.p}, q = {result.q}')
            solver1d.write()
    elapsed = time.time() - start_simulation
    print(f'time = {elapsed} s')


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()