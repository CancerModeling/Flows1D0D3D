import dolfin as df
from matplotlib import pyplot as plt
import time
import argparse
import os
import flows1d0d3d as f
import numpy as np
from _utils import read_mesh, open_input_pressures, AverageQuantityWriter
from implicit_pressure_solver import ImplicitPressureSolver


def cli():
    parser = argparse.ArgumentParser(description='Full 1d0d3d solver.')
    parser.add_argument("--data-folder", type=str, help="From which directory should get the input files?", default='../../data')
    parser.add_argument("--geometry-id", type=int, help="Which breast geometry should we use", default=1)

    args = parser.parse_args()
    return args


#list_taus = [1./2**k for k in [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]]
#list_taus = [1./2**k for k in [4, 3, 2]]
list_taus = [1./2**k for k in [16, 14, 12, 10, 8, 6, 4, 2]]


def run(num_steps, num_samples, tau):
    args = cli()

    data_folder = args.data_folder

    degree = 2
    gid = args.geometry_id
    t = 0
    t_preiter = 1

    mesh_3d_filename = os.path.join(data_folder, f'3d-meshes/test{gid}_full_1d0d3d_cm.xdmf')
    output_folder = 'tmp'

    list_runtimes_fixed_tau = []
    for i in range(num_samples):
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, f"1d-meshes/coarse-breast{gid}-geometry-with-extension.json"))
        solver1d.setup(degree, tau)

        t = solver1d.solve_flow(tau, t, int(t_preiter / tau))

        vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

        mesh = read_mesh(mesh_3d_filename)
        solver3d = ImplicitPressureSolver(mesh, output_folder, vessel_tip_pressures)
        df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
        df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
        # initialize the 3D solver:
        solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
        solver3d.solve()
        solver3d.write_solution(0.)
        solver3d.write_subdomains(output_folder)

        total_time = 0
        for i in range(num_steps):
            t = solver1d.solve_flow(tau, t, 1)
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)

            start = time.time()
            solver3d.solve()
            end = time.time()

            new_pressures = solver3d.get_pressures()
            solver1d.update_vessel_tip_pressures(new_pressures)
            total_time += (end - start)
        list_runtimes_fixed_tau.append(total_time)

    return np.mean(list_runtimes_fixed_tau)


def run_fixed_steps():
    num_steps = 50
    num_samples = 10

    list_runtimes = []

    average_runtimes_nl = np.min([0.0040298, 0.00153925, 0.00153967])

    for tau in list_taus:
        print(f'run tau = {tau}')
        list_runtimes.append(run(num_steps, num_samples, tau) / num_steps)
        print(f'elapsed: {list_runtimes}')

    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    def convert(ax):
        y1,y2 = ax.get_ylim()
        ax2.set_ylim(y1/average_runtimes_nl, y2/average_runtimes_nl)
        ax2.figure.canvas.draw()

    ax.callbacks.connect('ylim_changed', convert)

    ax.semilogx(list_taus, list_runtimes, label='3D')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel(r'absolute runtime per time step [s]')
    ax.hlines(np.min(average_runtimes_nl), list_taus[0], list_taus[-1], linestyle='dashed', label='nonlinear')
    ax.legend()
    ax.grid(True)

    ax2.set_ylabel(r'relative runtime w.r.t nonlinear scheme')

    plt.show()


if __name__ == '__main__':
    f.initialize()
    run_fixed_steps()
    f.finalize()