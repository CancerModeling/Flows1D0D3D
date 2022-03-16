from matplotlib import pyplot as plt
import time
import argparse
import os
import flows1d0d3d as f
import numpy as np


def cli():
    parser = argparse.ArgumentParser(description='Full 1d0d3d solver.')
    parser.add_argument('--t-end', type=float, help='When should the simulation stop', default=81)
    parser.add_argument("--data-folder", type=str, help="From which directory should get the input files?", default='../../data')
    parser.add_argument("--geometry-id", type=int, help="Which breast geometry should we use", default=1)

    args = parser.parse_args()
    return args


list_taus = [1./2**k for k in [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2]]


def run(num_steps, num_samples, tau):
    args = cli()

    data_folder = args.data_folder

    degree = 2
    gid = args.geometry_id
    t = 0

    list_runtimes_fixed_tau = []
    for i in range(num_samples):
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, f"1d-meshes/coarse-breast{gid}-geometry-with-extension.json"))
        solver1d.setup(degree, tau)
        start = time.time()
        t = solver1d.solve_flow(tau, t, num_steps)
        end = time.time()
        time_per_step = (end - start)
        list_runtimes_fixed_tau.append(time_per_step)

    return np.mean(list_runtimes_fixed_tau)


def run_fixed_steps():
    #num_steps = 100
    #num_samples = 20
    num_steps = 100
    num_samples = 3

    list_runtimes = []

    average_runtimes_nl = np.min([0.0040298, 0.00153925, 0.00153967])

    for tau in list_taus:
        list_runtimes.append(run(num_steps, num_samples, tau) / num_steps)

    print(f'elapsed: {list_runtimes}')

    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    def convert(ax):
        y1,y2 = ax.get_ylim()
        ax2.set_ylim(y1/average_runtimes_nl, y2/average_runtimes_nl)
        ax2.figure.canvas.draw()

    ax.callbacks.connect('ylim_changed', convert)

    ax.semilogx(list_taus, list_runtimes, label='linearized')
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