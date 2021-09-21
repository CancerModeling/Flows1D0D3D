import os
import json
import numpy as np
import argparse
from matplotlib import pyplot as plt
plt.style.use('seaborn-pastel')

parser = argparse.ArgumentParser(description='Animator for the vessel data.')
parser.add_argument('--filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
parser.add_argument('--vessel-by-edge-id', type=int, help='The edge id of the vessel tip to plot.', required=True)
parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
parser.add_argument('--dofs', type=int, nargs='+',  help='A list of dofs to observer.', default=[-1])

args = parser.parse_args()

directory_name = os.path.dirname(args.filepath)

with open(args.filepath) as f:
    metainfo = json.loads(f.read())

t = np.array(metainfo['times'])

start_index = sum(t < args.t_start)


def find_vessel_by_edge_id(edge_id):
    for v in metainfo['vertices']:
        if v['neighbor_edge_id'] == edge_id:
            return v
    return None


def reformat_data(d):
    if len(d.shape) == 1:
        d = d.reshape((len(d), 1))
    return d


def load_data(edge_id):
    v = find_vessel_by_edge_id(edge_id)
    d_list = []
    d_list.append(reformat_data(np.loadtxt(os.path.join(directory_name, v['filepaths']['filepath_p']))))
    if 'filepath_c' in v['filepaths']:
        d_list.append(reformat_data(np.loadtxt(os.path.join(directory_name, v['filepaths']['filepath_c']))))
    return d_list, v


data,vessel = load_data(args.vessel_by_edge_id)

num_plot_rows = len(data)
if vessel['outflow_type'] == 'rcl':
    num_plot_rows += 1

fig, axes = plt.subplots(num_plot_rows, len(args.dofs), squeeze=False, sharey='row')

indices = list(range(data[0].shape[1]))

num_dofs = vessel['num_dofs']
if vessel['outflow_type'] == 'rcl':
    num_p_dof = int(num_dofs / 2)
    num_q_dof = int(num_dofs / 2)
else:
    num_p_dof = num_dofs
    num_q_dof = 0

for i,dof in enumerate(args.dofs):
    ax = axes[0, i]
    ax.plot(t[start_index:], data[0][start_index:,dof] / 1.3333, label='{}'.format(indices[dof]), linewidth=3)
    ax.legend()
    if i == 0:
        ax.set_ylabel('$p_c$ [mmHg]')
    ax.set_xlabel('$t$')
    ax.grid(True)

    if vessel['outflow_type'] == 'rcl':
        ax = axes[1, i]
        ax.plot(t[start_index:], data[0][start_index:,dof+num_p_dof], label='{}'.format(indices[dof]), linewidth=3)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$q_c [cm^3 s^{-1}]$')
        ax.set_xlabel('$t$')
        ax.grid(True)

    if (len(data) > 1):
        ax = axes[2 if vessel['outflow_type'] == 'rcl' else 1, i]
        ax.plot(t[start_index:], data[1][start_index:, dof], label='{}'.format(indices[dof]), linewidth=3)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$c [mmol/cm^3]$')
        ax.set_xlabel('$t$')
        ax.grid(True)

plt.show()


# fig, axes = plt.subplots(1, 1, squeeze=False, sharey=True)
# ax = axes[0,0]
# p_A = data[start_index:,3] / 1.3333
# p_C = data[start_index:,4] / 1.3333
# p = p_C + vessel['R1'] / vessel['resistance'][0]  * (p_A - p_C)
# ax.plot(t[start_index:], p, label='{}'.format(indices[dof]), linewidth=3)
# ax.set_xlabel('$t$')
# ax.grid(True)
# plt.show()
