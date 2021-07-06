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

with open(args.filepath) as f:
    metainfo = json.loads(f.read())

t = np.array(metainfo['times'])

start_index = sum(t < args.t_start)


def find_vessel_by_edge_id(edge_id):
    for v in metainfo['vertices']:
        if v['neighbor_edge_id'] == edge_id:
            return v
    return None


def load_data(edge_id):
    v = find_vessel_by_edge_id(edge_id)
    print(v)
    d = np.loadtxt(v['filepath'])
    if len(d.shape) == 1:
        d = d.reshape((len(d), 1))
    return d, v


data,vessel = load_data(args.vessel_by_edge_id)

fig, axes = plt.subplots(1, len(args.dofs), squeeze=False, sharey=True)

indices = list(range(data.shape[1]))

for i,dof in enumerate(args.dofs):
    ax = axes[0,i]
    ax.plot(t[start_index:], data[start_index:,dof] / 1.3333, label='{}'.format(indices[dof]), linewidth=3)
    ax.legend()
    if i == 0:
        ax.set_ylabel('$p_c$ [mmHg]')
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

