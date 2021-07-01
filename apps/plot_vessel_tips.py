import json
import numpy as np
import argparse
import os
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
    return d


data = load_data(args.vessel_by_edge_id)

fig, axes = plt.subplots(1, len(args.dofs), squeeze=False)

indices = list(range(data.shape[1]))

for i,dof in enumerate(args.dofs):
    ax = axes[0,i]
    ax.plot(t[start_index:], data[start_index:,dof], label='{}'.format(indices[dof]))
    ax.legend()
    ax.set_ylabel('$p_c$')
    ax.set_xlabel('$t$')
    ax.grid(True)
plt.tight_layout()
plt.show()

