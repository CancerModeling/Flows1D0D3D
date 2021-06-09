import numpy as np
import argparse
from matplotlib import pyplot as plt
import os
import json

parser = argparse.ArgumentParser(description='Plots the windkessel.')
parser.add_argument('--neighbor-vessels', type=int, nargs='+', help='A list of ids of the neighboring vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--no-p-c', help='do not output p_c', action='store_true')
parser.add_argument('--no-q', help='do not output q', action='store_true')
parser.add_argument('--output-directory', type=str, default='output/')

args = parser.parse_args()

path = os.path.join(args.output_directory, 'windkessel_values.json')
print(path)

with open(path) as f:
    data = json.load(f)

start_index = int(sum(np.array(data['time']) < args.t_start))

num_plots = 2
if args.no_p_c:
    num_plots -= 1
if args.no_q:
    num_plots -= 1

fig = plt.figure()
fig.tight_layout()
axes = fig.subplots(num_plots, len(args.neighbor_vessels), sharey='row', sharex='col', squeeze=False)

vertex_data_to_show = [vertex for vertex in data['vertices'] if vertex['neighbor_edge_id'] in args.neighbor_vessels]

for i in range(len(vertex_data_to_show)):
    vertex_data = vertex_data_to_show[i]

    pidx = 0

    if not args.no_p_c:
        ax = axes[pidx, i]
        ax.plot(data['time'][start_index:], vertex_data['p_c'][start_index:], label=r'$p_c^{(' + str(vertex_data['neighbor_edge_id']) + r')}$')
        ax.legend()
        ax.grid(True)
        pidx += 1

    if not args.no_q:
        ax = axes[pidx, i]
        ax.plot(data['time'][start_index:], vertex_data['q'][start_index:], label=r'$q^{(' + str(vertex_data['neighbor_edge_id']) + r')}$')
        ax.legend()
        ax.grid(True)

plt.show()