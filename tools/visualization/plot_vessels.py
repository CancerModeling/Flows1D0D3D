import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--position', type=float, help='A list of ids of the vessels to plot.', default=0.5)
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--no-a', help='do not output A', action='store_true')
parser.add_argument('--no-q', help='do not output Q', action='store_true')
parser.add_argument('--no-p', help='do not output p', action='store_true')
parser.add_argument('--no-c', help='do not output c', action='store_true')
parser.add_argument('--no-legend', help='remove the legend from the plots', action='store_true')
parser.add_argument('--positive-q', help='reorients the flow q s.t. it is always positive', action='store_true')
parser.add_argument('--narrow', help='should a very narrow layout be used?', action='store_true')
parser.add_argument('--use-shifted-vessel-numbers', help='shifts the vessel index by one, to get intuitive labels for the CoW', action='store_true')

args = parser.parse_args()

directory = os.path.dirname(args.filepath)

if args.use_shifted_vessel_numbers:
    args.vessels = [vid-1 for vid in args.vessels]


with open(args.filepath) as f:
    meta = json.loads(f.read())


def find_vessel(vessel_id):
    for vessel in meta['vessels']:
        if vessel['edge_id'] == vessel_id:
            return vessel
    raise 'vessel with id {} not found'.format(vessel_id)


vessel_info = meta['vessels'][0]

num_rows = 0
if not args.no_a and 'a' in vessel_info['filepaths']:
    num_rows += 1
if not args.no_q:
    num_rows += 1
if not args.no_p:
    num_rows += 1
if not args.no_c and 'c' in vessel_info['filepaths']:
    num_rows += 1


fig = plt.figure(figsize=(4,2 * len(args.vessels))) if args.narrow else plt.figure()
fig.tight_layout()
axes = fig.subplots(num_rows, len(args.vessels), sharey='row', sharex='col', squeeze=False)

position = args.position 


for idx, vessel_id in enumerate(args.vessels):
    vessel_info = find_vessel(vessel_id)

    path = os.path.join(directory, vessel_info['filepaths']['q'])
    print('loading {}'.format(path))
    q = np.loadtxt(path, delimiter=',')
    q = q[:]
    q = q[:, int((q.shape[1]-1)*position)]
    
    if args.positive_q and q.mean() < 0:
        q *= -1

    if ('a' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['a'])
        print('loading {}'.format(path))
        a = np.loadtxt(path, delimiter=',')
        a = a[:]
        a = a[:, int((a.shape[1]-1) * position)]
    else:
        a = np.ones(q.shape) * vessel_info['A0']

    if ('p' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['p'])
        print('loading {}'.format(path))
        p = np.loadtxt(path, delimiter=',') / 1.333332
        p = p[:]
        p = p[:, int((p.shape[1]-2) * position)]
    else:
        p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1) / 1.33332

    if ('c' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['c'])
        print('loading {}'.format(path))
        c = np.loadtxt(path, delimiter=',')
        c = c[:]
        c = c[:, int((c.shape[1]-1)*position)]

    path = os.path.join(directory, meta['filepath_time'])
    print('loading {}'.format(path))
    t = np.loadtxt(path, delimiter=',')

    start_index = np.sum(t < args.t_start)
    end_index = np.sum(t < args.t_end)

    array_sizes = []
    if not args.no_a:
        array_sizes.append(len(a))
    if not args.no_q:
        array_sizes.append(len(q))
    if not args.no_p:
        array_sizes.append(len(p))
    if not args.no_c and 'c' in vessel_info['filepaths']:
        array_sizes.append(len(c))

    end_index = min(min(array_sizes), end_index)

    t = t[start_index:end_index]
    if not args.no_a:
        a = a[start_index:end_index]
    if not args.no_q:
        q = q[start_index:end_index]
    if not args.no_p:
        p = p[start_index:end_index]
    if not args.no_c and 'c' in vessel_info['filepaths']:
        c = c[start_index:end_index]

    row_idx = 0
    if not args.no_a and 'a' in vessel_info['filepaths']:
        ax = axes[row_idx, idx]
        label = None if args.no_legend else r'$A_{' + str(vessel_id) + '}$'
        ax.plot(t, a, label=label)
        if not args.no_legend:
            ax.legend()
        ax.grid(True)
        row_idx += 1
    if not args.no_p:
        ax = axes[row_idx, idx]
        label = None if args.no_legend else r'$p_{' + str(vessel_id) + '}$'
        ax.plot(t, p, label=label)
        if not args.no_legend:
            ax.legend()
        ax.grid(True)
        if idx == 0:
            ax.set_ylabel('p [mmHg]')
        row_idx += 1
    if not args.no_q:
        ax = axes[row_idx, idx]
        label = None if args.no_legend else r'$q_{' + str(vessel_id) + '}$'
        ax.plot(t, q, label=label)
        if not args.no_legend:
            ax.legend()
        ax.grid(True)
        if idx == 0:
            ax.set_ylabel('q $[cm^3/s]$')
        row_idx += 1
    if not args.no_c and 'c' in vessel_info['filepaths']:
        ax = axes[row_idx, idx]
        print(c/a)
        label = None if args.no_legend else r'$\Gamma_{' + str(vessel_id) + r'}/A_{' + str(vessel_id) + r'}$'
        ax.plot(t, c/a, label=label) 
        if not args.no_legend:
            ax.legend()
        ax.grid(True)
        row_idx += 1

for idx, vessel_id in enumerate(args.vessels):
    axes[-1, idx].set_xlabel('t [s]')

plt.tight_layout()
plt.show()

