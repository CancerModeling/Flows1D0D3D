import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--no-a', help='do not output A', action='store_true')
parser.add_argument('--no-q', help='do not output Q', action='store_true')
parser.add_argument('--no-p', help='do not output p', action='store_true')
parser.add_argument('--no-c', help='do not output c', action='store_true')
parser.add_argument('--use-shifted-vessel-numbers', help='shifts the vessel index by one, to get intuitive labels for the CoW', action='store_true')

args = parser.parse_args()

directory = os.path.dirname(args.filepath)

if args.use_shifted_vessel_numbers:
    args.vessels = [vid-1 for vid in args.vessels]


with open(args.filepath) as f:
    meta = json.loads(f.read())


def find_vessel(vessel_id):
    for vessel in  meta['vessels']:
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


fig = plt.figure()
fig.tight_layout()
axes = fig.subplots(num_rows, len(args.vessels), sharey='row', sharex='col', squeeze=False)


for idx, vessel_id in enumerate(args.vessels):
    vessel_info = find_vessel(vessel_id)

    if ('a' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['a'])
        print('loading {}'.format(path))
        a = np.loadtxt(path, delimiter=',')
        a = a[:]
        a = a[:, int(a.shape[1]/2)]

    path = os.path.join(directory, vessel_info['filepaths']['q'])
    print('loading {}'.format(path))
    q = np.loadtxt(path, delimiter=',')
    q = q[:]
    q = q[:, int(q.shape[1]/2)]

    if ('p' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['p'])
        print('loading {}'.format(path))
        p = np.loadtxt(path, delimiter=',') / 1.333332
        p = q[:]
        p = q[:, int(q.shape[1]/2)]
    else:
        p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1) / 1.33332

    if ('c' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['c'])
        print('loading {}'.format(path))
        c = np.loadtxt(path, delimiter=',')
        c = c[:]
        c = c[:, int(c.shape[1]/2)]

    path = os.path.join(directory, meta['filepath_time'])
    print('loading {}'.format(path))
    t = np.loadtxt(path, delimiter=',')

    start_index = np.sum(t < args.t_start)
    end_index = np.sum(t < args.t_end)

    array_sizes = []
    if not args.no_a and 'a' in vessel_info['filepaths']:
        array_sizes.append(len(a))
    if not args.no_q:
        array_sizes.append(len(q))
    if not args.no_p:
        array_sizes.append(len(p))
    if not args.no_c and 'c' in vessel_info['filepaths']:
        array_sizes.append(len(c))

    end_index = min(min(array_sizes), end_index)

    t = t[start_index:end_index]
    a = a[start_index:end_index]
    q = q[start_index:end_index]
    p = p[start_index:end_index]
    c = c[start_index:end_index]

    row_idx = 0
    if not args.no_a and 'a' in vessel_info['filepaths']:
        ax = axes[row_idx, idx]
        ax.plot(t, a, label='A_{}'.format(vessel_id))
        ax.legend()
        ax.grid(True)
        row_idx += 1
    if not args.no_q:
        ax = axes[row_idx, idx]
        ax.plot(t, q, label='q ({})'.format(vessel_id))
        ax.legend()
        ax.grid(True)
        row_idx += 1
    if not args.no_p:
        ax = axes[row_idx, idx]
        ax.plot(t, p, label='p ({})'.format(vessel_id))
        ax.legend()
        ax.grid(True)
        row_idx += 1
    if not args.no_c and 'c' in vessel_info['filepaths']:
        ax = axes[row_idx, idx]
        ax.plot(t, c/a, label='$\Gamma_{{{}}}/A_{{{}}}$'.format(vessel_id, vessel_id))
        ax.legend()
        ax.grid(True)
        row_idx += 1

plt.show()

