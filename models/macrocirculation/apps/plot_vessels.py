import numpy as np
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--no-a', help='do not output A', action='store_true')
parser.add_argument('--no-q', help='do not output Q', action='store_true')
parser.add_argument('--no-p', help='do not output p', action='store_true')
parser.add_argument('--no-c', help='do not output c', action='store_true')
parser.add_argument('--use-shifted-vessel-numbers', help='shifts the vessel index by one, to get intuitive labels for the CoW', action='store_true')

args = parser.parse_args()

if args.use_shifted_vessel_numbers:
    args.vessels = [vid-1 for vid in args.vessels]

num_rows = 0
if not args.no_a:
    num_rows += 1
if not args.no_q:
    num_rows += 1
if not args.no_p:
    num_rows += 1
if not args.no_c:
    num_rows += 1

fig = plt.figure()
fig.tight_layout()
axes = fig.subplots(num_rows, len(args.vessels), sharey='row', sharex='col')

for idx, vessel_id in enumerate(args.vessels):
    if len(args.vessels) == 1 and num_rows > 1:
        axidx = 0
        if not args.no_a:
            ax_a = axes[axidx]; axidx += 1
        if not args.no_q:
            ax_q = axes[axidx]; axidx += 1
        if not args.no_p:
            ax_p = axes[axidx]; axidx += 1
        if not args.no_c:
            ax_c = axes[axidx]; axidx += 1
    elif len(args.vessels) > 1 and num_rows > 1:
        axidx = 0
        if not args.no_a:
            ax_a = axes[axidx][idx]; axidx += 1
        if not args.no_q:
            ax_q = axes[axidx][idx]; axidx += 1
        if not args.no_p:
            ax_p = axes[axidx][idx]; axidx += 1
        if not args.no_c:
            ax_c = axes[axidx][idx]; axidx += 1
    elif len(args.vessels) > 1 and num_rows == 1:
        if not args.no_a:
            ax_a = axes[idx]
        if not args.no_q:
            ax_q = axes[idx]
        if not args.no_p:
            ax_p = axes[idx]
        if not args.no_c:
            ax_c = axes[idx]
    else:
        print('not implemented')

    a = np.loadtxt('output/data_A_vessel{:05d}.csv'.format(vessel_id), delimiter=',')
    a = a[1:]
    a = a[:, int(a.shape[1]/2)]

    q = np.loadtxt('output/data_Q_vessel{:05d}.csv'.format(vessel_id), delimiter=',')
    q = q[1:]
    q = q[:, int(q.shape[1]/2)]

    p = np.loadtxt('output/data_p_vessel{:05d}.csv'.format(vessel_id), delimiter=',')
    p = p[1:]
    p = p[:, int(p.shape[1]/2)]

    c = np.loadtxt('output/data_c_vessel{:05d}.csv'.format(vessel_id), delimiter=',')
    c = c[1:]
    c = c[:, int(c.shape[1]/2)]

    t = np.loadtxt('output/data_times.csv', delimiter=',')

    start_index = np.sum(t < args.t_start)
    end_index = np.sum(t < args.t_end)

    array_sizes = []
    if not args.no_a:
        array_sizes.append(len(a))
    if not args.no_q:
        array_sizes.append(len(q))
    if not args.no_p:
        array_sizes.append(len(p))
    if not args.no_c:
        array_sizes.append(len(c))

    end_index = min(min(array_sizes), end_index)
    print(end_index)

    t = t[start_index:end_index]
    a = a[start_index:end_index]
    q = q[start_index:end_index]
    p = p[start_index:end_index]
    c = c[start_index:end_index]

    if not args.no_a:
        ax_a.plot(t, a, label='A_{}'.format(vessel_id))
        ax_a.legend()
        ax_a.grid(True)
    if not args.no_q:
        ax_q.plot(t, q, label='q ({})'.format(vessel_id))
        ax_q.legend()
        ax_q.grid(True)
    if not args.no_p:
        ax_p.plot(t, p, label='p ({})'.format(vessel_id))
        ax_p.legend()
        ax_p.grid(True)
    if not args.no_c:
        ax_c.plot(t, c/a, label='$\Gamma_{{{}}}/A_{{{}}}$'.format(vessel_id, vessel_id))
        ax_c.legend()
        ax_c.grid(True)

plt.show()

