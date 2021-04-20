import numpy as np
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)

args = parser.parse_args()

fig = plt.figure()
fig.tight_layout()
axes = fig.subplots(4, len(args.vessels), sharey='row', sharex='col')

for idx, vessel_id in enumerate(args.vessels):
    if len(args.vessels) == 1:
        ax0 = axes[0]
        ax1 = axes[1]
        ax2 = axes[2]
        ax3 = axes[3]
    else:
        ax0 = axes[0][idx]
        ax1 = axes[1][idx]
        ax2 = axes[2][idx]
        ax3 = axes[3][idx]

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
    end_index = min([len(t), len(a), len(p), len(q), len(c)])

    #t = np.arange(0,len(p))

    t = t[start_index:end_index]
    a = a[start_index:end_index]
    p = p[start_index:end_index]
    q = q[start_index:end_index]
    c = c[start_index:end_index]

    ax0.plot(t, a, label='a ({})'.format(vessel_id))
    ax1.plot(t, q, label='q ({})'.format(vessel_id))
    ax2.plot(t, p, label='p ({})'.format(vessel_id))
    ax3.plot(t, c/a, label='c/a ({})'.format(vessel_id))
    ax0.legend()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax0.grid(True)
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)

plt.show()
