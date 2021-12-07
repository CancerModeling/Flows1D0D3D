import os
import numpy as np
import json
import matplotlib as mpl
from matplotlib import pyplot as plt
import argparse


font = {'family': 'normal', 'size': 22}
mpl.rc('font', **font)


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--filepaths', type=str, nargs='+', required=True)
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--labels', type=str, nargs='+', default=[])
parser.add_argument('--colors', type=str, nargs='+', default=[])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--output-folder', type=str, required=True)
parser.add_argument('--dataset-name', type=str, required=True)

args = parser.parse_args()

list_all_t = []
list_all_p = []
list_all_q = []

colors = args.colors

no_legend = True

labels = args.labels
if len(labels) != len(args.filepaths):
    labels = ['' for idx in range(args.filepaths)]


for filepath in args.filepaths:
    directory = os.path.dirname(filepath)

    with open(filepath) as f:
        meta = json.loads(f.read())

    def find_vessel(vessel_id):
        for vessel in  meta['vessels']:
            if vessel['edge_id'] == vessel_id:
                return vessel
        raise 'vessel with id {} not found'.format(vessel_id)


    vessel_info = meta['vessels'][0]

    list_t = []
    list_q = []
    list_p = []

    for idx, vessel_id in enumerate(args.vessels):
        vessel_info = find_vessel(vessel_id)

        path = os.path.join(directory, vessel_info['filepaths']['q'])
        print('loading {}'.format(path))
        q = np.loadtxt(path, delimiter=',')
        q = q[:]
        q = q[:, int(q.shape[1]/2)]

        if q.mean() < 0:
            q *= -1

        if ('a' in vessel_info['filepaths']):
            path = os.path.join(directory, vessel_info['filepaths']['a'])
            print('loading {}'.format(path))
            a = np.loadtxt(path, delimiter=',')
            a = a[:]
            a = a[:, int(a.shape[1]/2)]
        else:
            a = np.ones(q.shape) * vessel_info['A0']

        if ('p' in vessel_info['filepaths']):
            path = os.path.join(directory, vessel_info['filepaths']['p'])
            print('loading {}'.format(path))
            p = np.loadtxt(path, delimiter=',') / 1.333332
            print ('a', p)
            p = p[:]
            print ('b', p)
            p = p[:, int(p.shape[1]/2)]
            print ('c', p)
        else:
            p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1) / 1.33332

        path = os.path.join(directory, meta['filepath_time'])
        print('loading {}'.format(path))
        t = np.loadtxt(path, delimiter=',')

        start_index = np.sum(t < args.t_start)
        end_index = np.sum(t < args.t_end)

        list_p.append(p[start_index:end_index].tolist())
        list_q.append(q[start_index:end_index].tolist())
        list_t.append(t[start_index:end_index].tolist())

    list_all_p.append(list_p)
    list_all_q.append(list_q)
    list_all_t.append(list_t)


def flatten(t):
    if isinstance(t, list):
        return [item for sublist in t for item in flatten(sublist)]
    return [t]

max_p, min_p = np.array(flatten(list_all_p)).max(), np.array(flatten(list_all_p)).min()
max_q, min_q = np.array(flatten(list_all_q)).max(), np.array(flatten(list_all_q)).min()

if True:
    fig = plt.figure()

    for idx, vessel_id in enumerate(args.vessels):
        plt.clf()
        ax = fig.add_subplot(111) 
        for src_idx in range(len(list_all_p)):
            linestyle = '-' if src_idx % 2 == 0 else '--'
            ax.plot(list_all_t[src_idx][idx], list_all_p[src_idx][idx], label='$p_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=4-src_idx, linestyle=linestyle, color=colors[src_idx])
        if not no_legend:
            ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel('p [mmHg]')
        ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)
        plt.savefig(os.path.join(args.output_folder, '{}_p_{}.pdf'.format(args.dataset_name, vessel_id)), bbox_inches='tight')

    for idx, vessel_id in enumerate(args.vessels):
        plt.clf()
        ax = fig.add_subplot(111) 
        for src_idx in range(len(list_all_q)):
            ax.plot(list_all_t[src_idx][idx], list_all_q[src_idx][idx], label='$q_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=4-src_idx, linestyle=linestyle, color=colors[src_idx])
        if not no_legend:
            ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel(r'q [$cm^{3}/s$]')
        #ax.set_ylim(top=max_q+5, bottom=min_q-5)
        ax.grid(True)
        plt.savefig(os.path.join(args.output_folder, '{}_q_{}.pdf'.format(args.dataset_name, vessel_id)), bbox_inches='tight')

if True:
    fig, axes = plt.subplots(2, len(args.vessels), squeeze=False, sharey='row', sharex='col')

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[0, idx]
        for src_idx in range(len(list_all_p)):
            linestyle = '-' if src_idx % 2 == 0 else '--'
            ax.plot(list_all_t[src_idx][idx], list_all_p[src_idx][idx], label='$p_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=3-src_idx, linestyle=linestyle, color=colors[src_idx])
        if idx == 0:
            ax.set_ylabel('p [mmHg]')
        ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[1, idx]
        for src_idx in range(len(list_all_q)):
            linestyle = '-' if src_idx % 2 == 0 else '--'
            ax.plot(list_all_t[src_idx][idx], list_all_q[src_idx][idx], label='$q_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=3-src_idx, linestyle=linestyle, color=colors[src_idx])
        ax.set_xlabel('t [s]')
        if idx == 0:
            ax.set_ylabel(r'q [$cm^{3}/s$]')
        ax.grid(True)
    plt.show()

