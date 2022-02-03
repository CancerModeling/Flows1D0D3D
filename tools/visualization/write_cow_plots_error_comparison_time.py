import math
import os
import numpy as np
import json
import matplotlib
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--filepath-1', type=str, required=True)
parser.add_argument('--filepaths-2', type=str, nargs='+', required=True)
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--output-folder', type=str, required=True)
parser.add_argument('--dataset-name', type=str, required=True)

args = parser.parse_args()


def load_data(filepath, vessel_ids):
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

    for idx, vessel_id in enumerate(vessel_ids):
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

        list_p.append(p[start_index:end_index])
        list_q.append(q[start_index:end_index])
        list_t.append(t[start_index:end_index])

    list_p.append(list_p)
    list_q.append(list_q)
    list_t.append(list_t)

    return list_t, list_p, list_q



all_t_1, all_p_1, all_q_1 = load_data(args.filepath_1, args.vessels)

all_err_p = [] 
all_err_q = [] 
for filepath_2 in args.filepaths_2:
    all_t_2, all_p_2, all_q_2 = load_data(filepath_2, args.vessels)
    tmp_err_p = []
    tmp_err_q = []
    for vessel_idx in range(len(args.vessels)):
        tmp_err_p.append(np.abs(all_p_1[vessel_idx] - all_p_2[vessel_idx])/all_p_1[vessel_idx])
        tmp_err_q.append(np.abs(all_q_1[vessel_idx] - all_q_2[vessel_idx])/all_q_1[vessel_idx])
    all_err_p.append(tmp_err_p)
    all_err_q.append(tmp_err_q)


if True:
    fig, axes = plt.subplots(2, len(args.vessels), squeeze=False, sharey='row', sharex='col')

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[0, idx]
        for file_idx in range(len(args.filepaths_2)):
            ax.plot(all_t_1[idx], all_err_p[file_idx][idx], label=file_idx)
            print(all_err_p[file_idx][idx])
            ax.plot(all_t_1[idx], all_err_p[file_idx][idx].mean() * np.ones(len(all_t_1[idx])), '--', label=file_idx)
        if idx == 0:
            #ax.set_ylabel(r'$\frac{|p_{bi}-p_{one}|}{p_{bi}}$')
            ax.set_ylabel('relative error p')
        #ax.set_ylim(top=max_p+5, bottom=min_p-5)
        #ax.legend()
        ax.grid(True)

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[1, idx]
        for file_idx in range(len(args.filepaths_2)):
            ax.plot(all_t_1[idx], all_err_q[file_idx][idx], label=file_idx)
            ax.plot(all_t_1[idx], all_err_q[file_idx][idx].mean() * np.ones(len(all_t_1[idx])), '--', label=file_idx)
        ax.set_xlabel('t [s]')
        if idx == 0:
            #ax.set_ylabel(r'$\frac{|q_{bi}-q_{one}|}{q_{bi}}$')
            ax.set_ylabel('relative error q')
        #ax.legend()
        ax.grid(True)
    plt.show()


if True:
    fig, axes = plt.subplots(2, len(args.vessels), squeeze=False, sharey='row', sharex='col')

    #taus = [1.5625e-05 * 2**(i+1) for i in range(len(args.filepaths_2))]
    #taus = [i for i in range(len(args.filepaths_2))]
    #taus = [i for i in range(9)] + [math.log2(0.005/1.5625e-5)]
    tau_min = 1. / 2**16
    #taus = [1./2**k for k in [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6]]
    taus = [1./2**k for k in [16, 15, 14, 13, 12, 11, 10, 9, 8]]
    axis_values = [math.log(tau/tau_min, 2) for tau in taus]

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[0, idx]
        means_p = [all_err_p[file_idx][idx].mean() for file_idx in range(len(args.filepaths_2))]
        print(means_p)
        ax.plot(axis_values, means_p, '-')
        if idx == 0:
            #ax.set_ylabel(r'$\frac{|p_{bi}-p_{one}|}{p_{bi}}$')
            ax.set_ylabel('mean relative error p')
        #ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[1, idx]
        means_q = [all_err_q[file_idx][idx].mean() for file_idx in range(len(args.filepaths_2))]
        ax.plot(axis_values, means_q, '-')
        ax.set_xlabel(r'$\log_2(\frac{\tau}{\tau_{min}})$')
        if idx == 0:
            #ax.set_ylabel(r'$\frac{|q_{bi}-q_{one}|}{q_{bi}}$')
            ax.set_ylabel('mean relative error q')
        ax.grid(True)
    plt.show()
