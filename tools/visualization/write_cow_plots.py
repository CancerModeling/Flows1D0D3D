import os
import numpy as np
import json
import matplotlib as mpl
from matplotlib import pyplot as plt
import argparse


font = {'family': 'normal', 'size': 22}
mpl.rc('font', **font)


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--output-folder', type=str, required=True)
parser.add_argument('--dataset-name', type=str, required=True)

args = parser.parse_args()

directory = os.path.dirname(args.filepath)

with open(args.filepath) as f:
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

    list_p.append(p[start_index:end_index])
    list_q.append(q[start_index:end_index])
    list_t.append(t[start_index:end_index])

list_p = np.array(list_p)
list_q = np.array(list_q)

max_p, min_p = list_p.max(), list_p.min()
max_q, min_q = list_q.max(), list_q.min()


if True:
    fig = plt.figure()

    for idx, vessel_id in enumerate(args.vessels):
        plt.clf()
        ax = fig.add_subplot(111) 
        ax.plot(list_t[idx], list_p[idx], label='$p_{' + str(vessel_id+1) + '}$', linewidth=4)
        ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel('p [mmHg]')
        ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)
        plt.savefig(os.path.join(args.output_folder, '{}_p_{}.pdf'.format(args.dataset_name, vessel_id)), bbox_inches='tight')

    for idx, vessel_id in enumerate(args.vessels):
        plt.clf()
        ax = fig.add_subplot(111) 
        ax.plot(list_t[idx], list_q[idx], label='$q_{' + str(vessel_id+1) + '}$', linewidth=4)
        ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel(r'q [$cm^{3}/s$]')
        #ax.set_ylim(top=max_q+5, bottom=min_q-5)
        ax.grid(True)
        plt.savefig(os.path.join(args.output_folder, '{}_q_{}.pdf'.format(args.dataset_name, vessel_id)), bbox_inches='tight')

if True:
    plt.clf() 
    fig, axes = plt.subplots(2, len(args.vessels), squeeze=False, sharey='row', sharex='col')

    row_idx = 0
    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[row_idx,idx] 
        ax.plot(list_t[idx], list_p[idx], label='$p_{' + str(vessel_id+1) + '}$')
        #ax.legend()
        #ax.set_xlabel('t [s]')
        if idx == 0:
            ax.set_ylabel('p [mmHg]')
        ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)

    row_idx += 1

    for idx, vessel_id in enumerate(args.vessels):
        ax = axes[row_idx,idx] 
        ax.plot(list_t[idx], list_q[idx], label='$q_{' + str(vessel_id+1) + '}$')
        #ax.legend()
        ax.set_xlabel('t [s]')
        if idx == 0:
            ax.set_ylabel(r'q [$cm^{3}/s$]')
        #ax.set_ylim(top=max_q+5, bottom=min_q-5)
        ax.grid(True)
    plt.show()
