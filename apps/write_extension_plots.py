import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse


filepaths = [
    'output_full/heart_to_breast_1d_solution_nl.json',
    #'',
    'output_33vessel_extended/abstract_33_vessels.json'
]

labels = [
    'cow+extension+breast',
    #'extension+breast',
    'cow+extension',
]

vessels = [
    [35, 36, 37, 38, 39, 40, 41, 42],
    #[0, 1, 2, 3, 4, 5, 6, 7],
    [35, 36, 37, 38, 39, 40, 41, 42]
]

visible = [
    [True, True, True, True, True, True, True, True],
    #[True, True, True, True, True, True, True, True],
    [True, False, False, False, False, True, False, False],
]

positions = [
    [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    [0, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 0.5],
    [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
]

t_start = 0
t_end = 1

output_folder = 'tmp_for_plots'
dataset_name = 'extension_'

list_all_t = []
list_all_p = []
list_all_q = []

for filepath_idx, filepath in enumerate(filepaths):
    directory = os.path.dirname(filepath)

    with open(filepath) as f:
        meta = json.loads(f.read())

    def find_vessel(vessel_id):
        for vessel in meta['vessels']:
            if vessel['edge_id'] == vessel_id:
                return vessel
        raise 'vessel with id {} not found'.format(vessel_id)

    list_t = []
    list_q = []
    list_p = []

    for vessel_idx, vessel_id in enumerate(vessels[filepath_idx]):
        vessel_info = find_vessel(vessel_id)

        path = os.path.join(directory, meta['filepath_time'])
        print('loading {}'.format(path))
        t = np.loadtxt(path, delimiter=',')

        path = os.path.join(directory, vessel_info['filepaths']['q'])
        print('loading {}'.format(path))
        q = np.loadtxt(path, delimiter=',')
        q = q[:]
        q = q[:, int((q.shape[1]-1) * positions[filepath_idx][vessel_idx])]

        if q.mean() < 0:
            q *= -1

        if ('a' in vessel_info['filepaths']):
            path = os.path.join(directory, vessel_info['filepaths']['a'])
            print('loading {}'.format(path))
            a = np.loadtxt(path, delimiter=',')
            a = a[:]
            a = a[:, int((a.shape[1]-1) * positions[filepath_idx][vessel_idx])]
        else:
            a = np.ones(q.shape) * vessel_info['A0']

        if ('p' in vessel_info['filepaths']):
            path = os.path.join(directory, vessel_info['filepaths']['p'])
            print('loading {}'.format(path))
            p = np.loadtxt(path, delimiter=',') / 1.333332
            print ('a', p)
            p = p[:]
            print ('b', p)
            p = p[:, int((p.shape[1]-1) * positions[filepath_idx][vessel_idx])]
            print ('c', p)
        else:
            p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1) / 1.33332

        start_index = np.sum(t < t_start)
        end_index = np.sum(t < t_end)

        list_p.append(p[start_index:end_index].tolist())
        list_q.append(q[start_index:end_index].tolist())
        list_t.append(t[start_index:end_index].tolist())

    list_all_p.append(list_p)
    list_all_q.append(list_q)
    list_all_t.append(list_t)

max_p, min_p = np.array(list_all_p).max(), np.array(list_all_p).min()
max_q, min_q = np.array(list_all_q).max(), np.array(list_all_q).min()

if True:
    fig = plt.figure()

    for idx in range(len(list_all_p[0])):
        plt.clf()
        ax = fig.add_subplot(111) 
        for src_idx in range(len(list_all_p)):
            if not visible[src_idx][idx]:
                continue
            vessel_id = vessels[src_idx][idx]
            ax.plot(list_all_t[src_idx][idx], list_all_p[src_idx][idx], label='$p_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=4-src_idx)
        ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel('p [mmHg]')
        ax.set_ylim(top=max_p+5, bottom=min_p-5)
        ax.grid(True)
        plt.savefig(os.path.join(output_folder, '{}_p_{}.pdf'.format(dataset_name, vessel_id)))

    for idx in range(len(list_all_p[0])):
        plt.clf()
        ax = fig.add_subplot(111) 
        for src_idx in range(len(list_all_q)):
            if not visible[src_idx][idx]:
                continue
            vessel_id = vessels[src_idx][idx]
            ax.plot(list_all_t[src_idx][idx], list_all_q[src_idx][idx], label='$q_{' + str(vessel_id+1) + '}$ ' + labels[src_idx], linewidth=4-src_idx)
        ax.legend()
        ax.set_xlabel('t [s]')
        ax.set_ylabel(r'q [$cm^{3}/s$]')
        #ax.set_ylim(top=max_q+5, bottom=min_q-5)
        ax.grid(True)
        plt.savefig(os.path.join(output_folder, '{}_q_{}.pdf'.format(dataset_name, idx)))
