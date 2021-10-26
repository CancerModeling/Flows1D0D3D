import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Script for extracting pressures values from the input data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', required=True)
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--position', type=float, default=0.5)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--name-postfix', type=str, default='_in')
parser.add_argument('--output-filepath', type=str, required=True)
parser.add_argument('--periodic', help='adds a flag for periodicity', action='store_true')

args = parser.parse_args()

directory = os.path.dirname(args.filepath)

with open(args.filepath) as f:
    meta = json.loads(f.read())

def find_vessel(vessel_id):
    for vessel in meta['vessels']:
        if vessel['edge_id'] == vessel_id:
            return vessel
    raise 'vessel with id {} not found'.format(vessel_id)

data_list = []

for vessel_id in args.vessels:
    vessel_info = find_vessel(vessel_id)

    path = os.path.join(directory, meta['filepath_time'])
    print('loading {}'.format(path))
    t = np.loadtxt(path, delimiter=',')

    if ('a' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['a'])
        print('loading {}'.format(path))
        a = np.loadtxt(path, delimiter=',')
        a = a[:]
        a = a[:, int((a.shape[1]-1) * args.position)]
    else:
        a = np.ones(t.shape) * vessel_info['A0']

    if ('p' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['p'])
        print('loading {}'.format(path))
        p = np.loadtxt(path, delimiter=',')
        p = p[:]
        p = p[:, int((p.shape[1]-1) * args.position)]
    else:
        p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1)

    start_index = np.sum(t <= args.t_start-1e-3)
    end_index = np.sum(t <= args.t_end+1e-3)

    p = p[start_index:end_index].tolist()
    t = t[start_index:end_index].tolist()

    data = {}
    data['name'] = vessel_info['name'] + args.name_postfix
    data['p'] = p
    data['t'] = t
    data['periodic'] = args.periodic
    data_list.append(data)

serialized_data = json.dumps(data_list, indent=0)
print(serialized_data)

with open(args.output_filepath, 'w') as f:
    f.write(serialized_data)
    print('wrote to {}'.format(args.output_filepath))
