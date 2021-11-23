import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Script for extracting pressures values from the input data.')
parser.add_argument('--vertices', type=str, nargs='+', help='A list of vertex names for the vessels to check.', required=True)
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--t-end', type=float, default=10000)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--names', type=str, nargs='+', help='an explicit list of names', default=[])
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


def find_vessels_with_vertex(vertex_name):
    vessels = {}
    for vessel in meta['vessels']:
        if vessel["vertices"]["left"]["name"] == vertex_name or vessel["vertices"]["right"]["name"] == vertex_name:
            vessels[vessel["edge_id"]] = vessel
    return list(vessels.values())


def find_vessel_at_vertex(vertex_name):
    vessels = find_vessels_with_vertex(vertex_name)
    if len(vessels) != 1:
        raise RuntimeError("vertex {} has {} edge neighbors".format(vertex_name, len(vessels)))
    return vessels[0]


def get_position(vessel, vertex_name):
    if vessel["vertices"]["left"]["name"] == vertex_name:
        return 0.
    elif vessel["vertices"]["right"]["name"] == vertex_name:
        return 1.
    else:
        raise RuntimeError("vessel and vertex name ({}) do not match".format(vertex_name))


data_list = []

for idx, vertex_name in enumerate(args.vertices):
    vessel_info = find_vessel_at_vertex(vertex_name)
    position = get_position(vessel_info, vertex_name)

    path = os.path.join(directory, meta['filepath_time'])
    print('loading {}'.format(path))
    t = np.loadtxt(path, delimiter=',')

    if ('a' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['a'])
        print('loading {}'.format(path))
        a = np.loadtxt(path, delimiter=',')
        a = a[:]
        a = a[:, int((a.shape[1]-1) * position)]
    else:
        a = np.ones(t.shape) * vessel_info['A0']

    if ('p' in vessel_info['filepaths']):
        path = os.path.join(directory, vessel_info['filepaths']['p'])
        print('loading {}'.format(path))
        p = np.loadtxt(path, delimiter=',')
        p = p[:]
        p = p[:, int((p.shape[1]-1) * position)]
    else:
        p = vessel_info['G0'] * (np.sqrt(a/vessel_info['A0']) - 1)

    start_index = np.sum(t <= args.t_start-1e-3)
    end_index = np.sum(t <= args.t_end+1e-3)

    p = p[start_index:end_index].tolist()
    t = t[start_index:end_index].tolist()

    data = {}
    if len(args.names) == 0:
        data['name'] = vertex_name
    else:
        data['name'] = args.names[idx]
    data['p'] = p
    data['t'] = t
    data['periodic'] = args.periodic
    data_list.append(data)

serialized_data = json.dumps({'vertices': data_list}, indent=0)
print(serialized_data)

with open(args.output_filepath, 'w') as f:
    f.write(serialized_data)
    print('wrote to {}'.format(args.output_filepath))
