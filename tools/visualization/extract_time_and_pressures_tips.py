import os
import json
import numpy as np
import argparse


def reformat_data(d):
    if len(d.shape) == 1:
        d = d.reshape((len(d), 1))
    return d


def load_data_by_vessel(directory_name, v):
    d_list = {}
    if 'filepath_p' in v['filepaths']:
        d_list['p'] = reformat_data(np.loadtxt(os.path.join(directory_name, v['filepaths']['filepath_p'])))
    if 'filepath_c' in v['filepaths']:
        d_list['c'] = reformat_data(np.loadtxt(os.path.join(directory_name, v['filepaths']['filepath_c'])))
    if 'filepath_V' in v['filepaths']:
        d_list['V'] = reformat_data(np.loadtxt(os.path.join(directory_name, v['filepaths']['filepath_V'])))
    if 'filepath_p_out' in v:
        d_list['p_out'] = reformat_data(np.loadtxt(os.path.join(directory_name, v['filepath_p_out'])))
    return d_list, v


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Animator for the vessel data.')
    parser.add_argument('--filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
    parser.add_argument('--output-filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
    parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
    parser.add_argument('--t-stop', type=float, help='Start point when to plot', required=False)
    parser.add_argument('--dofs', type=int, nargs='+',  help='A list of dofs to observer.', default=[-1])

    args = parser.parse_args()

    directory_name = os.path.dirname(args.filepath)

    with open(args.filepath) as f:
        metainfo = json.loads(f.read())


    t = np.array(metainfo['times'])

    start_index = sum(t < args.t_start)

    if args.t_stop is not None:
        stop_index = sum(t <= args.t_stop)
    else:
        stop_index = len(t)

    data_list = []

    for vertex in metainfo['vertices']:
        data, vertex = load_data_by_vessel(directory_name, vertex)
        p = data['p'][start_index:stop_index, -1] * 1e3 # kg -> g conversion factor
        p_mean = p.mean()
        R = vertex['resistance'][-1] * 1e3 # kg -> g
        level = len(vertex['resistance'])
        vertex_id = vertex['vertex_id']
        coordinates = vertex['coordinates']
        data_list.append({
            'p': coordinates,
            'vertex_id': vertex_id,
            'pressure': p_mean,
            'concentration': 1,
            'R2': R,
            'radius': vertex['radii'][-1],
            'level': level
        })

    serialized_data = json.dumps({'vertices': data_list}, indent=0)
    print(serialized_data)

    with open(args.output_filepath, 'w') as f:
        f.write(serialized_data)
        print('wrote to {}'.format(args.output_filepath))
