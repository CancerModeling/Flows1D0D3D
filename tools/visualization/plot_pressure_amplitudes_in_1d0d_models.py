import numpy as np
import json
import os
import matplotlib as mpl 
from matplotlib import pyplot as plt
import argparse

import plot_vessel_tips as pvt


def load_1d_averaged_pressure_data(filepath, t_start, t_stop):
    directory_name = os.path.dirname(filepath)

    with open(filepath) as f:
        metainfo = json.loads(f.read())

    t = np.loadtxt(os.path.join(directory_name, metainfo['filepath_time']), delimiter=',')

    start_index = sum(t < t_start)

    if t_stop is not None:
        stop_index = sum(t <= t_stop)
    else:
        stop_index = len(t)

    radii_list = []
    mean_pressures_list = []
    min_pressures_list = []
    max_pressures_list = []

    for vessel in metainfo['vessels']:
        r = np.sqrt(vessel['A0']/np.pi)

        if 'p' in vessel['filepaths']:
            data = np.loadtxt(os.path.join(directory_name, vessel['filepaths']['p']), delimiter=',')
        elif 'a' in vessel['filepaths']:
            a = np.loadtxt(os.path.join(directory_name, vessel['filepaths']['a']), delimiter=',')
            data = vessel['G0'] * (np.sqrt(a / vessel['A0']) - 1)
        else:
            raise 'not implemented'

        num_time_points, num_spatial_points = data.shape

        pressures = data[start_index:stop_index,int(num_spatial_points/2)] / 1.3333
        pressure_min = pressures.min(axis=0)
        pressure_max = pressures.max(axis=0)
        pressure_mean = pressures.mean(axis=0)

        radii_list += [r]
        mean_pressures_list += [pressure_mean]
        max_pressures_list += [pressure_max]
        min_pressures_list += [pressure_min]

    radii_list = np.array(radii_list)
    permutation = np.argsort(radii_list)
    radii_list = radii_list[permutation]
    mean_pressures_list = np.array(mean_pressures_list)[permutation]
    max_pressures_list = np.array(max_pressures_list)[permutation]
    min_pressures_list = np.array(min_pressures_list)[permutation]

    return radii_list.tolist(), mean_pressures_list.tolist(), min_pressures_list.tolist(), max_pressures_list.tolist()


def load_0d_averaged_pressure_data(filepath, t_start, t_stop):
    directory_name = os.path.dirname(filepath)

    with open(filepath) as f:
        metainfo = json.loads(f.read())

    t = np.array(metainfo['times'])

    start_index = sum(t < t_start)

    if t_stop is not None:
        stop_index = sum(t <= t_stop)
    else:
        stop_index = len(t)

    radii_list = []
    mean_pressures_list = []
    min_pressures_list = []
    max_pressures_list = []

    for vessel in metainfo['vertices']:
        data,_ = pvt.load_data_by_vessel(directory_name, vessel)
        pressures = data['p'][start_index:stop_index,:] / 1.3333
        pressure_min = pressures.min(axis=0)
        pressure_max = pressures.max(axis=0)
        pressure_mean = pressures.mean(axis=0)
        radii = vessel['radii']

        radii_list += radii
        mean_pressures_list += pressure_mean.tolist()
        max_pressures_list += pressure_max.tolist()
        min_pressures_list += pressure_min.tolist()

    radii_list = np.array(radii_list)
    permutation = np.argsort(radii_list)
    radii_list = radii_list[permutation]
    mean_pressures_list = np.array(mean_pressures_list)[permutation]
    max_pressures_list = np.array(max_pressures_list)[permutation]
    min_pressures_list = np.array(min_pressures_list)[permutation]

    return radii_list.tolist(), mean_pressures_list.tolist(), min_pressures_list.tolist(), max_pressures_list.tolist()


def load_averaged_pressure_data(filepath, t_start, t_stop):
    with open(filepath) as f:
        metainfo = json.loads(f.read())

    if 'vertices' in metainfo:
        return load_0d_averaged_pressure_data(filepath, t_start, t_stop)
    elif 'vessels' in metainfo:
        return load_1d_averaged_pressure_data(filepath, t_start, t_stop)
    else:
        raise 'unknown'

color_max = ['#496380', '#6e95bf', '#92c6ff']
color_mean = ['#51815c', '#78bf87','#99f0aa'] 
color_min = ['#80504d', '#bf7774', '#ff9f9a']

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Animator for the vessel data.')
    parser.add_argument('--filepaths', type=str, help='Filepath to a file containing the pressures and flows.', nargs='+', required=True)
    parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
    parser.add_argument('--t-stop', type=float, help='Start point when to plot', required=False)

    args = parser.parse_args()

    show_data_points = True
    show_amplitude = True
    show_min_max = True

    if show_min_max:
        for idx, filepath in enumerate(args.filepaths):
            radii_list_, mean_p_list_, min_p_list_, max_p_list_ = load_averaged_pressure_data(filepath, args.t_start, args.t_stop)

            radii_list = np.array(radii_list_)
            permutation = np.argsort(radii_list_)
            radii_list = radii_list[permutation]
            mean_p_list = np.array(mean_p_list_)[permutation]
            max_p_list = np.array(max_p_list_)[permutation]
            min_p_list = np.array(min_p_list_)[permutation]

            label_max = 'maximum' if idx+1 == len(args.filepaths) else None
            label_mean = 'mean' if idx+1 == len(args.filepaths) else None
            label_min = 'minimum' if idx+1 == len(args.filepaths) else None

            plt.plot(radii_list_, max_p_list_, '.', label=label_max, color=color_max[idx])
            plt.plot(radii_list, mean_p_list, 'x', label=label_mean, color=color_mean[idx])
            plt.plot(radii_list, min_p_list, '+', label=label_min, color=color_min[idx])

        plt.gca().set_xscale('log', nonposx='clip')
        plt.gca().set_xlabel('r [cm]')
        plt.gca().set_ylabel('pressure [mmHg]')
        plt.gca().legend()
        plt.gca().invert_xaxis()
        plt.gca().grid(True)
        plt.show()

    if show_amplitude:
        for idx, filepath in enumerate(args.filepaths):
            radii_list_, mean_p_list_, min_p_list_, max_p_list_ = load_averaged_pressure_data(filepath, args.t_start, args.t_stop)

            radii_list = np.array(radii_list_)
            permutation = np.argsort(radii_list_)
            radii_list = radii_list[permutation]
            mean_p_list = np.array(mean_p_list_)[permutation]
            max_p_list = np.array(max_p_list_)[permutation]
            min_p_list = np.array(min_p_list_)[permutation]

            plt.plot(radii_list, (max_p_list-min_p_list), '.', color=color_max[idx])

        plt.gca().set_xscale('log', nonposx='clip')
        plt.gca().set_xlabel('r [cm]')
        plt.gca().set_ylabel('pressure amplitude [mmHg]')
        plt.gca().invert_xaxis()
        plt.gca().grid(True)
        plt.show()

    '''
        radii_list += radii_list_
        mean_p_list += mean_p_list_
        min_p_list += min_p_list_
        max_p_list += max_p_list_
        '''

    '''
    print('sorting')

    radii_list = np.array(radii_list)
    permutation = np.argsort(radii_list)
    radii_list = radii_list[permutation]
    mean_p_list = np.array(mean_p_list)[permutation]
    max_p_list = np.array(max_p_list)[permutation]
    min_p_list = np.array(min_p_list)[permutation]

    if show_data_points:
        plt.errorbar(radii_list, mean_p_list, yerr=(mean_p_list-min_p_list, max_p_list-mean_p_list), fmt='.')
        plt.gca().set_xscale('log', nonposx='clip')
        plt.gca().set_xlabel('r [cm]')
        plt.gca().set_ylabel('p [mmHg]')
        plt.gca().invert_xaxis()
        plt.gca().grid(True)
        plt.show()

    if show_amplitude:
        plt.plot(radii_list, (max_p_list-min_p_list)/2, '.')
        plt.gca().set_xscale('log', nonposx='clip')
        plt.gca().set_xlabel('r [cm]')
        plt.gca().set_ylabel('pressure amplitude [mmHg]')
        plt.gca().invert_xaxis()
        plt.gca().grid(True)
        plt.show()

    if show_min_max:
        plt.plot(radii_list, max_p_list, '.', label='maximum')
        plt.plot(radii_list, mean_p_list, 'x', label='mean')
        plt.plot(radii_list, min_p_list, '+', label='minimum')
        plt.gca().set_xscale('log', nonposx='clip')
        plt.gca().set_xlabel('r [cm]')
        plt.gca().set_ylabel('pressure [mmHg]')
        plt.gca().legend()
        plt.gca().invert_xaxis()
        plt.gca().grid(True)
        plt.show()
   '''
