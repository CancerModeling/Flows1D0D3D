import argparse
from matplotlib import pyplot as plt
import numpy as np
import json


parser = argparse.ArgumentParser(description='Plot the 3D averages')
parser.add_argument("--filepath", type=str, help="Filepath to the json file containing the average quantities.", required=True)
parser.add_argument('--t-start', type=float, help='Where should the plot start?', default=0)
parser.add_argument('--t-stop', type=float, help='Where should the plot be cut-off?', default=100)
args = parser.parse_args()

#filepath = 'tmp_transport/average_quantities.json'
#indices = [23, 60, 128, 145]
indices = [117, 23, 127, 145]
legend = ['1', '2', '3', '4']
#indices = range(31) 

with open(args.filepath, 'r') as file:
    data = json.loads(file.read())



plot_pressures = True 

t_start = args.t_start
t_stop = args.t_stop
#t_stop = 100


if True:
    fig = plt.figure()
    axes = fig.subplots(2, 1, sharey='row', sharex='col', squeeze=False)

    for vessel_idx, data_at_idx in data['quantities'].items():
        vessel_idx = int(vessel_idx)
        if vessel_idx not in indices:
            print ('not showing {}'.format(vessel_idx))
            #continue
        t = np.array(data['t'])
        start_index = np.sum(t < t_start)
        stop_index = np.sum(t < t_stop)
        t = t[start_index:stop_index]
        if plot_pressures:
            p_cap = np.array(data_at_idx['p_cap']) / 1333
            p_tis = np.array(data_at_idx['p_tis']) / 1333
            p_cap = p_cap[start_index:stop_index]
            p_tis = p_tis[start_index:stop_index]
            axes[0][0].plot(t, p_cap, label='{}'.format(vessel_idx))
            axes[1][0].plot(t, p_tis, label='{}'.format(vessel_idx))
            axes[0][0].set_ylabel(r'$\bar p_{cap}$ [mmHg]')
            axes[1][0].set_ylabel(r'$\bar p_{tis}$ [mmHg]')
        else:
            c_cap = data_at_idx['c_cap']
            c_tis = data_at_idx['c_tis']
            c_cap = c_cap[start_index:stop_index]
            c_tis = c_tis[start_index:stop_index]
            axes[0][0].plot(t, c_cap, label='{}'.format(vessel_idx))
            axes[1][0].plot(t, c_tis, label='{}'.format(vessel_idx))
            axes[0][0].set_ylabel(r'$\bar c_{cap}$ [$mmol/cm^3$]')
            axes[1][0].set_ylabel(r'$\bar c_{tis}$ [$mmol/cm^3$]')
    axes[1][0].set_xlabel('t [s]')
    axes[0][0].grid(True)
    axes[1][0].grid(True)
    plt.legend()
    plt.grid(True)
    plt.show()

if True:
    fig = plt.figure()
    axes = fig.subplots(2, 1, sharey='row', sharex='col', squeeze=False)

    for idx, vessel_idx in enumerate(indices):
        vessel_idx = int(vessel_idx)
        data_at_idx = data['quantities']['{}'.format(vessel_idx)]
        # if vessel_idx not in indices:
        #     print ('not showing {}'.format(vessel_idx))
        #     continue
        t = np.array(data['t'])
        start_index = np.sum(t < t_start)
        stop_index = np.sum(t < t_stop)
        t = t[start_index:stop_index]
        if plot_pressures:
            p_cap = np.array(data_at_idx['p_cap']) / 1333
            p_tis = np.array(data_at_idx['p_tis']) / 1333
            p_cap = p_cap[start_index:stop_index]
            p_tis = p_tis[start_index:stop_index]
            axes[0][0].plot(t, p_cap, label='{}'.format(vessel_idx))
            axes[1][0].plot(t, p_tis, label='{}'.format(vessel_idx))
            axes[0][0].set_ylabel(r'$\bar p_{cap}$ [mmHg]')
            axes[1][0].set_ylabel(r'$\bar p_{tis}$ [mmHg]')
        else:
            c_cap = data_at_idx['c_cap']
            c_tis = data_at_idx['c_tis']
            c_cap = c_cap[start_index:stop_index]
            c_tis = c_tis[start_index:stop_index]
            print(indices.index(vessel_idx), vessel_idx)
            print(c_cap)
            print(c_tis)
            axes[0][0].plot(t, c_cap, label=legend[indices.index(vessel_idx)])
            axes[1][0].plot(t, c_tis, label=legend[indices.index(vessel_idx)])
            axes[0][0].set_ylabel(r'$\bar c_{cap}$ [mmol/cm${}^3$]')
            axes[1][0].set_ylabel(r'$\bar c_{tis}$ [mmol/cm${}^3$]')
    axes[1][0].set_xlabel('t [s]')
    axes[0][0].grid(True)
    axes[1][0].grid(True)
    plt.legend()
    plt.grid(True)
    plt.show()
    
