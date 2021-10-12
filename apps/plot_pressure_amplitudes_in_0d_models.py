import numpy as np
import json
import os
import plot_vessel_tips as pvt
from matplotlib import pyplot as plt
import argparse

plt.style.use('seaborn-pastel')

parser = argparse.ArgumentParser(description='Animator for the vessel data.')
parser.add_argument('--filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
parser.add_argument('--t-stop', type=float, help='Start point when to plot', required=False)

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

radii_list = []
mean_pressures_list = []
min_pressures_list = []
max_pressures_list = []

for vessel in metainfo['vertices']:
    print(vessel)
    data,_ = pvt.load_data_by_vessel(directory_name, vessel)
    pressures = data['p'][start_index:stop_index,:] / 1.3333
    print(pressures.shape)
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

#plt.semilogx(radii_list, mean_pressures_list, 'o')
#plt.semilogx(radii_list, min_pressures_list)
#plt.semilogx(radii_list, max_pressures_list)
#plt.show()

show_data_points = True
show_amplitude = True
show_min_max = True

if show_data_points:
    plt.errorbar(radii_list, mean_pressures_list, yerr=(mean_pressures_list-min_pressures_list, max_pressures_list-mean_pressures_list), fmt='.')
    plt.gca().set_xscale('log', nonposx='clip')
    plt.gca().set_xlabel('r [cm]')
    plt.gca().set_ylabel('p [mmHg]')
    plt.gca().invert_xaxis()
    plt.gca().grid(True)
    plt.show()

if show_amplitude:
    plt.plot(radii_list, (max_pressures_list-min_pressures_list)/2, '.')
    plt.gca().set_xscale('log', nonposx='clip')
    plt.gca().set_xlabel('r [cm]')
    plt.gca().set_ylabel('pressure amplitude [mmHg]')
    plt.gca().invert_xaxis()
    plt.gca().grid(True)
    plt.show()

if show_min_max:
    plt.plot(radii_list, max_pressures_list, '.', label='maximum')
    plt.plot(radii_list, mean_pressures_list, 'x', label='mean')
    plt.plot(radii_list, min_pressures_list, '+', label='minimum')
    plt.gca().set_xscale('log', nonposx='clip')
    plt.gca().set_xlabel('r [cm]')
    plt.gca().set_ylabel('pressure [mmHg]')
    plt.gca().legend()
    plt.gca().invert_xaxis()
    plt.gca().grid(True)
    plt.show()


#print(metainfo)

'''
data, vessel = load_data(directory_name, metainfo, args.vessel_by_edge_id)
print(vessel)

num_plot_rows = 0
if 'p' in data and show_p:
    num_plot_rows += 1
if 'p' in data and show_Q:
    num_plot_rows += 1
if 'p' in data and 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_Q_total:
    num_plot_rows += 1
if 'V' in data and show_V:
    num_plot_rows += 1
if 'c' in data and show_c:
    num_plot_rows += 1
if 'c' in data and 'V' in data and show_N:
    num_plot_rows += 1
if 'c' in data and 'V' in data and 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_N_total:
    num_plot_rows += 1

fig, axes = plt.subplots(num_plot_rows, len(args.dofs), squeeze=False, sharey='row', sharex='col')

arbitrary_data_element = next(iter(data.values()))
indices = list(range(arbitrary_data_element.shape[1]))

num_dofs = vessel['num_dofs']
if vessel['outflow_type'] == 'rcl':
    num_p_dof = int(num_dofs / 2)
    num_q_dof = int(num_dofs / 2)
else:
    num_p_dof = num_dofs
    num_q_dof = 0

for i,dof in enumerate(args.dofs):
    next_ax_index = 0

    if 'p' in data:
        if show_p:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], data['p'][start_index:stop_index,dof] / 1.3333, label='{}'.format(indices[dof]), linewidth=3)
            ax.legend()
            if i == 0:
                ax.set_ylabel('$p_c$ [mmHg]')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

        if vessel['outflow_type'] == 'rcl':
            flows = data['p'][start_index:stop_index,dof+num_p_dof]
        else:
            pressures_here = data['p'][start_index:stop_index, dof]
            print('shape', data['p'].shape)
            if dof+1 < data['p'].shape[1]:
                pressures_next = data['p'][start_index:stop_index, dof+1]
            else:
                print(dof)
                pressures_next = vessel['p_out'] * np.ones(len(data['p'][start_index:stop_index, dof]))
            if 'resistance' in vessel:
                resistances = np.array(vessel['resistance'])
            else:
                print('r2')
                resistances = np.array([vessel['R2']])
            flows = (pressures_here - pressures_next) / resistances[dof]
            print(resistances[dof])

        if show_Q:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], flows, label='{}'.format(indices[dof]), linewidth=3)
            ax.legend()
            if i == 0:
                ax.set_ylabel('$q_c [cm^3 s^{-1}]$')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

        if 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_Q_total:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], 2**dof * flows, label='{}'.format(indices[dof]), linewidth=3)
            ax.legend()
            if i == 0:
                ax.set_ylabel('$q_{c,total} [cm^3 s^{-1}]$')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

    if 'c' in data and show_c:
        ax = axes[next_ax_index, i]; next_ax_index += 1
        ax.plot(t[start_index:stop_index], data['c'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$c [mmol cm^{-3}]$')
        if next_ax_index == num_plot_rows:
            ax.set_xlabel('$t$')
        ax.grid(True)

    if 'V' in data and show_V:
        ax = axes[next_ax_index, i]; next_ax_index += 1
        ax.plot(t[start_index:stop_index], data['V'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
        if show_V_from_p:
            ax.plot(t[start_index:stop_index], vessel['capacitances'][dof]*data['p'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=1)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$V [cm^3]$')
        if next_ax_index == num_plot_rows:
            ax.set_xlabel('$t$')
        ax.grid(True)

    if 'c' in data and 'V' in data and show_N:
        ax = axes[next_ax_index, i]; next_ax_index += 1
        ax.plot(t[start_index:stop_index], data['c'][start_index:stop_index, dof] * data['V'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$N [mmol]$')
        if next_ax_index == num_plot_rows:
            ax.set_xlabel('$t$')
        ax.grid(True)

    if 'c' in data and 'V' in data and 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_N_total:
        ax = axes[next_ax_index, i]; next_ax_index += 1
        ax.plot(t[start_index:stop_index], 2**(dof) * data['c'][start_index:stop_index, dof] * data['V'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
        ax.legend()
        if i == 0:
            ax.set_ylabel('$N_{total} [mmol]$')
        if next_ax_index == num_plot_rows:
            ax.set_xlabel('$t$')
        ax.grid(True)

plt.show()
'''
