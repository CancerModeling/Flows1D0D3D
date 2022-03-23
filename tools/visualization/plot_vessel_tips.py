import os
import json
import numpy as np
import argparse
from matplotlib import pyplot as plt


def find_vessel_by_edge_id(metainfo, edge_id):
    for v in metainfo['vertices']:
        if v['neighbor_edge_id'] == edge_id:
            return v
    return None


def reformat_data(d):
    if len(d.shape) == 1:
        d = d.reshape((len(d), 1))
    return d


def load_data_by_edge_id(directory_name, metainfo, edge_id):
    v = find_vessel_by_edge_id(metainfo, edge_id)
    return load_data_by_vessel(directory_name, v)


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
    plt.style.use('seaborn-pastel')

    parser = argparse.ArgumentParser(description='Animator for the vessel data.')
    parser.add_argument('--filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
    parser.add_argument('--vessel-by-edge-id', type=int, help='The edge id of the vessel tip to plot.', required=True)
    parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
    parser.add_argument('--t-stop', type=float, help='Start point when to plot', required=False)
    parser.add_argument('--show-p', help='plot the pressure p', action='store_true')
    parser.add_argument('--show-q', help='plot the flow q', action='store_true')
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

    show_Q_total = False
    show_N_total = False 
    show_p = args.show_p 
    show_Q = args.show_q 
    show_V = False 
    show_c = False 
    show_N = False 
    no_legend = True 

    show_V_from_p = False 

    data, vessel = load_data_by_edge_id(directory_name, metainfo, args.vessel_by_edge_id)
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
                if not no_legend:
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

                if dof+1 < data['p'].shape[1]:
                    pressures_next = data['p'][start_index:stop_index, dof+1]
                elif 'p_out' in data:
                    pressures_next = data['p_out'][start_index:stop_index].squeeze()
                else:
                    pressures_next = vessel['p_out'] * np.ones(len(data['p'][start_index:stop_index, dof]))

                if 'resistance' in vessel:
                    resistances = np.array(vessel['resistance'])
                else:
                    print('r2')
                    resistances = np.array([vessel['R2']])
                flows = (pressures_here - pressures_next) / resistances[dof]

            if show_Q:
                ax = axes[next_ax_index, i]; next_ax_index += 1
                ax.plot(t[start_index:stop_index], flows, label='{}'.format(indices[dof]), linewidth=3)
                if not no_legend:
                    ax.legend()
                if i == 0:
                    ax.set_ylabel('$q_c [cm^3 s^{-1}]$')
                if next_ax_index == num_plot_rows:
                    ax.set_xlabel('$t$')
                ax.grid(True)

            if 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_Q_total:
                ax = axes[next_ax_index, i]; next_ax_index += 1
                ax.plot(t[start_index:stop_index], 2**(dof+1) * flows, label='{}'.format(indices[dof]), linewidth=3)
                if not no_legend:
                    ax.legend()
                if i == 0:
                    ax.set_ylabel('$q_{c,total} [cm^3 s^{-1}]$')
                if next_ax_index == num_plot_rows:
                    ax.set_xlabel('$t$')
                ax.grid(True)

        if 'c' in data and show_c:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], data['c'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
            if not no_legend:
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
            if not no_legend:
                ax.legend()
            if i == 0:
                ax.set_ylabel('$V [cm^3]$')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

        if 'c' in data and 'V' in data and show_N:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], data['c'][start_index:stop_index, dof] * data['V'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
            if not no_legend:
                ax.legend()
            if i == 0:
                ax.set_ylabel('$N [mmol]$')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

        if 'c' in data and 'V' in data and 'furcation_number' in vessel and vessel['furcation_number'] > 1 and show_N_total:
            ax = axes[next_ax_index, i]; next_ax_index += 1
            ax.plot(t[start_index:stop_index], 2**(dof+1) * data['c'][start_index:stop_index, dof] * data['V'][start_index:stop_index, dof], label='{}'.format(indices[dof]), linewidth=3)
            if not no_legend:
                ax.legend()
            if i == 0:
                ax.set_ylabel('$N_{total} [mmol]$')
            if next_ax_index == num_plot_rows:
                ax.set_xlabel('$t$')
            ax.grid(True)

    plt.show()
