import json
import numpy as np
import argparse
import os
from matplotlib import pyplot as plt
plt.style.use('seaborn-pastel')

parser = argparse.ArgumentParser(description='Animator for the vessel data.')
parser.add_argument('--filepath', type=str, help='Filepath to a file containing the pressures and flows.', required=True)
parser.add_argument('--vessel-by-edge-id', type=int, help='The edge id of the vessel tip to plot.', required=True)
parser.add_argument('--t-start', type=float, help='Start point when to plot', default=0.)
parser.add_argument('--dofs', type=int, nargs='+',  help='A list of dofs to observer.', default=[-1])

args = parser.parse_args()

with open(args.filepath) as f:
    metainfo = json.loads(f.read())

t = np.array(metainfo['times'])

start_index = sum(t < args.t_start)


def find_vessel_by_edge_id(edge_id):
    for v in metainfo['vertices']:
        if v['neighbor_edge_id'] == edge_id:
            return v
    return None


def load_data(edge_id):
    v = find_vessel_by_edge_id(edge_id)
    print(json.dumps(v, indent=1))
    return np.loadtxt(v['filepath'])


data = load_data(args.vessel_by_edge_id)

fig, axes = plt.subplots(1, len(args.dofs), squeeze=False)

indices = list(range(data.shape[1]))

for i,dof in enumerate(args.dofs):
    ax = axes[0,i]
    ax.plot(t[start_index:], data[start_index:,dof], label='{}'.format(indices[dof]))
    ax.legend()
    ax.set_ylabel('$p_c$')
    ax.set_xlabel('$t$')
    ax.grid(True)
plt.tight_layout()
plt.show()

print(load_data(args.vessel_by_edge_id))

#print(metainfo)


'''

args = parser.parse_args()

data_sets = [load_data(args.output_directory, index, args.t_start) for index in list(args.vessels) ]

fig = plt.figure()
fig.tight_layout()
axes = fig.subplots(3, len(data_sets), sharey='row', sharex='col')

lines = []
t_index = 0
for dset_index in range(len(data_sets)):
    if len(data_sets) < 2:
        ax_Q = axes[0]
        ax_A = axes[1]
        ax_c = axes[2]
    else:
        ax_Q = axes[0, dset_index]
        ax_A = axes[1, dset_index]
        ax_c = axes[2, dset_index]
    data_Q = data_sets[dset_index]['Q'][t_index]
    data_A = data_sets[dset_index]['A'][t_index]
    data_c = data_sets[dset_index]['c'][t_index]
    print(data_c)
    grid = data_sets[dset_index]['grid']
    ax_Q.clear()
    l1, = ax_Q.plot(grid, data_Q, label='Q')
    ax_Q.legend()
    ax_Q.grid(True)
    ax_A.clear()
    l2, = ax_A.plot(grid, data_A, label='A')
    ax_A.legend()
    ax_A.grid(True)
    ax_c.clear()
    l3, = ax_c.plot(grid, data_c, label='c/a')
    ax_c.legend()
    ax_c.grid(True)
    lines.append([l1, l2, l3])
fig.suptitle('t={}'.format(data_sets[0]['t'][t_index]))

for ax, vid in zip(axes if len(data_sets) < 2 else axes[0], args.vessels):
    ax.set_title('vessel {}'.format(vid))

def animate(i):
    t_index = i % len(data_sets[0]['t'])
    for dset_index in range(len(data_sets)):
        if len(data_sets) < 2:
            ax_Q = axes[0]
            ax_A = axes[1]
            ax_c = axes[2]
        else:
            ax_Q = axes[0, dset_index]
            ax_A = axes[1, dset_index]
            ax_c = axes[2, dset_index]
        data_Q = data_sets[dset_index]['Q'][t_index]
        data_A = data_sets[dset_index]['A'][t_index]
        data_c = data_sets[dset_index]['c'][t_index]
        lQ = lines[dset_index][0]
        lA = lines[dset_index][1]
        lc = lines[dset_index][2]
        grid = data_sets[dset_index]['grid']
        lQ.set_ydata(data_Q)
        lA.set_ydata(data_A)
        lc.set_ydata(data_c)
        ax_Q.relim()
        ax_A.relim()
        ax_c.set_ylim([0, 1])
        ax_Q.autoscale_view()
        ax_A.autoscale_view()
        ax_c.autoscale_view()
    fig.suptitle('t={}'.format(data_sets[0]['t'][t_index]))

num_frames = len(data_sets[0]['t'])

anim = FuncAnimation(fig, animate, interval=20)

is_running = True
def toggle_animation(event):
    global is_running
    if is_running:
        anim.event_source.stop()
    else:
        anim.event_source.start()
    is_running = not is_running

fig.canvas.mpl_connect('button_press_event', toggle_animation)

plt.show()


anim.save('sine_wave.gif', writer='imagemagick')
'''