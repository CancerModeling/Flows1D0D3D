import json
import numpy as np
import argparse
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')


def load_data(vessel_id):
    filename_data = 'output/data_{}_vessel{:05d}.csv'
    filename_times = 'output/data_times.csv'

    data = {}
    with open(filename_times, 'r') as file:
        numeric_data = np.loadtxt(file)
        data['t'] = numeric_data
    for component in ['Q', 'A', 'p', 'c']:
        with open(filename_data.format(component, vessel_id), 'r') as file:
            numeric_data = np.loadtxt(file, delimiter=',')
            data['grid'] = numeric_data[0,:]
            data[component] = numeric_data[1:,:]
    data['c'] = data['c'] / data['A']
    return data


parser = argparse.ArgumentParser(description='Animator for the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[0])

args = parser.parse_args()

data_sets = [load_data(index) for index in list(args.vessels) ]

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

anim = FuncAnimation(fig, animate, interval=2)

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
