from matplotlib import pyplot as plt
import numpy as np
import json

filepath = 'tmp/average_quantities.json'
indices = [11, 145]
#indices = range(31) 

with open(filepath, 'r') as file:
    data = json.loads(file.read())


fig = plt.figure()
axes = fig.subplots(2, 1, sharey='row', sharex='col', squeeze=False)

plot_pressures = False 

t_start = 0 
t_stop = 40 


for vessel_idx, data_at_idx in data['quantities'].items():
    vessel_idx = int(vessel_idx)
    if vessel_idx not in indices:
        print ('not showing {}'.format(vessel_idx))
        continue
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
    