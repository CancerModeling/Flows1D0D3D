from matplotlib import pyplot as plt
import numpy as np
import json

filepath = 'tmp/average_quantities.json'
indices = [1, 2, 3]
indices = range(31) 

with open(filepath, 'r') as file:
    data = json.loads(file.read())


fig = plt.figure()
axes = fig.subplots(2, 1, sharey='row', sharex='col', squeeze=False)

plot_pressures = False 


for idx in indices:
    data_at_idx = None
    for d in data['quantities'].values():
        if d['idx'] == idx:
            data_at_idx = d
    print(data_at_idx)
    if plot_pressures:
        p_cap = np.array(data_at_idx['p_cap']) / 1333
        p_tis = np.array(data_at_idx['p_tis']) / 1333
        axes[0][0].plot(data['t'], p_cap, label='{}'.format(idx))
        axes[1][0].plot(data['t'], p_tis, label='{}'.format(idx))
        axes[0][0].set_ylabel(r'$\bar p_{cap}$ [mmHg]')
        axes[1][0].set_ylabel(r'$\bar p_{tis}$ [mmHg]')
    else:
        axes[0][0].plot(data['t'], data_at_idx['c_cap'], label='{}'.format(idx))
        axes[1][0].plot(data['t'], data_at_idx['c_tis'], label='{}'.format(idx))
        axes[0][0].set_ylabel(r'$\bar c_{cap}$ [$mmol/cm^3$]')
        axes[1][0].set_ylabel(r'$\bar c_{tis}$ [$mmol/cm^3$]')
axes[1][0].set_xlabel('t [s]')
axes[0][0].grid(True)
axes[1][0].grid(True)
plt.legend()
plt.grid(True)
plt.show()
    