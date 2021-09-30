import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.style.use('seaborn-pastel')

qoi = pd.read_csv('output_full_1d0d3d_pkj/qoi_3d.txt', delimiter=', ')
qoi = qoi.values
indices = list(range(qoi.shape[1]))

fig, axes = plt.subplots(3, 2)
labels = [r'$||p_{c}||$', r'$||p_{t}||$', \
          r'$||\phi_{\sigma c}||$', r'$||\phi_{\sigma t}||$', \
          r'$||\phi_T||$']
norm_labels = [r'$_{\infty}$', r'$_{1}$', r'$_{2}$']
titles = [r'Capillary pressure', r'Tissue pressure', \
          r'Capillary nutrients', r'Tissue nutrients', \
          r'Tumor']

sindex = 1; # start from second element as first element is time

for i in range(5):
    I = int(i/2)
    J = i%2
    ax = axes[I, J]
    for k in range(3):
        ax.plot(qoi[:, sindex + i*3 + k], \
                     label=labels[i] + norm_labels[k], \
                     linewidth = 3)
    ax.set_title(titles[i])
    ax.legend()
    ax.grid(True)

plt.show()
