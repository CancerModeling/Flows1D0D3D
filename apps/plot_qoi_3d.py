import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.style.use('seaborn-pastel')

qoi = pd.read_csv('output_full_1d0d3d_pkj/qoi_3d.txt', delimiter=', ')
qoi = qoi.values
indices = list(range(qoi.shape[1]))

fig, axes = plt.subplots(2, 2)
labels = [r'$||p_{c}||$', r'$||p_{t}||$', \
          r'$||\phi_{\sigma c}||$', r'$||\phi_{\sigma t}||$']
norm_labels = [r'$_{\infty}$', r'$_{1}$', r'$_{2}$']
titles = [r'Capillary pressure', r'Tissue pressure', \
          r'Capillary nutrients', r'Tissue nutrients']
for i in range(4):
    I = int(i/2)
    J = i%2
    ax = axes[I, J]
    for k in range(3):
        ax.plot(qoi[:, i*3 + k], \
                     label=labels[i] + norm_labels[k], \
                     linewidth = 3)
    ax.set_title(titles[i])
    ax.legend()
    ax.grid(True)

plt.show()