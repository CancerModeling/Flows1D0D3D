import os
import numpy as np
import json
from matplotlib import pyplot as plt
import argparse
import sympy as sp


def get_fourier_series(y, N):
    res = np.fft.rfft(y) / y.size * 2
    a0 = res[0]/2
    a = res[1:N].real
    b = -res[1:N].imag
    a0 = 0 if np.abs(a0)<1e-14 else a0
    a[np.abs(a)<1e-14] = 0
    b[np.abs(b)<1e-14] = 0

    t = sp.symbols('t')
    c = [sp.cos(2*sp.pi*k*t) for k in range(1,N)]
    s = [sp.sin(2*sp.pi*k*t) for k in range(1,N)]

    val = a0
    for cc,ss,aa,bb in zip(c,s,a,b):
        val += aa*cc + bb*ss

    return val, sp.lambdify(t, val)


parser = argparse.ArgumentParser(description='Plots the vessel data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of ids of the vessels to plot.', default=[15])
parser.add_argument('--position', type=float, help='A list of ids of the vessels to plot.', default=0.5)
parser.add_argument('--t-start', type=float, default=0)
parser.add_argument('--filepath', type=str, required=True)
parser.add_argument('--no-legend', help='remove the legend from the plots', action='store_true')
parser.add_argument('-N', type=int, help='Number of fourier coefficients', default=4)

args = parser.parse_args()

directory = os.path.dirname(args.filepath)


with open(args.filepath) as f:
    meta = json.loads(f.read())


def find_vessel(vessel_id):
    for vessel in meta['vessels']:
        if vessel['edge_id'] == vessel_id:
            return vessel
    raise 'vessel with id {} not found'.format(vessel_id)


vessel_info = meta['vessels'][0]


position = args.position 


for idx, vessel_id in enumerate(args.vessels):
    vessel_info = find_vessel(vessel_id)

    path = os.path.join(directory, vessel_info['filepaths']['q'])
    print('loading {}'.format(path))
    q = np.loadtxt(path, delimiter=',')
    q = q[:]
    q = q[:, int((q.shape[1]-1)*position)]

    path = os.path.join(directory, meta['filepath_time'])
    print('loading {}'.format(path))
    t = np.loadtxt(path, delimiter=',')
    
    start_index = np.sum(t <= args.t_start-1e-12)
    end_index = np.sum(t < args.t_start+1)

    end_index = min(len(q), end_index)

    t = t[start_index:end_index]
    print(t)
    q = q[start_index:end_index]

    fourier_val, fourier_fun = get_fourier_series(q, args.N)
    print(sp.ccode(fourier_val))


    plt.plot(t,q)
    plt.plot(t,fourier_fun(t))
    plt.tight_layout()
    plt.show()

