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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots the vessel data.')
    parser.add_argument('--filepath', type=str, required=True)
    parser.add_argument('--no-legend', help='remove the legend from the plots', action='store_true')
    parser.add_argument('-N', type=int, help='Number of fourier coefficients', default=4)
    parser.add_argument('--vessel-id', type=int, required=True, help='Vessel id for Fourier')
    parser.add_argument('--quantity', help='Which quantity should be approximated?', type=str, choices=['q', 'v', 'a'], default='v')

    args = parser.parse_args()

    directory = os.path.dirname(args.filepath)

    with open(args.filepath) as f:
        meta = json.loads(f.read())

    try:
        vessel_info = meta['vessel_data'][f'{args.vessel_id}']
    except:
        print(f'Vessel {args.vessel_id} not found. Only {meta["vessel_ids"]} available.'); exit()


    try:
        q = np.array(vessel_info[args.quantity])
    except:
        print(f'Quantity "{args.quantity}" not found'); exit()
    t = np.array(meta['time'])

    fourier_val, fourier_fun = get_fourier_series(q, args.N)
    print(sp.ccode(fourier_val))

    plt.plot(t,q)
    plt.plot(t,fourier_fun(t))
    plt.tight_layout()
    plt.show()

