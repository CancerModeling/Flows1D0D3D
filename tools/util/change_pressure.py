import numpy as np

from read_write_dgf import *
from read_write_vtk import write_vtk

def change_pressure(in_file, out_file, low_pres_old, high_pres_old, low_pres, high_pres):

    # read dgf file
    nodes, nodes_data, segments, segments_data = read_dgf(in_file)

    num_nodes = len(nodes)
    num_segments = len(segments)

    # in nodes_data, first element is pressure
    # search for old pressure and replace by new pressure
    for i in range(num_nodes):
        p = nodes_data[i][0]
        if abs(p - low_pres_old) < 1.e-5:
            p = low_pres
        elif abs(p - high_pres_old) < 1.e-5:
            p = high_pres

        nodes_data[i][0] = p

    # write to new dgf file
    write_dgf(out_file, nodes, nodes_data, segments, segments_data)

# low_pres_old = 20.
# high_pres_old = 25
# low_pres = 20000.
# high_pres = 50000.
# change_pressure("test_scaled_2.dgf", "diff_pres.dgf", low_pres_old, high_pres_old, low_pres, high_pres)

# run as follows
# python3 -B change_pressure.py
