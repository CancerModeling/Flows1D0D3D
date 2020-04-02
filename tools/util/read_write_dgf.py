import os
import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
import importlib

import logging


def read_dgf(in_file):

    # read dgf file
    inpf = open(in_file)

    # vertex and associated data
    nodes = []
    nodes_data = []

    # elements and associated data
    segments = []
    segments_data = []

    block_counter = 0

    while True:
        a = inpf.readline().split()
        a

        if a is None or len(a) == 0 or a is EOFError:
            break
        elif len(a) >= 3 and block_counter == 3:

            # we are reading node data

            # coordinates
            node_cord = []
            for i in range(3):
                node_cord.append(float(a[i]))

            # number of nodal data
            num_node_data = len(a) - 3
            node_data = []

            for i in range(num_node_data):
                node_data.append(float(a[3 + i]))

            # add to list
            nodes.append(node_cord)
            nodes_data.append(node_data)

        elif len(a) >= 3 and block_counter > 3:

            # we are reading segment data

            # connectivity
            segment = []
            for i in range(2):
                id = float(a[i])
                segment.append(int(id))

            # segment data
            num_seg_data = len(a) - 2
            segment_data = []

            for i in range(num_seg_data):
                segment_data.append(float(a[2 + i]))

            # add to list
            segments.append(segment)
            segments_data.append(segment_data)

        else:
            # increment counter
            block_counter += 1

    inpf.close()

    nodes = np.array(nodes)
    nodes_data = np.array(nodes_data)

    segments = np.array(segments)
    segments_data = np.array(segments_data)

    return nodes, nodes_data, segments, segments_data


def write_dgf(out_file, nodes, nodes_data, segments, segments_data):

    inpf = open(out_file,'w')

    num_nodes = len(nodes)
    num_segments = len(segments)

    # header
    inpf.write("DGF\n")

    # write coordinates of nodes
    inpf.write("Vertex\n")
    inpf.write("parameters 1\n")
    for i in range(num_nodes):
        inpf.write("{} {} {} {}\n".format(nodes[i][0], nodes[i][1], nodes[i][2], nodes_data[i][0]))

    # write segment connectivity
    inpf.write("#\n")
    inpf.write("SIMPLEX\n")
    inpf.write("parameters 2\n")
    for i in range(num_segments):
        inpf.write("{} {} {} {}\n".format(segments[i][0], segments[i][1], segments_data[i][0], segments_data[i][1]))

    # write boundary domain
    inpf.write("#\n")
    inpf.write("BOUNDARYDOMAIN\n")
    inpf.write("default 1\n")
    inpf.write("#\n")

    inpf.close()
