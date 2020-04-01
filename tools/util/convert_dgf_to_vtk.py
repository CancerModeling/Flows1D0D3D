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


def write_vtk(out_file, nodes, nodes_data, segments, segments_data):

    inpf = open(out_file,'w')

    # header
    inpf.write("# vtk DataFile Version 2.0\n")
    inpf.write("Network ETH\n")
    inpf.write("ASCII\n")
    inpf.write("DATASET POLYDATA\n")

    # write coordinates of nodes
    num_nodes = len(nodes)
    inpf.write("POINTS {} float\n".format(num_nodes))
    for i in range(num_nodes):
        inpf.write("{} {} {}\n".format(nodes[i][0], nodes[i][1], nodes[i][2]))

    # write segment connectivity
    num_segments = len(segments)
    inpf.write("\n")
    inpf.write("LINES {} {}\n".format(num_segments, 3 * num_segments))

    for i in range(num_segments):
        inpf.write("{} {} {}\n".format(2, segments[i][0], segments[i][1]))

    # write cell data
    num_seg_data = 0
    if len(segments_data) > 0:
        num_seg_data = len(segments_data[0])
        if num_seg_data > 2:
            num_seg_data = 2

    for j in range(num_seg_data):

        data_name = "radii"
        if j == 1:
            data_name = "segment_data_" + str(j)

        # header
        if j == 0:
            inpf.write("\n")
            inpf.write("CELL_DATA {}\n".format(num_segments))

        inpf.write("SCALARS {} float 1\n".format(data_name))
        inpf.write("LOOKUP_TABLE default\n")

        for i in range(num_segments):
            inpf.write("{}\n".format(segments_data[i][j]))

    # write node data
    num_node_data = 0
    if len(nodes_data) > 0:
        num_node_data = len(nodes_data[0])
        if num_node_data > 2:
            num_node_data = 2

    for j in range(num_node_data):

        data_name = "pressure"
        if j == 1:
            data_name = "node_data_" + str(j)

        # header
        if j == 0:
            inpf.write("\n")
            inpf.write("POINT_DATA {}\n".format(num_nodes))

        inpf.write("SCALARS {} float 1\n".format(data_name))
        inpf.write("LOOKUP_TABLE default\n")

        for i in range(num_nodes):
            inpf.write("{}\n".format(nodes_data[i][j]))


    inpf.close()



def convert_dgf_to_vtk(in_file, out_file):

    # read dgf file
    nodes, nodes_data, segments, segments_data = read_dgf(in_file)

    print("num nodes = {}".format(len(nodes)))
    # print(nodes)
    # print(nodes_data)

    print("num segments = {}".format(len(segments)))
    # print(segments)
    # print(segments_data)

    # write
    write_vtk(out_file, nodes, nodes_data, segments, segments_data)



convert_dgf_to_vtk("test.dgf", "test.vtk")
