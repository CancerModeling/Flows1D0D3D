import numpy as np

def read_vtk(in_file):

    print('VTK reader not implemented')


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
