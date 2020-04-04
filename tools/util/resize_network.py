import numpy as np

from read_write_dgf import *
from read_write_vtk import write_vtk

def get_min_max_radius(segments_data, rad_loc):
    vec = []
    for i in range(len(segments_data)):
        vec.append(segments_data[i][rad_loc])
    vec = np.array(vec)

    return np.max(vec), np.min(vec)

def get_min_max(nodes):

    # get min and max coordinates
    cord_min = [0., 0., 0.]
    cord_max = [0., 0., 0.]

    num_nodes = len(nodes)

    for j in range(3):

        vec = []
        for i in range(num_nodes):
            vec.append(nodes[i][j])

        vec = np.array(vec)
        cord_min[j] = np.min(vec)
        cord_max[j] = np.max(vec)

        # print(j)
        # print(vec)

    cord_max = np.array(cord_max)
    cord_min = np.array(cord_min)

    return cord_max, cord_min

def scale_network(in_file, out_file, scale, scale_radius):

    # scale_radius: scale the final radius by this value (usefule when needed to increase the radius so that we can consider coarser mesh in 3D)

    # read dgf file
    nodes, nodes_data, segments, segments_data = read_dgf(in_file)

    num_nodes = len(nodes)
    num_segments = len(segments)
    print("num nodes = {}".format(num_nodes))
    # print(nodes)
    # print(nodes_data)

    print("num segments = {}".format(num_segments))
    # print(segments)
    # print(segments_data)

    # scale the radius
    for i in range(num_segments):
        segments_data[i][0] = segments_data[i][0] * scale_radius

    cord_max, cord_min = get_min_max(nodes)

    print('\nMin-max coordinates')
    print(cord_max)
    print(cord_min)

    # we need to create a padding so that cylinders do not touch the boundary
    # padding size is dictated by maximum radius of cylinder
    max_rad, min_rad = get_min_max_radius(segments_data, 0)
    padding = 1.1 * max_rad
    print('\nPadding: {}'.format(padding))

    # now translate the cordinates of all nodes so that the minimum cord is
    # at (0,0,0)
    translation = np.multiply(cord_min, -1.)

    # add padding to translation
    for i in range(3):
        translation[i] += padding

    # translate all nodes
    for i in range(num_nodes):
        for j in range(3):
            nodes[i][j] += translation[j]

    # get min max and check
    cord_max_check, cord_min_check = get_min_max(nodes)

    print('\nMin-max coordinates after translation')
    print(cord_max_check)
    print(cord_min_check)

    # we reset min_cord to origin and max_cord to max_cord + 2 * padding
    for i in range(3):
        cord_min[i] = 0.
        cord_max[i] = cord_max[i] + 2. * padding

    # scaling method
    scale_vec = [1., 1., 1.]
    uniform_scaling = False
    if uniform_scaling:
        print('\nApplying uniform scaling')
        # below is uniform volumetric scaling method
        # now check the distance between min and max and scale this distance so
        # that new distance is equal to scale value
        distance = 0.
        for j in range(3):
            distance += (cord_max[j] - cord_min[j]) * (cord_max[j] - cord_min[j])
        distance = np.sqrt(distance)

        scale_factor = scale / distance

        # scaling vector
        for i in range(3):
            scale_vec[i] = scale_factor
    else:
        print('\nApplying non-uniform scaling')
        # get the scaling vector
        for j in range(3):
            scale_vec[j] = scale / cord_max[j]

    # scale all nodes
    for i in range(num_nodes):
        for j in range(3):
            nodes[i][j] *= scale_vec[j]

    # check of new max and min are within expected range
    cord_max_check, cord_min_check = get_min_max(nodes)
    print('\nMin-max coordinates after scaling')
    print(cord_max_check)
    print(cord_min_check)

    # scale the radius of segments_data
    # We assume radius is in the first
    scale_factor = np.min(scale_vec)
    for i in range(num_segments):
        segments_data[i][0] = segments_data[i][0] * scale_factor * scale_radius

    # output max and min radius
    max_rad, min_rad = get_min_max_radius(segments_data, 0)
    print('\nMax radius: {}, min radius: {}'.format(max_rad, min_rad))

    # write to new dgf file
    write_dgf(out_file, nodes, nodes_data, segments, segments_data)

    # also write to vtk
    write_vtk("test_scaled.vtk", nodes, nodes_data, segments, segments_data)


scale_network("test.dgf", "test_scaled_1.dgf", 2., 1.)

# run as follows
# python3 -B resize_network.py
