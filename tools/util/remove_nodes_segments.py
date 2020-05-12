import numpy as np

from read_write_dgf import *
from read_write_vtk import write_vtk

def add_unique(data, a):

    found = False
    for i in range(len(data)):
        if data[i] == a:
            found = True

    if found == False:
        data.append(a)


def get_segment_ids_neigh_nodes(nodes, segments, vertex):

    num_nodes = len(nodes)
    num_segments = len(segments)

    segs = []
    nds = []

    for i in range(num_segments):

        s = segments[i]

        for j in range(len(s)):

            if s[j] == vertex:
                segs.append(i)

                # add other node as neighbor
                if j == 0:
                    add_unique(nds, s[1])
                else:
                    add_unique(nds, s[0])



    return nds, segs


def remove_node(in_file, vertex):

    # read dgf file
    nodes_old, nodes_data_old, segments_old, segments_data_old = read_dgf(in_file)

    num_nodes_old = len(nodes_old)
    num_segments_old = len(segments_old)

    # get neighboring nodes and segments
    nds, segs = get_segment_ids_neigh_nodes(nodes_old, segments_old, vertex)
    print('nds: {}, segs: {}'.format(nds, segs))

    # vertex and associated data
    nodes = []
    nodes_data = []

    # elements and associated data
    segments = []
    segments_data = []

    for i in range(num_nodes_old):

        coords = nodes_old[i]

        if i != vertex:
            nodes.append(coords)
            nodes_data.append(nodes_data_old[i])


    for i in range(num_segments_old):

        s = segments_old[i]

        found = False
        for j in range(len(segs)):
            if segs[j] == i:
                found = True

        if found == False:
            snew = []
            if s[0] < vertex:
                snew.append(s[0])
            else:
                snew.append(s[0] - 1)

            if s[1] < vertex:
                snew.append(s[1])
            else:
                snew.append(s[1] - 1)

            segments.append(snew)
            segments_data.append(segments_data_old[i])


    # because of removal of vertex, another vertex may have 
    # become the boundary vertex
    for i in range(len(nodes)):

        for j in range(len(nds)):

            # account for substraction of 1 from vertex ids in new data
            search_nd = nds[j]
            if search_nd >= vertex:
                search_nd = search_nd - 1

            if search_nd == i:
                print('search_nd: {}, i: {}'.format(search_nd, i))
                # this node was neighbor to removed vertex

                # if the number of neighbors of this node is 1 then this is 
                # a boundary node
                nds_i, segs_i = get_segment_ids_neigh_nodes(nodes, segments, i)
                print('nds_i: {}, segs_i: {}'.format(nds_i, segs_i))
                if len(segs_i) == 1:
                    print('nd data:  {}'.format(nodes_data_old[vertex][0]))
                    # assign pressure that was assigned to removed vertex
                    nodes_data[i][0] = nodes_data_old[vertex][0]

    # return
    return nodes, nodes_data, segments, segments_data



nodes, nodes_data, segments, segments_data = remove_node("diff_pres.dgf", 14)
# write to new dgf file
write_dgf("diff_pres_rmv.dgf", nodes, nodes_data, segments, segments_data)
