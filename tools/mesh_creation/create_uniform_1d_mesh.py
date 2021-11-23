""" Simple script for converting Chengyue Wu's raw network data into a uniform mesh """

import numpy as np
from scipy.io import loadmat
import argparse
import json
import math


def remesh_vessel(reference_coordinates, reference_quantities, num_points):
    """
    Transforms the reference coordinates and reference quantities into a coarsened or refined representation consisting
    of the given number of points with linear interpolation.
    """
    num_quantities = len(reference_quantities)
    coordinates = np.array(reference_coordinates)
    quantities = np.array(reference_quantities)
    new_coordinates = np.linspace(0, reference_coordinates[-1], num_points).tolist()
    new_quantities = []
    for new_c in new_coordinates:
        idx = np.searchsorted(coordinates, new_c)
        if (idx >= len(coordinates)-1):
            idx = len(coordinates)-2
        left_c = coordinates[idx]
        right_c = coordinates[idx+1]
        left_q = quantities[:, idx]
        right_q = quantities[:, idx+1]
        new_q = right_q * (new_c - left_c)/(right_c - left_c) + left_q * (right_c - new_c)/(right_c - left_c)
        new_quantities.append(new_q.tolist())
    new_quantities = np.array(new_quantities).transpose()
    return new_coordinates, new_quantities.tolist()


def get_points(segment):
    points = []
    num_points = len(segment['YXZ'][0])
    for point_idx in range(num_points):
        point = [segment['YXZ'][1][point_idx], segment['YXZ'][0][point_idx], segment['YXZ'][2][point_idx]]
        points.append(point)
    return np.array(points)


def create_uniform_mesh(data, h):
    json_data = {}
    num_vertices = len(data['Network']['Nodes'][0][0][0])
    json_data['vertices'] = []
    for id in range (num_vertices):
        json_data['vertices'].append({'id': id, 'name': 'bg_{}'.format(id)})

    json_data['vessels'] = []
    num_segments = len(data["Network"]["Links"][0][0][0])
    for seg_id in range(num_segments):
        # extract vessel data
        segment = data["Network"]["Links"][0][0][0][seg_id]
        radii = segment['radii'][0].tolist()
        pv0 = segment['pv0'][0].tolist()

        # infer data
        points = get_points(segment)
        diff = points[1:] - points[0:-1]
        norms = np.linalg.norm(diff, axis=1)
        abstract_coordinates = [0] + np.cumsum(norms).tolist()
        vessel_length = abstract_coordinates[-1]
        print(vessel_length)

        if (h > vessel_length):
            raise "h ({}) > vessel_length ({}) not allowed".format(h, vessel_length)

        # remesh the vessel
        num_points = math.ceil(vessel_length / h) + 1
        quantities = [
            points.transpose()[0].tolist(),
            points.transpose()[1].tolist(),
            points.transpose()[2].tolist(),
            radii,
            pv0
        ]
        abstract_coordinates, quantities = remesh_vessel(abstract_coordinates, quantities, num_points)
        points = np.array(quantities[0:3]).transpose().tolist()
        radii = quantities[3]
        pv0 = quantities[4]

        # create vessel
        vessel = {
            'id': int(seg_id),
            'left_vertex_id': int(segment[0][0][0]-1),
            'right_vertex_id': int(segment[1][0][0]-1),
            'name': 'bg_{}'.format(seg_id),
            'embedded_coordinates': points,
            'number_edges': len(points)-1,
            'vessel_length': vessel_length,
            'wall_thickness': 0.005,      # completely artificial
            'elastic_modulus': 1.3e6,     # completely artificial
            'gamma': 2,                   # Poisseuile flow in the breast geometry
            'radii': radii,
        }
        json_data['vessels'].append(vessel)

    return json_data


def reduce_network(data, max_num_vessels):
    vertices_to_edges = []
    for i in range(len(data['vertices'])):
        vertices_to_edges.append([])

    for vessel in data['vessels']:
        vertices_to_edges[vessel['left_vertex_id']].append(vessel['id'])
        vertices_to_edges[vessel['right_vertex_id']].append(vessel['id'])

    #start_vertex = 135
    start_vertex = 0 

    added_edges = []
    active_vertices = set([start_vertex])
    finished_vertices = set()

    num_edges_to_add = max_num_vessels

    while (len(active_vertices) > 0):
        vidx = active_vertices.pop()
        finished_vertices.add(vidx)
        print(vertices_to_edges, vidx, len(vertices_to_edges))
        for eidx in vertices_to_edges[vidx]:
            if len(added_edges) < num_edges_to_add and not eidx in set(added_edges):
                print('added', eidx)
                added_edges.append(eidx)
            lvid = data['vessels'][eidx]['left_vertex_id']
            rvid = data['vessels'][eidx]['right_vertex_id']
            if not lvid in finished_vertices and not lvid in active_vertices:
                active_vertices.add(lvid)
            if not rvid in finished_vertices and not rvid in active_vertices:
                active_vertices.add(rvid)

    added_vertices = [start_vertex]
    for eidx in added_edges:
        lvid = data['vessels'][eidx]['left_vertex_id']
        rvid = data['vessels'][eidx]['right_vertex_id']
        if not lvid in set(added_vertices):
            added_vertices.append(lvid)
        if not rvid in set(added_vertices):
            added_vertices.append(rvid)

    vertex_renumbering = {}
    for vid in added_vertices:
        print(vid)
        vertex_renumbering[vid] = len(vertex_renumbering)

    edges_list = []
    for idx, eidx in enumerate(added_edges):
        edge = data['vessels'][eidx]
        edge['id'] = idx
        edge['left_vertex_id'] = vertex_renumbering[edge['left_vertex_id']]
        edge['right_vertex_id'] = vertex_renumbering[edge['right_vertex_id']]
        edges_list.append(edge)

    vertices_list = []
    for vid in added_vertices:
        v = data['vertices'][vid]
        v['id'] = vertex_renumbering[vid]
        vertices_list.append(v)

    data['vessels'] = edges_list
    data['vertices'] = vertices_list

    print(added_edges)
    print(added_vertices)
    print(vertex_renumbering)

    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-filepath', type=str, help='Path to the mat file containing the data.')
    parser.add_argument('--output-filepath', type=str, help='Path to the json file containing the refined mesh.')
    parser.add_argument('--mesh-width', type=float, help='The mesh width, which we want to achieve. Should not be larger than 0.18!')
    parser.add_argument('--json-indent', action='store_true', help='Should the to json output file contain indents?')
    parser.add_argument('--max-num-vessels', type=int, help='The maximum number of vessels allowed', default=100000)

    args = parser.parse_args()

    data = loadmat(args.input_filepath)

    with open(args.output_filepath, 'w') as f:
        uniform_mesh = create_uniform_mesh(data, args.mesh_width)
        uniform_mesh = reduce_network(uniform_mesh, args.max_num_vessels)
        # min_vessel_length = np.array([d['vessel_length'] for d in json_data['vessels']]).min()
        # print('min_vessel_length = ', min_vessel_length)
        json_string = json.dumps(uniform_mesh, indent=args.json_indent)
        f.write(json_string)
