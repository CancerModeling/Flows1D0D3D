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
    json_data['num_vertices'] = len(data['Network']['Nodes'][0][0][0])
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
            'embedded_coordinates': points,
            'abstract_coordinates': abstract_coordinates,
            'vessel_length': vessel_length,
            'micro_edge_length': vessel_length/(num_points-1),
            'radii': radii,
            'pv0': pv0
        }
        json_data['vessels'].append(vessel)

    return json_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-filepath', type=str, help='Path to the mat file containing the data.')
    parser.add_argument('--output-filepath', type=str, help='Path to the json file containing the refined mesh.')
    parser.add_argument('--mesh-width', type=float, help='The mesh width, which we want to achieve. Should not be larger than 0.18!')
    parser.add_argument('--json-indent', action='store_true', help='Should the to json output file contain indents?')

    args = parser.parse_args()

    data = loadmat(args.input_filepath)

    with open(args.output_filepath, 'w') as f:
        uniform_mesh = create_uniform_mesh(data, args.mesh_width)
        # min_vessel_length = np.array([d['vessel_length'] for d in json_data['vessels']]).min()
        # print('min_vessel_length = ', min_vessel_length)
        json_string = json.dumps(uniform_mesh, indent=args.json_indent)
        f.write(json_string)