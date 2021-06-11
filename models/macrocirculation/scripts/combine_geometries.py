"""
Simple script which combines two input geometries by adding artificial vessels between them.
"""
import json
import numpy as np


class ConnectionPair:
    def __init__(self, name0, name1, length, number_macro_edges, number_micro_edges, start_radius):
        self.name0 = name0
        self.name1 = name1
        self.length = length
        self.number_macro_edges = number_macro_edges
        self.number_micro_edges = number_micro_edges
        self.start_radius = start_radius


def find_vertex_by_name(geo, name):
    candidates = [v for v in geo['vertices'] if v.get('name') == name]
    if len(candidates) > 1:
        raise Exception('candidate for connection not unique')
    if len(candidates) < 1:
        raise Exception('candidate for connection not found')
    return candidates[0]


def find_edges_neighbors_by_vertex_id(geo, v):
    return [e for e in geo['vessels'] if e['left_vertex_id'] == v['id'] or e['right_vertex_id'] == v['id']]


geometry1_path = '../data/network-33-vessels-with-connection.json'
geometry2_path = '../data/coarse-network-geometry.json'
output_path = '../data/combined-network-geometry.json'

to_connect = [
    ConnectionPair('cw_con_1', 'bg_132', 10., 6, 240, start_radius=0.24),
    ConnectionPair('cw_con_2', 'bg_135', 10., 6, 240, start_radius=0.24)
]

with open(geometry1_path) as f:
    geo = json.load(f)

with open(geometry2_path) as f:
    geo2 = json.load(f)

max_vertex_id = max([v["id"] for v in geo['vertices']])
max_vessel_id = max([e["id"] for e in geo['vessels']])

for v in geo2['vertices']:
    v["id"] += max_vertex_id+1

for e in geo2['vessels']:
    e["id"] += max_vessel_id+1
    e["left_vertex_id"] += max_vertex_id+1
    e["right_vertex_id"] += max_vertex_id+1

geo['vertices'] += geo2['vertices']
geo['vessels'] += geo2['vessels']

max_vertex_id = max([v["id"] for v in geo['vertices']])
max_vessel_id = max([e["id"] for e in geo['vessels']])

for cpair in to_connect:
    v0 = find_vertex_by_name(geo, cpair.name0)
    v1 = find_vertex_by_name(geo, cpair.name1)

    neighbors_v0 = find_edges_neighbors_by_vertex_id(geo, v0)
    neighbors_v1 = find_edges_neighbors_by_vertex_id(geo, v1)

    elastic_modulus_v0 = np.mean([e['elastic_modulus'] for e in neighbors_v0])
    wall_thickness_v0 = np.mean([e['wall_thickness'] for e in neighbors_v0])
    radius_v0 = cpair.start_radius
    elastic_modulus_v1 = np.mean([e['elastic_modulus'] for e in neighbors_v1])
    wall_thickness_v1 = np.mean([e['wall_thickness'] for e in neighbors_v1])
    radius_v1 = np.mean([e.get('radius') or np.mean(e['radii']) for e in neighbors_v1])

    wall_thickness_interp = np.linspace(wall_thickness_v0, wall_thickness_v1, cpair.number_macro_edges+2)[1:-1]
    elastic_modulus_interp = np.linspace(elastic_modulus_v0, elastic_modulus_v1, cpair.number_macro_edges+2)[1:-1]
    radius_interp = np.linspace(radius_v0, radius_v1, cpair.number_macro_edges+2)[1:-1]

    num_micro_edges = max(2, int(np.floor(cpair.number_micro_edges / cpair.number_macro_edges)))

    vessel_length = cpair.length / cpair.number_macro_edges

    # create new microvertices
    new_microvertices = [{'id': i+max_vertex_id+1 } for i in range(cpair.number_macro_edges-1)]
    geo['vertices'] += new_microvertices
    max_vertex_id += len(new_microvertices)

    # we collect all the microvertices in one list
    all_microvertices = [v0] + new_microvertices + [v1]

    new_microedges = []
    for i in range(cpair.number_macro_edges):
        left_vertex = all_microvertices[i]
        right_vertex = all_microvertices[i+1]
        d = {
            'id': max_vessel_id + i + 1,
            'name': 'con_vertex({})_vertex({})_({})'.format(cpair.name0, cpair.name1, i),
            'vessel_length': vessel_length,
            'radius': radius_interp[i],
            'wall_thickness': wall_thickness_interp[i],
            'elastic_modulus': elastic_modulus_interp[i],
            'number_edges': num_micro_edges,
            'left_vertex_id': left_vertex['id'],
            'right_vertex_id': right_vertex['id'],
        }
        new_microedges.append(d)
    geo['vessels'] += new_microedges
    max_vessel_id += len(new_microedges)

with open(output_path, 'w') as f:
    f.write(json.dumps(geo, indent=2))
