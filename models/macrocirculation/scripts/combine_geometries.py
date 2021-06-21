"""
Simple script which combines two input geometries by adding artificial vessels between them.
"""
import json
import numpy as np


class ConnectionPair:
    def __init__(self, name0, name1, name2, length, number_micro_edges, radius, elastic_modulus, wall_thickness):
        self.name0 = name0
        self.name1 = name1
        self.name2 = name2
        self.length = length
        self.number_micro_edges = number_micro_edges
        self.radius = radius
        self.elastic_modulus = elastic_modulus
        self.wall_thickness = wall_thickness


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
    ConnectionPair('cw_con_1', 'bg_132', 'bg_141', [5, 2.5, 5], [120, 60, 120], radius=0.06, elastic_modulus=1300000.0, wall_thickness=0.005),
    ConnectionPair('cw_con_2', 'bg_135', 'bg_119', [5, 2.5, 7], [120, 60, 168], radius=0.06, elastic_modulus=1300000.0, wall_thickness=0.005)
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

max_vessel_id = max([e["id"] for e in geo['vessels']])
max_vertex_id = max([v["id"] for v in geo['vertices']])

for cpair in to_connect:
    v0 = find_vertex_by_name(geo, cpair.name0)
    v1 = find_vertex_by_name(geo, cpair.name1)
    v2 = find_vertex_by_name(geo, cpair.name2)

    new_vertex = {'id': max_vertex_id+1, 'name': 'mp_({})_({})_({})'.format(v0['name'], v1['name'], v2['name'])}
    max_vertex_id += 1

    geo['vertices'] += [new_vertex]

    new_microedges = []
    new_microedges.append({
        'id': max_vessel_id + 1,
        'name': 'con_vertex({})_vertex({})'.format(v0['name'], new_vertex['name']),
        'vessel_length': cpair.length[0],
        'radius': cpair.radius,
        'wall_thickness': cpair.wall_thickness,
        'elastic_modulus': cpair.elastic_modulus,
        'gamma': 2,
        'number_edges': cpair.number_micro_edges[0],
        'left_vertex_id': v0['id'],
        'right_vertex_id': new_vertex['id'],
    })
    max_vessel_id += 1
    new_microedges.append({
        'id': max_vessel_id + 1,
        'name': 'con_vertex[{}]_vertex[{}]'.format(new_vertex['name'], v1['name']),
        'vessel_length': cpair.length[1],
        'radius': cpair.radius,
        'wall_thickness': cpair.wall_thickness,
        'elastic_modulus': cpair.elastic_modulus,
        'gamma': 2,
        'number_edges': cpair.number_micro_edges[1],
        'left_vertex_id': new_vertex['id'],
        'right_vertex_id': v1['id'],
    })
    max_vessel_id += 1
    new_microedges.append({
        'id': max_vessel_id + 1,
        'name': 'con_vertex[{}]_vertex[{}]'.format(new_vertex['name'], v2['name']),
        'vessel_length': cpair.length[2],
        'radius': cpair.radius,
        'wall_thickness': cpair.wall_thickness,
        'elastic_modulus': cpair.elastic_modulus,
        'gamma': 2,
        'number_edges': cpair.number_micro_edges[2],
        'left_vertex_id': new_vertex['id'],
        'right_vertex_id': v0['id'],
    })
    max_vessel_id += 1

    geo['vessels'] += new_microedges

with open(output_path, 'w') as f:
    f.write(json.dumps(geo, indent=2))
