"""
Simple script which combines two input geometries by adding artificial vessels between them.
"""
import json


def find_vertex_by_name(geo, name):
    candidates = [v for v in geo['vertices'] if v.get('name') == name]
    if len(candidates) > 1:
        raise Exception('candidate for connection not unique')
    if len(candidates) < 1:
        raise Exception('candidate for connection not found')
    return candidates[0]


def find_edges_neighbors_by_vertex_id(geo, v):
    return [e for e in geo['vessels'] if e['left_vertex_id'] == v['id'] or e['right_vertex_id'] == v['id']]


input_path = '../data/boundary-combined-network-geometry.json'
output_path = '../data/boundary-combined-network-geometry-v2.json'

with open(input_path) as f:
    data = json.load(f)

for d in data['vertices']:
    if d['name'].startswith('bg_'):
        d['peripheral_resistance'] = d['peripheral_resistance']/4.

with open(output_path, 'w') as f:
    f.write(json.dumps(data, indent=2))
