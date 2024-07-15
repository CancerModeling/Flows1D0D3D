import json 
import math
import pathlib
import os
import argparse


def _find_edge(graph, v):
	return [e for e in graph['vessels'] if e['left_vertex_id'] == v['id'] or e['right_vertex_id'] == v['id']][0]

def _find_vertex(graph, vid):
	return [v for v in graph['vertices'] if v['id'] == vid][0]


def distribute(graph, r, c, predicate=lambda x: True):
	outlets = [v for v in graph['vertices'] if 'Outflow' in v['name'] and predicate(v)]
	initial_areas = {v['id']: _find_edge(graph, v)['radius']**2 * math.pi for v in outlets}
	total_area = sum([a for (v, a) in initial_areas.items()])
	for (vid, a) in initial_areas.items():
		v = _find_vertex(graph, vid)
		ratio = a / total_area
		v['peripheral_resistance'] = r / ratio 
		v['peripheral_compliance'] = c * ratio 
	return graph


if __name__ == '__main__':
	#parser = argparse.ArgumentParser()
	#parser.add_argument('--filepath', type=str, help='Path to the json file specifying the mesh.', required=True)
	#parser.add_argument('--r', type=float, help='Total resitance to distribute', required=True)
	#parser.add_argument('--c', type=float, help='Total capacitance to distribute', required=True)
	#parser.add_argument('--omit', action='store_true', help='Should the connecting vessel be omitted?')
	#args = parser.parse_args()
	#filepath = args.filepath 
	#r = args.r
	#c = args.c

	filepath = '../data/1d-meshes/Graph0.json'

	with open(filepath, 'r') as file:
		data = file.read()
		graph = json.loads(data)
		graph = distribute(graph, r=59.7, c=0.0116, predicate=lambda v: v['id'] > 24)
		graph = distribute(graph, r=84.8, c=0.0082, predicate=lambda v: v['id'] <= 24)
		p = pathlib.Path(filepath)
		folder = pathlib.Path(filepath).parent.resolve()
		with open(os.path.join(folder, p.stem + '-rcr' + p.suffix),'w') as file:
			file.write(json.dumps(graph, indent=2))
	