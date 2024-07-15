import json 
import pathlib
import os
import argparse
import copy


def subdivide(graph, name, name_vertex=None):
	last_vertex_id = max([v['id'] for v in graph['vertices']])
	last_edge_id = max([e['id'] for e in graph['vessels']])
	old_edge = [e for e in graph['vessels'] if e['name'] == name][0]
	new_vertex = {'id': last_vertex_id+1 } 
	if name_vertex is not None:
		new_vertex['name'] = name_vertex
	graph['vertices'].append(new_vertex)
	new_edge = dict(old_edge)
	new_edge['id'] = last_edge_id+1
	graph['vessels'].append(new_edge)
	old_edge['right_vertex_id'] = new_vertex['id']
	new_edge['left_vertex_id'] = new_vertex['id']
	old_edge['vessel_length'] *= 0.5
	new_edge['vessel_length'] *= 0.5
	num_coordinates = len(old_edge['embedded_coordinates'])
	assert len(old_edge['embedded_coordinates']) == old_edge['number_edges']+1
	assert len(old_edge['embedded_coordinates']) > 2 
	old_edge['embedded_coordinates'] = old_edge['embedded_coordinates'][0:int(num_coordinates/2)+1]
	old_edge['number_edges'] = len(old_edge['embedded_coordinates'])-1
	new_edge['embedded_coordinates'] = new_edge['embedded_coordinates'][int(num_coordinates/2):-1]
	new_edge['number_edges'] = len(new_edge['embedded_coordinates'])-1

def _get_next_vertex_id(graph):
	return max([v['id'] for v in graph['vertices']]) + 1

def _get_next_edge_id(graph):
	return max([v['id'] for v in graph['vessels']]) + 1

def _find_edge_by_name(graph, name):
	return [e for e in graph['vessels'] if e['name'] == name][0]

def _find_edge_by_id(graph, id):
	return [e for e in graph['vessels'] if e['id'] == id][0]

def _find_vertex_by_id(graph, id):
	return [v for v in graph['vertices'] if v['id'] == id][0]

def _find_left_vertex(graph, edge):
	candidates = [v for v in graph['vertices'] if v['id'] == edge['left_vertex_id']]
	assert len(candidates) == 1
	return candidates[0]

def _find_right_vertex(graph, edge):
	candidates = [v for v in graph['vertices'] if v['id'] == edge['right_vertex_id']]
	assert len(candidates) == 1
	return candidates[0]

def _find_common_vertex(graph, edges):
	vertex_set = { edges[0]['left_vertex_id'], edges[0]['right_vertex_id'] }
	for e in edges:
		vertex_set &= { e['left_vertex_id'], e['right_vertex_id'] }
	assert len(vertex_set) == 1
	vid = vertex_set.pop()
	return _find_vertex_by_id(graph, vid)
	
def _find_other_vertex(graph, edge, vertex):
	if edge['left_vertex_id'] == vertex['id']:
		return _find_right_vertex(graph, edge)
	elif edge['right_vertex_id'] == vertex['id']:
		return _find_left_vertex(graph, edge)
	else:
		raise RuntimeError('The vertex is not connected to the edge?')

def _delete_vertex(graph, vertex):
	vid = vertex['id']
	graph['vertices'] = [v for v in graph['vertices'] if v['id'] != vid]
	for v in graph['vertices']:
		assert v['id'] != vid
		if v['id'] > vid:
			v['id'] -= 1
	for e in graph['vessels']:
		assert e['left_vertex_id'] != vid 
		assert e['right_vertex_id'] != vid 
		if e['left_vertex_id'] > vid:
			e['left_vertex_id'] -= 1
		if e['right_vertex_id'] > vid:
			e['right_vertex_id'] -= 1

def _delete_edge(graph, edge):
	eid = edge['id']
	graph['vessels'] = [e for e in graph['vessels'] if e['id'] != eid]
	for e in graph['vessels']:
		assert e['id'] != eid
		if e['id'] > eid:
			e['id'] -= 1

def _add_vertex(graph):
	vertex = {'id': _get_next_vertex_id(graph) }
	graph['vertices'].append(vertex)
	return vertex

def _add_edge(graph, edge_blueprint):
	edge = edge_blueprint.copy()
	edge['id'] = _get_next_edge_id(graph)
	graph['vessels'].append(edge)
	return edge 

def _split_edge(graph, edge):
	v1 = _find_left_vertex(graph, edge)
	v2 = _add_vertex(graph)
	v3 = _add_vertex(graph)
	v4 = _find_right_vertex(graph, edge)
	e2 = _add_edge(graph, edge)
	e3 = _add_edge(graph, edge)
	#edge['name'] += ' (1)'
	#e2['name'] += ' (2)'
	#e3['name'] += ' (3)'
	edge['right_vertex_id'] = v2['id'] 
	e2['left_vertex_id'] = v2['id'] 
	e2['right_vertex_id'] = v3['id'] 
	e3['left_vertex_id'] = v3['id'] 
	e3['right_vertex_id'] = v4['id'] 
	coords = edge['embedded_coordinates'] 
	for idx, e in enumerate([edge, e2, e3]):
		e['embedded_coordinates'] = coords[int((len(coords))/3*idx):int((len(coords))/3*(idx+1))+1]
		e['number_edges'] = len(e['embedded_coordinates'])-1
	return [v1,v2,v3,v4], [edge, e2, e3]

def _get_edge_neighbors_of_vertex(graph, vertex):
	return [e for e in graph['vessels'] if e['right_vertex_id'] == vertex['id'] or e['left_vertex_id'] == vertex['id']]

def _get_adjacent_vertex_ids_from_edges(edges):
	connected_vertices = [e['left_vertex_id'] for e in edges] + [e['right_vertex_id'] for e in edges]
	connected_vertices = list(set(connected_vertices)) 
	return connected_vertices

def prepare_bifurcation(graph, edges):
	common_vertex = _find_common_vertex(graph, edges)
	print(common_vertex)

	for idx,edge in enumerate(edges):
		[v1,v2,v3,v4], [e1,e2,e3] = _split_edge(graph, edge)
		if v4['id'] == common_vertex['id']:
			e3['name'] = 'gap'
			v3['name'] = f'coupling_{idx}_outer'
			v2['name'] = f'coupling_{idx}_inner'
		elif v1['id'] == common_vertex['id']:
			e1['name'] = 'gap'
			v2['name'] = f'coupling_{idx}_outer'
			v3['name'] = f'coupling_{idx}_inner'
		else:
			raise RuntimeError('Not possible.')
		e2['name'] = f'vessel_coupling_{idx}'
	
	graph_with_gap = copy.deepcopy(graph)
	for edge in _get_edge_neighbors_of_vertex(graph_with_gap, common_vertex):
		_delete_edge(graph_with_gap, edge)
	_delete_vertex(graph_with_gap, common_vertex)

	graph_gap = copy.deepcopy(graph)
	edges = _get_edge_neighbors_of_vertex(graph_gap, common_vertex)
	connected_vertices = _get_adjacent_vertex_ids_from_edges(edges)
	edges = [e for vid in connected_vertices for e in _get_edge_neighbors_of_vertex(graph_gap, _find_vertex_by_id(graph_gap, vid))]
	edge_ids = list(set([e['id'] for e in edges]))
	edges_to_delete = [e for e in graph_gap['vessels'] if e['id'] not in edge_ids]
	for e in edges_to_delete:
		_delete_edge(graph_gap, e)
	connected_vertices = _get_adjacent_vertex_ids_from_edges(graph_gap['vessels'])
	vertices_to_delete = [v for v in graph_gap['vertices'] if v['id'] not in connected_vertices]
	for v in vertices_to_delete:
		_delete_vertex(graph_gap, v)
	for v in graph_gap['vertices']:
		name = v['name']
		# outer and inner switch roles
		if '_inner' in name:
			v['name'] = v['name'].replace('_inner', '_outer')
		if '_outer' in name:
			v['name'] = v['name'].replace('_outer', '_inner')
	
	return graph, graph_with_gap, graph_gap

def subdivide_for_coupling(graph, name, omit_connecting_edge=True):
	last_vertex_id = max([v['id'] for v in graph['vertices']])
	last_edge_id = max([e['id'] for e in graph['vessels']])
	old_edge = [e for e in graph['vessels'] if e['name'] == name][0]
	v1 = [v for v in graph['vertices'] if v['id'] == old_edge['left_vertex_id']][0]
	v4 = [v for v in graph['vertices'] if v['id'] == old_edge['right_vertex_id']][0]
	# create new vertices
	v2 = {'id': last_vertex_id+1}
	v3 = {'id': last_vertex_id+2}
	graph['vertices'].append(v2)
	graph['vertices'].append(v3)
	# create edges
	edge2 = dict(old_edge)
	edge2['id'] = last_edge_id + 2 # choose the largest edge index of edge 2, since it is not always connected
	if not omit_connecting_edge:
		graph['vessels'].append(edge2)
	edge3 = dict(old_edge)
	edge3['id'] = last_edge_id + 1
	graph['vessels'].append(edge3)
	# fix connectivity
	old_edge['right_vertex_id'] = v2['id']
	edge2['left_vertex_id'] = v2['id']
	edge2['right_vertex_id'] = v3['id']
	edge3['left_vertex_id'] = v3['id']
	edge3['right_vertex_id'] = v4['id']
	# fix edge names
	old_edge['name'] = 'vessel_coupling_1'
	edge2['name'] = 'gap'
	edge3['name'] = 'vessel_coupling_2'
	# rename vertices
	v1['name'] = 'coupling_1_inner'
	v2['name'] = 'coupling_1_outer'
	v3['name'] = 'coupling_2_outer'
	v4['name'] = 'coupling_2_inner'
	# fix length
	old_edge['vessel_length'] /= 3
	edge2['vessel_length'] /= 3
	edge3['vessel_length'] /= 3
	# fix embedded coordinates 
	assert len(old_edge['embedded_coordinates']) == old_edge['number_edges']+1
	assert len(old_edge['embedded_coordinates']) > 3 
	num_coordinates = len(old_edge['embedded_coordinates'])
	old_edge['embedded_coordinates'] = old_edge['embedded_coordinates'][0:int(num_coordinates/3)+1]
	old_edge['number_edges'] = len(old_edge['embedded_coordinates'])-1
	edge2['embedded_coordinates'] = edge2['embedded_coordinates'][int(num_coordinates/3):int(num_coordinates/3*2)+1]
	edge2['number_edges'] = len(edge2['embedded_coordinates'])-1
	edge3['embedded_coordinates'] = edge3['embedded_coordinates'][int(num_coordinates/3*2):-1]
	edge3['number_edges'] = len(edge3['embedded_coordinates'])-1
	# generate the coupled graph for testing
	new_graph = {}
	old_edge,edge2,edge3 = dict(old_edge), dict(edge2), dict(edge3)
	v1,v2,v3,v4 = dict(v1),dict(v2),dict(v3),dict(v4)
	old_edge['id'] = 0
	edge2['id'] = 1
	edge3['id'] = 2
	v1['id'] = 0
	v2['id'] = 1
	v3['id'] = 2
	v4['id'] = 3
	v1['name'] = 'coupling_1_outer'
	v2['name'] = 'coupling_1_inner'
	v3['name'] = 'coupling_2_inner'
	v4['name'] = 'coupling_2_outer'
	old_edge['left_vertex_id'] = 0
	old_edge['right_vertex_id'] = 1
	edge2['left_vertex_id'] = 1
	edge2['right_vertex_id'] = 2
	edge3['left_vertex_id'] = 2
	edge3['right_vertex_id'] = 3
	new_graph['vessels'] = [old_edge,edge2,edge3]
	new_graph['vertices'] = [v1,v2,v3,v4]
	return graph, new_graph


def subdivide_bifurcation(graph, names):
	pass


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--filepath', type=str, help='Path to the json file specifying the mesh.', required=False)
	parser.add_argument('--names', type=str, help='Names of the vessels', nargs='+')
	parser.add_argument('--ids', type=int, help='Names of the vessels', nargs='+')
	args = parser.parse_args()

	#filepath = '/home/wagneran/Flows1D0D3D_aneurysm/data/1d-meshes/bifurcation.json'
	#name = 'belongs to vmtk branch 21'
	#filepath = args.filepath 

	filepath = args.filepath
	names = args.names
	ids = args.ids

	with open(filepath, 'r') as file:
		data = file.read()
		graph = json.loads(data)

		cli_edges = []
		if args.names:
			cli_edges += [_find_edge_by_name(graph, name) for name in args.names] 
		if args.ids:
			cli_edges += [_find_edge_by_id(graph, idx) for idx in args.ids]
		
		# cli_edges += [_find_edge_by_id(graph, idx) for idx in [0,1,2]]

		graph, graph_with_gap, graph_gap = prepare_bifurcation(graph, cli_edges)

		p = pathlib.Path(filepath)
		folder = pathlib.Path(filepath).parent.resolve()
		with open(os.path.join(folder, p.stem + '-with-gap' + p.suffix),'w') as file:
			file.write(json.dumps(graph_with_gap, indent=2))
		with open(os.path.join(folder, p.stem + '-gap' + p.suffix), 'w') as file:
			file.write(json.dumps(graph_gap, indent=2))

		print('graph')
		print(json.dumps(graph, indent=1))
		print('graph with gap')
		print(json.dumps(graph_with_gap, indent=1))
		print('graph gap')
		print(json.dumps(graph_gap, indent=1))



