import json 
import pathlib
import os


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


def subdivide_for_coupling(graph, name):
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
	edge3 = dict(old_edge)
	edge2['id'] = last_edge_id + 1
	edge3['id'] = last_edge_id + 2
	graph['vessels'].append(edge2)
	graph['vessels'].append(edge3)
	# fix connectivity
	old_edge['right_vertex_id'] = v2['id']
	edge2['left_vertex_id'] = v2['id']
	edge2['right_vertex_id'] = v3['id']
	edge3['left_vertex_id'] = v3['id']
	edge3['right_vertex_id'] = v4['id']
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


filepath = '/home/andreas/Flows1D0D3D/data/1d-meshes/Graph0.json'
name = 'belongs to vmtk branch 21'
with open(filepath, 'r') as file:
	data = file.read()
	graph = json.loads(data)
	graph, new_graph = subdivide_for_coupling(graph, name)
	p = pathlib.Path(filepath)
	folder = pathlib.Path(filepath).parent.resolve()
	with open(os.path.join(folder, p.stem + '_with_gap' + p.suffix),'w') as file:
		file.write(json.dumps(graph, indent=2))
	with open(os.path.join(folder, p.stem + '_gap' + p.suffix), 'w') as file:
		file.write(json.dumps(new_graph, indent=2))
  