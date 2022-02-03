import json


filepath_tips = '../../python/examples/tmp/heart_to_breast_1d_solution_tips_li.json'
filepath_avg = '../../python/examples/tmp/average_quantities.json'

with open(filepath_tips, 'r') as file:
    data_tips = json.loads(file.read())

with open(filepath_avg, 'r') as file:
    data_avg = json.loads(file.read())

zone_ids = [21, 4, 24, 30]

for zid in zone_ids:
    for vdata in data_tips['vertices']:
        eid = int(vdata['neighbor_edge_id'])
        vid = vdata['vertex_id']
        pzid = adata = data_avg['quantities']['{}'.format(vid)]['idx']
        if pzid == zid:
            print('edge = {} with vertex = {} and zone = {}'.format(eid, vid, pzid))