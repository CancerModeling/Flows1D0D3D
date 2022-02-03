from matplotlib import pyplot as plt
import numpy as np
import json
import os


def reformat_data(d):
    if len(d.shape) == 1:
        d = d.reshape((len(d), 1))
    return d


def midpoint(d):
    return 0.5 * (d[1:] + d[:-1])


class TipDataLoader:
    def __init__(self, filepath, t_start=0, t_end=100) -> None:
        self.directory = os.path.dirname(filepath)
        with open(filepath) as file:
            self.meta = json.loads(file.read())
        self.t_full = self.t = np.array(self.meta['times'])
        self.start_index = np.sum(self.t < t_start)
        self.end_index = np.sum(self.t < t_end)
        self.t = self.t[self.start_index: self.end_index] 
        self.t_mid = 0.5 * (self.t_full[1:] + self.t_full[:-1])
        self.t_mid = self.t_mid[self.start_index: self.end_index-1]
        self.tau = self.t_full[1:] - self.t_full[:-1]
        self.tau = self.tau[self.start_index: self.end_index-1]

    def find_vessel_by_edge_id(self, edge_id):
        for v in self.meta['vertices']:
            if v['neighbor_edge_id'] == edge_id:
                return v
        return None
    
    def _load_path(self, filepath):
        return reformat_data(np.loadtxt(os.path.join(self.directory, filepath)))[self.start_index: self.end_index, :]
    
    def load_data_of_vessel(self, v):
        data = {}
        if 'filepath_p' in v['filepaths']:
            data['p'] = self._load_path(v['filepaths']['filepath_p'])
        if 'filepath_c' in v['filepaths']:
            data['c'] = self._load_path(v['filepaths']['filepath_c'])
        if 'filepath_V' in v['filepaths']:
            data['V'] = self._load_path(v['filepaths']['filepath_V'])
        if 'filepath_p_out' in v:
            data['p_out'] = self._load_path(v['filepath_p_out'])

        q_list = [] 
        q_total_list = [] 
        for dof in range(v['num_dofs']):
            pressures_here = data['p'][:, dof]

            if dof+1 < data['p'].shape[1]:
                pressures_next = data['p'][:, dof+1]
            elif 'p_out' in data:
                pressures_next = data['p_out'].squeeze()
            else:
                pressures_next = vessel['p_out'] * np.ones(len(data['p'][:, dof]))

            if 'resistance' in vessel:
                resistances = np.array(vessel['resistance'])
            else:
                resistances = np.array([vessel['R2']])
            q = ((pressures_here - pressures_next) / resistances[dof])
            q_total = q * 2 ** (dof + 1)
            q_list.append(q.tolist())
            q_total_list.append( q_total.tolist() )
        
        data['q'] = np.array(q_list).transpose()
        data['q_total'] = np.array(q_total_list).transpose()

        return data
    
    def load_data_by_edge_id(self, edge_id):
        v = self.find_vessel_by_edge_id(edge_id)
        return self.load_data_by_vessel(v)




filepath = 'tmp_transport/heart_to_breast_1d_solution_tips_li.json'
t_start = 24
t_end = 25 
loader = TipDataLoader(filepath, t_start, t_end)

flows = []
for vessel in loader.meta['vertices']:
    d = loader.load_data_of_vessel(vessel)
    flow_last = (midpoint(d['q_total'][:,-1]) * loader.tau / (t_end - t_start)).sum()
    flow_blast = (midpoint(d['q_total'][:,-2]) * loader.tau / (t_end - t_start)).sum()
    print(flow_blast, flow_last)
    flows.append(flow_last)
flows = np.array(flows)

print('sum {}, mean {}, min = {}, max = {}'.format(flows.sum(), flows.mean(), flows.min(), flows.max()))