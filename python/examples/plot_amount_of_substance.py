from matplotlib import pyplot as plt
import numpy as np
import json
import os


class DataLoader:
    def __init__(self, filepath, t_start=0, t_end=100) -> None:
        self.directory = os.path.dirname(filepath)
        with open(filepath) as file:
            self.meta = json.loads(file.read())
        path_time = os.path.join(self.directory, self.meta['filepath_time'])
        self.t_full = self.t = np.loadtxt(path_time, delimiter=',')
        self.start_index = np.sum(self.t < t_start)
        self.end_index = np.sum(self.t < t_end)
        self.t = self.t[self.start_index: self.end_index] 
        self.t_mid = 0.5 * (self.t_full[1:] + self.t_full[:-1])
        self.t_mid = self.t_mid[self.start_index: self.end_index-1]
        self.tau = self.t_full[1:] - self.t_full[:-1]
        self.tau = self.tau[self.start_index: self.end_index-1]
    
    def restrict(self, quantity):
        return quantity.take(indices=range(self.start_index, self.end_index), axis=0)

    def find_vessel(self, vessel_id):
        for vessel in self.meta['vessels']:
            if vessel['edge_id'] == vessel_id:
                return vessel
        return RuntimeError("vessel with id {} not found".format(vessel_id))
    
    def get_quantity(self, vessel_id, quantity_name):
        vessel = self.find_vessel(vessel_id)
        path = os.path.join(self.directory, vessel['filepaths'][quantity_name])
        quantity = np.loadtxt(path, delimiter=',')
        return self.restrict(quantity)

    def get_quantity_at(self, vessel_id, quantity_name, position):
        quantity = self.get_quantity(vessel_id, quantity_name)
        quantity = quantity[:, int((quantity.shape[1]-1) * position)]
        return quantity

    def get_q_at(self, vessel_id, position):
        return self.get_quantity_at(vessel_id, 'q', position)

    def get_total_flow_at(self, vessel_ids, position):
        q_total = 0
        for vessel_id in vessel_ids:
            q = self.get_quantity_at(vessel_id, 'q', position)
            q = (q[1:] + q[:-1]) * 0.5
            q_total += np.sum(q*self.tau)
        return q_total

    def get_c_at(self, vessel_id, position):
        return self.get_quantity_at(vessel_id, 'c', position)
    
    def get_a_at(self, vessel_id, position):
        vessel = self.find_vessel(vessel_id)
        if ('a' in vessel['filepaths']):
            return self.get_quantity_at(vessel_id, 'a', position)
        else:
            return np.ones(self.t.shape) * vessel['A0']
    

def get_amount_of_substance(loader, vessel_ids, position):
    n_total = 0
    for vessel_id in vessel_ids:
        c = loader.get_c_at(vessel_id, position)
        a = loader.get_a_at(vessel_id, position)
        q = loader.get_q_at(vessel_id, position)
        tau = loader.tau
        n = c/a*q
        n = (n[1:] + n[:-1]) * 0.5
        n_total += np.cumsum(n*tau)
    return n_total


def get_q(loader, vessel_ids, position):
    t = loader.t
    q_total = np.zeros(t.shape) 
    for vessel_id in vessel_ids:
        q_total += loader.get_q_at(vessel_id, position)
    return q_total


filepath = 'tmp_transport/heart_to_breast_1d_solution_nl.json'
t_start = 23
t_end = 25 
loader = DataLoader(filepath, t_start, t_end)

head_vessel_ids = [13,16,4,5]
arm_vessel_ids = [14,15]
breast_vessel_ids = [35, 36]
torso_vessel_ids = [7]
position = 0.5

n_head = get_amount_of_substance(loader, head_vessel_ids, position)
n_arm = get_amount_of_substance(loader, arm_vessel_ids, position)
n_breast = get_amount_of_substance(loader, breast_vessel_ids, position)
n_torso = get_amount_of_substance(loader, torso_vessel_ids, position)

if False:
    plt.plot(loader.t_mid, n_head, label='head')
    plt.plot(loader.t_mid, n_arm, label='arm')
    plt.plot(loader.t_mid, n_breast, label='breast')
    plt.plot(loader.t_mid, n_torso, label='torso')
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel(r'$\int n_i dt$')
    plt.grid(True)
    plt.show()

if False:
    n_total = n_head + n_arm + n_breast + n_torso
    plt.plot(loader.t_mid, n_head/n_total, label='head')
    plt.plot(loader.t_mid, n_arm/n_total, label='arm')
    plt.plot(loader.t_mid, n_breast/n_total, label='breast')
    plt.plot(loader.t_mid, n_torso/n_total, label='torso')
    print(n_breast[-1]/n_total)
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel(r'$\frac{\int n_i dt}{\int n_{total} dt}$')
    plt.grid(True)
    plt.show()

q_head = get_q(loader, head_vessel_ids, position)
q_arm = get_q(loader, arm_vessel_ids, position)
q_breast = get_q(loader, breast_vessel_ids, position)
q_torso = get_q(loader, torso_vessel_ids, position)

q_total = q_head + q_arm + q_breast + q_torso

flow_head = loader.get_total_flow_at(head_vessel_ids, position)
flow_arm = loader.get_total_flow_at(arm_vessel_ids, position)
flow_breast = loader.get_total_flow_at(breast_vessel_ids, position)
flow_torso = loader.get_total_flow_at(torso_vessel_ids, position)
flow_total = flow_head + flow_arm + flow_breast + flow_torso
print('flow into head = {}\t % {}'.format(flow_head, flow_head/flow_total))
print('flow into arm = {}\t % {}'.format(flow_arm, flow_arm/flow_total))
print('flow into breast = {}\t % {}'.format(flow_breast, flow_breast/flow_total))
print('flow into torso = {}\t % {}'.format(flow_torso, flow_torso/flow_total))
print(flow_arm/flow_total)
print(flow_breast/flow_total)
print(flow_torso/flow_total)

if True:
    plt.plot(loader.t, q_torso, label='torso', color='tab:blue')
    plt.plot(loader.t, np.ones(loader.t.shape) * q_torso.mean(), '-.', label='torso', color='tab:blue')
    plt.plot(loader.t, q_head, label='head', color='tab:orange')
    plt.plot(loader.t, np.ones(loader.t.shape) * q_head.mean(), '-.', color='tab:orange')
    plt.plot(loader.t, q_arm, label='arm', color='tab:red')
    plt.plot(loader.t, np.ones(loader.t.shape) * q_arm.mean(), '-.', color='tab:red')
    plt.plot(loader.t, q_breast, label='breast', color='tab:green')
    plt.plot(loader.t, np.ones(loader.t.shape) * q_breast.mean(), '-.', color='tab:green')
    plt.legend()
    plt.ylim(-40, 100)
    plt.xlabel('t [s]')
    plt.ylabel(r'')
    plt.grid(True)
    plt.show()