import time
import json
from dataclasses import dataclass
import numpy as np
import dolfin as df


def read_mesh(mesh_filename):
    mesh = df.Mesh()
    f = df.XDMFFile(df.MPI.comm_world, mesh_filename)
    f.read(mesh)
    return mesh


def open_input_pressures(filepath):
    with open(filepath) as file:
        tip_pressures = json.loads(file.read())
        tip_pressures = {int(k):int(v) for k,v in tip_pressures.items()}
    return tip_pressures


def setup_subdomains(mesh, points, weights=None):
    subdomains = df.MeshFunction('size_t', mesh, mesh.topology().dim())
    subdomains.set_all(0)

    if weights is None:
        weights = np.ones(len(points))

    for cell in df.cells(mesh):
        mp = cell.midpoint()
        mp_array = mp.array().reshape((1, 3))
        diff = mp_array - points
        dist = np.linalg.norm(diff, axis=1)

        index = np.argmin(dist * weights)
        subdomains[cell] = index

    dx = df.Measure('dx', domain=mesh, subdomain_data=subdomains)

    return subdomains, dx


@dataclass
class Point:
    """Point mock data."""
    x: float
    y: float
    z: float


@dataclass
class VesselTipData:
    """Mock vessel tip coupling data."""
    p: Point
    vertex_id: int
    pressure: float
    concentration: float
    R2: float
    radius_first: float
    radius_last: float
    level: int


def vessel_tip_coupling_data_to_str(data_list):
    """A list of vessel tip data elements is converted into a string."""
    s = []
    for v in data_list:
        s.append('VesselTipData(')
        s.append('  p = Point(x={}, y={}, z={}),'.format(v.p.x, v.p.y, v.p.z))
        s.append('  vertex_id = {},'.format(v.vertex_id))
        s.append('  pressure = {},'.format(v.pressure))
        s.append('  concentration = {},'.format(v.concentration))
        s.append('  R2 = {},'.format(v.R2))
        s.append('  radius_first = {},'.format(v.radius_first))
        s.append('  radius_last = {},'.format(v.radius_last))
        s.append('  level = {}'.format(v.level))
        s.append('),')
    return '\n'.join(s)


class AverageQuantityWriter:
    def __init__(self, coupling_data):
        num_outlets = len(coupling_data)
        self.point_to_vertex_id = np.zeros(num_outlets)
        self.quantities = {}
        self.t = []
        for idx, vessel_tip in enumerate(coupling_data):
            self.point_to_vertex_id[idx] = vessel_tip.vertex_id
            self.quantities[vessel_tip.vertex_id] = { 'idx': idx, 'c_cap': [], 'c_tis': [], 'p_cap': [], 'p_tis': []} 

    def update(self, t, p_cap, p_tis, c_cap=None, c_tis=None):
        self.t.append(t)
        for vertex_id in p_cap.keys():
            self.quantities[vertex_id]['p_cap'].append(p_cap[vertex_id])
        for vertex_id in p_tis.keys():
            self.quantities[vertex_id]['p_tis'].append(p_tis[vertex_id])
        if c_cap is not None:
            for vertex_id in c_cap.keys():
                self.quantities[vertex_id]['c_cap'].append(c_cap[vertex_id])
        if c_tis is not None:
            for vertex_id in c_tis.keys():
                self.quantities[vertex_id]['c_tis'].append(c_tis[vertex_id])

    def write(self, filepath):
        with open(filepath, 'w') as file:
            file.write(json.dumps({'t': self.t, 'quantities': self.quantities}))


def viscosity_bloodplasma(r: float):
    from math import pow, exp

    # radius to diameter
    d = 2 * r

    # dimensionless diameter:  d / 1 micrometer
    d_tilde = (1e-2 * d) / 1e-6

    mu_0p45 = 6.0 * exp(-0.085 * d_tilde) + 3.2 - 2.44 * exp(-0.06 * pow(d_tilde, 0.645))

    # viscosity of blood plasma [Pa s]
    mu_p = 1e-3

    # the blood viscosity in [Pa s]
    # Note: we set the hematocrit to H = 0.45
    #       thus ((1 - H)^C -1)/ ((1 - 0.45)^C - 1) simplifies to 1
    mu_bl = mu_p * (1 + (mu_0p45 - 1) * pow(d_tilde / (d_tilde - 1.1), 2)) * pow(d_tilde / (d_tilde - 1.1), 2)

    # we convert to the cm units:
     # [Pa s] = 10 [Ba s]
    mu_bl_cm = mu_bl * 10

    return mu_bl_cm


class StopWatch:
    def __init__(self):
        self.start_time = None
        self.elapsed_times = []
        self.num_steps = []

    def start(self):
        assert self.start_time is None
        self.start_time = time.time()

    def end(self, num_steps=1):
        assert self.start_time is not None
        elapsed = (time.time() - self.start_time) / num_steps
        self.elapsed_times.append(elapsed)
        self.num_steps.append(num_steps)
        self.start_time = None
        return elapsed

    def average(self):
        return np.array(self.elapsed_times).mean()

    def total(self):
        return np.sum(np.array(self.elapsed_times) * np.array(self.num_steps))

    def __str__(self):
        return f'avg = {self.average()}, total = {self.total()}, number of data points = {len(self.elapsed_times)}'
