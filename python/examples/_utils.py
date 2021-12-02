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
