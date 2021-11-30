import json
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
