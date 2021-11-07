import os
import flows1d0d3d as f
import dolfin as df
import numpy as np


# constants:
rho = 997 * 1e-3
K = 1e-13 * 1e4
mu = 4e-3 * 10


class SimpleCapillaryPressureSolver:
    def __init__(self, mesh_filename, output_folder, vessel_tip_pressures):
        self._setup_mesh(mesh_filename)
        self._setup_vessel_tip_pressures(vessel_tip_pressures)
        self._setup_subdomains()
        self._setup_problem()
        self.file_p_cap = df.File(os.path.join(output_folder, "p_cap.pvd"), "compressed")

    def _setup_mesh(self, mesh_filename):
        self.mesh = df.Mesh()
        f = df.XDMFFile(df.MPI.comm_world, mesh_filename)
        f.read(self.mesh)

    def _setup_vessel_tip_pressures(self, list_vessel_tip_pressures):
        num_outlets = len(list_vessel_tip_pressures)
        self.points = np.zeros((num_outlets, 3))
        self.point_to_vertex_id = np.zeros(num_outlets)
        self.pressures = np.zeros(num_outlets)
        self.R2 = np.zeros(num_outlets)
        self.level = np.zeros(num_outlets)
        for idx, vessel_tip in enumerate(list_vessel_tip_pressures):
            p = vessel_tip.p
            self.points[idx,:] = np.array([p.x, p.y, p.z])
            self.point_to_vertex_id[idx] = vessel_tip.vertex_id
            self.pressures[idx] = vessel_tip.pressure
            self.R2[idx] = vessel_tip.R2
            self.level[idx] = vessel_tip.level

    def _setup_subdomains(self):
        self.subdomains = df.MeshFunction('size_t', self.mesh, self.mesh.topology().dim())
        self.subdomains.set_all(0)

        for cell in df.cells(self.mesh):
            mp = cell.midpoint()
            mp_array = mp.array().reshape((1,3))
            diff = mp_array - self.points
            dist = np.linalg.norm(diff, axis=1)
            index = np.argmin(dist)
            self.subdomains[cell] = index

        self.dx = df.Measure('dx', domain=self.mesh, subdomain_data=self.subdomains)

    def _setup_problem(self):
        V = df.FunctionSpace(self.mesh, 'P', 1)
        phi_c = df.TrialFunction(V)
        psi_c = df.TestFunction(V)

        J = 0
        J += df.inner( df.Constant(rho * K / mu) * df.grad(phi_c), df.grad(psi_c) ) * df.dx

        volumes = []
        for k in range(len(self.pressures)):
            volume = df.assemble(df.Constant(1) * self.dx(k))
            volumes.append(volume)

        R2_ven = 1e6
        p_ven = 10 * 1333

        # rhs contributions
        for k in range(len(self.pressures)):
            J -= df.Constant( 2**(self.level[k]-1) / self.R2[k] / volumes[k] ) * df.Constant(self.pressures[k]) * psi_c * self.dx(k)
            J -= df.Constant( 1. / R2_ven / volumes[k] ) * df.Constant(p_ven) * psi_c * self.dx(k)

        # lhs contributions
        for k in range(len(self.pressures)):
            J += df.Constant( 2**(self.level[k]-1) / self.R2[k] / volumes[k] ) * phi_c * psi_c * self.dx(k)
            J += df.Constant( 1. / R2_ven / volumes[k] ) * phi_c * psi_c * self.dx(k)

        self.J = J
        self.bcs = []
        self.current = df.Function(V, name='u')

    def solve(self):
        df.solve(df.lhs(self.J) == df.rhs(self.J), self.current, self.bcs)

    def write_subdomains(self, output_folder):
        f_sd = df.XDMFFile(df.MPI.comm_world, os.path.join(output_folder, 'subdomains.xdmf'))
        f_sd.write(self.subdomains)

    def write_solution(self, t):
        self.file_p_cap << (self.current, t)

    def update_vessel_tip_pressures(self, list_vessel_tip_pressures):
        for idx, vessel_tip in enumerate(list_vessel_tip_pressures):
            assert vessel_tip.vertex_id == self.point_to_vertex_id[idx]
            self.pressures[idx] = vessel_tip.pressure
        # we have to reassemble the system:
        self._setup_problem()

    def get_pressures(self):
        pass


def run():
    data_folder = '../../data'
    output_folder = './tmp'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau = 1. / 2**16
    tau_out = 1. / 2**6
    t_end = 10.
    t = 0

    s = f.HeartToBreast1DSolver()
    s.set_output_folder('./tmp')
    #s.set_path_inflow_pressures(os.path.join(data_folder, ""));
    s.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"));
    s.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"));
    s.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"));
    s.setup(degree, tau)

    vessel_tip_pressures = s.get_vessel_tip_pressures()

    # new_pressures = {}
    # for tip_pressures in vessel_tip_pressures:
    #     new_pressures[tip_pressures.vertex_id] = 22. * 1333

    # s.update_vessel_tip_pressures(new_pressures)

    solver3d = SimpleCapillaryPressureSolver(mesh_3d_filename, output_folder, vessel_tip_pressures)
    solver3d.write_solution(0.)

    for i in range(int(t_end / tau_out)):
        print ('iter = {}, t = {}'.format(i, t))
        t = s.solve_flow(tau, t, int(tau_out / tau))
        s.write_output(t)
        vessel_tip_pressures = s.get_vessel_tip_pressures()
        solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
        solver3d.solve()
        solver3d.write_solution(t)
        max_value = solver3d.current.vector().max()
        min_value = solver3d.current.vector().min()
        if df.MPI.rank(df.MPI.comm_world) == 0:
            print('max = {}, min = {}'.format(max_value, min_value))
        #for tip_pressure in vessel_tip_pressures:
        #    print(tip_pressure)


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
