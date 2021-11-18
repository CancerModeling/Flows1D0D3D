import os
import flows1d0d3d as f
import dolfin as df
import numpy as np


# constants:
rho_c = 997 * 1e-3
rho_t = 1060 * 1e-3
K_c = 1e-13 * 1e4 # open?
K_t = 1e-18 * 1e4 # Tobias paper homogenisierung
mu_c = 0.0805438
mu_t = 0.79722 * 1e-3 * 10
L_cv = 1e-8  # Ottesen, Olufsen, Larsen, p. 153
p_ven = 10 * 1333

L_ct = 1e-6 #?
S_ct = 1.   #?

sigma = 0.  #?
pi_int = 6.6e3
pi_bl = 3.3e4

L_tl = 1e-8 #?
p_lym = 1333.2

D_c = 1. #?
D_t = 1. #?
lambda_t = 1e-8 #??


def read_mesh(mesh_filename):
    mesh = df.Mesh()
    f = df.XDMFFile(df.MPI.comm_world, mesh_filename)
    f.read(mesh)
    return mesh


def setup_subdomains(mesh, points):
    subdomains = df.MeshFunction('size_t', mesh, mesh.topology().dim())
    subdomains.set_all(0)

    for cell in df.cells(mesh):
        mp = cell.midpoint()
        mp_array = mp.array().reshape((1,3))
        diff = mp_array - points
        dist = np.linalg.norm(diff, axis=1)
        index = np.argmin(dist)
        subdomains[cell] = index

    dx = df.Measure('dx', domain=mesh, subdomain_data=subdomains)

    return subdomains, dx


class PressureSolver:
    def __init__(self, mesh, output_folder, vessel_tip_pressures):
        self.mesh = mesh
        self._setup_vessel_tip_pressures(vessel_tip_pressures)
        self._setup_subdomains()
        self._setup_problem()
        self._setup_solver()
        self.file_p_cap = df.File(os.path.join(output_folder, "p_cap.pvd"), "compressed")
        self.file_p_tis = df.File(os.path.join(output_folder, "p_tis.pvd"), "compressed")

    def _setup_vessel_tip_pressures(self, list_vessel_tip_pressures):
        num_outlets = len(list_vessel_tip_pressures)
        self.points = np.zeros((num_outlets, 3))
        self.point_to_vertex_id = np.zeros(num_outlets)
        self.pressures = np.zeros(num_outlets)
        self.average_pressures = np.zeros(num_outlets)
        self.R2 = np.zeros(num_outlets)
        self.level = np.zeros(num_outlets)
        for idx, vessel_tip in enumerate(list_vessel_tip_pressures):
            p = vessel_tip.p
            self.points[idx, :] = np.array([p.x, p.y, p.z])
            self.point_to_vertex_id[idx] = vessel_tip.vertex_id
            self.pressures[idx] = vessel_tip.pressure
            self.R2[idx] = vessel_tip.R2
            self.level[idx] = vessel_tip.level

    def _setup_subdomains(self):
        self.subdomains, self.dx = setup_subdomains(self.mesh, self.points)

    def _setup_problem(self):
        P1El = df.FiniteElement('P', self.mesh.ufl_cell(), 1)
        MixedEl = df.MixedElement([P1El, P1El])
        V = df.FunctionSpace(self.mesh, MixedEl)

        phi_c, phi_t = df.TrialFunctions(V)
        psi_c, psi_t = df.TestFunctions(V)

        self.current = df.Function(V, name='u') if not hasattr(self, 'current') else self.current

        volumes = []
        for k in range(len(self.pressures)):
            volume = df.assemble(df.Constant(1) * self.dx(k))
            volumes.append(volume)
        self.volumes = volumes

        self.update_average_pressures()

        J = 0
        J += df.inner( df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c) ) * df.dx
        J += df.inner( df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t) ) * df.dx

        # q_cv
        J -= df.Constant( rho_c * L_cv ) * (df.Constant(p_ven) - phi_c) * psi_c * df.dx

        # q_ca
        print(self.average_pressures)
        for k in range(len(self.pressures)):
            J -= df.Constant(rho_c * 2**(self.level[k]-1) / self.R2[k] / volumes[k]) * (df.Constant(self.pressures[k]) - phi_c) * psi_c * self.dx(k)
            #J -= df.Constant(rho_c * 2**(self.level[k]-1) / self.R2[k] / volumes[k]) * (df.Constant(self.pressures[k]) - df.Constant(self.average_pressures[k])) * psi_c * self.dx(k)

        # q_ct
        J -= df.Constant(rho_t * L_ct * S_ct ) * (phi_t-phi_c) * psi_c * df.dx
        J -= df.Constant(rho_t * L_ct * S_ct * sigma * (pi_int-pi_bl)) * psi_c * df.dx

        # q_tc
        J -= df.Constant(rho_t * L_ct * S_ct ) * (phi_c-phi_t) * psi_t * df.dx
        J -= df.Constant(rho_t * L_ct * S_ct * sigma * (pi_bl-pi_int)) * psi_t * df.dx

        # q_tl
        #J -= df.Constant(rho_t * L_tl) * (p_lym - phi_t) * psi_t * df.dx
        J -= df.Constant(rho_t * L_tl) * (p_lym - phi_t) * psi_t * df.dx

        # setup diagonal preconditioner
        P = 0
        P += df.inner( df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        P += df.inner( df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx
        for k in range(len(self.pressures)):
            P -= df.Constant( rho_c * 2**(self.level[k]-1) / self.R2[k] / volumes[k]) * (- phi_c) * psi_c * self.dx(k)
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_c) * psi_c * df.dx
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_t) * psi_t * df.dx
        P += df.Constant(0) * psi_t * df.dx

        self.J = J
        self.P = P
        self.bcs = []

    def _setup_solver(self):
        A, b = df.assemble_system(df.lhs(self.J), df.rhs(self.J), self.bcs)
        P, _ = df.assemble_system(df.lhs(self.P), df.rhs(self.P), self.bcs)
        self.solver = df.KrylovSolver('gmres', 'amg')
        self.solver.set_operators(A, P)

    def solve(self):
        A, b = df.assemble_system(df.lhs(self.J), df.rhs(self.J), self.bcs)
        self.solver.solve(self.current.vector(), b)

    def write_subdomains(self, output_folder):
        f_sd = df.XDMFFile(df.MPI.comm_world, os.path.join(output_folder, 'subdomains.xdmf'))
        f_sd.write(self.subdomains)

    def write_solution(self, t):
        self.file_p_cap << (self.current.split()[0], t)
        self.file_p_tis << (self.current.split()[1], t)

    def update_vessel_tip_pressures(self, list_vessel_tip_pressures):
        for idx, vessel_tip in enumerate(list_vessel_tip_pressures):
            assert vessel_tip.vertex_id == self.point_to_vertex_id[idx]
            self.pressures[idx] = vessel_tip.pressure
        # we have to reassemble the system:
        self._setup_problem()

    def get_pressures(self):
        new_data = {}
        for k in range(len(self.pressures)):
            integrated_pressure = df.assemble(self.current.split()[0] * self.dx(k))
            average_pressure = integrated_pressure / self.volumes[k]
            new_data[int(self.point_to_vertex_id[k])] = average_pressure
        return new_data

    def update_average_pressures(self):
        for k in range(len(self.pressures)):
            integrated_pressure = df.assemble(self.current.split()[0] * self.dx(k))
            self.average_pressures[k] = integrated_pressure / self.volumes[k]


def run():
    data_folder = '../../data'
    output_folder = './tmp'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau = 1. / 2**16
    tau_out = 1. / 2**6
    t_end = 20.
    t = 0
    t_coup_start = 0.
    t_coup_interval_transport = 1000

    solver1d = f.HeartToBreast1DSolver()
    solver1d.set_output_folder('./tmp')
    #s.set_path_inflow_pressures(os.path.join(data_folder, ""));
    solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"));
    solver1d.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"));
    solver1d.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"));
    solver1d.setup(degree, tau)

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    mesh = read_mesh(mesh_3d_filename)
    solver3d = PressureSolver(mesh, output_folder, vessel_tip_pressures)
    #solver3d.current.interpolate(df.Constant((33. * 1333, 33. * 1333)))
    solver3d.write_solution(0.)

    solver3d_transport = TransportSolver(mesh, output_folder, solver3d.current, vessel_tip_pressures)
    solver3d_transport.write_solution(0.)

    for i in range(int(t_end / tau_out)):
        print ('iter = {}, t = {}'.format(i, t))
        t = solver1d.solve_flow(tau, t, int(tau_out / tau))
        solver1d.write_output(t)
        if t > t_coup_start:
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            print('start solving pressures')
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
            solver3d.solve()
            solver3d.write_solution(t)
            print('endsolving pressures')
            new_pressures = solver3d.get_pressures()

            solver1d.update_vessel_tip_pressures(new_pressures)
            solver3d_transport.update_average_pressures(new_pressures)

            # solve the transport:
            if i % t_coup_interval_transport == 0:
                print('start solving transport')
                # REMOVE ME START
                for x in vessel_tip_pressures:
                    x.concentration = 1.
                # REMOVE ME END
                solver3d_transport.update_vessel_tip_concentrations(vessel_tip_pressures)

                solver3d_transport.solve()
                solver3d_transport.write_solution(t)
                print('end solving transport')

            max_value = solver3d.current.vector().max()
            min_value = solver3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('max = {}, min = {}'.format(max_value, min_value))


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
