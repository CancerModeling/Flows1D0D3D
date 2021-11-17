import math
import os
import flows1d0d3d as f
import dolfin as df
import numpy as np
from dataclasses import dataclass
from typing import Tuple

# constants:
# constants:
rho_c = 997 * 1e-3
rho_t = 1060 * 1e-3
K_c = 1e-13 * 1e4  # open?
# K_c = 1e-12 * 1e4  # wildest guess!!!
K_t = 1e-18 * 1e4  # Tobias paper homogenisierung
mu_c = 0.0805438
mu_t = 0.79722 * 1e-3 * 10
# L_cv = 1e-7  # guess
# L_cv = 0.667 / 1333  # Ottesen, Olufsen, Larsen, p. 153
# L_cv = 4.641e-7  # Liverpool!
L_cv = 1e-8  # guess
#L_cv = 2.2e-3 * 10-4  # p.142 carlo
p_ven = 10 * 1333

#L_ct = 1e-6  # ?
L_ct = 1e-11  # -> L_cap Tobias paper
L_ct = 1e-9  # guess

r_c_bar = 3.375e-4
l_c_bar = 0.06
S_ct = 460 * r_c_bar * l_c_bar * 2 * np.pi / 1e-3  # ?

sigma = 0.  # ?
pi_int = 6.6e3
pi_bl = 3.3e4

L_tl = 1e-8  # ?
p_lym = 1333.2

D_c = 1.  # ?
D_t = 1.  # ?
lambda_t = 1.  # ??

# Liverpool values
#K_c = 4.28e-4 * 1e-6
#mu_c = 1.


def read_mesh(mesh_filename):
    mesh = df.Mesh()
    f = df.XDMFFile(df.MPI.comm_world, mesh_filename)
    f.read(mesh)
    return mesh


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


class PressureSolver:
    def __init__(self, mesh, output_folder, vessel_tip_pressures):
        self.mesh = mesh
        self._setup_vessel_tip_pressures(vessel_tip_pressures)
        self._setup_subdomains()
        self._setup_function_spaces()
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
        self.radii = np.zeros(num_outlets)
        self.R2 = np.zeros(num_outlets)
        self.level = np.zeros(num_outlets)
        self.total_flows = np.zeros(num_outlets)
        for idx, vessel_tip in enumerate(list_vessel_tip_pressures):
            p = vessel_tip.p
            self.points[idx, :] = np.array([p.x, p.y, p.z])
            self.point_to_vertex_id[idx] = vessel_tip.vertex_id
            self.pressures[idx] = vessel_tip.pressure
            self.R2[idx] = vessel_tip.R2
            self.radii[idx] = vessel_tip.radius
            self.level[idx] = vessel_tip.level
            self.total_flows[idx] = 2**(vessel_tip.level-1) * (vessel_tip.pressure - 30 * 1333) / vessel_tip.R2

    def _setup_subdomains(self):
        weights = 1./self.radii**2
        weights = np.ones(self.radii.shape)
        self.subdomains, self.dx = setup_subdomains(self.mesh, self.points, weights)

    def _setup_function_spaces(self):
        PEl = df.FiniteElement('P', self.mesh.ufl_cell(), 1)
        REl = df.FiniteElement('R', self.mesh.ufl_cell(), 0)
        SpatialEl = [PEl, PEl]
        elements = SpatialEl + [REl for idx in range(len(self.points))]
        self.first_multiplier_index = len(SpatialEl)
        MEl = df.MixedElement(elements)
        self.V = df.FunctionSpace(self.mesh, MEl)

    def _setup_problem(self):
        trial_functions = df.TrialFunctions(self.V)
        test_functions = df.TestFunctions(self.V)

        phi_c = trial_functions[0]
        phi_t = trial_functions[1]
        psi_c = test_functions[0]
        psi_t = test_functions[1]

        llambda = trial_functions[self.first_multiplier_index:]
        mu = test_functions[self.first_multiplier_index:]

        self.current = df.Function(self.V, name='u') if not hasattr(self, 'current') else self.current

        #phi_t = df.Constant(20 * 1333)

        volumes = []
        for k in range(len(self.pressures)):
            volume = df.assemble(df.Constant(1) * self.dx(k))
            volumes.append(volume)
        self.volumes = volumes

        self.update_average_pressures()

        J = 0
        J += df.inner(df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        J += df.inner(df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx

        coeff_ca = []
        coeff_cv = []
        list_L_ct = []
        list_L_tl = []
        for k in range(len(self.pressures)):
            coeff_ca.append(rho_c * 2 ** (self.level[k] - 1) / self.R2[k] / volumes[k])
            #coeff_cv.append(rho_c * L_cv)
            #L_cv = 0.95 * self.total_flows[k] / (self.volumes[k] * p_ven)
            #L_ct = 0.05 * self.total_flows[k] / (self.volumes[k] * 16 * 1333 * S_ct)
            #L_tl = 0.05 * self.total_flows[k] / (self.volumes[k] * p_lym)
            coeff_cv.append(rho_c * L_cv)
            list_L_ct.append(L_ct)
            list_L_tl.append(L_tl)
        # mean value:
        coeff_ca_mean = np.array(coeff_ca).mean()
        coeff_ca = np.ones(len(coeff_ca)) * coeff_ca_mean
        '''
        total_volume = np.sum(volumes)
        for k in range(len(self.pressures)):
            coeff_ca.append(rho_c * 2 ** (self.level[k] - 1) / self.R2[k] / total_volume)
            coeff_cv.append(rho_c * L_cv)
        coeff_ca_mean = np.array(coeff_ca).mean()
        coeff_ca = np.ones(len(coeff_ca)) * coeff_ca_mean
        '''

        # llambda[k] = - 1/Omega_k int_Omega[k] p_c dx
        for k in range(len(self.pressures)):
            alpha = df.Constant(coeff_ca[k] + coeff_cv[k])
            J += - alpha * llambda[k] * mu[k] * self.dx(k) + alpha * phi_c * mu[k] * self.dx(k)

        # q_cv
        for k in range(len(self.pressures)):
            #J -= df.Constant(coeff_cv[k]) * (df.Constant(p_ven) - llambda[k]) * psi_c * self.dx(k)
            J -= df.Constant(coeff_cv[k]) * (df.Constant(p_ven) - phi_c) * psi_c * self.dx(k)

        # q_ca
        for k in range(len(self.pressures)):
            J -= df.Constant(coeff_ca[k]) * (df.Constant(self.pressures[k]) - llambda[k]) * psi_c * self.dx(k)

        # q_ct
        for k in range(len(self.pressures)):
            J -= df.Constant(rho_t * list_L_ct[k] * S_ct) * (phi_t - phi_c) * psi_c * self.dx(k)
            J -= df.Constant(rho_t * list_L_ct[k] * S_ct * sigma * (pi_bl - pi_int)) * psi_c * self.dx(k)

        # q_tc
        for k in range(len(self.pressures)):
            J -= df.Constant(rho_t * list_L_ct[k] * S_ct) * (phi_c - phi_t) * psi_t * self.dx(k)
            J -= df.Constant(rho_t * list_L_ct[k] * S_ct * sigma * (pi_int - pi_bl)) * psi_t * self.dx(k)

        # q_tl
        for k in range(len(self.pressures)):
            J -= df.Constant(rho_t * list_L_tl[k]) * (p_lym - phi_t) * psi_t * self.dx(k)

        # setup diagonal preconditioner
        P = 0
        P += df.inner(df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        P += df.inner(df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx
        for k in range(len(self.pressures)):
            P -= df.Constant(rho_t * list_L_ct[k] * S_ct) * (-phi_c) * psi_c * self.dx(k)
            P -= df.Constant(rho_t * list_L_ct[k] * S_ct) * (-phi_t) * psi_t * self.dx(k)
            P -= df.Constant(rho_t * list_L_tl[k]) * (- phi_t) * psi_t * self.dx(k)
        P += df.Constant(0) * psi_t * df.dx
        for k in range(len(self.pressures)):
            alpha = df.Constant(coeff_ca[k] + coeff_cv[k])
            P += - alpha * llambda[k] * mu[k] * self.dx(k)

        self.J = J
        self.P = P
        self.bcs = []

    def _setup_solver(self):
        self.A, self.b = df.assemble_system(df.lhs(self.J), df.rhs(self.J), self.bcs)
        self.P, _ = df.assemble_system(df.lhs(self.P), df.rhs(self.P), self.bcs)
        # self.solver = df.KrylovSolver('gmres', 'jacobi')
        self.solver = df.KrylovSolver('gmres', 'amg')
        # print(self.solver.parameters.keys())
        self.solver.parameters['monitor_convergence'] = True
        self.solver.parameters['nonzero_initial_guess'] = True
        # self.solver.parameters['absolute_tolerance'] = 1e-12
        # self.solver.parameters['relative_tolerance'] = 1e-12
        # self.solver.parameters['maximum_iterations'] = 10000
        # self.solver = df.KrylovSolver('gmres', 'jacobi')
        self.solver.set_operators(self.A, self.P)
        #self.solver.set_operator(A)

    def solve(self):
        system_assembler = df.SystemAssembler(df.lhs(self.J), df.rhs(self.J), self.bcs)
        system_assembler.assemble(self.b)
        self.solver.solve(self.current.vector(), self.b)

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

    def get_pressures_from_multiplier(self):
        new_data = {}
        for k in range(len(self.pressures)):
            average_pressure = self.current.split()[self.first_multiplier_index+k].vector().get_local()[0]
            new_data[int(self.point_to_vertex_id[k])] = average_pressure
        return new_data

    def update_average_pressures(self):
        for k in range(len(self.pressures)):
            integrated_pressure = df.assemble(self.current.split()[0] * self.dx(k))
            self.average_pressures[k] = integrated_pressure / self.volumes[k]


@dataclass
class MockPoint:
    x: float
    y: float
    z: float


@dataclass
class MockVesselTipPressures:
    p: Tuple[float, float, float]
    vertex_id: int
    pressure: float
    concentration: float
    R2: float
    level: float


def run():
    data_folder = '../../data'
    output_folder = './tmp'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau_out = 1. / 2 ** 3
    #tau_out = 1.
    #tau_coup = 2.
    tau_coup = tau_out
    t_end = 80.
    t = 0
    t_coup_start = 2.
    use_fully_coupled = False

    if use_fully_coupled:
        tau = 1. / 2 ** 16
        solver1d = f.FullyCoupledHeartToBreast1DSolver()
        solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"))
        solver1d.set_path_linear_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))
        solver1d.set_path_coupling_conditions(os.path.join(data_folder, "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"))
    else:
        tau = 1. / 2 ** 4
        solver1d = f.LinearizedHeartToBreast1DSolver()
        solver1d.set_path_inflow_pressures(os.path.join(data_folder, "1d-input-pressures/from-33-vessels-with-small-extension.json"))
        solver1d.set_path_geometry(os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"))
    solver1d.set_output_folder(output_folder)
    solver1d.setup(degree, tau)

    coupling_interval_0d3d = int(round(tau_coup / tau_out))

    t = solver1d.solve_flow(tau, t, int(6 / tau))

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    mesh = read_mesh(mesh_3d_filename)
    solver3d = PressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    solver3d.write_solution(0.)
    solver3d.write_subdomains(output_folder)

    for i in range(0, int(t_end / tau_out)):
        print('iter = {}, t = {}'.format(i, t))
        t = solver1d.solve_flow(tau, t, int(tau_out / tau))
        solver1d.write_output(t)
        if (t > t_coup_start) and (i % coupling_interval_0d3d == 0):
            vessel_tip_pressures = solver1d.get_vessel_tip_pressures()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('start solving 3D pressures')
            solver3d.update_vessel_tip_pressures(vessel_tip_pressures)
            solver3d.solve()
            solver3d.write_solution(t)
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('end solving 3D pressures')
            new_pressures = solver3d.get_pressures()
            # new_pressures_comp = solver3d.get_pressures_from_multiplier()
            # print(new_pressures)
            # print(new_pressures_comp)
            solver1d.update_vessel_tip_pressures(new_pressures)

            max_value = solver3d.current.vector().max()
            min_value = solver3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))


if __name__ == '__main__':
    f.initialize()
    run()
    f.finalize()
