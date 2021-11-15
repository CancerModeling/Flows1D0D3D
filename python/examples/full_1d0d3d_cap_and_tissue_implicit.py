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


def setup_subdomains(mesh, points):
    subdomains = df.MeshFunction('size_t', mesh, mesh.topology().dim())
    subdomains.set_all(0)

    for cell in df.cells(mesh):
        mp = cell.midpoint()
        mp_array = mp.array().reshape((1, 3))
        diff = mp_array - points
        dist = np.linalg.norm(diff, axis=1)

        weight = np.ones(len(dist)) # propto r
        np.argmin(dist * weight)

        index = np.argmin(dist)
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
        for k in range(len(self.pressures)):
            coeff_ca.append(rho_c * 2 ** (self.level[k] - 1) / self.R2[k] / volumes[k])
            coeff_cv.append(rho_c * L_cv)
        # mean value:
        #coeff_ca_mean = np.array(coeff_ca).mean()
        #coeff_ca = np.ones(len(coeff_ca)) * coeff_ca_mean

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
        J -= df.Constant(rho_t * L_ct * S_ct) * (phi_t - phi_c) * psi_c * df.dx
        J -= df.Constant(rho_t * L_ct * S_ct * sigma * (pi_bl - pi_int)) * psi_c * df.dx

        # q_tc
        J -= df.Constant(rho_t * L_ct * S_ct) * (phi_c - phi_t) * psi_t * df.dx
        J -= df.Constant(rho_t * L_ct * S_ct * sigma * (pi_int - pi_bl)) * psi_t * df.dx

        # q_tl
        J -= df.Constant(rho_t * L_tl) * (p_lym - phi_t) * psi_t * df.dx

        # setup diagonal preconditioner
        P = 0
        P += df.inner(df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        P += df.inner(df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_c) * psi_c * df.dx
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_t) * psi_t * df.dx
        P -= df.Constant(rho_t * L_tl) * (- phi_t) * psi_t * df.dx
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
        #self.solver = df.KrylovSolver('gmres', 'jacobi')
        self.solver = df.KrylovSolver('gmres', 'amg')
        #print(self.solver.parameters.keys())
        self.solver.parameters['monitor_convergence'] = True
        #self.solver.parameters['absolute_tolerance'] = 1e-12
        #self.solver.parameters['relative_tolerance'] = 1e-12
        #self.solver.parameters['maximum_iterations'] = 10000
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
    tau = 1. / 2 ** 16
    tau_out = 1. / 2 ** 3
    t_end = 20.
    t = 0
    t_coup_start = 2.

    solver1d = f.HeartToBreast1DSolver()
    solver1d.set_output_folder('./tmp')
    # s.set_path_inflow_pressures(os.path.join(data_folder, ""));
    solver1d.set_path_nonlinear_geometry(os.path.join(data_folder, "1d-meshes/33-vessels-with-small-extension.json"));
    solver1d.set_path_linear_geometry(
        os.path.join(data_folder, "1d-meshes/coarse-breast-geometry-with-extension.json"));
    solver1d.set_path_coupling_conditions(os.path.join(data_folder,
                                                       "1d-coupling/couple-33-vessels-with-small-extension-to-coarse-breast-geometry-with-extension.json"));
    solver1d.setup(degree, tau)

    vessel_tip_pressures = solver1d.get_vessel_tip_pressures()

    # print('[')
    # for v in vessel_tip_pressures:
    #     print(
    #         "MockVesselTipPressures(p=MockPoint(x={:.3}, y={:.3}, z={:.3}), vertex_id={}, pressure={}, concentration={}, R2={}, level={}),".format(
    #             v.p.x, v.p.y, v.p.z, v.vertex_id, 33 * 1333, 0, v.R2, v.level
    #         ))
    # print(']')
    # return 0

    mesh = read_mesh(mesh_3d_filename)
    solver3d = PressureSolver(mesh, output_folder, vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    solver3d.write_solution(0.)

    for i in range(int(t_end / tau_out)):
        print('iter = {}, t = {}'.format(i, t))
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

            max_value = solver3d.current.vector().max()
            min_value = solver3d.current.vector().min()
            if df.MPI.rank(df.MPI.comm_world) == 0:
                print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))

in_p = 33 * 1333

mock_vessel_tip_pressures = [
    MockVesselTipPressures(p=MockPoint(x=3.57, y=7.65, z=4.37), vertex_id=2, pressure=in_p, concentration=0,
                           R2=120495648.5176453, level=12),
    MockVesselTipPressures(p=MockPoint(x=2.21, y=8.29, z=5.57), vertex_id=3, pressure=in_p, concentration=0,
                           R2=129807914.06789602, level=12),
    MockVesselTipPressures(p=MockPoint(x=4.93, y=9.81, z=5.49), vertex_id=7, pressure=in_p, concentration=0,
                           R2=122967826.3673249, level=12),
    MockVesselTipPressures(p=MockPoint(x=3.33, y=8.37, z=4.61), vertex_id=11, pressure=in_p, concentration=0,
                           R2=160414934.95677227, level=11),
    MockVesselTipPressures(p=MockPoint(x=2.21, y=9.81, z=14.7), vertex_id=23, pressure=in_p, concentration=0,
                           R2=125046304.98609628, level=12),
    MockVesselTipPressures(p=MockPoint(x=1.01, y=9.89, z=11.9), vertex_id=27, pressure=in_p, concentration=0,
                           R2=137279441.42421895, level=12),
    MockVesselTipPressures(p=MockPoint(x=1.01, y=10.1, z=13.8), vertex_id=28, pressure=in_p, concentration=0,
                           R2=146632792.38299593, level=12),
    MockVesselTipPressures(p=MockPoint(x=7.25, y=9.89, z=13.6), vertex_id=39, pressure=in_p, concentration=0,
                           R2=137941990.30448535, level=11),
    MockVesselTipPressures(p=MockPoint(x=4.29, y=9.57, z=15.5), vertex_id=40, pressure=in_p, concentration=0,
                           R2=150775781.71597233, level=11),
    MockVesselTipPressures(p=MockPoint(x=2.29, y=10.1, z=0.693), vertex_id=56, pressure=in_p, concentration=0,
                           R2=124923507.01379998, level=11),
    MockVesselTipPressures(p=MockPoint(x=0.773, y=8.93, z=3.57), vertex_id=60, pressure=in_p, concentration=0,
                           R2=176278128.67142257, level=12),
    MockVesselTipPressures(p=MockPoint(x=5.81, y=4.69, z=9.49), vertex_id=76, pressure=in_p, concentration=0,
                           R2=163053931.97067818, level=11),
    MockVesselTipPressures(p=MockPoint(x=2.77, y=9.01, z=11.7), vertex_id=87, pressure=in_p, concentration=0,
                           R2=154092211.81762356, level=13),
    MockVesselTipPressures(p=MockPoint(x=2.45, y=9.89, z=11.8), vertex_id=92, pressure=in_p, concentration=0,
                           R2=163470061.43666664, level=11),
    MockVesselTipPressures(p=MockPoint(x=5.89, y=8.77, z=14.1), vertex_id=100, pressure=in_p, concentration=0,
                           R2=119275602.28663285, level=11),
    MockVesselTipPressures(p=MockPoint(x=2.77, y=10.1, z=0.613), vertex_id=102, pressure=in_p, concentration=0,
                           R2=196577356.68226808, level=12),
    MockVesselTipPressures(p=MockPoint(x=3.65, y=2.69, z=2.13), vertex_id=104, pressure=in_p, concentration=0,
                           R2=163470061.43666664, level=11),
    MockVesselTipPressures(p=MockPoint(x=7.01, y=9.25, z=3.33), vertex_id=106, pressure=in_p, concentration=0,
                           R2=168610289.46967137, level=12),
    MockVesselTipPressures(p=MockPoint(x=4.69, y=1.49, z=5.49), vertex_id=111, pressure=in_p, concentration=0,
                           R2=162225664.565375, level=11),
    MockVesselTipPressures(p=MockPoint(x=8.21, y=3.81, z=7.89), vertex_id=112, pressure=in_p, concentration=0,
                           R2=131139462.43191028, level=11),
    MockVesselTipPressures(p=MockPoint(x=5.09, y=5.65, z=10.7), vertex_id=116, pressure=in_p, concentration=0,
                           R2=171210843.470033, level=12),
    MockVesselTipPressures(p=MockPoint(x=7.41, y=9.25, z=11.7), vertex_id=117, pressure=in_p, concentration=0,
                           R2=163470061.43666664, level=11),
    MockVesselTipPressures(p=MockPoint(x=8.29, y=9.97, z=12.5), vertex_id=120, pressure=in_p, concentration=0,
                           R2=148851838.98346952, level=12),
    MockVesselTipPressures(p=MockPoint(x=9.17, y=7.81, z=4.21), vertex_id=126, pressure=in_p, concentration=0,
                           R2=147730112.62774938, level=11),
    MockVesselTipPressures(p=MockPoint(x=4.61, y=1.73, z=4.93), vertex_id=127, pressure=in_p, concentration=0,
                           R2=159467496.38934076, level=11),
    MockVesselTipPressures(p=MockPoint(x=8.69, y=5.49, z=6.05), vertex_id=128, pressure=in_p, concentration=0,
                           R2=134281999.0028481, level=11),
    MockVesselTipPressures(p=MockPoint(x=9.09, y=7.73, z=6.77), vertex_id=129, pressure=in_p, concentration=0,
                           R2=129818181.5848586, level=11),
    MockVesselTipPressures(p=MockPoint(x=5.25, y=6.37, z=0.853), vertex_id=130, pressure=in_p, concentration=0,
                           R2=153966783.72596282, level=12),
    MockVesselTipPressures(p=MockPoint(x=6.21, y=8.37, z=1.49), vertex_id=141, pressure=in_p, concentration=0,
                           R2=163470061.43666655, level=11),
    MockVesselTipPressures(p=MockPoint(x=5.97, y=9.09, z=0.933), vertex_id=143, pressure=in_p, concentration=0,
                           R2=163470061.43666673, level=11),
    MockVesselTipPressures(p=MockPoint(x=7.81, y=8.53, z=0.133), vertex_id=145, pressure=in_p, concentration=0,
                           R2=163470061.43666705, level=11),
]


def run2():
    data_folder = '../../data'
    output_folder = './tmp'
    mesh_3d_filename = '../../data/3d-meshes/test_full_1d0d3d_cm.xdmf'

    os.makedirs(output_folder, exist_ok=True)

    degree = 2
    tau = 1. / 2 ** 16
    tau_out = 1. / 2 ** 6
    t_end = 20.
    t = 0
    t_coup_start = 1.

    in_p = 33. * 1333


    mesh = read_mesh(mesh_3d_filename)
    solver3d = PressureSolver(mesh, output_folder, mock_vessel_tip_pressures)
    df.assign(solver3d.current.sub(0), df.interpolate(df.Constant(33 * 1333), solver3d.V.sub(0).collapse()))
    df.assign(solver3d.current.sub(1), df.interpolate(df.Constant(13 * 1333), solver3d.V.sub(1).collapse()))
    solver3d.write_solution(0.)

    solver3d.update_vessel_tip_pressures(mock_vessel_tip_pressures)
    print('start solving')
    solver3d.solve()
    print('end solving')
    solver3d.write_solution(1.)
    print('endsolving pressures')

    max_value = solver3d.current.vector().max()
    min_value = solver3d.current.vector().min()
    if df.MPI.rank(df.MPI.comm_world) == 0:
        print('max = {}, min = {}'.format(max_value / 1333, min_value / 1333))

    average_pressures = solver3d.get_pressures()

    for idx, key in enumerate(average_pressures.keys()):
        v = solver3d.current.split(True)[solver3d.first_multiplier_index+idx].vector()[:]
        print('id = {}, lambda = {}, avg_p = {}'.format(key, v, average_pressures[key]))


if __name__ == '__main__':
    f.initialize()
    #run2()
    run()
    f.finalize()
