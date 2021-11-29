import os
import numpy as np
import dolfin as df
from _utils import setup_subdomains
from parameters import FlowModelParameters


class ImplicitPressureSolver:
    def __init__(self, mesh, output_folder, vessel_tip_pressures, flow_config=FlowModelParameters()):
        self.mesh = mesh
        self.flow_config = flow_config
        self._setup_vessel_tip_pressures(vessel_tip_pressures)
        self._setup_subdomains()
        self._setup_volumes()
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
            self.radii[idx] = vessel_tip.radius_first
            self.level[idx] = vessel_tip.level
            self.total_flows[idx] = 2**(vessel_tip.level-1) * (vessel_tip.pressure - 30 * 1333) / vessel_tip.R2

    def _setup_subdomains(self):
        weights = np.ones(self.radii.shape)
        self.subdomains, self.dx = setup_subdomains(self.mesh, self.points, weights)
    
    def _setup_volumes(self):
        volumes = []
        for k in range(len(self.pressures)):
            volume = df.assemble(df.Constant(1) * self.dx(k))
            volumes.append(volume)
        self.volumes = volumes
        self.volume = df.assemble(df.Constant(1) * self.dx)

    def _setup_function_spaces(self):
        PEl = df.FiniteElement('P', self.mesh.ufl_cell(), 1)
        REl = df.FiniteElement('R', self.mesh.ufl_cell(), 0)
        SpatialEl = [PEl, PEl]
        elements = SpatialEl + [REl for idx in range(len(self.points))]
        self.first_multiplier_index = len(SpatialEl)
        MEl = df.MixedElement(elements)
        self.V = df.FunctionSpace(self.mesh, MEl)

    def _get_coeff_ca(self):
        coeff_ca = []
        for k in range(len(self.pressures)):
            coeff_ca.append(self.flow_config.rho_c * 2 ** (self.level[k] - 1) / self.R2[k] / self.volumes[k])
        # mean value:
        coeff_ca_mean = np.array(coeff_ca).mean()
        coeff_ca = np.ones(len(coeff_ca)) * coeff_ca_mean
        return coeff_ca

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

        rho_c = self.flow_config.rho_c
        rho_t = self.flow_config.rho_t
        K_c = self.flow_config.K_c
        K_t = self.flow_config.K_t
        mu_c = self.flow_config.mu_c
        mu_t = self.flow_config.mu_t
        L_tl = self.flow_config.L_cv
        L_cv = self.flow_config.L_cv
        L_ct = self.flow_config.L_ct
        S_ct = self.flow_config.S_ct

        self.update_average_pressures()

        J = 0
        J += df.inner(df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        J += df.inner(df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx

        coeff_ca = []
        for k in range(len(self.pressures)):
            coeff_ca.append(rho_c * 2 ** (self.level[k] - 1) / self.R2[k] / self.volumes[k])
        # mean value:
        coeff_ca_mean = np.array(coeff_ca).mean()
        coeff_ca = np.ones(len(coeff_ca)) * coeff_ca_mean

        # llambda[k] = - 1/Omega_k int_Omega[k] p_c dx
        for k in range(len(self.pressures)):
            alpha = df.Constant(coeff_ca[k] + rho_c * L_cv)
            J += - alpha * llambda[k] * mu[k] * self.dx(k) + alpha * phi_c * mu[k] * self.dx(k)

        # q_cv
        J -= self.get_q_cv(phi_c) * psi_c * self.dx

        # q_ca
        for k in range(len(self.pressures)):
            J -= self.get_q_ca(k, llambda[k]) * psi_c * self.dx(k)

        # q_ct
        J -= self.get_q_ct(phi_c, phi_t) * psi_c * self.dx

        # q_tc
        J -= self.get_q_tc(phi_c, phi_t) * psi_t * self.dx

        # q_tl
        J -= self.get_q_tl(phi_t) * psi_t * self.dx

        # setup diagonal preconditioner
        P = 0
        P += df.inner(df.Constant(rho_c * K_c / mu_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        P += df.inner(df.Constant(rho_t * K_t / mu_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_c) * psi_c * self.dx
        P -= df.Constant(rho_t * L_ct * S_ct) * (-phi_t) * psi_t * self.dx
        P -= df.Constant(rho_t * L_tl) * (- phi_t) * psi_t * self.dx
        P += df.Constant(0) * psi_t * df.dx
        for k in range(len(self.pressures)):
            alpha = df.Constant(coeff_ca[k] + rho_c * L_cv)
            P += - alpha * llambda[k] * mu[k] * self.dx(k)

        self.J = J
        self.P = P
        self.bcs = []
    
    def get_q_ca(self, k, llambda):
        coeff_ca = self._get_coeff_ca()
        return df.Constant(coeff_ca[k]) * (df.Constant(self.pressures[k]) - llambda) 
    
    def get_q_cv(self, phi_c):
        # parameters:
        rho_c = self.flow_config.rho_c
        L_cv = self.flow_config.L_cv
        p_ven = self.flow_config.p_ven
        # term:
        return df.Constant(rho_c * L_cv) * (df.Constant(p_ven) - phi_c) 
    
    def get_q_ct(self, phi_c, phi_t):
        # parameters:
        rho_t = self.flow_config.rho_t
        L_ct = self.flow_config.L_ct
        S_ct = self.flow_config.S_ct
        sigma = self.flow_config.sigma
        pi_bl = self.flow_config.pi_bl
        pi_int = self.flow_config.pi_int
        # term:
        q_ct = df.Constant(rho_t * L_ct * S_ct) * (phi_t - phi_c)
        q_ct += df.Constant(rho_t * L_ct * S_ct * sigma * (pi_bl - pi_int))
        return q_ct

    def get_q_tc(self, phi_c, phi_t):
        return - self.get_q_ct(phi_c, phi_t)
    
    def get_q_tl(self, phi_t):
        # parameters:
        rho_t = self.flow_config.rho_t
        L_tl = self.flow_config.L_tl
        p_lym = self.flow_config.p_lym
        # term:
        return df.Constant(rho_t * L_tl) * (p_lym - phi_t)

    def _setup_solver(self):
        self.A, self.b = df.assemble_system(df.lhs(self.J), df.rhs(self.J), self.bcs)
        self.P, _ = df.assemble_system(df.lhs(self.P), df.rhs(self.P), self.bcs)
        self.solver = df.KrylovSolver('gmres', 'amg')
        self.solver.parameters['monitor_convergence'] = True
        self.solver.parameters['nonzero_initial_guess'] = True
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

    def get_pressures_from_multiplier(self):
        new_data = {}
        # trial_functions = df.TrialFunctions(self.V)
        # llambda = trial_functions[self.first_multiplier_index:]
        for k in range(len(self.pressures)):
            #average_pressure = self.current.split()[self.first_multiplier_index+k].vector().get_local()[0]
            average_pressure = df.assemble(self.current.split()[self.first_multiplier_index+k] * self.dx)/self.volume
            new_data[int(self.point_to_vertex_id[k])] = average_pressure
        return new_data

    def update_average_pressures(self):
        for k in range(len(self.pressures)):
            average_pressure = df.assemble(self.current.split()[0] * self.dx(k)) / self.volumes[k]
            self.average_pressures[k] = average_pressure
