import os
import numpy as np
import dolfin as df
from _utils import setup_subdomains
from parameters import TransportModelParameters


class TransportSolver:
    def __init__(self, mesh, output_folder, flow_solver, list_vessel_tips, transport_config=TransportModelParameters()):
        self.mesh = mesh
        self.transport_config = transport_config
        self.flow_solver = flow_solver
        self.list_vessel_tips = list_vessel_tips
        self._setup_vessel_tip_pressures(list_vessel_tips)
        self._setup_subdomains()
        self._setup_space()
        self._setup_problem()
        self.file_c_cap = df.File(os.path.join(output_folder, "c_cap.pvd"), "compressed")
        self.file_c_tis = df.File(os.path.join(output_folder, "c_tis.pvd"), "compressed")

    def _setup_vessel_tip_pressures(self, list_vessel_tips):
        num_outlets = len(list_vessel_tips)
        self.points = np.zeros((num_outlets, 3))
        self.point_to_vertex_id = np.zeros(num_outlets)
        self.pressures = np.zeros(num_outlets)
        self.concentrations = np.zeros(num_outlets)
        self.average_pressures = np.zeros(num_outlets)
        self.R2 = np.zeros(num_outlets)
        self.level = np.zeros(num_outlets)
        for idx, vessel_tip in enumerate(list_vessel_tips):
            p = vessel_tip.p
            self.points[idx, :] = np.array([p.x, p.y, p.z])
            self.point_to_vertex_id[idx] = vessel_tip.vertex_id
            self.pressures[idx] = vessel_tip.pressure
            self.concentrations[idx] = vessel_tip.concentration
            self.R2[idx] = vessel_tip.R2
            self.level[idx] = vessel_tip.level

    def _setup_subdomains(self):
        self.subdomains, self.dx = setup_subdomains(self.mesh, self.points)

    def _setup_space(self):
        P1El = df.FiniteElement('P', self.mesh.ufl_cell(), 1)
        MixedEl = df.MixedElement([P1El, P1El])
        self.V = df.FunctionSpace(self.mesh, MixedEl)

    def _setup_problem(self):
        phi_c, phi_t = df.TrialFunctions(self.V)
        psi_c, psi_t = df.TestFunctions(self.V)

        volumes = []
        for k in range(len(self.pressures)):
            volume = df.assemble(df.Constant(1) * self.dx(k))
            volumes.append(volume)
        self.volumes = volumes

        p_c = self.flow_solver.current.split()[0]
        p_t = self.flow_solver.current.split()[1]

        config = self.flow_solver.flow_config

        K_c = config.K_c
        K_t = config.K_t
        mu_c = config.mu_c
        mu_t = config.mu_t
        rho_c = config.rho_c

        config = self.transport_config
        D_c = config.D_c
        D_t = config.D_t
        lambda_t = config.lambda_t

        v_c = -K_c/mu_c * df.grad(p_c)
        v_t = -K_t/mu_t * df.grad(p_t)

        J = 0
        J += - df.inner(v_c * phi_c - df.Constant(D_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        J += - df.inner(v_t * phi_t - df.Constant(D_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx
        #J += - df.inner(- df.Constant(D_c) * df.grad(phi_c), df.grad(psi_c)) * df.dx
        #J += - df.inner(- df.Constant(D_t) * df.grad(phi_t), df.grad(psi_t)) * df.dx

        # t_cv
        q_cv = self.flow_solver.get_q_cv(p_c)
        from_c_to_v = df.conditional(df.le(q_cv, 0), df.Constant(1), df.Constant(0))
        from_v_to_c = df.conditional(df.ge(q_cv, 0), df.Constant(1), df.Constant(0))
        J -= q_cv * from_c_to_v * phi_c * psi_c * df.dx
        J -= q_cv * from_v_to_c * df.Constant(0.0) * psi_c * df.dx
        # t_ct
        q_ct = self.flow_solver.get_q_ct(p_c, p_t) 
        from_c_to_t = df.conditional(df.le(q_ct, 0), df.Constant(1), df.Constant(0))
        from_t_to_c = df.conditional(df.ge(q_ct, 0), df.Constant(1), df.Constant(0))
        J -= q_ct * from_c_to_t * phi_c * psi_c * df.dx
        J -= q_ct * from_t_to_c * phi_t * psi_c * df.dx
        # t_ca
        for k in range(len(self.pressures)):
            p_bar = self.average_pressures[k]
            q_ca = self.flow_solver.get_q_ca(k, p_bar)
            from_c_to_a = df.conditional(df.le(q_ca, 0), df.Constant(1), df.Constant(0))
            from_a_to_c = df.conditional(df.ge(q_ca, 0), df.Constant(1), df.Constant(0))
            J -= q_ca * from_a_to_c * df.Constant(max(self.concentrations[k], 0)) * psi_c * self.dx(k)
            J -= q_ca * from_c_to_a * df.Constant(0) * psi_c * self.dx(k)

        # t_tc
        q_tc = self.flow_solver.get_q_tc(p_c, p_t) 
        J -= q_tc * from_c_to_t * phi_c * psi_t * df.dx
        J -= q_tc * from_t_to_c * phi_t * psi_t * df.dx

        # t_tl
        q_tl = self.flow_solver.get_q_tl(p_t) 
        from_t_to_l = df.conditional(df.le(q_tl, 0), df.Constant(1), df.Constant(0))
        from_l_to_t = df.conditional(df.ge(q_tl, 0), df.Constant(1), df.Constant(0))
        J -= q_tl * from_t_to_l * phi_t * psi_t * df.dx
        J -= q_tl * from_l_to_t * df.Constant(0) * psi_t * df.dx
        #J -= q_tl * phi_t * psi_t * df.dx

        # r_t
        J -= - lambda_t * phi_t *psi_t * df.dx

        self.J = J
        self.bcs = []
        self.current = df.Function(self.V, name='u') if not hasattr(self, 'current') else self.current

    def update_average_pressures(self, average_pressures):
        for idx in range(len(self.points)):
            self.average_pressures[idx] = average_pressures[self.point_to_vertex_id[idx]]
        # we have to reassemble the system:
        self._setup_problem()

    def update_vessel_tip_concentrations(self, list_vessel_tip_concentrations):
        for idx, vessel_tip in enumerate(list_vessel_tip_concentrations):
            assert vessel_tip.vertex_id == self.point_to_vertex_id[idx]
            self.concentrations[idx] = vessel_tip.concentration
        # we have to reassemble the system:
        self._setup_problem()

    def solve(self):
        A, b = df.assemble_system(df.lhs(self.J), df.rhs(self.J), self.bcs)
        solver = df.KrylovSolver('gmres', 'jacobi')
        solver.parameters['monitor_convergence'] = True
        solver.parameters['nonzero_initial_guess'] = True
        solver.parameters['absolute_tolerance'] = 1e-12
        solver.parameters['relative_tolerance'] = 1e-12
        solver.set_operator(A)
        solver.solve(self.current.vector(), b)

    def write_solution(self, t):
        self.file_c_cap << (self.current.split()[0], t)
        self.file_c_tis << (self.current.split()[1], t)