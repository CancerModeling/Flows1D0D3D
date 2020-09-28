# import os
import numpy as np
import sys
import copy
import os
import pathlib

def bool_to_string(val):

    if val:
        return 'true'
    else:
        return 'false'

###########
# class-def
###########

class DefaultSimParams:
    'Class to store default parameters'

    def __init__(self):
        self.pp_tag = 't1'
        self.model_name = 'NetFVFE'
        self.test_name = 'none'
        self.output_debug_info = True
        self.advection_active = False
        self.network_decouple_nutrients = False
        self.coupled_1d3d = True
        self.solve_ecm = False
        self.solve_pres_with_net_update = True
        self.assembly_method = 2
        self.dimension = 3
        self.n_mpi = 1
        self.run_screen = False
        
        ## domain, mesh, and time
        self.L = 2.
        self.num_elems = 16
        self.final_time = 1.
        self.delta_t = 0.05
        self.total_outputs = 4
        self.dt_output = int(np.floor(self.final_time / self.delta_t) / self.total_outputs)
        if self.dt_output < 1:
            self.dt_output = 1

        ## solver
        self.linear_solver_max_iter = 250
        self.linear_solver_tol = 1.e-8
        self.nonlinear_solver_max_iter = 50
        self.nonlinear_solver_tol = 1.e-7

        ## 3D tumor parameters
        self.chi_c = 0.
        self.chi_h = 0.
        self.ic_nutrient_value = 0.5
        self.ECM_ic_val = 1.
        self.MDE_ic_val = 0.
        self.D_sigma = 3.
        self.D_TAF = 5.0e-1
        self.D_MDE = 1.
        self.bar_phi_ECM_P = 0.5

        ## lambda and sigma
        self.lambda_P = 5.
        self.lambda_A = 0.005
        self.lambda_Ph = 0.5
        self.lambda_PH = 1.
        self.lambda_HP = 1.
        self.lambda_HN = 1.
        self.lambda_TAF = 1.e+1
        self.lambda_ECM_D = 0.
        self.lambda_ECM_P = 0.
        self.lambda_MDE_D = 0.2
        self.lambda_MDE_P = 0.4
        self.sigma_PH  = 0.55
        self.sigma_HP = 0.65
        self.sigma_HN = 0.44

        ## cahn-hilliard
        self.mobility_P = 50.
        self.mobility_H = 0.5 * self.mobility_P
        self.E_phi_T = 0.045
        self.epsilon_P = 5.e-3
        self.epsilon_H = 5.e-3

        ## tumor ic
        self.tumor_ic_type = [1]
        self.tumor_ic_center = [[0.5 * self.L, 0.4 * self.L, 0.5 * self.L]]
        self.tumor_ic_radius = [[0.15 * self.L, 0.15 * self.L, 0.15 * self.L]]
        self.tumor_ic_file = 'tum_ic_data_' + self.pp_tag + '.csv'

        ## 3D flow parameters
        self.tissue_flow_visc = 1.
        self.tissue_flow_coef = 1.e-9
        self.tissue_flow_density = 1.
        self.tissue_pressure_ic_val = 0.
        
        ## 1D vessel parameters
        self.vessel_D_sigma = 1.e-1
        self.osmotic_coeff = 0.95
        self.scenario = 'two_vessels'
        self.L_p = 1.e-7
        self.L_s = 4.5
        self.vessel_in_nutrient = 1.
        self.vessel_in_nutrient_vein = 0.
        self.vessel_blood_density = 1.
        self.vessel_viscosity = 1.

        self.network_init_file = 'two_vessels'
        self.create_init_vessel = True
        self.vessel_num_refinement = 3
        self.discrete_cyl_length = 20
        self.discrete_cyl_angle = 20
        self.network_coupling_theta = 1.
        self.vessel_pressures = [3000., 2000., 1100., 1600.]
        self.vessel_radius = [0.046875, 0.0625]
        self.vessel_line_1 = [0.1875, 0.1875, 0., 0.1875, 0.1875, 2.]
        self.vessel_line_2 = [1.8125, 1.8125, 0., 1.8125, 1.8125, 2.]        
        self.identify_vein_pressure = 0.99 * self.vessel_pressures[1]
        self.identify_artery_radius = self.vessel_radius[1]
        self.coupling_3d1d_integration_method = 2
        self.disable_remove_redundant_vessel = False

        ## exprimental
        self.outlet_apply_neumann = False
        self.outlet_neumann_val = 0.
        self.inlet_apply_neumann = False
        self.inlet_neumann_val = 0.

        ## growth params
        self.network_active = True
        self.network_update_interval = 2
        self.network_update_TAF_th = 1.e-4
        self.log_normal_mean = 1.
        self.log_normal_std_dev = 0.2
        self.network_radius_exponent = 3.
        self.network_bifurc_prob = 0.94
        self.network_min_radius = 9.e-3
        self.network_sprout_prob = 0.93
        self.min_length_for_sprouting = 0.13
        self.seed = 100
        self.network_update = True

        self.run_path = ''


    def set_time(self, time, dt, total_out):
        self.final_time = time
        self.delta_t = dt
        if total_out > 0:
            self.total_outputs = total_out
            self.dt_output = int(np.floor(self.final_time / self.delta_t) / self.total_outputs)
        else:
            self.total_outputs = int(self.final_time / self.delta_t)
            self.dt_output = 1

        if self.dt_output < 1:
            self.dt_output = 1

    def set_pressures(self, pressures):
        self.vessel_pressures = pressures
        self.identify_vein_pressure = 0.999*self.vessel_pressures[1]

    def set_vessel_radius(self, radius):
        self.vessel_radius = radius
        self.identify_artery_radius = radius[1]
        

    # static method
    def write_str(self, to_inp = True):

        strs = ''
        if to_inp == False:
            strs += '\n---       Parameter values       ---\n\n'
        
        space = 30 + 2
        strs += '{}= {}\n'.format('{0:{space}}'.format('model_name', space=space), self.model_name)
        strs += '{}= {}\n'.format('{0:{space}}'.format('pp_tag', space=space), self.pp_tag)
        strs += '{}= {}\n'.format('{0:{space}}'.format('test_name', space=space), self.test_name)
        strs += '{}= {}\n'.format('{0:{space}}'.format('advection_active', space=space), bool_to_string(self.advection_active))
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_decouple_nutrients', space=space), bool_to_string(self.network_decouple_nutrients))
        strs += '{}= {}\n'.format('{0:{space}}'.format('coupled_1d3d', space=space), bool_to_string(self.coupled_1d3d))
        strs += '{}= {}\n'.format('{0:{space}}'.format('solve_ecm', space=space), bool_to_string(self.solve_ecm))
        strs += '{}= {}\n'.format('{0:{space}}'.format('solve_pres_with_net_update', space=space), bool_to_string(self.solve_pres_with_net_update))
        strs += '{}= {}\n'.format('{0:{space}}'.format('assembly_method', space=space), self.assembly_method)
        strs += '{}= {}\n'.format('{0:{space}}'.format('dimension', space=space), self.dimension)
        
        ## domain and time
        strs += '\n'
        if to_inp == False:
            strs += '{}= {}\n'.format('{0:{space}}'.format('L', space=space), self.L)
        else:
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_xmin', space=space), 0.)
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_xmax', space=space), self.L)
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_ymin', space=space), 0.)
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_ymax', space=space), self.L)
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_zmin', space=space), 0.)
            strs += '{}= {}\n'.format('{0:{space}}'.format('domain_zmax', space=space), self.L)
        strs += '{}= {}\n'.format('{0:{space}}'.format('mesh_n_elements', space=space), self.num_elems)
        strs += '{}= {}\n'.format('{0:{space}}'.format('final_time', space=space), self.final_time)
        strs += '{}= {}\n'.format('{0:{space}}'.format('time_step', space=space), self.delta_t)
        strs += '{}= {}\n'.format('{0:{space}}'.format('perform_output', space=space), 'true')
        strs += '{}= {}\n'.format('{0:{space}}'.format('output_interval', space=space), self.dt_output)
        strs += '{}= {}\n'.format('{0:{space}}'.format('output_tag', space=space), self.pp_tag)
        strs += '{}= {}\n'.format('{0:{space}}'.format('output_debug_info', space=space), bool_to_string(self.output_debug_info))

        ## solver
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('linear_solver_max_iter', space=space), self.linear_solver_max_iter)
        strs += '{}= {}\n'.format('{0:{space}}'.format('linear_solver_tol', space=space), self.linear_solver_tol)
        strs += '{}= {}\n'.format('{0:{space}}'.format('nonlinear_solver_max_iter', space=space), self.nonlinear_solver_max_iter)
        strs += '{}= {}\n'.format('{0:{space}}'.format('nonlinear_solver_tol', space=space), self.nonlinear_solver_tol)

        ## 3D tumor params
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('chi_c', space=space), self.chi_c)
        strs += '{}= {}\n'.format('{0:{space}}'.format('chi_h', space=space), self.chi_h)
        strs += '{}= {}\n'.format('{0:{space}}'.format('ic_nutrient_value', space=space), self.ic_nutrient_value)
        strs += '{}= {}\n'.format('{0:{space}}'.format('ECM_ic_val', space=space), self.ECM_ic_val)
        strs += '{}= {}\n'.format('{0:{space}}'.format('MDE_ic_val', space=space), self.MDE_ic_val)
        strs += '{}= {}\n'.format('{0:{space}}'.format('D_sigma', space=space), self.D_sigma)
        strs += '{}= {}\n'.format('{0:{space}}'.format('D_TAF', space=space), self.D_TAF)
        strs += '{}= {}\n'.format('{0:{space}}'.format('D_MDE', space=space), self.D_MDE)
        strs += '{}= {}\n'.format('{0:{space}}'.format('bar_phi_ECM_P', space=space), self.bar_phi_ECM_P)

        ## lambdas and sigmas
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_P', space=space), self.lambda_P)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_A', space=space), self.lambda_A)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_Ph', space=space), self.lambda_Ph)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_PH', space=space), self.lambda_PH)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_HP', space=space), self.lambda_HP)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_HN', space=space), self.lambda_HN)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_TAF', space=space), self.lambda_TAF)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_ECM_D', space=space), self.lambda_ECM_D)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_ECM_P', space=space), self.lambda_ECM_P)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_MDE_D', space=space), self.lambda_MDE_D)
        strs += '{}= {}\n'.format('{0:{space}}'.format('lambda_MDE_P', space=space), self.lambda_MDE_P)
        strs += '{}= {}\n'.format('{0:{space}}'.format('sigma_PH', space=space), self.sigma_PH)
        strs += '{}= {}\n'.format('{0:{space}}'.format('sigma_HP', space=space), self.sigma_HP)
        strs += '{}= {}\n'.format('{0:{space}}'.format('sigma_HN', space=space), self.sigma_HN)

        ## cahn-hilliard
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('bar_M_P', space=space), self.mobility_P)
        strs += '{}= {}\n'.format('{0:{space}}'.format('bar_M_H', space=space), self.mobility_H)
        strs += '{}= {}\n'.format('{0:{space}}'.format('bar_E_phi_T', space=space), self.E_phi_T)
        strs += '{}= {}\n'.format('{0:{space}}'.format('epsilon_P', space=space), self.epsilon_P)
        strs += '{}= {}\n'.format('{0:{space}}'.format('epsilon_H', space=space), self.epsilon_H)

        ## tumor ic (this data can also be provided in seperate file easier when we have multiple cores)
        strs += '\n'
        if to_inp == False:
            strs += '{}= {}\n'.format('{0:{space}}'.format('ic_tumor_type', space=space), self.tumor_ic_type)
            strs += '{}= {}\n'.format('{0:{space}}'.format('tumor_ic_center', space=space), self.tumor_ic_center)
            strs += '{}= {}\n'.format('{0:{space}}'.format('tumor_ic_radius', space=space), self.tumor_ic_radius)
        else:
            self.tumor_ic_file = 'tum_ic_data_' + self.pp_tag + '.csv'
            strs += '{}= {}\n'.format('{0:{space}}'.format('ic_tumor_file', space=space), self.tumor_ic_file)

        ## 3D flow
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_flow_viscosity', space=space), self.tissue_flow_visc)
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_flow_K', space=space), self.tissue_flow_coef)
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_flow_density', space=space), self.tissue_flow_density)

        if self.create_init_vessel:
            self.tissue_pressure_ic_val = self.vessel_pressures[1]
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_pressure_ic_val', space=space), self.tissue_pressure_ic_val)

        ## 1D vessel params 
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_D_sigma', space=space), self.vessel_D_sigma)
        strs += '{}= {}\n'.format('{0:{space}}'.format('osmotic_reflection_coeff', space=space), self.osmotic_coeff)
        strs += '{}= {}\n'.format('{0:{space}}'.format('scenario', space=space), self.scenario)
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_flow_L_p', space=space), self.L_p)
        strs += '{}= {}\n'.format('{0:{space}}'.format('tissue_nut_L_s', space=space), self.L_s)
        strs += '{}= {}\n'.format('{0:{space}}'.format('assembly_factor_p_t', space=space), 1./self.L_p)
        strs += '{}= {}\n'.format('{0:{space}}'.format('assembly_factor_c_t', space=space), 1./self.L_s)
        strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_in_nutrient', space=space), self.vessel_in_nutrient)
        strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_in_nutrient_vein', space=space), self.vessel_in_nutrient_vein)
        strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_blood_density', space=space), self.vessel_blood_density)
        strs += '{}= {}\n'.format('{0:{space}}'.format('init_vessel_viscosity', space=space), self.vessel_viscosity)

        strs += '{}= {}\n'.format('{0:{space}}'.format('outlet_apply_neumann', space=space), bool_to_string(self.outlet_apply_neumann))
        strs += '{}= {}\n'.format('{0:{space}}'.format('outlet_neumann_val', space=space), self.outlet_neumann_val)
        strs += '{}= {}\n'.format('{0:{space}}'.format('inlet_apply_neumann', space=space), bool_to_string(self.inlet_apply_neumann))
        strs += '{}= {}\n'.format('{0:{space}}'.format('inlet_neumann_val', space=space), self.inlet_neumann_val)

        strs += '\n'
        if to_inp:
            if self.create_init_vessel:
                self.network_init_file = 'two_vessels_' + self.pp_tag + '.dgf'
            strs += '{}= {}\n'.format('{0:{space}}'.format('network_init_file', space=space), self.network_init_file)

        strs += '{}= {}\n'.format('{0:{space}}'.format('create_init_vessel', space=space), self.create_init_vessel)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_init_refinement', space=space), self.vessel_num_refinement)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_discret_cyl_length', space=space), self.discrete_cyl_length)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_discret_cyl_angle', space=space), self.discrete_cyl_angle)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_coupling_method_theta', space=space), self.network_coupling_theta)
        strs += '{}= {}\n'.format('{0:{space}}'.format('coupling_3d1d_integration_method', space=space), self.coupling_3d1d_integration_method)
        strs += '{}= {}\n'.format('{0:{space}}'.format('disable_remove_redundant_vessel', space=space), bool_to_string(self.disable_remove_redundant_vessel))
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_compute_elem_weights', space=space), "true")
        
        if to_inp:
            if self.create_init_vessel:
                self.identify_vein_pressure = 0.999*self.vessel_pressures[1]
        strs += '{}= {}\n'.format('{0:{space}}'.format('identify_vein_pressure', space=space), self.identify_vein_pressure)

        if self.create_init_vessel:
            self.identify_artery_radius = self.vessel_radius[1]
        strs += '{}= {}\n'.format('{0:{space}}'.format('identify_artery_radius', space=space), self.identify_artery_radius)
        
        if to_inp == False:
            strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_pressures', space=space), self.vessel_pressures)
            strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_radius', space=space), self.vessel_radius)
            strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_line_1', space=space), self.vessel_line_1)
            strs += '{}= {}\n'.format('{0:{space}}'.format('vessel_line_2', space=space), self.vessel_line_2)
        

        ## growth params
        strs += '\n'
        strs += '{}= {}\n'.format('{0:{space}}'.format('is_network_active', space=space), bool_to_string(self.network_active))
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_update_interval', space=space), self.network_update_interval)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_update_taf_threshold', space=space), self.network_update_TAF_th)
        strs += '{}= {}\n'.format('{0:{space}}'.format('log_normal_mean', space=space), self.log_normal_mean)
        strs += '{}= {}\n'.format('{0:{space}}'.format('log_normal_std_dev', space=space), self.log_normal_std_dev)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_radius_exponent_gamma', space=space), self.network_radius_exponent)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_bifurcate_probability', space=space), self.network_bifurc_prob)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_min_radius', space=space), self.network_min_radius)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_sprouting_prob', space=space), self.network_sprout_prob)
        strs += '{}= {}\n'.format('{0:{space}}'.format('min_length_for_sprouting', space=space), self.min_length_for_sprouting)
        strs += '{}= {}\n'.format('{0:{space}}'.format('seed', space=space), self.seed)
        strs += '{}= {}\n'.format('{0:{space}}'.format('network_update', space=space), self.network_update)

        if to_inp == False:
            strs += '\n------------------------------------\n'

        return strs

    def show(self):
        print(self.write_str(False))


###########
# class-def
###########

class SimParams(DefaultSimParams):
    'Class to hold simulation parameters'

    def __init__(self):
        super(SimParams, self).__init__()


###########
# problem setup
###########

def gen_tumor_ic_file(file_dir, dp):

    # write to file
    inpf = open(file_dir + dp.tumor_ic_file, 'w')
    inpf.write("type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz\n")

    # type = 1 -- spherical tumor core
    # type = 2 -- elliptical tumor core (sharp)
    # type = 3 -- spherical tumor core and then spherical hypoxic core
    # type = 5 -- spherical tumor core (sharp)
    ic_type = dp.tumor_ic_type
    center = dp.tumor_ic_center
    r = dp.tumor_ic_radius
    for i in range(len(center)):
        inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(ic_type[i], center[i][0], center[i][1], center[i][2], r[i][0], r[i][1], r[i][2], r[i][0], r[i][1], r[i][2]))

    inpf.close()

def gen_init_network_file(file_dir, dp):

    L = dp.L
    R = dp.vessel_radius
    pressures = dp.vessel_pressures
    l1 = dp.vessel_line_1
    l2 = dp.vessel_line_2

    # write to file
    inpf = open(file_dir + dp.network_init_file, 'w')
    inpf.write("DGF\n")

    # vertex
    inpf.write("Vertex\n")
    inpf.write("parameters {}\n".format(1))
    inpf.write("{} {} {} {}\n".format(l1[0], l1[1], l1[2], pressures[0]))
    inpf.write("{} {} {} {}\n".format(l1[3], l1[4], l1[5], pressures[1]))
    inpf.write("{} {} {} {}\n".format(l2[0], l2[1], l2[2], pressures[2]))
    inpf.write("{} {} {} {}\n".format(l2[3], l2[4], l2[5], pressures[3]))

    # segments
    inpf.write("#\n")
    inpf.write("SIMPLEX\n")
    inpf.write("parameters {}\n".format(2))
    inpf.write("{} {} {} {}\n".format(0, 1, R[0], 0.0075))
    inpf.write("{} {} {} {}\n".format(2, 3, R[1], 0.0075))
    # 
    inpf.write("#\n")
    inpf.write("BOUNDARYDOMAIN\n")
    inpf.write("default {}\n".format(1))

    # 
    inpf.write("#")
    inpf.close()

def setup():

  n = 1
  dps = []
  for i in range(n):

    # create sim params
    dp = SimParams()
    dp.model_name = 'NetFVFE'
    dp.test_name = 'none'
    dp.coupled_1d3d = False
    dp.pp_tag = 'angio_' + str(i)
    dp.run_path = 'tests_sim/' + dp.model_name + '/' + dp.pp_tag + '/' 

    dps.append(dp)

  return dps

def run():

  dps = setup()
  run_screen = 1 # 1 - true, 0 - false

  for dp in dps:
    print('setting up input parameters ... ')

    pathlib.Path('tests_sim').mkdir(parents=True, exist_ok=True)
    pathlib.Path('tests_sim/' + dp.model_name).mkdir(parents=True, exist_ok=True)
    pathlib.Path(dp.run_path).mkdir(parents=True, exist_ok=True) 

    # create input file
    fo = open(dp.run_path + 'input.in', 'w')
    fo.write(dp.write_str())
    fo.close()

    # create tumor ic file
    gen_tumor_ic_file(dp.run_path, dp)

    # create network file
    if dp.create_init_vessel:
      gen_init_network_file(dp.run_path, dp)
    else:
      if dp.run_path != './':
        os.system('cp ' + dp.network_init_file + ' ' + dp.run_path)

    # run
    print('calling executible ... ')
    os.system('./run.sh ' + dp.model_name + ' ' + dp.run_path + ' ' + str(dp.n_mpi) + ' ' + str(run_screen) + ' ' + dp.pp_tag)

#### 
if len(sys.argv) > 1:
  if str(sys.argv[1]) == 'run':
    run()
else:
    print('To run sim, run script with argument "run"')
