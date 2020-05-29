import os
import numpy as np
import sys

def add(param_index, param_val, index, val):

    param_index.append(index)
    param_val.append(val)

def bool_to_string(val):

    if val:
        return 'true'
    else:
        return 'false'


def gen_tumor_ic_file(filename, type, center, r):

    # write to file
    inpf = open(filename,'w')
    inpf.write("type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz\n")

    # type = 1 -- spherical tumor core
    # type = 2 -- elliptical tumor core (sharp)
    # type = 3 -- spherical tumor core and then spherical hypoxic core
    # type = 5 -- spherical tumor core (sharp)
    inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(type, center[0], center[1], center[2], r[0], r[1], r[2], r[0], r[1], r[2]))

    inpf.close()


def gen_init_network_file(L, R, pressures, filename):

    # write to file
    inpf = open(filename,'w')
    inpf.write("DGF\n")

    # vertex
    inpf.write("Vertex\n")
    inpf.write("parameters {}\n".format(1))
    inpf.write("{} {} {} {}\n".format(L - 2.*R, L - 2.*R, 0.0, pressures[0]))
    inpf.write("{} {} {} {}\n".format(L - 2.*R, L - 2.*R, L, pressures[1]))
    inpf.write("{} {} {} {}\n".format(2.*R, 2.*R, 0.0, pressures[2]))
    inpf.write("{} {} {} {}\n".format(2.*R, 2.*R, L, pressures[3]))

    # segments
    inpf.write("#\n")
    inpf.write("SIMPLEX\n")
    inpf.write("parameters {}\n".format(2))
    inpf.write("{} {} {} {}\n".format(0, 1, R, 0.0075))
    inpf.write("{} {} {} {}\n".format(2, 3, R, 0.0075))
    # 
    inpf.write("#\n")
    inpf.write("BOUNDARYDOMAIN\n")
    inpf.write("default {}\n".format(1))

    # 
    inpf.write("#")
    inpf.close()


def network_input(L, pp_tag, vessel_radius, pressures, network_file, file_dir, num_refinement, param_index, param_val):

    add(param_index, param_val, 'is_network_active', 'true')
    
    # network file
    add(param_index, param_val, 'network_init_file', network_file)
    add(param_index, param_val, 'network_init_refinement', num_refinement)

    # set below to reasonable value such as 1, 4, 10 if want to grow network
    add(param_index, param_val, 'network_update_interval', 100000)

    # control parameters for growth algorithm
    add(param_index, param_val, 'vessel_lambda_g', 0.5)
    add(param_index, param_val, 'vessel_R_factor', 1.)
    add(param_index, param_val, 'log_normal_mean', 0.01)
    add(param_index, param_val, 'log_normal_std_dev', 0.1)
    add(param_index, param_val, 'network_radius_exponent_gamma', 2.)
    add(param_index, param_val, 'network_no_branch_dist', 10)
    add(param_index, param_val, 'network_new_veesel_max_angle', 0.4)
    add(param_index, param_val, 'network_branch_angle', 0.5)
    add(param_index, param_val, 'network_update_taf_threshold', 0.1)
    add(param_index, param_val, 'network_vessel_no_taf_dist', 0)
    add(param_index, param_val, 'network_nonlocal_search_num_points', 3)
    add(param_index, param_val, 'network_nonlocal_search_length_factor', 5.)
    add(param_index, param_val, 'network_local_search', 'false')
    add(param_index, param_val, 'network_no_new_node_search_factor', 0.25)

    # to identify veins so that we can apply correct bc
    ident_pres = 0.9999 * pressures[1]
    add(param_index, param_val, 'identify_vein_pressure', ident_pres)
    
    # generate network file
    gen_init_network_file(L, vessel_radius, pressures, file_dir + network_file)


def input(L, file_dir, pp_tag, model_name, test_name, output_debug_info, advection_active, decouple_nutrients, num_elems, final_time, delta_t, dt_output, chi_c, D_sigma, nutrient_ic_val, tissue_flow_coef, vessel_D_sigma, L_p, L_s, vessel_pressures, vessel_radius, tumor_ic_type, tumor_ic_center, tumor_ic_radius, vessel_num_refinement, lambda_P, mobility_P):

    # Info:
    # Creates input.in for model
    # Creates initial network file
    # Creates tumor initial condition file

    param_index = []
    param_val = []
    break_points = []
    break_msg = []

    # files
    # input_file = 'input_' + str(pp_tag) + '.in'
    input_file = 'input.in'
    tumor_ic_file = 'tum_ic_data_' + str(pp_tag) + '.csv'
    network_file = 'two_vessels_' + str(pp_tag) + '.dgf'
    
    ## model
    break_points.append(len(param_val))
    break_msg.append('# model')
    
    # specify model such as NetFVFE, NetFVFE, NetFC, AvaFV
    add(param_index, param_val, 'model_name', model_name)

    # specify test (if any) which solves sub-system
    # disable line below if running full system or specify empty string ''
    add(param_index, param_val, 'test_name', test_name)

    # domain
    add(param_index, param_val, 'dimension', 3)
    add(param_index, param_val, 'domain_xmin', 0.)
    add(param_index, param_val, 'domain_xmax', L)
    add(param_index, param_val, 'domain_ymin', 0.)
    add(param_index, param_val, 'domain_ymax', L)
    add(param_index, param_val, 'domain_zmin', 0.)
    add(param_index, param_val, 'domain_zmax', L)    

    # there are various ways to assemble the source terms
    # 1 - use implicit whenever possible
    # 2 - first project species to [0,1] and follow 1
    # 3 - project to [0,1] and handle all source terms explicitly
    # In experience, 2 is more robust
    add(param_index, param_val, 'assembly_method', 2)

    # simplification of computation
    add(param_index, param_val, 'advection_active', bool_to_string(advection_active))
    add(param_index, param_val, 'network_decouple_nutrients', bool_to_string(decouple_nutrients))

    # control parameters for 1d-3d coupling
    add(param_index, param_val, 'network_discret_cyl_length', 20)
    add(param_index, param_val, 'network_discret_cyl_angle', 20)
    add(param_index, param_val, 'network_compute_elem_weights', 'true')
    add(param_index, param_val, 'network_coupling_method_theta', 1.0)

    
    ## restart info
    break_points.append(len(param_val))
    break_msg.append('\n# restart')
    add(param_index, param_val, 'restart', 'false')
    # add(param_index, param_val, 'mesh_restart_file', 'mesh.e')
    # add(param_index, param_val, 'solution_restart_file', 'restart_220.e')

    ## mesh
    break_points.append(len(param_val))
    break_msg.append('\n# mesh')
    add(param_index, param_val, 'mesh_n_elements', num_elems)

    ## time
    break_points.append(len(param_val))
    break_msg.append('\n# time')
    init_t = 0.
    add(param_index, param_val, 'time_step', delta_t)
    add(param_index, param_val, 'initial_time', init_t)
    add(param_index, param_val, 'initial_step', 0)
    add(param_index, param_val, 'final_time', final_time)

    ## output
    break_points.append(len(param_val))
    break_msg.append('\n# output')
    if dt_output < 1:
        dt_output = 1
    add(param_index, param_val, 'perform_output', 'true')
    add(param_index, param_val, 'output_interval', dt_output)
    add(param_index, param_val, 'restart_save', 'false')
    add(param_index, param_val, 'restart_save_interval', 1)
    add(param_index, param_val, 'output_tag', str(pp_tag))
    add(param_index, param_val, 'output_debug_info', bool_to_string(output_debug_info))


    ## solver
    break_points.append(len(param_val))
    break_msg.append('\n# solver')
    add(param_index, param_val, 'linear_solver_max_iter', 250)
    add(param_index, param_val, 'linear_solver_tol', 1.0e-8)
    add(param_index, param_val, 'nonlinear_solver_max_iter', 50)
    add(param_index, param_val, 'nonlinear_solver_tol', 1.0e-7)

    ## nutrient
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient')
    add(param_index, param_val, 'lambda_P', lambda_P)
    add(param_index, param_val, 'lambda_A', 0.005)
    add(param_index, param_val, 'lambda_Ph', 0.5)
    add(param_index, param_val, 'D_sigma', D_sigma)
    add(param_index, param_val, 'delta_sigma', 1.)
    add(param_index, param_val, 'chi_c', chi_c)

    ## tumor
    break_points.append(len(param_val))
    break_msg.append('\n# tumor')
    add(param_index, param_val, 'bar_M_P', mobility_P)
    add(param_index, param_val, 'bar_E_phi_T', 0.045)
    add(param_index, param_val, 'epsilon_T', 5.0e-3)

    ## hypoxic
    break_points.append(len(param_val))
    break_msg.append('\n# hypoxic')
    add(param_index, param_val, 'bar_M_H', 0.5 * mobility_P)
    add(param_index, param_val, 'lambda_PH', 1.)
    add(param_index, param_val, 'lambda_HP', 1.)
    add(param_index, param_val, 'lambda_HN', 1.)
    add(param_index, param_val, 'sigma_PH', 0.55)
    add(param_index, param_val, 'sigma_HP', 0.65)
    add(param_index, param_val, 'sigma_HN', 0.44)

    ## necrotic
    break_points.append(len(param_val))
    break_msg.append('\n# necrotic')
    add(param_index, param_val, 'bar_M_N', 0.)

    ## TAF
    break_points.append(len(param_val))
    break_msg.append('\n# TAF')
    add(param_index, param_val, 'D_TAF', 1.0e-2)
    add(param_index, param_val, 'delta_TAF', 1.0)
    add(param_index, param_val, 'lambda_TAF', 1.e+1)

    ## ECM
    break_points.append(len(param_val))
    break_msg.append('\n# ECM')
    add(param_index, param_val, 'lambda_ECM_D', 0.)
    add(param_index, param_val, 'lambda_ECM_P', 0.)
    add(param_index, param_val, 'bar_phi_ECM_P', 0.5)
    add(param_index, param_val, 'chi_h', 0.0)

    ## MDE
    break_points.append(len(param_val))
    break_msg.append('\n# MDE')
    add(param_index, param_val, 'D_MDE', 10.)
    add(param_index, param_val, 'delta_MDE', 0.1)
    add(param_index, param_val, 'lambda_MDE_D', 0.3)
    add(param_index, param_val, 'lambda_MDE_P', 0.4)

    ## flow 1D 
    break_points.append(len(param_val))
    break_msg.append('\n# flow 1D')
    add(param_index, param_val, 'init_vessel_viscosity', 1.)
    add(param_index, param_val, 'vessel_in_pressure', 1.)
    add(param_index, param_val, 'vessel_in_nutrient', 1.)
    add(param_index, param_val, 'vessel_in_nutrient_vein', 0.)
    add(param_index, param_val, 'vessel_blood_density', 1.)
    add(param_index, param_val, 'vessel_D_sigma', vessel_D_sigma)
    add(param_index, param_val, 'osmotic_reflection_coeff', 1.)
    
    ## flow 2D/3D 
    break_points.append(len(param_val))
    break_msg.append('\n# flow 21/3DD')
    add(param_index, param_val, 'tissue_flow_viscosity', 1.)
    add(param_index, param_val, 'tissue_flow_K', tissue_flow_coef)
    add(param_index, param_val, 'tissue_flow_density', 1.)
    
    # coupling strength between 1d-3d pressure
    add(param_index, param_val, 'tissue_flow_L_p', L_p)

    # coupling strength between 1d-3d nutrients
    add(param_index, param_val, 'tissue_nut_L_s', L_s)

    # below is the factor for pressure and nutrient equation. This factor is 
    # multiplied to both sides of equation to improve the coupling.
    add(param_index, param_val, 'assembly_factor_p_t', 1./L_p)
    add(param_index, param_val, 'assembly_factor_c_t', 1./L_s)

    add(param_index, param_val, 'tissue_pressure_bc_val', 0.)
    add(param_index, param_val, 'bc_tissue_pressure_north', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_south', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_east', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_west', 'false')

    # pressure ic (set it equal to lower pressure in artery)
    add(param_index, param_val, 'tissue_pressure_ic_val', vessel_pressures[1])

    ## nutrient ic
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient ic')
    add(param_index, param_val, 'ic_nutrient_value', nutrient_ic_val)

    ## tumor ic
    break_points.append(len(param_val))
    break_msg.append('\n# tumor ic')
    add(param_index, param_val, 'ic_tumor_file', tumor_ic_file)
    gen_tumor_ic_file(file_dir + tumor_ic_file, tumor_ic_type, tumor_ic_center, tumor_ic_radius)

    ## ECM ic
    break_points.append(len(param_val))
    break_msg.append('\n# ECM ic')
    add(param_index, param_val, 'ECM_ic_val', 1.)

    ## MDE ic
    break_points.append(len(param_val))
    break_msg.append('\n# MDE ic')
    add(param_index, param_val, 'MDE_ic_val', 0.0)

    ## nutrient bc
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient bc')
    add(param_index, param_val, 'bc_nutrient_north', 'false')
    add(param_index, param_val, 'bc_nutrient_south', 'false')
    add(param_index, param_val, 'bc_nutrient_east', 'false')
    add(param_index, param_val, 'bc_nutrient_west', 'false')

    ## network 
    break_points.append(len(param_val))
    break_msg.append('\n# network')
    network_input(L, str(pp_tag), vessel_radius, vessel_pressures, network_file, file_dir, vessel_num_refinement, param_index, param_val)


    break_points_new = [-1 for x in range(len(param_val))]
    break_msg_new = ["" for x in range(len(param_val))]
    for i in range(len(break_points)):
        ii = break_points[i]
        break_points_new[ii] = 1
        break_msg_new[ii] = break_msg[i]


    # get longest index name
    length = max(len(x) for x in param_index)
    space = length + 5

    # write to file
    inpf = open(file_dir + input_file,'w')

    # Model
    for i in range(len(param_index)):
        if break_points_new[i] == 1:
            inpf.write("{}\n".format(break_msg_new[i]))

        str_prt = "{0:{space}} = ".format(param_index[i], space = space)
        str_prt += "{0}\n".format(param_val[i])
        inpf.write(str_prt)
    
    inpf.close()


## run sim
L = 2.
pp_tag = 't1'
model_name = 'NetFVFE'
test_name = 'test_nut'
# test_name = 'test_net_tum_2'
output_debug_info = True
advection_active = True
decouple_nutrients = False 
num_elems = 20
final_time = 5.
delta_t = 0.05
dt_output = 10
chi_c = 0.
D_sigma = 1.
nutrient_ic_val = 0.
if test_name != 'test_nut':
    nutrient_ic_val = 0.6
tissue_flow_coef = 1.e-7
vessel_D_sigma = 5.e-3
L_p = 1.e-5
L_s = 1.
vessel_pressures = [100000., 50000., 10000., 20000.]
vessel_radius = 0.05 * L
tumor_ic_type = 1
tumor_ic_center = [0.5 * L, 0.5 * L, 0.5 * L]
tumor_ic_radius = [0.15 * L, 0.15 * L, 0.15 * L]
vessel_num_refinement = 4
lambda_P = 5.
mobility_P = 50.

# file path to run sim
file_dir = pp_tag + '/'
os.system("./run.sh " + model_name + ' ' + pp_tag + ' 0')


# create input file
input(L, file_dir, pp_tag, model_name, test_name, output_debug_info, advection_active, decouple_nutrients, num_elems, final_time, delta_t, dt_output, chi_c, D_sigma, nutrient_ic_val, tissue_flow_coef, vessel_D_sigma, L_p, L_s, vessel_pressures, vessel_radius, tumor_ic_type, tumor_ic_center, tumor_ic_radius, vessel_num_refinement, lambda_P, mobility_P)

# run
os.system("./run.sh " + model_name + ' ' + pp_tag + ' 1')
