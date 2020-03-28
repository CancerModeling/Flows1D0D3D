import os
import numpy as np
import pandas as pd
import csv
import sys
import matplotlib.pyplot as plt
import importlib

import logging

def add(param_index, param_val, index, val):

    param_index.append(index)
    param_val.append(val)


def gen_tumor_ic_file(filename):

    # write to file
    inpf = open(filename,'w')
    inpf.write("type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz\n")

    inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(3, 0.5, 1.0, 0.0, 0.15, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0))
    inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(3, 1.5, 1.0, 0.0, 0.15, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0))

    inpf.close()


def gen_init_network_file(filename):

    # write to file
    inpf = open(filename,'w')
    inpf.write("DGF\n")

    # vertex
    inpf.write("Vertex\n")
    inpf.write("parameters {}\n".format(1))
    inpf.write("{} {} {} {}\n".format(300.0e-6, 300.0e-6, 0.0, 10.0))
    inpf.write("{} {} {} {}\n".format(300.0e-6, 300.0e-6, 385.0e-6, 8.0))
    inpf.write("{} {} {} {}\n".format(85.0e-6, 85.0e-6, 0.0, 3.0))
    inpf.write("{} {} {} {}\n".format(85.0e-6, 85.0e-6, 385.0e-6, 4.0))

    # segments
    inpf.write("#\n")
    inpf.write("SIMPLEX\n")
    inpf.write("parameters {}\n".format(2))
    inpf.write("{} {} {} {}\n".format(0, 1, 10.0e-06, 0.0075))
    inpf.write("{} {} {} {}\n".format(2, 3, 20.0e-06, 0.0075))
    # 
    inpf.write("#\n")
    inpf.write("BOUNDARYDOMAIN\n")
    inpf.write("default {}\n".format(1))

    # 
    inpf.write("#")
    inpf.close()


def network_input(param_index, param_val):

    add(param_index, param_val, 'is_network_active', 'true')
    init_file = 'test_single_line.dgf'
    add(param_index, param_val, 'network_init_file', init_file)
    add(param_index, param_val, 'network_init_refinement', 5)
    add(param_index, param_val, 'vessel_lambda_g', 0.5)
    add(param_index, param_val, 'vessel_R_factor', 1.)
    add(param_index, param_val, 'network_update_interval', 1)
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
    add(param_index, param_val, 'network_discret_cyl_length', 10)
    add(param_index, param_val, 'network_discret_cyl_angle', 10)
    add(param_index, param_val, 'network_coupling_method_theta', 0.0)
    add(param_index, param_val, 'identify_vein_pressure', 8.0)    

    gen_init_network_file(init_file)


def input():

    # Info:
    # Creates input.in for model
    # Creates initial network file
    # Creates tumor initial condition file

    param_index = []
    param_val = []
    break_points = []
    break_msg = []
    
    # model
    L = 385.0e-6
    break_points.append(len(param_val))
    break_msg.append('# model')
    add(param_index, param_val, 'model_name', 'NetFV')
    add(param_index, param_val, 'dimension', 3)
    add(param_index, param_val, 'domain_xmin', 0.)
    add(param_index, param_val, 'domain_xmax', L)
    add(param_index, param_val, 'domain_ymin', 0.)
    add(param_index, param_val, 'domain_ymax', L)
    add(param_index, param_val, 'domain_zmin', 0.)
    add(param_index, param_val, 'domain_zmax', L)    
    add(param_index, param_val, 'assembly_method', 1)
    

    break_points.append(len(param_val))
    break_msg.append('\n# restart')
    add(param_index, param_val, 'restart', 'false')
    # add(param_index, param_val, 'mesh_restart_file', 'mesh.e')
    # add(param_index, param_val, 'solution_restart_file', 'restart_220.e')

    # mesh
    break_points.append(len(param_val))
    break_msg.append('\n# mesh')
    num_elems = 50
    add(param_index, param_val, 'mesh_n_elements', num_elems)

    # time
    break_points.append(len(param_val))
    break_msg.append('\n# time')
    final_t = 1.0
    init_t = 0.
    delta_t = final_t / 10000
    add(param_index, param_val, 'time_step', delta_t)
    add(param_index, param_val, 'initial_time', init_t)
    add(param_index, param_val, 'initial_step', 0)
    add(param_index, param_val, 'final_time', final_t)

    # output
    break_points.append(len(param_val))
    break_msg.append('\n# output')
    dt_output = 50
    add(param_index, param_val, 'perform_output', 'true')
    add(param_index, param_val, 'output_interval', dt_output)
    add(param_index, param_val, 'restart_save', 'false')
    add(param_index, param_val, 'restart_save_interval', 1)

    # solver
    break_points.append(len(param_val))
    break_msg.append('\n# solver')
    add(param_index, param_val, 'linear_solver_max_iter', 250)
    add(param_index, param_val, 'linear_solver_tol', 1.0e-8)
    add(param_index, param_val, 'nonlinear_solver_max_iter', 50)
    add(param_index, param_val, 'nonlinear_solver_tol', 1.0e-8)

    # nutrient
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient')
    add(param_index, param_val, 'lambda_P', 2.)
    add(param_index, param_val, 'lambda_A', 0.005)
    add(param_index, param_val, 'lambda_Ph', 0.5)
    add(param_index, param_val, 'D_sigma', 5.0e-7)
    add(param_index, param_val, 'delta_sigma', 0.01)
    add(param_index, param_val, 'chi_c', 0.)

    # tumor
    break_points.append(len(param_val))
    break_msg.append('\n# tumor')
    add(param_index, param_val, 'bar_M_P', 2.)
    add(param_index, param_val, 'bar_E_phi_T', 0.045)
    add(param_index, param_val, 'epsilon_T', 0.005)

    # hypoxic
    break_points.append(len(param_val))
    break_msg.append('\n# hypoxic')
    add(param_index, param_val, 'bar_M_H', 1.)
    add(param_index, param_val, 'lambda_PH', 1.)
    add(param_index, param_val, 'lambda_HP', 1.)
    add(param_index, param_val, 'lambda_HN', 1.)
    add(param_index, param_val, 'sigma_PH', 0.55)
    add(param_index, param_val, 'sigma_HP', 0.65)
    add(param_index, param_val, 'sigma_HN', 0.44)

    # necrotic
    break_points.append(len(param_val))
    break_msg.append('\n# necrotic')
    add(param_index, param_val, 'bar_M_N', 0.)

    # TAF
    break_points.append(len(param_val))
    break_msg.append('\n# TAF')
    add(param_index, param_val, 'D_TAF', 1.)
    add(param_index, param_val, 'delta_TAF', 0.1)
    add(param_index, param_val, 'lambda_TAF', 10.)

    # ECM
    break_points.append(len(param_val))
    break_msg.append('\n# ECM')
    add(param_index, param_val, 'lambda_ECM_D', 1.)
    add(param_index, param_val, 'lambda_ECM_P', 1.)
    add(param_index, param_val, 'bar_phi_ECM_P', 0.5)
    add(param_index, param_val, 'chi_h', 0.035)

    # MDE
    break_points.append(len(param_val))
    break_msg.append('\n# MDE')
    add(param_index, param_val, 'D_MDE', 10.)
    add(param_index, param_val, 'delta_MDE', 0.1)
    add(param_index, param_val, 'lambda_MDE_D', 0.3)
    add(param_index, param_val, 'lambda_MDE_P', 0.4)

    # flow 1D 
    break_points.append(len(param_val))
    break_msg.append('\n# flow 1D')
    add(param_index, param_val, 'init_vessel_viscosity', 3.5e-3)
    add(param_index, param_val, 'vessel_in_pressure', 1.)
    add(param_index, param_val, 'vessel_in_nutrient', 1.)
    add(param_index, param_val, 'vessel_in_nutrient_vein', 0.)
    add(param_index, param_val, 'vessel_blood_density', 1.)
    add(param_index, param_val, 'vessel_D_sigma', 5.0e-7)
    add(param_index, param_val, 'osmotic_reflection_coeff', 1.)
    

    # flow 2D/3D 
    break_points.append(len(param_val))
    break_msg.append('\n# flow 21/3DD')
    add(param_index, param_val, 'tissue_flow_viscosity', 1.0e-3)
    add(param_index, param_val, 'tissue_flow_K', 1.0e-18)
    add(param_index, param_val, 'tissue_flow_density', 1.)
    add(param_index, param_val, 'tissue_flow_L_p', 1.0e-12)
    add(param_index, param_val, 'tissue_nut_L_s', 3.5e-5)

    add(param_index, param_val, 'coupling_factor_p_t', 1.0e+12)

    add(param_index, param_val, 'tissue_pressure_bc_val', 0.)
    add(param_index, param_val, 'tissue_pressure_ic_val', 0.)
    add(param_index, param_val, 'bc_tissue_pressure_north', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_south', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_east', 'false')
    add(param_index, param_val, 'bc_tissue_pressure_west', 'false')

    # nutrient ic
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient ic')
    add(param_index, param_val, 'ic_nutrient_value', 0.0)

    # tumor ic
    break_points.append(len(param_val))
    break_msg.append('\n# tumor ic')
    add(param_index, param_val, 'ic_tumor_file', 'tum_ic_data.csv')
    gen_tumor_ic_file('tum_ic_data.csv')

    # ECM ic
    break_points.append(len(param_val))
    break_msg.append('\n# ECM ic')
    add(param_index, param_val, 'ECM_ic_domain_type', '')
    add(param_index, param_val, 'ECM_ic_val', 1.)

    # MDE ic
    break_points.append(len(param_val))
    break_msg.append('\n# MDE ic')
    add(param_index, param_val, 'MDE_ic_val', 0.5)

    # nutrient bc
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient bc')
    add(param_index, param_val, 'bc_nutrient_north', 'false')
    add(param_index, param_val, 'bc_nutrient_south', 'false')
    add(param_index, param_val, 'bc_nutrient_east', 'false')
    add(param_index, param_val, 'bc_nutrient_west', 'false')

    # network 
    break_points.append(len(param_val))
    break_msg.append('\n# network')
    network_input(param_index, param_val)


    break_points_new = [-1 for x in xrange(len(param_val))]
    break_msg_new = ["" for x in xrange(len(param_val))]
    for i in xrange(len(break_points)):
        ii = break_points[i]
        break_points_new[ii] = 1
        break_msg_new[ii] = break_msg[i]


    # get longest index name
    length = max(len(x) for x in param_index)
    space = length + 5

    # write to file
    inpf = open('input.in','w')

    # Model
    for i in xrange(len(param_index)):
        if break_points_new[i] == 1:
            inpf.write("{}\n".format(break_msg_new[i]))

        str_prt = "{0:{space}} = ".format(param_index[i], space = space)
        str_prt += "{0}\n".format(param_val[i])
        inpf.write(str_prt)
    
    inpf.close()


if len(sys.argv) > 1:
    select_method = str(sys.argv[1])
    if select_method == 'input':
        input()
