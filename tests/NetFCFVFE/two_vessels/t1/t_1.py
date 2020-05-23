import numpy as np
import sys

def add(param_index, param_val, index, val):

    param_index.append(index)
    param_val.append(val)


def gen_tumor_ic_file(L, filename):

    # write to file
    inpf = open(filename,'w')
    inpf.write("type, cx, cy, cz, tum_rx, tum_ry, tum_rz, hyp_rx, hyp_ry, hyp_rz\n")

    # type = 1 -- spherical tumor core
    # type = 3 -- spherical tumor core and then spherical hypoxic core
    # type = 5 -- spherical tumor core (sharp)
    tum_ic_type = 1
    inpf.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(tum_ic_type, 0.5*L, 0.5*L, 0.5*L, 0.15*L, 0.15*L, 0.15*L, 0.15*L, 0.15*L, 0.15*L))

    inpf.close()


def get_pressure_in_vessel():

    P_1 = 100000.
    P_2 = 50000.
    P_3 = 10000.
    P_4 = 20000.

    return np.array([P_1, P_2, P_3, P_4])


def gen_init_network_file(L, filename):

    # write to file
    inpf = open(filename,'w')
    inpf.write("DGF\n")

    R = 0.05 * L
    
    pressures = get_pressure_in_vessel()

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

    return pressures[1]


def network_input(L, param_index, param_val):

    add(param_index, param_val, 'is_network_active', 'true')
    
    # network file
    init_file = 'two_vessels.dgf'
    add(param_index, param_val, 'network_init_file', init_file)
    add(param_index, param_val, 'network_init_refinement', 4)

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
    
    # generate network file
    P_2 = 0.9999 * gen_init_network_file(L, init_file)

    # to identify veins so that we can apply correct bc
    add(param_index, param_val, 'identify_vein_pressure', P_2)


def input():

    # Info:
    # Creates input.in for model
    # Creates initial network file
    # Creates tumor initial condition file

    param_index = []
    param_val = []
    break_points = []
    break_msg = []
    
    ## model
    L = 2.
    break_points.append(len(param_val))
    break_msg.append('# model')
    
    # specify model such as NetFVFE, NetFVFE, NetFC, AvaFV
    add(param_index, param_val, 'model_name', 'NetFCFVFE')

    # specify test (if any) which solves sub-system
    # disable line below if running full system or specify empty string ''
    add(param_index, param_val, 'test_name', 'test_net_tum_2')

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
    add(param_index, param_val, 'advection_active', 'true')
    add(param_index, param_val, 'network_decouple_nutrients', 'true')

    # control parameters for 1d-3d coupling
    add(param_index, param_val, 'network_discret_cyl_length', 20)
    add(param_index, param_val, 'network_discret_cyl_angle', 20)
    add(param_index, param_val, 'network_compute_elem_weights', 'true')
    add(param_index, param_val, 'network_coupling_method_theta', 1.0)

    # set below to reasonable value such as 1, 4, 10 if want to grow network
    add(param_index, param_val, 'network_update_interval', 100000)
    
    ## restart info
    break_points.append(len(param_val))
    break_msg.append('\n# restart')
    add(param_index, param_val, 'restart', 'false')
    # add(param_index, param_val, 'mesh_restart_file', 'mesh.e')
    # add(param_index, param_val, 'solution_restart_file', 'restart_220.e')

    ## mesh
    break_points.append(len(param_val))
    break_msg.append('\n# mesh')
    num_elems = 20
    add(param_index, param_val, 'mesh_n_elements', num_elems)

    ## time
    break_points.append(len(param_val))
    break_msg.append('\n# time')
    final_t = 5.0
    init_t = 0.
    delta_t = 0.05
    add(param_index, param_val, 'time_step', delta_t)
    add(param_index, param_val, 'initial_time', init_t)
    add(param_index, param_val, 'initial_step', 0)
    add(param_index, param_val, 'final_time', final_t)

    ## output
    break_points.append(len(param_val))
    break_msg.append('\n# output')
    total_outputs = 20
    dt_output = int(np.floor(final_t / delta_t) / total_outputs)
    if dt_output < 1:
        dt_output = 1
    add(param_index, param_val, 'perform_output', 'true')
    add(param_index, param_val, 'output_interval', dt_output)
    add(param_index, param_val, 'restart_save', 'false')
    add(param_index, param_val, 'restart_save_interval', 1)

    ## solver
    break_points.append(len(param_val))
    break_msg.append('\n# solver')
    add(param_index, param_val, 'linear_solver_max_iter', 250)
    add(param_index, param_val, 'linear_solver_tol', 1.0e-7)
    add(param_index, param_val, 'nonlinear_solver_max_iter', 50)
    add(param_index, param_val, 'nonlinear_solver_tol', 1.0e-7)

    ## nutrient
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient')
    add(param_index, param_val, 'lambda_P', 5.)
    add(param_index, param_val, 'lambda_A', 0.005)
    add(param_index, param_val, 'lambda_Ph', 0.5)
    add(param_index, param_val, 'D_sigma', 1.)
    add(param_index, param_val, 'delta_sigma', 1.)
    add(param_index, param_val, 'chi_c', 0.)

    ## tumor
    break_points.append(len(param_val))
    break_msg.append('\n# tumor')
    add(param_index, param_val, 'bar_M_P', 50.)
    add(param_index, param_val, 'bar_E_phi_T', 0.045)
    add(param_index, param_val, 'epsilon_T', 5.0e-3)

    ## hypoxic
    break_points.append(len(param_val))
    break_msg.append('\n# hypoxic')
    add(param_index, param_val, 'bar_M_H', 25.)
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
    add(param_index, param_val, 'vessel_D_sigma', 5.0e-3)
    add(param_index, param_val, 'osmotic_reflection_coeff', 1.)
    
    ## flow 2D/3D 
    break_points.append(len(param_val))
    break_msg.append('\n# flow 21/3DD')
    add(param_index, param_val, 'tissue_flow_viscosity', 1.)
    add(param_index, param_val, 'tissue_flow_K', 1.e-7)
    add(param_index, param_val, 'tissue_flow_density', 1.)
    
    # coupling strength between 1d-3d pressure
    L_p = 1.0e-5
    add(param_index, param_val, 'tissue_flow_L_p', L_p)

    # coupling strength between 1d-3d nutrients
    L_s = 1.
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
    pressures = get_pressure_in_vessel()
    add(param_index, param_val, 'tissue_pressure_ic_val', pressures[1])

    ## nutrient ic
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient ic')
    add(param_index, param_val, 'ic_nutrient_value', 0.6)

    ## tumor ic
    break_points.append(len(param_val))
    break_msg.append('\n# tumor ic')
    add(param_index, param_val, 'ic_tumor_file', 'tum_ic_data.csv')
    gen_tumor_ic_file(L, 'tum_ic_data.csv')

    ## ECM ic
    break_points.append(len(param_val))
    break_msg.append('\n# ECM ic')
    add(param_index, param_val, 'ECM_ic_domain_type', '')
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
    network_input(L, param_index, param_val)


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
    inpf = open('input.in','w')

    # Model
    for i in range(len(param_index)):
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
