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



def input(L, file_dir, pp_tag, model_name, output_debug_info, num_elems, final_time, delta_t, dt_output, chi_c, D_sigma, nutrient_ic_val,  tumor_ic_type, tumor_ic_center, tumor_ic_radius, lambda_P, mobility_P):

    # Info:
    # Creates input.in for model
    # Creates tumor initial condition file

    param_index = []
    param_val = []
    break_points = []
    break_msg = []

    # files
    # input_file = 'input_' + str(pp_tag) + '.in'
    input_file = 'input.in'
    tumor_ic_file = 'tum_ic_data_' + str(pp_tag) + '.csv'
    
    ## model
    break_points.append(len(param_val))
    break_msg.append('# model')
    
    # specify model such as NetFVFE, NetFVFE, NetFC, AvaFV, TwoSpecies
    add(param_index, param_val, 'model_name', model_name)

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
    add(param_index, param_val, 'D_sigma', D_sigma)
    add(param_index, param_val, 'chi_c', chi_c)

    ## tumor
    break_points.append(len(param_val))
    break_msg.append('\n# tumor')
    add(param_index, param_val, 'bar_M_P', mobility_P)
    add(param_index, param_val, 'bar_E_phi_T', 0.045)
    add(param_index, param_val, 'epsilon_T', 5.0e-3)

    ## nutrient ic
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient ic')
    add(param_index, param_val, 'ic_nutrient_value', nutrient_ic_val)

    ## tumor ic
    break_points.append(len(param_val))
    break_msg.append('\n# tumor ic')
    add(param_index, param_val, 'ic_tumor_file', tumor_ic_file)
    gen_tumor_ic_file(file_dir + tumor_ic_file, tumor_ic_type, tumor_ic_center, tumor_ic_radius)

    ## nutrient bc
    break_points.append(len(param_val))
    break_msg.append('\n# nutrient bc')
    add(param_index, param_val, 'bc_nutrient_north', 'false')
    add(param_index, param_val, 'bc_nutrient_south', 'false')
    add(param_index, param_val, 'bc_nutrient_east', 'false')
    add(param_index, param_val, 'bc_nutrient_west', 'false')


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
pp_tag = 'test'
model_name = 'TwoSpecies'
output_debug_info = True
num_elems = 20
final_time = 5.
delta_t = 0.05
dt_output = 10
chi_c = 0.
D_sigma = 1.
nutrient_ic_val = 0.5
tumor_ic_type = 1
tumor_ic_center = [0.5 * L, 0.5 * L, 0.5 * L]
tumor_ic_radius = [0.15 * L, 0.15 * L, 0.15 * L]
lambda_P = 5.
mobility_P = 50.

# file path to run sim
file_dir = pp_tag + '/'
os.system("./run.sh " + model_name + ' ' + pp_tag + ' 0')


# create input file
input(L, file_dir, pp_tag, model_name, output_debug_info, num_elems, final_time, delta_t, dt_output, chi_c, D_sigma, nutrient_ic_val,  tumor_ic_type, tumor_ic_center, tumor_ic_radius, lambda_P, mobility_P)

# run
os.system("./run.sh " + model_name + ' ' + pp_tag + ' 1')
