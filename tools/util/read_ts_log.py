import os
import sys
import csv
import numpy as np
import pandas as pd

def read_ts(in_file, get_num_sys):

    # read dgf file
    inpf = open(in_file)

    setup_t = 0.
    
    num_sys = get_num_sys
    sys_names = []

    num_step_data = 6 + num_sys
    step_data_names = ['time', 'solve time', 'pres solve time', 'net update time', 'nonlinear iter', 'pres nonlinear iter']

    step_data = []
    
    block_counter = 0

    while True:
        a = inpf.readline().split()
        a

        if a is None or len(a) == 0 or a is EOFError:
            break
        elif a[0] == '#':
            block_counter += 1
        elif block_counter == 1:
            # read setup time
            setup_t = float(a[1])
            block_counter += 1
        elif block_counter == 3:
            # read sys names
            n = int(a[1])

            if num_sys == 0:
                return n, setup_t, sys_names, step_data_names, step_data

            block_counter += 1
        elif block_counter == 4:
            
            # read system name
            sys_names.append(str(a[0]))
            step_data_names.append(str(a[0]))
        elif block_counter == 5:
            # skip sys assembly time data
            block_counter += 2
        elif block_counter == 8:

            n = int(a[1])
            if n != num_step_data:
                print('Error: StepLog has {} data while calculation suggest it should have {} data'.format(n, num_step_data))
            block_counter += 1
        elif block_counter == 9:
            # step log header (skip reading)
            block_counter += 1
        elif block_counter == 10:
            data = []
            for i in range(num_step_data):
                data.append(float(a[i]))

            step_data.append(data)

    sys_names = np.array(sys_names)
    step_data = np.array(step_data)
    inpf.close()

    return num_sys, setup_t, sys_names, step_data_names, step_data


def get_n(in_file):

    num_sys, setup_t, sys_names, step_data_names, step_data = read_ts(in_file, 0)
    return num_sys


def test():
    in_file = 'test.txt'

    # first read num sys
    num_sys = get_n(in_file)

    # now read all data
    n, setup_t, sys_names, step_data_names, step_data = read_ts(in_file, num_sys)

    print('setup_t: {} ms = {} mili-sec'.format(setup_t, setup_t / 1.e3))
    print('num_sys: {}'.format(num_sys))
    s = ''
    for i in range(num_sys):
        s += sys_names[i] + ' '
    print(s)
    s= ''
    print('num step data names: {}'.format(len(step_data_names)))
    for i in range(len(step_data_names)):
        s = step_data_names[i] + ' '
    print(s)
    print('num step data: {}'.format(len(step_data)))
    for i in range(len(step_data)):
        s = ' '
        d = step_data[i]
        for j in range(len(d)):
            s += str(d[j]) + ' '

        print(s)

test()
