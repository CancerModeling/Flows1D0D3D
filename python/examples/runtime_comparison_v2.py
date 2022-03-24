import numpy as np
import os
import re
import subprocess
from matplotlib import pyplot as plt


folder_bin = '/home/wagneran/Flows1D0D3D_v2/cmake-build-release/bin/macrocirculation/'
folder_python_examples = '/home/wagneran/Flows1D0D3D_v2/python/examples'
folder_mesh = '/home/wagneran/Flows1D0D3D_v2/data/1d-meshes/'


def run_bin(num_processors, folder, exe, parameters):
    os.environ['OMP_NUM_THREADS'] = '1'
    params = [
        'mpirun', '-np', f'{num_processors}',
        os.path.join(folder, exe),
    ] + parameters
    print(params)
    output = subprocess.check_output(params, cwd=folder)

    output = output.decode('utf-8')

    regex = r"^time = (.*) s$"

    match = re.search(regex, output, re.MULTILINE)

    if match == None:
        print('Error: Incomplete output')
        return float('nan')
    else:
        return float(match.group(1))


def run_py(num_processors, folder, exe, parameters):
    os.environ['OMP_NUM_THREADS'] = '1'
    params = [
                 'mpirun', '-np', f'{num_processors}',
                 'python',
                 os.path.join(folder, exe),
             ] + parameters
    print(params)
    output = subprocess.check_output(params, cwd=folder)

    output = output.decode('utf-8')

    regex = r"^time = (.*) s$"

    match = re.search(regex, output, re.MULTILINE)

    if match == None:
        print('Error: Incomplete output')
        return float('nan')
    else:
        return float(match.group(1))


def runtimes_breast_geometry(num_processors, num_repetitions):
    elapsed_times = []
    num_processors = [1, 2, 4, 6]

    exe = 'MacrocirculationLinearFlowBreastGeometry'
    params = ['--t-end', '0.1', '--tau-out', '1']

    for n in num_processors:
        elapsed = []
        for _ in range(num_repetitions):
            elapsed.append(run_bin(n, folder_bin, exe, params))
        elapsed_times.append(np.min(elapsed))

    return elapsed_times


def runtimes_cow_geometry(num_processors, num_repetitions):
    elapsed_times = []
    num_processors = [1, 2, 4, 6]

    exe = 'MacrocirculationNonlinear1DSolver'
    mesh_file = os.path.join(folder_mesh, '33-vessels-refined.json')
    params = ['--t-end', '0.05',
              '--tau-out', '1',
              '--tau', f'{1/2**17}',
              '--update-interval-transport', '10000000',
              '--mesh-file', mesh_file
              ]

    for n in num_processors:
        elapsed = []
        for _ in range(num_repetitions):
            elapsed.append(run_bin(n, folder_bin, exe, params))
        elapsed_times.append(np.min(elapsed))

    return elapsed_times


def runtimes_3d_geometry(num_processors, num_repetitions):
    elapsed_times = []
    num_processors = [1, 2, 4, 6]

    exe = 'run_implicit_pressure_solver.py'
    params = ['--tip-pressures-input-file', 'tip_pressures1.json', '--refinements', '0']

    for n in num_processors:
        elapsed = []
        for _ in range(num_repetitions):
            elapsed.append(run_py(n, folder_python_examples, exe, params))
        elapsed_times.append(np.min(elapsed))
        print(elapsed_times)

    return elapsed_times


if __name__ == '__main__':
    num_processors = [1, 2, 4, 6]

    #elapsed_times = runtimes_breast_geometry(num_processors, 3)
    #elapsed_times = [44.7567, 25.6128, 11.8786, 8.68186]

    #elapsed_times = runtimes_cow_geometry(num_processors, 3)
    #elapsed_times = [16.8023, 9.04702, 5.28342, 4.05778]

    elapsed_times = runtimes_3d_geometry(num_processors, 1)
    #elapsed_times = [4.5778210163116455, 3.0047760009765625, 2.2704617977142334, 2.009453535079956]

    print(num_processors)
    print(elapsed_times)

    plt.plot(num_processors, elapsed_times[0] / np.array(elapsed_times), 'o-')
    #plt.legend()
    plt.grid(True)
    plt.xlabel('number of processors')
    plt.ylabel('speed up')
    plt.show()

