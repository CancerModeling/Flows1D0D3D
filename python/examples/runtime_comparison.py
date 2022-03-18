from matplotlib import pyplot as plt
import numpy as np
import os
import re
import subprocess


def run(log_tau, use_fully_coupled):
    try:
        params = [
            'mpirun', '-np', '1',
            'python', './main_flow.py',
            '--geometry-id', '1',
            '--tau-coup-log2', '3',
            '--tau-log2', f'{log_tau}',
            '--output-folder', 'tmp',
            '--disable-output',
            '--t-3dcoup-start', '0',
            '--t-preiter', '1',
            '--t-end', f'{1.}'
        ]
        if use_fully_coupled:
            params.append('--use-fully-coupled')
        print(params)
        output = subprocess.check_output(params)
    except:
        print('some error :(')
        return float('nan')

    output = output.decode('utf-8')

    regex = r"^elapsed simulation time: (.*)$"

    match = re.search(regex, output, re.MULTILINE)

    if match == None:
        print('Error: Incomplete output')
        return float('nan')
    else:
        return float(match.group(1))


if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '1'

    num_samples = 2

    # log_taus = [6, 5, 4, 3]

    #log_taus = [16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3]
    log_taus = [16, 14, 12, 10, 8, 6, 4, 3]

    taus = [1/2**k for k in log_taus]

    elapsed_time_fully_coupled = np.min([run(16, True) for i in range(num_samples)])
    print(f'elapsed_time_fully_coupled {elapsed_time_fully_coupled}')

    # elapsed_times = [135.19850397109985, 102.90438175201416, 95.01434254646301, 91.54877495765686, 87.79965162277222, 84.6084852218628, 21.64975118637085, 10.549214601516724]

    elapsed_times = []
    for log_tau in log_taus:
        elapsed_time = np.min([run(log_tau, False) for i in range(num_samples)])
        elapsed_times.append(elapsed_time)
        print(elapsed_time)
    print(elapsed_times)

    #elapsed_time_fully_coupled = 37
    #elapsed_times = [17.576059579849243, 8.549134492874146, 3.9048962593078613, 1.300276279449463]

    plt.semilogx(taus, elapsed_time_fully_coupled / np.array(elapsed_times))
    plt.hlines(1, taus[0], taus[-1])
    plt.ylabel('speedup')
    plt.xlabel(r'$\tau$')
    plt.grid(True)
    plt.show()
