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

    log_taus = [16, 14, 12, 10, 8, 6, 4, 3]

    taus = [1/2**k for k in log_taus]

    #elapsed_time_fully_coupled = np.min([run(16, True) for i in range(num_samples)])
    #print(f'elapsed_time_fully_coupled {elapsed_time_fully_coupled}')

    # elapsed_times = []
    # for log_tau in log_taus:
    #     elapsed_time = np.min([run(log_tau, False) for i in range(num_samples)])
    #     elapsed_times.append(elapsed_time)
    #     print(elapsed_time)
    # print(elapsed_times)

    elapsed_time_fully_coupled = 264.72977447509766
    elapsed_times = [149.2980239391327, 59.18277359008789, 29.591812133789062, 19.72035551071167, 14.324291229248047, 11.722833633422852, 10.837181091308594, 10.617775201797485]

    plt.semilogx(taus, elapsed_time_fully_coupled / np.array(elapsed_times), label='breast model')
    #plt.hlines(1, taus[0], taus[-1], linestyle='--', label='fully coupled model')
    plt.legend()
    plt.ylabel('speed up fully coupled model')
    plt.xlabel(r'$\tau$')
    plt.grid(True)
    plt.show()
