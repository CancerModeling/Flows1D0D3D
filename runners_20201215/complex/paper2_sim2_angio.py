import sys

from seeds import seeds
from config._default_sim_params  import RatbrainSimParams as DefaultSimParams
from config._setup_directories import setup_directories 
from config._call_executibles import call_executibles 


class SimParams(DefaultSimParams):
    def __init__(self):
        super(SimParams, self).__init__()
        # Play with parameters here... 


def setup_params():
  num_elems = 55
  dps = []
  log_normal_mean = 1.5
  network_update_interval = 10
  scale = 0.0025
  lower_bound = 0.1
  upper_bound = 0.9
  num_eigenfunctions = 17
  
  for i in range(0, 10):
    # create sim params
    dp = SimParams()
    dp.seed = seeds[i] 
    dp.network_update_interval = network_update_interval
    dp.E_phi_T = 0.03
    
    dp.num_elems = num_elems
    dp.log_normal_mean = log_normal_mean
    dp.pp_tag = 'complex_with_angio_explicit_0p03_psi_0p0025_noise_{}_{}'.format(num_elems, i)
    dp.lambda_TAF_deg = 0 
    dp.sigma_HTAF = 0 
    dp.delta_t = 0.0025
    dp.dt_output = 20
    dp.final_time = 8.
    dp.run_path = 'sim/' + dp.pp_tag + '/' 

    dp.hyp_noise_num_eigenfunctions = 17
    dp.hyp_noise_seed = dp.seed
    dp.hyp_noise_scale = scale
    dp.hyp_noise_lower_bound = lower_bound
    dp.hyp_noise_upper_bound = upper_bound

    dp.pro_noise_num_eigenfunctions = 17 
    dp.pro_noise_seed = dp.seed + 1
    dp.pro_noise_scale = scale
    dp.pro_noise_lower_bound = lower_bound
    dp.pro_noise_upper_bound = upper_bound

    dp.scheme_name = 'solve_explicit'
    
    dp.coupled_1d3d = True
    dp.n_mpi = 0
    
    # dp.coupled_1d3d = False
    # dp.n_mpi = 12

    dps.append(dp)
  return dps


if len(sys.argv) > 1:
  if str(sys.argv[1]) == 'run':
    dps = setup_params()
    setup_directories(dps)
    call_executibles(dps, run_screen=False)
  elif str(sys.argv[1]) == 'run_with_screen':
    dps = setup_params()
    setup_directories(dps)
    call_executibles(dps, run_screen=True)
  elif str(sys.argv[1]) == 'prepare':
    dps = setup_params()
    setup_directories(dps)
  else:
    print('?')
else:
    print('To run sim, run script with argument "run"')

