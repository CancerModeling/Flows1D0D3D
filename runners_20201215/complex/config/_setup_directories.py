import copy
import os
import pathlib

from config._gen_init_network_file  import gen_init_network_file
from config._gen_tumor_ic_file import gen_tumor_ic_file


def setup_directories(dps):
  for dp in dps:
    print('setting up directories ... ')

    pathlib.Path('sim').mkdir(parents=True, exist_ok=True)
    pathlib.Path('sim/' + dp.model_name).mkdir(parents=True, exist_ok=True)
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
        os.system('cp ' + 'data/' + dp.network_init_file + ' ' + dp.run_path)
