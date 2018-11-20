#!/bin/python
from __future__ import print_function
import subprocess
import numpy as np
import os
import re
import sys

import kinetics_stats_by_temperature as kstats
import edit_config
#import create_configs
import grid_params

def get_KDs_from_dict(KDs_dict, type0, type1):
    for key, KDs in KDs_dict.items():
        if re.search(' {} - {}$'.format(type0, type1), key):
            return KDs
        if re.search(' {} - {}$'.format(type1, type0), key):
            return KDs
    raise RuntimeError("Couldn't find KD in output")


if __name__ != "__main__":
    sys.exit(-1)

imp_folder= os.environ['HOME'] + '/imp_git/fast/'
imp= imp_folder + '/setup_environment.sh'
fg_sim= imp + " " + imp_folder + "/bin/fg_simulation"
sim_time_factor=0.25
config_file= sys.argv[1]
site_site_KD= float(sys.argv[2])
params=grid_params.GridParams()
kap_range=5.5
Ks= grid_params.compute_k_from_r_KD_and_quadratic_model\
    ( kap_range,
      site_site_KD,
      params.QuadraticLogKDParams )
assert(len(Ks)>0)
kap_k= Ks[0]
edit_config.do_edit(config_file,'config.pb', kap_k=kap_k)
initial_cmd= fg_sim + " --configuration config.pb --output output.pb --short_sim_factor 0.000001"
subprocess.check_call(initial_cmd, shell=True)
round=1
while True:
    print("Round {}".format(round))
    print("Sim time factor {}".format(sim_time_factor))
    os.rename('output.pb', 'prev_output.pb')
    restart_cmd= fg_sim + " --restart prev_output.pb --output output.pb --short_sim_factor {sim_time_factor}"\
        .format(sim_time_factor= sim_time_factor)
    subprocess.check_call(restart_cmd, shell=True)
    try:
        kstats.do_all_stats(['output.pb'], 0.5e-6, verbose=False)
        [KDs_dict, _] =         kstats.do_all_stats(['output.pb'], 0.0001e-6, verbose=False)
        if KDs_dict is not None:
            KDs= get_KDs_from_dict(KDs_dict,'fg0','kap20')
            print(KDs_dict)
            pretty_KDs= [kstats.pretty_molarity(KD) for KD in KDs]
            print("KDs round #{}: {}  ( {} - {} )".format(round, *pretty_KDs))
            LKD= np.log10(KDs)
            print("Difference in log space {:.3f} {:.3f}".format(LKD[0]-LKD[1],LKD[2]-LKD[0]))
            if LKD[0]-LKD[1]<np.log10(1.2):
                break
    except RuntimeError as e:
        print("Runtime error exception during KD optimization: '{}'".format(e))
        print("Warning couldn't get stats at round {}".format(round))
    sim_time_factor= sim_time_factor*1.25
    round= round+1

print("Finished succesfully")
