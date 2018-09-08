#!/bin/python
from __future__ import print_function
import subprocess
import numpy as np
import os
import re
import sys

import glob
import kinetics_stats_by_temperature as kstats
import edit_config
import create_configs
import grid_params
import pandas as pd

from IMP.npctransport import *

AVOGADRO=6.0221409E+23 # molecules/mole
L_per_A3=1E-27 # Liters per A^3

def get_concentration(n, side_A):
    volume_L = side_A**3 * L_per_A3
    return n / (volume_L * AVOGADRO)

def get_KDs_from_dict(KDs_dict, type0, type1):
    for key, KDs in KDs_dict.items():
        if re.search(' {} - {}$'.format(type0, type1), key):
            return KDs
        if re.search(' {} - {}$'.format(type1, type0), key):
            return KDs
    assert(False)

def get_stats_entry_for_output_file(output_file, type0='fg0', type1='kap20'):
    assert(re.match('fg', type0) and re.match('kap', type1))
    [KDs_dicts, KDs_dicts_new, outputs]= kstats.do_all_stats([output_file],
                                                             0.1e-6,
                                                             verbose=False,
                                                             return_outputs=True)
    KDs_dict= get_KDs_from_dict(KDs_dicts, type0, type1)
    KDs_dict_new= get_KDs_from_dict(KDs_dicts_new, type0, type1)
    output= outputs[0]
    assign= output.assignment
    box_side= assign.box_side.value
    assert( assign.fgs[0].type == type0 and
            assign.floaters[0].type == type1 )
    n_kaps= assign.floaters[0].number.value
    n_fgs= assign.fgs[0].number.value
    C_floaters= get_concentration(n_kaps, box_side)
    C_fgs= get_concentration(n_fgs, box_side)
    kap_range= None
    kap_k= None
    for interaction in assign.interactions:
        if interaction.type0==type0 and interaction.type1==type1:
            kap_range= assign.interactions[0].interaction_range.value
            kap_k= assign.interactions[0].interaction_k.value
    assert(kap_range is not None and kap_k is not None)
    params=grid_params.GridParams()
    input_site_KD= create_configs.compute_KD_from_k_r_and_quadratic_model(kap_k,
                                                                          kap_range,
                                                                          params.QuadraticLogKDParams)
    entry= { 'box_side_A': box_side,
             'n_kaps': n_kaps,
             'n_fgs': n_fgs,
             'fg_seq': fg_seq,
             'kap_valency': assign.floaters[0].interactions.value,
             'KD_sites_M_from_konoff': KDs_dict[0],
             'fbound_sites_fg': KDs_dict_new['fbound_sitesA'],
             'fbound_sites_kap': KDs_dict_new['fbound_sitesB'],
             'KD_sites_M_from_fbound': KDs_dict_new['KD_sites'],
             'input_site_KD': input_site_KD,
             'k_on_per_ns_per_missing_ss_contact': KDs_dict_new['k_on_per_ns_per_missing_ss_contact'],
             'k_off_per_ns': KDs_dict_new['k_off_per_ns_per_ss_contact'],
             'fbound_chains_fg': KDs_dict_new['fbound_chainsA'],
             'fbound_chains_kap': KDs_dict_new['fbound_chainsB'],
             'KD_chains_M_from_fbound': KDs_dict_new['KD_chains'],
             'energy_kcal_per_mole': KDs_dicts_new['energy'][0]

         }
    return entry


def put_data_in_csv_file(data, filename):
    ''' append if file exists already '''
    if(len(data)==0):
        return
    df= pd.DataFrame(data=data)
    df['kap_C_M']= get_concentration(df['n_kaps'], df['box_side_A'])
    df['fg_C_M']= get_concentration(df['n_fgs'], df['box_side_A'])
    df['k_on_per_ns_per_M']= df['k_on_per_ns_per_missing_ss_contact'] * AVOGADRO * df['box_side_A']**3 * L_per_A3
    if os.path.exists(filename):
        open_mode='a'
        is_header= False
    else:
        open_mode='w'
        is_header= True
    with open(filename, open_mode) as f:
        df.to_csv(path_or_buf=f,
                  sep=' ',
                  header= is_header)
#    print(df.head().to_string())


if __name__ != "__main__":
    sys.exit(-1)
if len(sys.argv)>1:
    stats_csv_file=sys.argv[1]
else:
    stats_csv_file='stats.csv'
if os.path.exists(stats_csv_file):
    os.remove(stats_csv_file)
data= []
MAX_DATA_ENTRIES= 500
for fg_seq in ['FSSSSS', 'FFSSSS', 'FFFSSS', 'FFFFSS', 'FFFFFF']:
    output_files=glob.glob('Output/{}/*.pb'.format(fg_seq))
    print("Processing {:d} output files".format(len(output_files)))
    for output_file in output_files:
        try:
            stats_entry= get_stats_entry_for_output_file(output_file, 'fg0', 'kap20')
            data.append(stats_entry)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            print("Skipped {}".format(output_file))
            raise
        print("Processed {}".format(output_file))
        if len(data)>MAX_DATA_ENTRIES:
            put_data_in_csv_file(data, stats_csv_file)
            data=[]
put_data_in_csv_file(data, stats_csv_file)

print("Finished succesfully")
