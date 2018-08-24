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
    output= Output()
    params=grid_params.GridParams()
    with open(output_file, "rb") as f:
        output.ParseFromString(f.read())
    assign= output.assignment
    box_side= assign.box_side.value
    assert( assign.fgs[0].type=type0 &&
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
    input_site_KD= create_configs.compute_KD_from_k_r_and_quadratic_model(kap_k,
                                                                          kap_range,
                                                                          params.QuadraticLogKDParams)
    [KDs_dicts, KDs_dicts_new]= kstats.do_all_stats([output_file], 0.1e-6, verbose=False)
    KDs_dict= get_KDs_from_dict(KDs_dicts, 'fg0', 'kap20')
    KDs_dict_new= get_KDs_from_dict(KDs_dicts_new, 'fg0', 'kap20')
    print(KDs_dict_new)
    entry= { 'box_side_A': box_side,
             'n_kaps': n_kaps,
             'n_fgs': n_fgs,
             'fg_seq': fg_seq,
             'kap_valency': assign.floaters[0].interactions.value,
             'KD': KDs_dict[0],
             'fbound_sites_A': KDs_dict_new['fboundA'],
             'fbound_sites_B': KDs_dict_new['fboundB'],
             'KD_sites': KDs_dict_new['KD'],
             'input_site_KD': input_site_KD,
             'k_on_per_ns_per_missing_ss_contact': KDs_dict_new['k_on_per_ns_per_missing_ss_contact'],
             'k_off_per_ns': KDs_dict_new['k_off_per_ns_per_ss_contact']
         }
    return entry


if __name__ != "__main__":
    sys.exit(-1)
for fg_seq in ['FFFFSS', 'FSSSSS']:
    output_files=glob.glob('Output/{}/*.pb'.format(fg_seq))
    print("processing {:d} output files".format(len(output_files)))
    data= []
    for output_file in output_files:
        try:
            stats_entry= get_stats_entry_for_output_file(output_file, 'fg0', 'kap20')
            data.append(stats_entry)
        except:
            print("Skipping file {}".format(output_file))
            raise
        print("Processed {}".format(output_file))
        break
    columns= ['box_side_A', 'n_kaps', 'n_fgs', 'fg_seq', 'kap_valency',
              'KD', 'fbound_sites_A', 'fbound_sites_B', 'KD_sites']
    df= pd.DataFrame(columns=columns)
    df= df.append(data, ignore_index=True)
    df['kap_C']= get_concentration(df['n_kaps'], df['box_side_A'])
    df['fg_C']= get_concentration(df['n_fgs'], df['box_side_A'])
    df['k_on_per_ns_per_M']:
    print(df.head().to_string())

print("Finished succesfully")
