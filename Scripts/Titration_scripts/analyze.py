#!/bin/python
from __future__ import print_function
import subprocess
import numpy as np
import os
import os.path
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

def get_stats_entry_for_output_file(output_file, type0='fg0', type1='kap20',
                                    stats_min_time_sec= 0.1e-6,
                                    is_new_tf_stats= False):
    assert(re.match('fg', type0) and re.match('kap', type1))
    [KDs_dicts, KDs_dicts_new, outputs]= kstats.do_all_stats([output_file],
                                                             stats_min_time_sec,
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
    kap_valency= assign.floaters[0].interactions.value

    # Get back to new for all if using corrected runs, in old runs,
    # wrong stats whenever more than 1 interaction sites on second
    # interactor (which is kap in our case):
    if is_new_tf_stats:
        SITES_B_VERSION= 'NEW'
    else:
        SITES_B_VERSION= ('OLD' if kap_valency>=2 else 'NEW')

    entry= { 'box_side_A':                     box_side,
             'n_kaps':                         n_kaps,
             'n_fgs':                          n_fgs,
             'fg_seq':                         fg_seq,
             'kap_valency':                    kap_valency,
             'KD_sites_M_from_konoff':         KDs_dict[0],
             'KD_sites_M_from_konoff_lbound':  KDs_dict[1],
             'KD_sites_M_from_konoff_ubound':  KDs_dict[2],
             'fbound_sites_fg':                KDs_dict_new['fbound_sitesA_NEW'],
             'fbound_sites_fg_lbound':         KDs_dict_new['fbound_sitesA_lbound_NEW'],
             'fbound_sites_fg_ubound':         KDs_dict_new['fbound_sitesA_ubound_NEW'],
             'fbound_sites_kap':               KDs_dict_new['fbound_sitesB_{}'.format(SITES_B_VERSION)],
             'fbound_sites_kap_lbound':        KDs_dict_new['fbound_sitesB_lbound_{}'.format(SITES_B_VERSION)],
             'fbound_sites_kap_ubound':        KDs_dict_new['fbound_sitesB_ubound_{}'.format(SITES_B_VERSION)],
             'KD_sites_M_from_fbound':         KDs_dict_new['KD_sites_{}'.format(SITES_B_VERSION)],
             'KD_sites_M_from_fbound_lbound':  KDs_dict_new['KD_sites_lbound_{}'.format(SITES_B_VERSION)],
             'KD_sites_M_from_fbound_ubound':  KDs_dict_new['KD_sites_ubound_{}'.format(SITES_B_VERSION)],
             'input_site_KD':                  input_site_KD,
             'k_on_per_ns_per_missing_ss_contact':        KDs_dict_new['k_on_per_ns_per_missing_ss_contact'],
             'k_on_per_ns_per_missing_ss_contact_lbound': KDs_dict_new['k_on_per_ns_per_missing_ss_contact_lbound'],
             'k_on_per_ns_per_missing_ss_contact_ubound': KDs_dict_new['k_on_per_ns_per_missing_ss_contact_ubound'],
             'k_off_per_ns':                   KDs_dict_new['k_off_per_ns_per_ss_contact'],
             'k_off_per_ns_lbound':            KDs_dict_new['k_off_per_ns_per_ss_contact_lbound'],
             'k_off_per_ns_ubound':            KDs_dict_new['k_off_per_ns_per_ss_contact_ubound'],
             'fbound_chains_fg':               KDs_dict_new['fbound_chainsA'],
             'fbound_chains_fg_lbound':        KDs_dict_new['fbound_chainsA_lbound'],
             'fbound_chains_fg_ubound':        KDs_dict_new['fbound_chainsA_ubound'],
             'fbound_chains_kap':              KDs_dict_new['fbound_chainsB'],
             'fbound_chains_kap_lbound':       KDs_dict_new['fbound_chainsB_lbound'],
             'fbound_chains_kap_ubound':       KDs_dict_new['fbound_chainsB_ubound'],
             'fbound_chains_kap_NEW':              KDs_dict_new['fbound_chainsB_NEW'],
             'fbound_chains_kap_NEW_lbound':       KDs_dict_new['fbound_chainsB_lbound_NEW'],
             'fbound_chains_kap_NEW_ubound':       KDs_dict_new['fbound_chainsB_ubound_NEW'],
             'KD_chains_M_from_fbound':        KDs_dict_new['KD_chains'],
             'KD_chains_M_from_fbound_lbound': KDs_dict_new['KD_chains_lbound'],
             'KD_chains_M_from_fbound_ubound': KDs_dict_new['KD_chains_ubound'],
             'energy_kcal_per_mole':           KDs_dicts_new['energy'][0],
             'energy_kcal_per_mole_lbound':    KDs_dicts_new['energy'][0] - 1.96*KDs_dicts_new['energy'][1],
             'energy_kcal_per_mole_ubound':    KDs_dicts_new['energy'][0] + 1.96*KDs_dicts_new['energy'][1],
             'output_file':                    output_file
         }
    for key, value in KDs_dicts_new.iteritems():
        if re.match("Rg", key):
            entry[key]= value
        if re.match("dmax", key):
            entry[key]= value
    return entry


def append_data_to_csv_file(data, filename):
    ''' append if file exists already '''
    if(len(data)==0):
        return
    df= pd.DataFrame(data=data)
    df['kap_C_M']= get_concentration(df['n_kaps'], df['box_side_A'])
    df['fg_C_M']= get_concentration(df['n_fgs'], df['box_side_A'])
    df['k_on_per_ns_per_M']= df['k_on_per_ns_per_missing_ss_contact'] \
        * AVOGADRO * df['box_side_A']**3 * L_per_A3
    df['k_on_per_ns_per_M_lbound']= df['k_on_per_ns_per_missing_ss_contact_lbound'] \
        * AVOGADRO * df['box_side_A']**3 * L_per_A3
    df['k_on_per_ns_per_M_ubound']= df['k_on_per_ns_per_missing_ss_contact_ubound'] \
        * AVOGADRO * df['box_side_A']**3 * L_per_A3
    if os.path.exists(filename):
        open_mode='a'
        is_header= False
    else:
        open_mode='w'
        is_header= True
    with open(filename, open_mode) as f:
        print("Appending to {} with header={} open_mode {}".format(filename, is_header, open_mode))
        df.to_csv(path_or_buf=f,
                  sep=' ',
                  mode= open_mode,
                  header= is_header,
                  encoding= 'utf-8')
#    print(df.head().to_string())

def get_processed_files(stats_filename):
    if os.path.exists(stats_filename):
        with open(stats_filename, 'r') as f:
            data= pd.read_csv(stats_filename, sep=' ')
        return list(data['output_file'])
    else:
        return []

def check_output_file(output_file, verbose):
    '''
    Check if file shoudl be processed (not processed beore, exists, a file,
    and not empty).
    '''
    if output_file in processed_files:
        if verbose:
            print("Skipping " + output_file + " because it was processed")
        return False
    if not os.path.isfile(output_file):
        if verbose:
            print("Warning: skipping {} - not exists or not a regular file" \
                  .format(output_file))
        return False
    if(os.stat(output_file).st_size == 0):
        if verbose:
            print("Warning: skipping {} - empty file" \
                  .format(output_file))
        return False
    return True


if __name__ != "__main__":
    sys.exit(-1)
if len(sys.argv)>1:
    stats_csv_file=sys.argv[1]
else:
    stats_csv_file='stats.csv'
if len(sys.argv)>2:
    stats_min_time_sec= float(sys.argv[2])
else:
    stats_min_time_sec= 0.1e-6
if len(sys.argv)>3:
    ''' if true, always uses new version of kap/tf site statistics '''
    is_new_tf_stats= bool(sys.argv[3])
else:
    is_new_tf_stats= False
processed_files= get_processed_files(stats_csv_file)
#if os.path.exists(stats_csv_file):
#    os.remove(stats_csv_file)
data= []
MAX_DATA_ENTRIES= 500
for fg_seq in ['F', 'FSSSSS', 'FFSSSS', 'FFFSSS', 'FFFFSS', 'FFFFFF']:
    output_files=glob.glob('Output/{}/*.pb'.format(fg_seq))
    print("Processing {:d} output files".format(len(output_files)))
    for output_file in output_files:
        try:
            if not check_output_file(output_file, verbose= True):
                continue
            stats_entry= get_stats_entry_for_output_file(output_file, 'fg0', 'kap20',
                                                         stats_min_time_sec,
                                                         is_new_tf_stats)
            data.append(stats_entry)
            processed_files.append(output_file)
        except (KeyboardInterrupt, SystemExit):
            raise
        except RuntimeError as e:
            print("Skipped {} due to a Runtime Error exception - consider uncommenting 'raise' from analyze.py script".format(output_file))
            print(e)
#            raise
        print("Processed {}".format(output_file))
        if len(data)>MAX_DATA_ENTRIES:
            append_data_to_csv_file(data, stats_csv_file)
            data=[]
append_data_to_csv_file(data, stats_csv_file)

print("Finished succesfully")
