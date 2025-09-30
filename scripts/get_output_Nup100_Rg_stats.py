import IMP
import IMP.algebra
import IMP.container
import IMP.core
from IMP.npctransport import *
import glob
import math
import numpy as np
import pandas as pd
import tempfile
import time


def get_Nup100_GLFG_distances(init_output_fname):
    ''' 
    Compute distances between Nup100 GLFGs in positions 1-29, 1-15 15-29,  12-24
    corresponding to residues 10-570, 10-290, 290-570, 230-470
    
    :param init_output_fname: path to output file from which to compute statistics
                              note the file will be rewritten so better work on a copy
    :return: .
        '''
    # temporarily look the output file for writing just in case (using os)
    sd = IMP.npctransport.SimulationData(init_output_fname, False)
    m = sd.get_model()
    chains = sd.get_fg_chain_roots()
    df = pd.DataFrame(columns=['d1_29', 'd1_15', 'd15_29', 'd12_24'])
    for chain in chains:
        type = repr(chain)[1:-2]
        if type=='Nup100':
                children = chain.get_children()
                xyz1 = IMP.core.XYZ(children[-1])
                xyz15 = IMP.core.XYZ(children[-15])
                xyz29 = IMP.core.XYZ(children[-29])
                xyz12 = IMP.core.XYZ(children[-12])
                xyz24 = IMP.core.XYZ(children[-24])
                d1_29 = IMP.core.get_distance(xyz1, xyz29)
                d1_15 = IMP.core.get_distance(xyz1, xyz15)
                d15_29 = IMP.core.get_distance(xyz15, xyz29)
                d12_24 = IMP.core.get_distance(xyz12, xyz24)
                df.loc[len(df)] = [d1_29, d1_15, d15_29, d12_24]
    return df

def call_get_Nup100_GLFG_distances_safely(init_output_fname):
    '''
    Call get_Nup100_GLFG_distances() safely, by copying the input
    file to a temporary copy, work on the copy, and delete it when done
    '''
    try:
        with tempfile.NamedTemporaryFile(suffix=".pb") as f:
            f.write(open(init_output_fname, "rb").read())
            f.flush()
            df = get_Nup100_GLFG_distances(f.name)
        return df
    except:
        print("Skipping file due to error: " + init_output_fname)
        return None

def print_df_stats(df):
    print(f"number of entries: {len(df.index)}")
    print("Mean:")
    print(df.mean())
    print("Standard-deviation:")
    print(df.std())
    print("Std-error of the mean")
    print(df.sem())

if __name__ == '__main__':
    fnames= glob.glob(sys.argv[1]+"*.pb")
    dfs = []
    for fname in fnames:
        print("Handling " + fname)
        df = call_get_Nup100_GLFG_distances_safely(fname)
        dfs.append(df)
        if len(dfs)%20==0:
            df = pd.concat(dfs)
            print_df_stats(df)
            dfs = [df]
    df = pd.concat(dfs)
    print_df_stats(df)
    del dfs