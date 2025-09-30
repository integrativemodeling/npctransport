#!/usr/bin/python

import IMP.algebra
from IMP.em import *
import glob
import re
import numpy as np
import scipy.ndimage.interpolation
import multiprocessing
import pandas as pd
import sys
M=110

def load_nup_from_file(nup_filename):
    '''
    Load nup file from space delimited txt file nup_filename
    (from files produced by get_xyz_stats_v*.py)
    Apply 8-fold rotation to symmetrize
    '''
    global M
    a = pd.read_csv(file_path, sep=' ', header=None).to_numpy()
    print(a.shape)
    assert(len(a.shape)==3 and M==len(a[0])
           and len(a[0])==len(a[1]) and len(a[1])==len(a[2]))
    # Apply 8-fold symmetry
    ret= a*0
    for i in range(8):
        ret= ret+scipy.ndimage.interpolation.rotate(a,
                                                    angle=i*45,
                                                    reshape=False)
    return ret

def add_nup_to_dmaps(dmaps, nup_name, weight=1.0):
    '''
    Add all files with nup type beginning with nup_name
    to dmap with specified weight (i.e. S*.<nup_name>*.txt)
    '''
    global M
    nup_filenames=glob.glob("S*.{}*.txt".format(nup_name))
    for nup_filename in nup_filenames:
        if re.search('anchor', nup_filename):
            print("Skipping anchor {}".format(nup_filename))
            continue
        if re.search(nup_name+"[0-9]", nup_filename):
            continue # wrong nup if prefix + digit
        print("Processing {}".format(nup_filename))
        nup_dmap= DensityMap(dh)
        a=load_nup_from_file(nup_filename)
        nup_dmap.set_void_map(a.shape[0], a.shape[1], a.shape[2])
        for xi in range(M):
            for yi in range(M):
                for zi in range(M):
                    for dmap in dmaps:
                        index=dmap.xyz_ind2voxel(xi,yi,zi)
                        cur_val=dmap.get_value(index)
                        nup_dmap.set_value(index, weight*a[xi][yi][zi])
                        dmap.set_value(index, cur_val + weight*a[xi][yi][zi])
        write_map(nup_dmap,
                  nup_filename.replace("txt", "mrc"),
                  MRCReaderWriter())

def get_dh(matrix_3D_shape, voxel_size_A=10.0):
    box_dimensions = matrix_3D_shape * voxel_size_A
    bb= IMP.algebra.BoundingBox3D( IMP.algebra.Vector3D(-box_dimensions/2),
                                   IMP.algebra.Vector3D(+box_dimensions/2) );
    dh=create_density_header(bb, 10.0)
    print("Origin", dh.get_xorigin())
    print("Resolution", dh.get_resolution())
    return dh

def write_diffuser(name_prefix,
                     R_MIN= 14,
                     R_MAX= 28,
                   WINDOW= 8):
    '''
    Convert diffusers of different radii
    with name prefix name_prefix, e.g. kap14..kap28
    '''
    global M
    global VOXEL_SIZE_A
    dh = get_dh(np.array([M,M,M]), VOXEL_SIZE_A)
    dmap= DensityMap(dh)
    dmap.set_void_map(M, M, M)
    R_BASE= R_MIN
    dmap_cur= DensityMap(dh)
    dmap_cur.set_void_map(M, M, M)
    for R in range(R_MIN,R_MAX+1,2):
        diffuser_name= "{}{:d}".format(name_prefix, R)
        add_nup_to_dmaps([dmap, dmap_cur], diffuser_name, weight=R**3)
        if(R-R_BASE == WINDOW-2):
            cur_file= "{}s{:d}-{:d}.mrc".format(name_prefix,
                                                R_BASE,
                                                R)
            write_map(dmap_cur,
                      cur_file,
                      MRCReaderWriter())
            print("Wrote {}".format(cur_file))
            R_BASE= R+2
        dmap_cur= DensityMap(dh)
        dmap_cur.set_void_map(M, M, M)
    write_map(dmap, "{}s.mrc".format(name_prefix) ,
    MRCReaderWriter())
    return dmap

def write_fg(nup, dmap):
    global M
    dh = get_dh(np.array([M,M,M]), VOXEL_SIZE_A)
    dmap_cur= DensityMap(dh)
    dmap_cur.set_void_map(M, M, M)
    add_nup_to_dmaps([dmap, dmap_cur], nup)
    write_map(dmap_cur, "%s.mrc" % nup, MRCReaderWriter())

def write_fgs(fg_nup_list):
    '''
    Convert FGs to MRC
    '''

    # FG repeats
    global M
    dh = get_dh(np.array([M,M,M]), VOXEL_SIZE_A)
    dmap=DensityMap(dh)
    dmap.set_void_map(M, M, M)
    for nup in fg_nup_list:
        write_fg(nup, dmap)
    write_map(dmap, "FGs.mrc", MRCReaderWriter())

def write_obstacles():
    # Obstacles
    global M
    dh = get_dh(np.array([M,M,M]), VOXEL_SIZE_A)
    dmap=DensityMap(dh)
    dmap.set_void_map(110,110,110)
    add_nup_to_dmaps([dmap], "obstaclesInflated")
    write_map(dmap, "obstacles.mrc", MRCReaderWriter())


# kaps and inerts
M = int(sys.argv[1]) if len(sys.argv)>1 else 110 # expected box dimensions in voxels
VOXEL_SIZE_A = float(sys.argv[2]) if len(sys.argv)>2 else 10.0
write_diffuser("kap", 14, 34, WINDOW=10)
write_diffuser("inert", 14, 34, WINDOW=10)
fg_nup_list= ["Nsp1","Nup100","Nup116","Nup42","Nup159","Nup49","Nup57","Nup1","Nup60","Nup145", "fg0"]
write_fgs(fg_nup_list)
#write_obstacles()
