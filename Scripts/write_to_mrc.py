#!/usr/bin/python

import IMP.algebra
from IMP.em import *
import scipy.ndimage.interpolation
import numpy as np
import glob
M=110

def load_nup_from_file(nup_filename):
    '''
    Load nup file from space delimited txt file nup_filename
    (from files produced by get_xyz_stats_v*.py)
    Apply 8-fold rotation to symmetrize
    '''
    global M
    a=np.zeros([M,M,M])
    f= open(nup_filename, 'r')
    xi= 0
    yi= 0
    zi= 0
    for line in f:
        fields=[float(x) for x in line.split()]
        for zi,value in enumerate(fields):
            a[xi][yi][zi]= value
        N=len(fields)
        assert(N==110)
        xi=(xi+1) % N
        if(xi==0):
            yi=yi+1
    # Apply 8-fold symmetry
    ret= a*0
    for i in range(8):
        ret= ret+scipy.ndimage.interpolation.rotate(a,
                                                    angle=i*45,
                                                    reshape=False)
    return ret

def add_nup_to_dmap(dmap, nup_name, weight=1.0):
    '''
    Add all files with nup type beginning with nup_name
    to dmap with specified weight (i.e. S*.<nup_name>*.txt)
    '''
    global M
    nup_filenames=glob.glob("S*.{}*.txt".format(nup_name))
    for nup_filename in nup_filenames:
        a=load_nup_from_file(nup_filename)
        for xi in range(M):
            for yi in range(M):
                for zi in range(M):
                    index=dmap.xyz_ind2voxel(xi,yi,zi)
                    cur_val=dmap.get_value(index)
                    dmap.set_value(index, cur_val + weight*a[xi][yi][zi])

bb= IMP.algebra.BoundingBox3D( IMP.algebra.Vector3D([-550,-550,-550]),
                               IMP.algebra.Vector3D([+550,+550,+550]) );
dh=create_density_header(bb, 10.0)
print "Origin", dh.get_xorigin()
print "Resolution", dh.get_resolution()

def write_diffuser(name_prefix,
                     R_MIN=14,
                     R_MAX=28):
    '''
    Convert diffusers of different radii
    with name prefix name_prefix, e.g. kap14..kap28
    '''
    dmap=DensityMap(dh)
    dmap.set_void_map(110,110,110)
    WINDOW=8
    R_BASE=R_MIN
    dmap_cur=DensityMap(dh)
    dmap_cur.set_void_map(110,110,110)
    for R in range(R_MIN,R_MAX+1,2):
        diffuser_name= "{}{:d}".format(name_prefix, R)
        add_nup_to_dmap(dmap, diffuser_name, weight=R**3)
        add_nup_to_dmap(dmap_cur, diffuser_name, weight=R**3)
        if(R-R_BASE == WINDOW-2):
            write_map(dmap_cur,
                      "{}s{:d}-{:d}.mrc".format(diffuser_name,
                                                R_BASE,
                                                R),
                      MRCReaderWriter())
        R_BASE=R+2
        dmap_cur=DensityMap(dh)
        dmap_cur.set_void_map(110,110,110)
    write_map(dmap, "{}s.mrc".format(diffuser_name) ,
    MRCReaderWriter())
    return dmap

def write_fgs(fg_nup_list):
    '''
    Convert FGs to MRC
    '''

    # FG repeats
    dmap=DensityMap(dh)
    dmap.set_void_map(110,110,110)
    for nup in fg_nup_list:
        add_nup_to_dmap(dmap, nup)
        dmap_cur=DensityMap(dh)
        dmap_cur.set_void_map(110,110,110)
        add_nup_to_dmap(dmap_cur, nup)
        write_map(dmap_cur, "%s.mrc" % nup, MRCReaderWriter())
    write_map(dmap, "FGs.mrc", MRCReaderWriter())

def write_obstacles():
    # Obstacles
    dmap=DensityMap(dh)
    dmap.set_void_map(110,110,110)
    add_nup_to_dmap(dmap, "obstaclesInflated")
    write_map(dmap, "obstacles.mrc", MRCReaderWriter())


# kaps and inerts
write_diffuser("kap")
write_diffuser("inert")
fg_nup_list= ["Nsp1","Nup100","Nup116","Nup42","Nup159","Nup49","Nup57","Nup1","Nup60","Nup145"]
write_fgs(fg_nup_list)
write_obstacles()
