#!/usr/bin/python

import IMP.algebra
from IMP.em import *
import scipy.ndimage.interpolation
import numpy as np
import glob

M=110

def load_nup(nup_name, folders_prefix="../STATS_StepN"):
    global M
    a=np.zeros([M,M,M])
    folders= glob.glob(folders_prefix+"*/")
    n_opened= 0
    for folder in folders:
        try:
            f= open(folder+"/S0.%s.txt" % nup_name, 'r')
            n_opened= n_opened+1
        except:
            print("Couldn't open %s in folder %s" % (nup_name, folder))
            continue
        xi= 0
        yi= 0
        zi= 0
        for line in f:
            fields=[float(x) for x in line.split()]
            for zi,value in enumerate(fields):
                a[xi][yi][zi]= a[xi][yi][zi]+value
            N=len(fields)
            assert(N==110)
            xi=(xi+1) % N
            if(xi==0):
                yi=yi+1
    assert(n_opened>0)
    # Apply 8-fold symmetry
    ret= a*0
    for i in range(8):
        ret= ret+scipy.ndimage.interpolation.rotate(a,
                                                    angle=i*45,
                                                    reshape=False)
    return ret

def add_nup_to_dmap(dmap, nup_name, weight=1.0):
    global M
    a=load_nup(nup_name)
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

# kaps
dmap=DensityMap(dh)
dmap.set_void_map(110,110,110)
WINDOW=8
R_MIN=14
R_MAX=28
R_BASE=R_MIN
dmap_cur=DensityMap(dh)
dmap_cur.set_void_map(110,110,110)
for R in range(R_MIN,R_MAX+1,2):
    kap= "kap%d" % R
    add_nup_to_dmap(dmap, kap, weight=R**3)
    add_nup_to_dmap(dmap_cur, kap, weight=R**3)
    if(R-R_BASE == WINDOW-2):
        write_map(dmap_cur, "Kaps%d-%d.mrc" % (R_BASE,R), MRCReaderWriter())
        R_BASE=R+2
        dmap_cur=DensityMap(dh)
        dmap_cur.set_void_map(110,110,110)
write_map(dmap, "Kaps.mrc" , MRCReaderWriter())

# FG repeats
dmap=DensityMap(dh)
dmap.set_void_map(110,110,110)
for nup in ["Nsp1","Nup100","Nup116","Nup42","Nup159","Nup49","Nup57","Nup1","Nup60","Nup145"]:
    add_nup_to_dmap(dmap, nup)
    dmap_cur=DensityMap(dh)
    dmap_cur.set_void_map(110,110,110)
    add_nup_to_dmap(dmap_cur, nup)
    write_map(dmap_cur, "%s.mrc" % nup, MRCReaderWriter())
write_map(dmap, "ALL.mrc", MRCReaderWriter())

# Obstacles
dmap=DensityMap(dh)
dmap.set_void_map(110,110,110)
add_nup_to_dmap(dmap, "obstaclesInflated")
write_map(dmap, "obstacles.mrc", MRCReaderWriter())
