from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from numpy.testing import assert_approx_equal
import scipy.ndimage.filters
import math
import my_util
import sys
#import seaborn as sns

def get_zr_from_files(fnames):
    '''
    Load txt files of xyz voxel densities, sum them, and convert them to a matrix of rxz grid.
    The input files are assumed to be dxd rows of d columns, readable by np.loadtxt, and they are reshaped
    to a dxdxd matrix of xyz densities. It is assumed that the origin of each axis is halfway between
    the d/2 and d/2+1 coordinates. The units of the axes are arbitraty.

    Returns a dxd matrix, representing densities of as a function of R and Z. The densities are normalized
    to reflect densities in the original x,y,z space - in other words, they
    are divided by R to prevent overcounting. The first coordinate is R=0.0-0.5 in the same units of x,y,z, the second
    coordinate R=0.5-1.5, etc.
    '''
    assert(len(fnames)>=1)
    xyz= np.loadtxt(fnames[0], dtype=int)
    for fname in fnames[1:]:
        xyz= xyz+np.loadtxt(fname, dtype=int)
    dold= xyz.shape
    d= int(math.sqrt(dold[0]))
    assert_approx_equal(d, math.sqrt(dold[0]))
    assert_approx_equal(d, dold[1])
    xyz=np.reshape(xyz, [d,d,d])
    print(xyz.shape)
    R= 0.5*d
    zr= np.zeros([d,d])
    for xi in range(d):
        x= xi-R+0.5 # +0.5 to reflect the middle of each grid point
        for yi in range(d):
            y= yi-R+0.5
            r= math.sqrt(x**2+y**2)
            assert(r>0)
            ri= int(round(r))
            zr[:,ri]=zr[:,ri]+xyz[xi,yi,:]*(R/r) # normalize by perimeter length (which disfavors small r otherwise)
    return zr

def print_up_to_maximal_R(zr, max_R_nm):
    dim_z= zr.shape[0]
    mid_zi = int(round(0.5*dim_z))
    print("Sum over center Rs up to {} nm".format(max_R_nm))
    count_by_z= np.sum(zr[:,0:max_R_nm],axis=1)
    print("From z={:.0f} to z={:.0f} nm".format(0, -R))
    for zi in my_util.range_inclusive(mid_zi-1, 0, step=-1):
        print("{:04.0f}".format(count_by_z[zi]), end=" ")
        if( (mid_zi - zi)  % 10 == 0):
            print("    ", end="")
    print()
    print("From z={:.0f} to z={:.0f} nm".format(0, R))
    for zi in range(mid_zi, dim_z):
        print("{:04.0f}".format(count_by_z[zi]), end=" ")
        if((zi - mid_zi)  % 10 == 9):
            print("    ", end="")
    print()
    print("--")


zr= get_zr_from_files(sys.argv[1:])
print(zr.shape)
R=zr.shape[1]/2+0.5
for zi in range(zr.shape[1]):
    z= zi-R
#    print("z=",z,"nm")
#    print(b[:,zi])
axes=["Z","R"]
for iaxis, axis in enumerate(axes):
    print("Sum over {}:".format(axis))
    for x in np.sum(zr, axis=iaxis):
        print("{:.1f}".format(x), end=" ")
    print()
for max_R_nm in my_util.range_inclusive(5,55,5):
    print_up_to_maximal_R(zr, max_R_nm)

plt.figure()
zr= scipy.ndimage.filters.gaussian_filter(zr, 5)
xticklabels=[]
for i in range(zr.shape[0]):
    if i % 10 == 0:
        xticklabels.append(str(i))
    else:
        xticklabels.append('')
yticklabels=[]
for i in range(zr.shape[1]):
    if int(i-0.5*zr.shape[1]) % 10 == 0:
        yticklabels.append(str(i-R+1))
    else:
        yticklabels.append('')
ax= sns.heatmap(zr,xticklabels=xticklabels, yticklabels=yticklabels)
ax.set_xlim([0,100.0])
ax.invert_yaxis()
plt.show()
