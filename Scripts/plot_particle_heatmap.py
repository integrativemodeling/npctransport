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

def heatmap_to_xy(heatmap):
    IS_LOG= False
    if IS_LOG:
        H=np.log10(heatmap*10.0+0.1) - np.log10(0.1)
    else:
        max_ratio= 5.0
        min_element= np.max(np.max(heatmap)) / max_ratio
        assert(min_element > 0.0)
        H= heatmap/min_element
    x=[]
    y=[]
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            n= int(round(H[i,j]))
            x.extend([float(i)]*n)
            y.extend([float(j)]*n)
    return x,y

def plot_zr_contourf(zr, ax):
    R=zr.shape[1]/2+0.5
#    zr= scipy.ndimage.filters.gaussian_filter(zr, 3)
    ZZ,RR= np.meshgrid(range(zr.shape[0]), range(zr.shape[1]))
    NN= zr[ZZ, RR]
    ax.contourf(RR, ZZ, NN, cmap=plt.cm.get_cmap('Blues', 9))
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
#    ax.invert_yaxis()
    xticks=range(0,111,10)
    yticks= [s+0.5 for s in range(int(R)-100,int(R)+101,10)]
    yticklabels= [ytick-R for ytick in yticks]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(["{:d}".format(xx) for xx in xticks])
    ax.set_yticklabels(yticklabels)
    ax.set_xlim([0,55.0])
    ax.set_ylim([R-50.0,R+50.0])

def plot_zr_heatmap(zr):
    R=zr.shape[1]/2+0.5
    xticks=range(0,111,10)
    yticks= [s+0.5 for s in range(int(R)-100,int(R)+101,10)]
    yticklabels= [ytick-R for ytick in yticks]
    ax = sns.heatmap(zr, cmap=plt.cm.get_cmap('Blues', 9),
                     robust=False, square=True)
    print(xticks)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(["{:d}".format(xx) for xx in xticks])
    ax.set_yticklabels(yticklabels)
    ax.set_xlim([0,55.0])
    ax.set_ylim([R-50.0,R+50.0])
    return ax

def do_plots(zr):
    fig, ax= plt.subplots(1, 2, sharey=True)
    plot_zr_contourf(zr, ax[0])
    ax[0].set_aspect(1.0, adjustable='box-forced')
    zr= scipy.ndimage.filters.gaussian_filter(zr, 3)
    plt.axes(ax[1])
    assert(ax[1] == plot_zr_heatmap(np.flipud(zr)))
#    plt.axes().set_aspect('equal', 'box')
#    ax[0].invert_yaxis()
#    ax[1].invert_yaxis()


# Prepare:
zr= get_zr_from_files(sys.argv[1:])
# Print stuff:
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
# Plot:
zr= scipy.ndimage.filters.gaussian_filter(zr, 3)
do_plots(zr)
Lzr= np.log(zr+0.1) - np.log(0.1)
do_plots(Lzr)
plt.show()
