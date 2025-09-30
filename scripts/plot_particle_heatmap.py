from __future__ import print_function
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from numpy.testing import assert_approx_equal
import scipy.ndimage
import scipy.ndimage.filters
import math
import my_util
import sys
import os.path
#import seaborn as sns
import argparse
import IMP.em

def parse_commandline():
    parser = argparse.ArgumentParser \
             ( description='plot particle heatmap based on npctransport simulation,' +
               ' taking as input 3D spatial distribution matrices in txt form' +
               ' produced by e.g. get_xyz_from_files.py',
               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_files', metavar='filename', type=str, nargs='+',
                        help='input files with 3D spatial distribution matrices')
    parser.add_argument('--z_levels', metavar='z_level', type=int, nargs='*',
                        default=[0],
                        help='z levels for which to show distributions, in nm')
#    parser.add_argument('--split', default=24, type=int, dest='split')
    parser.add_argument('--show_zr', action='store_true', help='show plot in ZR space')
    parser.add_argument('--write_to_mrc', metavar='mrc_filename', type=str, default='',
                        help='optional mrc file to which XYZ density will be saved')
    args = parser.parse_args()
    return args


def get_xyz_from_files(fnames):
    '''
    Load txt files of xyz voxel densities and sum them
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
    return xyz

def apply_8fold_sym(xyz):
    ''' apply 8-fold symmetry along z-axis to an x,y,z matrix '''
    print("Original xyz max {0:.3f} min {1:.3f}".format(xyz.max(),
                                                        xyz.min()))
#    xyz= scipy.ndimage.shift(xyz,-0.5)
#    print("Shifted+positived xyz max {0:.3f} min {1:.3f}".format(xyz.max(),
#                                                       xyz.min()))
#    print("Shifted:", xyz.shape)
    d= xyz.shape[0]
    R= d*0.5
    ret_value= 0.0*xyz
    for i in range(8):
        rot_i= scipy.ndimage.rotate(xyz,angle=i*45.0)
#        print(rot_i.shape)
        R_i= rot_i.shape[0]*0.5
        offset= int(R_i-R)
#        print(R_i, R, offset)
        rot_i= rot_i[offset:d+offset,offset:d+offset,:]
        print("Rotated voxels max {0:.3f} min {1:.3f} angle {2:.1f}".format(rot_i.max(),
                                                                            rot_i.min(),
                                                                            i*45.0))
        ret_value= ret_value+rot_i
        print("Symmetrized voxels max {0:.3f} min {1:.3f}".format(ret_value.max(),
                                                                  ret_value.min()))
    ret_value[ret_value<0.0]=0.0
    print("Positived Symmetrized voxels max {0:.3f} min {1:.3f}".format(ret_value.max(),
                                                                    ret_value.min()))
    return ret_value

def get_zr_from_xyz(xyz):
    '''
    Convert an x x y x z 3D matrix of voxel densities to a 2-d matrix of r x z pixels
    It is assumed that the origin of each axis is halfway between
    the d/2 and d/2+1 coordinates. The units of the axes are arbitraty.

    Returns a dxd matrix, representing densities of as a function of R and Z. The densities are normalized
    to reflect densities in the original x,y,z space - in other words, they
    are divided by R to prevent overcounting. The first coordinate is R=0.0-0.5 in the same units of x,y,z, the second
    coordinate R=0.5-1.5, etc.
    '''
    d= xyz.shape[0]
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
    xyz= get_xyz_from_files(fnames)
    return get_zr_from_xyz(xyz)

def write_xyz_to_mrc(xyz, mrc_filename):
    VOXEL_SIZE_A= 10.0
    M= xyz.shape[0]
    assert(M==xyz.shape[1])
    assert(M==xyz.shape[2])
    M_A= M*VOXEL_SIZE_A
    bb= IMP.algebra.BoundingBox3D( IMP.algebra.Vector3D([-M_A/2, -M_A/2, -M_A/2]),
                                   IMP.algebra.Vector3D([+M_A/2, +M_A/2, +M_A/2]) )
    dh= IMP.em.create_density_header(bb, int(VOXEL_SIZE_A))
    print(dh)
    print( "MRC origin {0:}".format(dh.get_xorigin()) )
    print( "MRC resolution {0:.1f} A".format(dh.get_resolution()) )
    print( "MRC spacing {0:.1f} A".format(dh.get_spacing()) )
    dmap= IMP.em.DensityMap(dh)
    dmap.set_void_map(M,M,M)
    for xi in range(M):
        for yi in range(M):
            for zi in range(M):
                index=dmap.xyz_ind2voxel(xi,yi,zi)
                dmap.set_value(index, xyz[xi][yi][zi])
    IMP.em.write_map(dmap, mrc_filename, IMP.em.MRCReaderWriter())
    print("Wrote {0}, hopefully".format(mrc_filename))


def print_up_to_maximal_R(zr, max_R_nm):
    dim_z= zr.shape[0]
    mid_zi = int(round(0.5*dim_z))
    R= mid_zi
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
    ax.contourf(RR, ZZ, NN, cmap=plt.cm.get_cmap('Blues', 10))
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

def plot_zr_heatmap(zr, zr_obstacles, ax):
    R=zr.shape[1]/2+0.5
    xticks=range(0,111,10)
    yticks= [s+0.5 for s in range(int(R)-100,int(R)+101,10)]
    yticklabels= [ytick-R for ytick in yticks]
    sns.heatmap(zr, cmap=plt.cm.get_cmap('Blues', 10),
                robust=False, square=True,
                mask=(zr_obstacles>=0.1),
                ax=ax)
    sns.heatmap(zr_obstacles, cmap=plt.cm.get_cmap('Reds', 10),
                robust=False, square=True,
                mask=(zr_obstacles<0.1),
                ax=ax )
    print(xticks)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(["{:d}".format(xx) for xx in xticks])
    ax.set_yticklabels(yticklabels)
    ax.set_xlim([0,55.0])
    ax.set_ylim([R-50.0,R+50.0])
    return ax

def do_plots(zr, zr_obstacles, title):
    fig, ax= plt.subplots(1, 2, sharey=True)
    plot_zr_contourf(zr, ax[0])
    ax[0].set_aspect(1.0, adjustable='box-forced')
    zr= scipy.ndimage.filters.gaussian_filter(zr, 3)
    plt.axes(ax[1])
    assert(ax[1] == plot_zr_heatmap( np.flipud(zr),
                                     np.flipud(zr_obstacles),
                                     ax[1] ))
    for axis in ax:
        axis.set_title(title)
#    plt.axes().set_aspect('equal', 'box')
#    ax[0].invert_yaxis()
#    ax[1].invert_yaxis()
    return ax

def analyze_zr(zr):
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


# *************************
# ********* MAIN: *********
# *************************
args=parse_commandline()
# Prepare:
xyz= get_xyz_from_files(args.input_files)
xyz= apply_8fold_sym(xyz)
xyz_obstacles= get_xyz_from_files(['/salilab/park4/barak/Tmp/S0.obstacles.txt'])
#analyze_zr(zr)
title= os.path.splitext(os.path.basename(sys.argv[1]))[0]
if(len(args.input_files)>1):
    title= "{0:} and {1:d} more".format(title, len(args.input_files)-1)

if(len(args.write_to_mrc)>0):
    write_xyz_to_mrc(xyz, args.write_to_mrc)

# Plot zr:
if(args.show_zr):
    zr=get_zr_from_xyz(xyz)
    zr_obstacles= get_zr_from_xyz(xyz_obstacles)
    #zr= scipy.ndimage.filters.gaussian_filter(zr, 3)
    do_plots(zr, zr_obstacles, title)
    Lzr= np.log(zr+0.1) - np.log(0.1)
    do_plots(Lzr, zr_obstacles, title)
    plt.show()
    sys.exit()

# Plot xy
#xyz= scipy.ndimage.filters.gaussian_filter(xyz, 0.025)
xyz_obstacles= apply_8fold_sym(xyz_obstacles)
xyz_obstacles= scipy.ndimage.filters.gaussian_filter(xyz_obstacles, 0.15)
xyz_obstacles[xyz_obstacles<0.0]= 0.0
for z_nm in args.z_levels:
    d= xyz.shape[0]
    zi= d/2+z_nm
    assert(zi==int(zi) and zi>=0 and zi<d)
    xy0= xyz[:,:,d/2+z_nm]
    xy0_obstacles= xyz_obstacles[:,:,d/2+z_nm]
    #do_plots(xy0, title)
    Lxy0= np.log(xy0+0.1) - np.log(0.1)
    plt.figure()
    R=xyz.shape[0]/2+0.5
    ticks= [s+0.5 for s in range(int(R)-100,int(R)+101,10)]
    ticklabels= [tick-R for tick in ticks]
    ax=sns.heatmap(xy0_obstacles, cmap=plt.cm.get_cmap('Reds', 10),
                   robust=False, square=True,
                   mask=(xy0_obstacles==0) #,ax=ax
    )
    ax = sns.heatmap(xy0, cmap=plt.cm.get_cmap('Blues', 10),
                     robust=False, square=True,
                     mask=(xy0_obstacles>0.1), ax=ax)
    print("xy0 ob mins {0:.2f} max {1:.2f}".format(xy0_obstacles.min(),
                                                   xy0_obstacles.max()))
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticklabels)
    ax.set_xlim([R-50.0,R+50.0])
    ax.set_ylim([R-50.0,R+50.0])
    ztitle = (title + " z={0:d} nm").format(z_nm)
    ax.set_title(ztitle)
plt.show()
