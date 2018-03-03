#!/usr/bin/bash /Users/barak/imp/fast/setup_environment.sh /anaconda/bin/python

#########################################
# Example: imppy Scripts/write_obstacles_voxel_map.py InputData/wholeNPC_0.rmf3  > obstacles.txt
#
# Author: Barak Raveh barak@salilab.org
# Data written: Oct 27, 2016
# Date last udpated: Oct 27, 2016
#########################################
from IMP.npctransport import *
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.algebra
import RMF
import math
import re
import sys
import numpy as np

input_rmf= sys.argv[1]
IS_REMOVE_GLE1_AND_NUP42=True
VOXEL_SIZE_A= 10
VOXELS_N= 110
VOXELS_LENGTH_A= VOXEL_SIZE_A * VOXELS_N
VOXELS_MIN_A= -0.5 * VOXELS_LENGTH_A
VOXELS_MAX_A= VOXELS_MIN_A + VOXELS_LENGTH_A

fgs_anchors_regexp= {}
fgs_anchors_regexp['Nsp1.*_601-636_']= 'Nsp1'
fgs_anchors_regexp['Nup1.*_301-350_']= 'Nup1'
if True or not IS_REMOVE_GLE1_AND_NUP42:
    fgs_anchors_regexp['Nup42.*_364-413_'] ='Nup42'
fgs_anchors_regexp['Nup49.*_201-269_'] ='Nup49'
fgs_anchors_regexp['Nup57.*_201-286_'] ='Nup57'
fgs_anchors_regexp['Nup60.*_351-398_'] ='Nup60'
fgs_anchors_regexp['Nup100.*_551-575_'] ='Nup100'
fgs_anchors_regexp['Nup116.*_751-775_'] ='Nup116'
fgs_anchors_regexp['Nup145.*_201-225_'] ='Nup145'
fgs_anchors_regexp['Nup159.*_1082-1116_'] ='Nup159'

Z_TRANSFORM=0 #-75.5
TUNNEL_RADIUS=350
NUCLEAR_ENVELOPE_WIDTH=300
OBSTACLE_SCALE_FACTOR=1.0 # inflate obstacles a bit to prevent artificial cavities
obstacles={}

def do_stats(XYZ):
  for YZ in XYZ:
    for Z in YZ:
      for n_xyz in Z:
        print "%2.1f " % n_xyz,
      print ""
            # print separator line of -1s
            #            if Z is not None:
            #                for n_xyz in Z:
            #                    print >>F, -1,
            #                print >>F, ""

def add_obstacles(XYZ, coords_list, radius_A):
    '''
    XYZ- 3D voxels map
    coords_list - list of 3D coordinate vectors
    radius - radius of all obstacles of this type
    '''
    radius_NM= 0.1*radius_A
    volume_NM3= (4.0/3.0)*math.pi*radius_NM**3
    for coords in coords_list:
        x=coords[0]
        y=coords[1]
        z=coords[2] + Z_TRANSFORM
        if not (VOXELS_MIN_A < x < VOXELS_MAX_A and
                VOXELS_MIN_A < y < VOXELS_MAX_A and
                VOXELS_MIN_A < z < VOXELS_MAX_A):
            continue
        nx= math.floor( (x - VOXELS_MIN_A) / VOXEL_SIZE_A)
        ny= math.floor( (y - VOXELS_MIN_A) / VOXEL_SIZE_A)
        nz= math.floor( (z - VOXELS_MIN_A) / VOXEL_SIZE_A)
        XYZ[nx][ny][nz]= XYZ[nx][ny][nz]+ 1 #volume_NM3
        if radius_A  >VOXEL_SIZE_A:
            add_obstacles(XYZ,
                [ [x + VOXEL_SIZE_A, y, z],
                    [x - VOXEL_SIZE_A, y, z],
                    [x, y + VOXEL_SIZE_A, z],
                    [x, y - VOXEL_SIZE_A, z],
                    [x, y, z + VOXEL_SIZE_A],
                    [x, y, z - VOXEL_SIZE_A] ],
                radius_A - VOXEL_SIZE_A)
    return obstacles



def handle_xyz_children(parent):
    for p in parent.get_children():
        xyzr=IMP.core.XYZR(p)
        coords=xyzr.get_coordinates()
        radius=round(xyzr.get_radius(),1)
        for i in range(8):
            R=IMP.algebra.get_rotation_about_normalized_axis([0,0,1],
                                                             i*math.pi/4.0)
#            print parent.get_parent().get_name(), parent.get_name(), p.get_name(),
            coords_i=R*coords
            distance=math.sqrt(coords_i[0]**2+coords_i[1]**2)
            if(distance-radius>TUNNEL_RADIUS and
                    abs(coords_i[2])+radius < 0.5*NUCLEAR_ENVELOPE_WIDTH): # Filter
                continue
#            print i, coords_i, radius, distance,
            is_anchor=False
            for anchor_name,fg_name in fgs_anchors_regexp.iteritems():
                if re.search(anchor_name,p.get_name()):
#                    print "FG-Anchor_"+fgs_anchors_regexp[anchor_name]
                    is_anchor=True
                    print "Anchor", fg_name, coords_i
                    break
            if not is_anchor:
#                print "Obstacle"
                if radius in obstacles:
                    obstacles[radius].append(coords_i)
                else:
                    obstacles[radius]=[coords_i]

def handle_representation(r):
#    print r.get_name()
    if r.get_name()=="Beads":
        print "Beads", r.get_parent().get_name(), ";", r.get_name()
        handle_xyz_children(r)
    if re.search("Res:1$", r.get_name()): # what to do abuut res10 vs res1?
        print "Res:1", r.get_parent().get_name(), ";", r.get_name()
        for rr in r.get_children():
            handle_xyz_children(rr)

def get_coarse_grained_obstacles(obstacles):
  # bin obstacles by distance from central axis
  in_obstacles={}
  for (radius, coords) in obstacles.iteritems():
    for coord in coords:
      s=IMP.algebra.Sphere3D(coord,radius)
      dXY=math.sqrt(coord[0]**2+coord[1]**2)
      dXY=math.ceil(dXY/25.0)*25
      if dXY in in_obstacles:
          in_obstacles[dXY].append(s)
      else:
          in_obstacles[dXY]=[s]
  # coarse-grain each bin using an increasing resolution from center outward
  intermediate_obstacles=[]
  for dXY,spheres in in_obstacles.iteritems():
    resolution=min(20,  (dXY/150.0)**2)
    intermediate_obstacles=intermediate_obstacles+ \
        IMP.algebra.get_simplified_from_volume(spheres, resolution)
#    print "IN: ", spheres
#  print "INTERMEDIATE: ", dXY, " - ",  intermediate_obstacles
  out_obstacles={}
  for sphere in intermediate_obstacles:
     R=max(4.0, round(sphere.get_radius()))
     if not R in out_obstacles:
         out_obstacles[R]=[]
     out_obstacles[R].append(sphere.get_center())
#  print out_obstacles
  return out_obstacles




# ********* MAIN: *********
# Load information from RMF:
f=RMF.open_rmf_file_read_only(input_rmf)
m=IMP.Model()
h=IMP.rmf.create_hierarchies(f,m)
IMP.rmf.load_frame(f,0)
print("RMF", input_rmf, "loaded")
for nup in h[0].get_children():
    nup_name=nup.get_name()
    is_single_spoke= not re.search("@", nup_name) or re.search("@11$", nup_name)
    is_gle1_or_nup42= re.match("Gle1", nup_name) or (False and re.match("Nup42", nup_name))
    if is_single_spoke and not (is_gle1_or_nup42 and IS_REMOVE_GLE1_AND_NUP42):
        for r in nup.get_children():
            handle_representation(r)
# Add particles to configuration file
obstacles=get_coarse_grained_obstacles(obstacles)
XYZ=np.zeros([VOXELS_N, VOXELS_N, VOXELS_N])
for (radius, coords) in obstacles.iteritems():
   add_obstacles(XYZ,
       coords,
       radius*OBSTACLE_SCALE_FACTOR) # inflate obstacles a bit to prevent artificial holes
do_stats(XYZ)
