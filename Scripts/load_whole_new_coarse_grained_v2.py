#!/usr/bin/bash /Users/barak/imp/fast/setup_environment.sh /anaconda/bin/python

#########################################
# Example: imppy Scripts/load_whole.py config.pb InputData/wholeNPC_0.rmf3  > config.txt
#
# TODO: for each FG, map each beads has interaction sites and their strength factor (due to e.g. density)
# TODO: we might want to give the anchor a radius of its own
#
# Author: Barak Raveh barak@salilab.org
# Data written: Aug 16, 2016
# Date last udpated: Aug 22, 2016 or later
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

IS_REMOVE_GLE1_AND_NUP42=True
outfile = sys.argv[1]
input_rmf = sys.argv[2] #'wholeNPC_0.rmf3')

fgs_anchors_regexp={}
fgs_anchors_regexp['Nsp1.*_601-636_']='Nsp1'
fgs_anchors_regexp['Nup1.*_301-350_']='Nup1'
if True or not IS_REMOVE_GLE1_AND_NUP42:
    fgs_anchors_regexp['Nup42.*_364-413_']='Nup42'
fgs_anchors_regexp['Nup49.*_201-269_']='Nup49'
fgs_anchors_regexp['Nup57.*_201-286_']='Nup57'
fgs_anchors_regexp['Nup60.*_351-398_']='Nup60'
fgs_anchors_regexp['Nup100.*_551-575_']='Nup100'
fgs_anchors_regexp['Nup116.*_751-775_']='Nup116'
fgs_anchors_regexp['Nup145.*_201-225_']='Nup145'
fgs_anchors_regexp['Nup159.*_1082-1116_']='Nup159'
# TODO: adjust motif density in different nups, rg?
fgs_nres={}
fgs_nres['Nsp1']=600
fgs_nres['Nup1']=700 # double check
if True or not IS_REMOVE_GLE1_AND_NUP42:
    fgs_nres['Nup42']=363
fgs_nres['Nup49']=200
fgs_nres['Nup57']=200
fgs_nres['Nup60']=110 # double check
fgs_nres['Nup100']=550
fgs_nres['Nup116']=750
fgs_nres['Nup145']=200
fgs_nres['Nup159']=600 # double check, what to do with C'?

Z_TRANSFORM=0 #-75.5
FG_BEADS_PER_RES=20 #10
FG_RADIUS_PER_BEAD=8.5*math.sqrt(FG_BEADS_PER_RES/15.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.5
obstacles={}
fgs={}
TUNNEL_RADIUS=350
NUCLEAR_ENVELOPE_WIDTH=300
kap_k=4.75
kap_range=4.5
sigma0_deg=45
sigma1_deg=45
kap_interaction_sites=4

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=10
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=1
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=5
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=20 # TODO: change to 10
    # non-specific attraction
    config.nonspecific_range.lower= 8.0
    config.nonspecific_k.lower= 0.01
    config.slack.lower = 30
    config.number_of_trials=1
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water
    config.statistics_interval_ns=10.0
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
#config.dump_interval_ns=0.01
#config.simulation_time_ns=5
    config.dump_interval_ns=1000
    config.simulation_time_ns=1000
    config.box_is_on.lower=1
    config.box_side.lower=2000
    config.slab_is_on.lower=1
    config.slab_thickness.lower=NUCLEAR_ENVELOPE_WIDTH
    config.tunnel_radius.lower=TUNNEL_RADIUS
    config.is_xyz_hist_stats=1
    return config


def add_fgs(config, type_name, nres, anchor_coordinates):
    nbeads = 1 + int(math.ceil(nres / FG_BEADS_PER_RES)) # +1 for anchor
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= nbeads,
                                      number=len(anchor_coordinates),
                                      radius=FG_RADIUS_PER_BEAD,
                                      interactions= FG_INTERACTIONS_PER_BEAD,
                                      rest_length_factor = FG_REST_LENGTH_FACTOR)
    for coordinates in anchor_coordinates:
        pos =fgs.anchor_coordinates.add()
        pos.x= coordinates[0]
        pos.y= coordinates[1]
        pos.z= coordinates[2] + Z_TRANSFORM
    return fgs

def add_obstacles(config, type_name, coords_list, radius):
    '''
    config - config file
    type_name - type name
    coords_list - list of 3D coordinate vectors
    radius - radius of all obstacles of this type

    Retruns obstacles protobuf object or None if len(coords)==0
    '''
    if len(coords_list)==0:
        return None
    obstacle = IMP.npctransport.add_obstacle_type \
        (config, type_name=type_name, R =radius)
    for coords in coords_list:
        pos=obstacle.xyzs.add()
        pos.x=coords[0]
        pos.y=coords[1]
        pos.z=coords[2] + Z_TRANSFORM
    return obstacles



def handle_xyz_children(config, parent):
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
                    if fg_name in fgs:
                        fgs[fg_name].append(coords_i)
                    else:
                        fgs[fg_name]=[coords_i]
                    break
            if not is_anchor:
#                print "Obstacle"
                if radius in obstacles:
                    obstacles[radius].append(coords_i)
                else:
                    obstacles[radius]=[coords_i]

def handle_representation(config, r):
#    print r.get_name()
    if r.get_name()=="Beads":
        print "Beads", r.get_parent().get_name(), ";", r.get_name()
        handle_xyz_children(config,r)
    if re.search("Res:1$", r.get_name()): # what to do abuut res10 vs res1?
        print "Res:1", r.get_parent().get_name(), ";", r.get_name()
        for rr in r.get_children():
            handle_xyz_children(config, rr)

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
config= get_basic_config()
# Load information from RMF:
f=RMF.open_rmf_file_read_only(input_rmf)
m=IMP.Model()
h=IMP.rmf.create_hierarchies(f,m)
IMP.rmf.load_frame(f,0)
for nup in h[0].get_children():
    nup_name=nup.get_name()
    is_single_spoke= not re.search("@", nup_name) or re.search("@11$", nup_name)
    is_gle1_or_nup42= re.match("Gle1", nup_name) or (False and re.match("Nup42", nup_name))
    if is_single_spoke and not (is_gle1_or_nup42 and IS_REMOVE_GLE1_AND_NUP42):
        for r in nup.get_children():
            handle_representation(config, r)
# Add particles to configuration file
obstacles=get_coarse_grained_obstacles(obstacles)
for (radius, coords) in obstacles.iteritems():
    add_obstacles(config, "obstacles %.1f" % radius, coords, radius)
for (name,coords) in fgs.iteritems():
    add_fgs(config, name, fgs_nres[name], coords)
# Add fg-fg interactions:
for name0 in fgs.iterkeys():
    for name1 in fgs.iterkeys():
        interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                           name0= name0,
                                                           name1= name1,
                                                           interaction_k= 0.01,
                                                           interaction_range= 6)
# Add kaps and inerts:
nonspecifics={}
kaps={}
rrange=range(14,30,2)
for radius in rrange:
    inert_name="R%d" % radius
    nonspecifics[radius]= IMP.npctransport.add_float_type(config,
                                                  number=200,
                                                  radius=radius,
                                                  type_name=inert_name,
                                                  interactions=0)
    kap_name="kap%d" % radius
    kaps[radius]= IMP.npctransport.add_float_type(config,
                                             number=200,
                                             radius=radius,
                                             type_name=kap_name,
                                             interactions=kap_interaction_sites)
    for fg_name in fgs.iterkeys():
        IMP.npctransport.add_interaction(config,
                                         name0=fg_name,
                                         name1=inert_name,
                                         interaction_k=0,
                                         interaction_range=0)
        IMP.npctransport.add_interaction(config,
                                         name0=fg_name,
                                         name1=kap_name,
                                         interaction_k=kap_k,
                                         interaction_range=kap_range,
                                         range_sigma0_deg=sigma0_deg,
                                         range_sigma1_deg=sigma1_deg,
        )


# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
print config
