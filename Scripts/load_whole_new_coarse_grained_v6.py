#########################################
# Example: imppy Scripts/load_whole_new_coarse_grained_v6.py config.pb InputData/wholeNPC_0.rmf3  > config.txt
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
import pandas as pd
from FGParams import *

def get_bead_nup_and_range(bead_name):
    result= re.search('^([A-Za-z0-9]*)_([0-9]*)-([0-9]*)_bead', bead_name)
    if result is None:
        return False
    assert(len(result.groups(0))==3)
    nup= result.group(1)
    res_from= int(result.group(2))
    res_to= int(result.group(3))
    return (nup, res_from, res_to)

def is_bead_in_fg_domain(bead_name, fg_params):
    [bead_nup, bead_from, bead_to]= get_bead_nup_and_range(bead_name)
    if bead_nup not in fg_params.keys():
        return False
    fg_froms=[fgp.res_from for fgp in fg_params[bead_nup].itervalues()]
    fg_tos=[fgp.res_to for fgp in fg_params[bead_nup].itervalues()]
    fg_from= min(fg_froms)
    fg_to= max(fg_tos)
    return (min(bead_to,fg_to)-max(bead_from,fg_from)>=0)

def is_anchor_bead(bead_name, fg_params):
    [bead_nup, bead_from, bead_to]= get_bead_nup_and_range(bead_name)
    if bead_nup not in fg_params.keys():
        return False
    anchor_fgp= fg_params[bead_nup]['anchor']
    return  \
        bead_from==anchor_fgp.res_from and \
        bead_to  ==anchor_fgp.res_to


def get_fg_params(IS_REMOVE_GLE1_AND_NUP42):
    # TODO: adjust motif density in different nups, rg?
    fg_params={}
    nan= float('nan')
    default_fgp=FGParams( \
                          res_from=nan,
                          res_to= nan,
                          self_range= 6.00,
                          kap_k= 3.58,
                          kap_range= 4.95,
                          nonspec_range= 5.00,
                          backbone_k= 0.0075,
                          backbone_tau= 50 )
    fg_params['Nsp1']= {
        'Nsp1N': FGParams_from_default(
            default_fgp, 1 ,180, 1.55, 0.13),
        'Nsp1C': FGParams_from_default(
            default_fgp, 181, 550, 1.50, 0.01),
        'Nsp1s': FGParams( \
            res_from= 551,
            res_to= 600,
            self_k= 1.5,
            self_range= 6.00,
            kap_k= nan,
            kap_range= nan,
            nonspec_k= 0.13,
            nonspec_range= 5.00,
            backbone_k= 0.0075,
            backbone_tau= 50),
        'anchor': FGParams(res_from= 601,
                           res_to= 636)
    }
    fg_params['Nup100']= {
        'Nup100': FGParams_from_default(
            default_fgp, 1, 550, 1.55, 0.13),
        'Nup100s': FGParams(
            res_from= 551,
            res_to= 800,
            self_k= 1.50,
            self_range= 6.00,
            kap_k= nan,
            kap_range= nan,
            nonspec_k= 0.07,
            nonspec_range= 5.00,
            backbone_k= 0.0075,
            backbone_tau= 50),
        'anchor': FGParams(res_from= 801,
                           res_to= 815)
    }
    # Nup100s?
    fg_params['Nup116']= {
        'Nup116': FGParams_from_default(
            default_fgp, 1, 750, 1.55, 0.13),
        'Nup116s': FGParams(
            res_from= 751,
            res_to= 950,
            self_k= 1.50,
            self_range= 6.00,
            kap_k= nan,
            kap_range= nan,
            nonspec_k= 0.07,
            nonspec_range= 5.00,
            backbone_k= 0.0075,
            backbone_tau= 50),
        'anchor': FGParams(res_from= 951,
                           res_to= 965)
    }
    fg_params['Nup159']= {
        'Nup159': FGParams_from_default(
            default_fgp, 441, 881, 1.50, 0.07),
        'Nup159s': FGParams(
            res_from= 881,
            res_to= 1081,
            self_k= 1.50,
            self_range= 6.00,
            kap_k= nan,
            kap_range= nan,
            nonspec_k= 0.07,
            nonspec_range= 5.00,
            backbone_k= 0.0075,
            backbone_tau= 50),
        'anchor': FGParams(res_from= 1082,
                           res_to= 1116)
    }
    if not IS_REMOVE_GLE1_AND_NUP42:
        fg_params['Nup42']= {
            'Nup42': FGParams_from_default(
                default_fgp, 1, 363, 1.55, 0.13),
            'anchor': FGParams(res_from= 364,
                               res_to= 413)
        }
    fg_params['Nup49']= {
        'Nup49': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13), # TODO: there may be an FG motif in 233-236
        'anchor': FGParams(res_from= 201,
                           res_to= 269)
    }
    fg_params['Nup57']= {
        'Nup57': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13),
        'anchor': FGParams(res_from= 201,
                           res_to= 286)
    }
    fg_params['Nup145']= {
        'Nup145': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13),
        'anchor': FGParams(res_from= 201,
                           res_to= 225)
        # Nup145Ns?
    }
    fg_params['Nup1']= {
        'anchor': FGParams(res_from= 151,
                           res_to= 200), # SJ's anchor is at 301-351, so need some tweaking
        'Nup1s': FGParams_from_default(
            default_fgp, 201, 325, 1.6, 0.15), # NEED TO BE PARAMS OF NO FGs
        'Nup1m': FGParams_from_default( # TRIPLE CHECK - this region has more charges (Timney et al.)
            default_fgp, 325, 797, 1.55, 0.04),
        'Nup1c': FGParams_from_default( # TRIPLE CHECK - this region has less charges (Timney et al.)
            default_fgp, 798, 1076, 1.6, 0.15),
    }
    fg_params['Nup60']= {
        'anchor': FGParams(res_from= 351,
                           res_to= 398),
        'Nup60': FGParams_from_default(
            default_fgp, 399, 539, 1.55, 0.09)
    }
    return fg_params


def test_is_bead_in_fg_domain():
    fg_params= get_fg_params(IS_REMOVE_GLE1_AND_NUP42=False)
    assert(get_bead_nup_and_range('Nsp1_1-50_bead_floppy')==('Nsp1',1,50))
    assert(is_bead_in_fg_domain('Nsp1_1-50_bead_floppy', fg_params))
    assert(is_bead_in_fg_domain('Nsp1_550-600_bead_floppy', fg_params))
    assert(is_bead_in_fg_domain('Nsp1_550-636_bead_floppy', fg_params))
    assert(is_bead_in_fg_domain('Nsp1_1-750_bead_floppy', fg_params))
    assert(is_bead_in_fg_domain('Nsp1_636-750_bead_floppy', fg_params))
    assert(not is_bead_in_fg_domain('Nsp1_637-750_bead_floppy', fg_params))
    assert(is_anchor_bead('Nup60_351-398_bead_floppy', fg_params))
    assert(not is_anchor_bead('Nup60_350-398_bead_floppy', fg_params))
    assert(not is_anchor_bead('Nup60_351-397_bead_floppy', fg_params))
    assert(not is_anchor_bead('Nup60_1-7_bead_floppy', fg_params))

test_is_bead_in_fg_domain()
sys.exit()


IS_REMOVE_GLE1_AND_NUP42=True
outfile = sys.argv[1]
input_rmf = sys.argv[2] #'wholeNPC_0.rmf3')
COARSE_GRAINED_OBSTACLES= True
Z_TRANSFORM=0 #-75.5
FG_BEADS_PER_RES=20 #10
FG_RADIUS_PER_BEAD=6.0*(FG_BEADS_PER_RES/15.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.5*math.sqrt(15.0/FG_BEADS_PER_RES) # make sure that radius*rest_length scale with sqrt(# beads)
obstacles={}
fgs={}
TUNNEL_RADIUS=375
NUCLEAR_ENVELOPE_WIDTH=250
OBSTACLE_SCALE_FACTOR=2.0 # inflate obstacles a bit to prevent artificial cavities
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
    config.backbone_k.lower=0.2
    #config.time_step_factor.lower=0.3
    config.time_step_factor.lower=4
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=10 # TODO: change to 10
    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
    config.slack.lower = 30
    config.number_of_trials=1
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water
    config.statistics_interval_ns=10.0
    config.output_statistics_interval_ns=100.0
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


def add_fgs(config, type_name, nup_fg_params, anchor_coordinates):
    fg_froms= [fgp.res_from for fgp in nup_fg_params.itervalues()].sort()
    fg_tos=[fgp.res_to for fgp in nup_fg_params.itervalues()].sort()
    fg_names=
    for i,fg_from in enumerate(fg_froms):
        if i>0:
            assert(fg_from==fg_tos[i-1]+1)
    nres= fg_tos[-1]-fg_froms[0]+1
    nbeads = 1 + int(math.ceil((nres+0.0) / FG_BEADS_PER_RES)) # +1 for anchor
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
`   # TODO: add suffixes according to information in nup_fg_params
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



def handle_xyz_children(config, parent, fg_params):
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
            if is_anchor_bead(p.get_name(), fg_params):
                [fg_name,,]= get_bead_nup_and_range(p.get_name())
                if fg_name in fgs:
                    fgs[fg_name].append(coords_i)
                else:
                    fgs[fg_name]=[coords_i]
                break
            else:
                if is_bead_in_fg_domain(p.get_name(), fg_params):
                    continue
#                print "Obstacle"
                if radius in obstacles:
                    obstacles[radius].append(coords_i)
                else:
                    obstacles[radius]=[coords_i]

def handle_representation(config, r, fg_params):
#    print r.get_name()
    if r.get_name()=="Beads":
        print "Beads", r.get_parent().get_name(), ";", r.get_name()
        handle_xyz_children(config,r, fg_params)
    if re.search("Res:1$", r.get_name()): # what to do abuut res10 vs res1?
        print "Res:1", r.get_parent().get_name(), ";", r.get_name()
        for rr in r.get_children():
            handle_xyz_children(config, rr, fg_params)

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
    resolution=min(20,  (dXY/110.0)**1.35)
    intermediate_obstacles=intermediate_obstacles+ \
        IMP.algebra.get_simplified_from_volume(spheres, resolution)
#    print "IN: ", spheres
#  print "INTERMEDIATE: ", dXY, " - ",  intermediate_obstacles
  out_obstacles={}
  for sphere in intermediate_obstacles:
     R=max(4.0, round(sphere.get_radius()))
     if not (R in out_obstacles):
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
fg_params= get_fg_params(IS_REMOVE_GLE1_AND_NUP42)
for nup in h[0].get_children():
    nup_name=nup.get_name()
    is_single_spoke= not re.search("@", nup_name) or re.search("@11$", nup_name)
    is_gle1_or_nup42= re.match("Gle1", nup_name) or (False and re.match("Nup42", nup_name))
    if is_single_spoke and not (is_gle1_or_nup42 and IS_REMOVE_GLE1_AND_NUP42):
        for r in nup.get_children():
            handle_representation(config, r, fg_params)
# Add particles to configuration file
if(COARSE_GRAINED_OBSTACLES):
    obstacles=get_coarse_grained_obstacles(obstacles)
for (radius, coords) in obstacles.iteritems():
    add_obstacles(config,
        "obstacles %.1f" % radius,
        coords,
        radius*OBSTACLE_SCALE_FACTOR) # inflate obstacles a bit to prevent artificial holes
for (name,coords) in fgs.iteritems():
    add_fgs(config, name, fgs_params[name], coords)
# TODO: fix add_fgs so that it will be based on fg_params
# Add fg-fg interactions:
# TODO: add interaction for all FG subtypes based on fg_params
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
