#########################################
# Example: imppy Scripts/load_whole_new_coarse_grained_v7.py config.pb InputData/wholeNPC_0.rmf3  > config.txt
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
from collections import OrderedDict
from FGParams import *

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
fgs_to_anchor_coords={}
TUNNEL_RADIUS=375
NUCLEAR_ENVELOPE_WIDTH=250
OBSTACLE_SCALE_FACTOR=2.0 # inflate obstacles a bit to prevent artificial cavities
kap_k=4.75
kap_range=4.5
sigma0_deg=45
sigma1_deg=45
kap_interaction_sites=4


def get_bead_nup_and_range(bead_name):
    '''
    Parses bead name to find nup name (first word)
    and residue range (second and third words)

    Returns nup name, first res and last res
    '''
    result= re.search('^([A-Za-z0-9]*)\.*[0-9]*_([0-9]*)-([0-9]*)_bead', bead_name)
#    print(bead_name)
#    print result.groups()
    if result is None:
        return False # not a bead leaf (which has a standard form e.g. Nup1_1-50_bead...)
#    assert(len(result.groups(0))==3)
    nup= result.group(1)
    res_from= int(result.group(2))
    res_to= int(result.group(3))
    return (nup, res_from, res_to)

def is_bead_in_fg_domain(bead_name, fg_params):
    try:
        [bead_nup, bead_from, bead_to]= get_bead_nup_and_range(bead_name)
    except:
#        print(bead_name, " is not an anchor bead")
        return False
    if bead_nup not in fg_params.keys():
        return False
    fg_froms=[fgp.res_from for fgp in fg_params[bead_nup].itervalues()]
    fg_tos=[fgp.res_to for fgp in fg_params[bead_nup].itervalues()]
    fg_from= min(fg_froms)
    fg_to= max(fg_tos)
    return (min(bead_to,fg_to)-max(bead_from,fg_from)>=0)

def is_anchor_bead(bead_name, fgs_regions_to_params):
    try:
        [bead_nup, bead_from, bead_to]= get_bead_nup_and_range(bead_name)
    except:
#        print(bead_name, " is not an anchor bead")
        return False
    if bead_nup not in fgs_regions_to_params.keys():
        return False
    anchor_fgp= fgs_regions_to_params[bead_nup]['anchor']
    return  \
        bead_from==anchor_fgp.res_from and \
        bead_to  ==anchor_fgp.res_to

def get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42):
    # TODO: adjust motif density in different nups, rg?
    fgs_regions_to_params={}
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
    fgs_regions_to_params['Nsp1']= {
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
    fgs_regions_to_params['Nup100']= {
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
    fgs_regions_to_params['Nup116']= {
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
    fgs_regions_to_params['Nup159']= {
        'Nup159': FGParams_from_default(
            default_fgp, 442, 881, 1.50, 0.07),
        'Nup159s': FGParams(
            res_from= 882,
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
        fgs_regions_to_params['Nup42']= {
            'Nup42': FGParams_from_default(
                default_fgp, 1, 363, 1.55, 0.13),
            'anchor': FGParams(res_from= 364,
                               res_to= 413)
        }
    fgs_regions_to_params['Nup49']= {
        'Nup49': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13), # TODO: there may be an FG motif in 233-236
        'anchor': FGParams(res_from= 201,
                           res_to= 269)
    }
    fgs_regions_to_params['Nup57']= {
        'Nup57': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13),
        'anchor': FGParams(res_from= 201,
                           res_to= 286)
    }
    fgs_regions_to_params['Nup145']= {
        'Nup145': FGParams_from_default(
            default_fgp, 1, 200, 1.55, 0.13),
        'anchor': FGParams(res_from= 201,
                           res_to= 225)
        # Nup145Ns?
    }
    fgs_regions_to_params['Nup1']= {
        'anchor': FGParams(res_from= 151,
                           res_to= 200), # SJ's anchor is at 301-351, so need some tweaking
        'Nup1s': FGParams_from_default(
            default_fgp, 201, 325, 1.6, 0.15), # NEED TO BE PARAMS OF NO FGs
        'Nup1m': FGParams_from_default( # TRIPLE CHECK - this region has more charges (Timney et al.)
            default_fgp, 326, 797, 1.55, 0.04),
        'Nup1c': FGParams_from_default( # TRIPLE CHECK - this region has less charges (Timney et al.)
            default_fgp, 798, 1076, 1.6, 0.15),
    }
    fgs_regions_to_params['Nup60']= {
        'anchor': FGParams(res_from= 351,
                           res_to= 398),
        'Nup60': FGParams_from_default(
            default_fgp, 399, 539, 1.55, 0.09)
    }
    return fgs_regions_to_params


def test_is_bead_in_fg_domain():
    fgs_regions_to_params= get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42)
    assert(get_bead_nup_and_range('Nsp1_1-50_bead_floppy')==('Nsp1',1,50))
    assert(is_bead_in_fg_domain('Nsp1_1-50_bead_floppy', fgs_regions_to_params))
    assert(is_bead_in_fg_domain('Nsp1_550-600_bead_floppy', fgs_regions_to_params))
    assert(is_bead_in_fg_domain('Nsp1_550-636_bead_floppy', fgs_regions_to_params))
    assert(is_bead_in_fg_domain('Nsp1_1-750_bead_floppy', fgs_regions_to_params))
    assert(is_bead_in_fg_domain('Nsp1_636-750_bead_floppy', fgs_regions_to_params))
    assert(not is_bead_in_fg_domain('Nsp1_637-750_bead_floppy', fgs_regions_to_params))
    assert(is_anchor_bead('Nup60_351-398_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60_350-398_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60_351-397_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60_1-7_bead_floppy', fgs_regions_to_params))
    assert(is_anchor_bead('Nup60.1_351-398_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60.1_350-398_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60.1_351-397_bead_floppy', fgs_regions_to_params))
    assert(not is_anchor_bead('Nup60.1_1-7_bead_floppy', fgs_regions_to_params))

def get_fg_regions_to_params_ordered(fg_regions_to_params):
    ''' return an OrderedDict object mapping from type to FGParams,
        with keys sorted from either the anchor (if it exists) or from
        the type corresponding to the first residue. If anchor is C',
        res_from and res_to are reversed to reflect that
    '''
    assert(len(fg_regions_to_params) >= 1)
    IS_ANCHOR= ('anchor' in fg_regions_to_params)
    ordered_list= [(key,value) for key,value in \
                   sorted(fg_regions_to_params.iteritems(), \
                          key=lambda (k,v): (v.res_from,k)) ]
    # Verify continuous seq:
    for i,(region_name,fg_params) in enumerate(ordered_list):
        if i>0:
            assert(fg_params.res_from == last_res_to+1)
        last_res_to= fg_params.res_to
    if (IS_ANCHOR):
        IS_REVERSE = ordered_list[-1][0]=='anchor' and \
                     len(ordered_list)>1
        assert(IS_REVERSE or (ordered_list[0][0] == 'anchor'))
    if IS_REVERSE: # reverse = C to N instead of N to C
        ordered_list.reverse()
        for (k, v) in ordered_list:
            tmp= v.res_from
            v.res_from= v.res_to
            v.res_to= tmp
    fg_regions_to_params_ordered= OrderedDict(ordered_list)
    return fg_regions_to_params_ordered

def get_bead_suffixes_and_is_anchor(fg_regions_to_params, fg_beads_per_res):
    # get a list of bead names
    fg_regions_to_params_ordered \
        = get_fg_regions_to_params_ordered(fg_regions_to_params)
    bead_suffixes= []
    total_nres= 0
    is_anchor= False
    for region, fg_params in fg_regions_to_params_ordered.iteritems():
        if(region=='anchor'):
            is_anchor= True
            assert(len(bead_suffixes)==0)
            bead_suffixes.append(region)
        else:
            nbeads= len(bead_suffixes)
            total_nres_beads= nbeads * fg_beads_per_res
            redundant_nres = total_nres_beads - total_nres
            region_nres= abs(fg_params.res_to - fg_params.res_from + 1)
            region_nres_to_add = region_nres - redundant_nres
            nbeads = int(math.ceil((region_nres_to_add+0.0) / FG_BEADS_PER_RES))
            for i in range(nbeads):
                bead_suffixes.append(region)
            total_nres= total_nres + region_nres
    return bead_suffixes, is_anchor

def add_fgs(config, type_name, fg_params, anchor_coordinates):
    bead_suffixes, is_anchor \
        = get_bead_suffixes_and_is_anchor(fg_params, FG_BEADS_PER_RES)
    nbeads = len(bead_suffixes)
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= nbeads,
                                      number=len(anchor_coordinates),
                                      radius=FG_RADIUS_PER_BEAD,
                                      interactions= FG_INTERACTIONS_PER_BEAD,
                                      rest_length_factor = FG_REST_LENGTH_FACTOR)
    if 'anchor' in bead_suffixes:
        assert(bead_suffixes[0]=='anchor' and 'anchor' not in bead_suffixes[1:])
        for coordinates in anchor_coordinates:
#            print(type_name, "Coords:", coordinates)
            pos =fgs.anchor_coordinates.add()
            pos.x= coordinates[0]
            pos.y= coordinates[1]
            pos.z= coordinates[2] + Z_TRANSFORM
    fgs.type_suffix_list.extend(bead_suffixes)
#    for bead_suffix in bead_suffixes:
#        suffixes= fgs.type_suffix_list.append(bead_suffix)
    return fgs




def handle_xyz_children(config, parent, fg_params):
    ''' handle a node whose children have xyz cooridnates '''
    for child in parent.get_children():
        xyzr=IMP.core.XYZR(child)
        coords=xyzr.get_coordinates()
        radius=round(xyzr.get_radius(),1)
        # Apply 8-fold symmetry
        for i in range(8):
            R=IMP.algebra.get_rotation_about_normalized_axis([0,0,1],
                                                             i*math.pi/4.0)
#            print parent.get_parent().get_name(), parent.get_name(), child.get_name(),
            coords_i=R*coords
            distance=math.sqrt(coords_i[0]**2+coords_i[1]**2)
            if(distance-radius>TUNNEL_RADIUS and
                    abs(coords_i[2])+radius < 0.5*NUCLEAR_ENVELOPE_WIDTH): # Filter
                continue
#            print i, coords_i, radius, distance,
            if is_anchor_bead(child.get_name(), fg_params):
#                print("Anchor bead found", child.get_name())
                [fg_name,res_from,res_to]= get_bead_nup_and_range(child.get_name())
                if fg_name in fgs_to_anchor_coords:
                    fgs_to_anchor_coords[fg_name].append(coords_i)
                else:
                    fgs_to_anchor_coords[fg_name]=[coords_i]
            else:
                if is_bead_in_fg_domain(child.get_name(), fg_params):
                    print "# Removing obstacle bead that is in fg_domain", child.get_name()
                    continue
#                print "Obstacle"
                if radius in obstacles:
                    obstacles[radius].append(coords_i)
                else:
                    obstacles[radius]=[coords_i]

def handle_nup_node(config, nup_node, nup_name, fg_params):
    node_name= nup_node.get_name()
#    print "Handling ", node_name
    if node_name=="Beads":
#        print "Beads", nup_name, ";", node_name
        handle_xyz_children(config, nup_node, fg_params)
    if re.search("Res:1$", node_name): # TODO: what to do about res10 vs res1?
        # node is root of res1 - grandparents are individual residues
#        print "Res:1", nup_name, ";", node_name
        for res1_node in nup_node.get_children():
            handle_xyz_children(config, res1_node, fg_params)

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

def add_obstacle_type(config, type_name, coords_list, radius):
    '''
    Add obstacles of specified type and radius to config
    using specified coordinates

    config - config file
    type_name - type name
    coords_list - list of 3D coordinate vectors
    radius - radius of all obstacles of this type
    Do nothing if len(coords)==0
    '''
    if len(coords_list)==0:
        return
    obstacle = IMP.npctransport.add_obstacle_type \
        (config, type_name=type_name, R =radius)
    for coords in coords_list:
        pos=obstacle.xyzs.add()
        pos.x=coords[0]
        pos.y=coords[1]
        pos.z=coords[2] + Z_TRANSFORM


def add_obstacles_from_rmf(config, input_rmf, fg_params):
    global obstacles
    # Add all obstacles based on RMF file input_rmf
    f=RMF.open_rmf_file_read_only(input_rmf)
    m=IMP.Model()
    h=IMP.rmf.create_hierarchies(f,m)
#    print("Hierarchy", h)
#    print("Hierarchy root", h[0].get_name())
    IMP.rmf.load_frame(f,0)
    for nup_root in h[0].get_children():
        nup_name=nup_root.get_name()
#        print("Nup root", nup_name)
        is_single_spoke= not re.search("@", nup_name) or re.search("@11$", nup_name)
        is_gle1_or_nup42= re.match("Gle1", nup_name) or (False and re.match("Nup42", nup_name))
        if is_single_spoke and is_gle1_or_nup42:
            print("Gle1/Nup42 detected in main spoke:", nup_name)
        is_remove= (is_gle1_or_nup42 and IS_REMOVE_GLE1_AND_NUP42)
        if is_single_spoke and not is_remove:
            for nup_node in nup_root.get_children():
                handle_nup_node(config, nup_node, nup_name, fg_params)
        # Add particles to configuration file
    if(COARSE_GRAINED_OBSTACLES):
        obstacles=get_coarse_grained_obstacles(obstacles)
    for (radius, coords) in obstacles.iteritems():
        add_obstacle_type(config,
                      "obstacles %.1f" % radius,
                      coords,
                      radius*OBSTACLE_SCALE_FACTOR) # inflate obstacles a bit to prevent artificial holes

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


# *************************
# ********* MAIN: *********
# *************************
test_is_bead_in_fg_domain() # just to test it works correctly
config= get_basic_config()
# Add FGs:
fgs_regions_to_params= get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42)
add_obstacles_from_rmf(config, input_rmf, fgs_regions_to_params)
#print("FGs to anchor keys", fgs_to_anchor_coords.keys())
for (name,coords) in fgs_regions_to_params.iteritems(): # TODO: get coords right
    add_fgs(config,
            name,
            fgs_regions_to_params[name],
            fgs_to_anchor_coords[name])
# Add fg-fg interactions:
for fg0, regions_to_params0 in fgs_regions_to_params.iteritems():
    for fg1, regions_to_params1 in fgs_regions_to_params.iteritems():
        for region0, params0 in regions_to_params0.iteritems():
            for region1, params1 in regions_to_params1.iteritems():
#                print(fg0, region0, fg1, region1)
                if params0.self_k is None or params1.self_k is None: # (e.g. anchor)
                    continue
                self_k = 0.5*(params0.self_k+params1.self_k)
                self_range = 0.5*(params0.self_range+params1.self_range)
                interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                                   name0= fg0+region0,
                                                                   name1= fg1+region1,
                                                                   interaction_k= self_k,
                                                                   interaction_range= self_range)
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
    # Add interactions:
    for fg_name, regions_to_params in fgs_regions_to_params.iteritems():
        for region, params in regions_to_params.iteritems():
            IMP.npctransport.add_interaction(config,
                                             name0=fg_name + region,
                                             name1=inert_name,
                                             interaction_k=0,
                                             interaction_range=0)
            IMP.npctransport.add_interaction(config,
                                             name0=fg_name + region,
                                             name1=kap_name,
                                             interaction_k= params.kap_k,
                                             interaction_range= params.kap_range,
                                             range_sigma0_deg= sigma0_deg,
                                             range_sigma1_deg= sigma1_deg,
        )

# Add obstacles


# dump to file
f=open(outfile, "wb")
f.write(config.SerializeToString())
print config
