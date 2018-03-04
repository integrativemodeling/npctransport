#########################################
# Example: imppy Scripts/load_whole_new_coarse_grained_v7.py config.pb InputData/wholeNPC_0.rmf3  > config.txt
#
# TODO: for each FG, map each beads has interaction sites and their strength factor (due to e.g. density)
# TODO: we might want to give the anchor a radius of its own
#
# Author: Barak Raveh barak@salilab.org
# Data written: Aug 16, 2016
# Date last udpated: Mar 3, 2018 or later
#########################################
from __future__ import print_function
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
import argparse
from collections import OrderedDict
from collections import defaultdict
from FGParamsFactory import *


IS_REMOVE_GLE1_AND_NUP42=True
#outfile = sys.argv[1]
#input_rmf = sys.argv[2] #'wholeNPC_0.rmf3')
COARSE_GRAINED_OBSTACLES= True
Z_TRANSFORM=0 #-75.5
FG_BEADS_PER_RES=20 #10
FG_RADIUS_PER_BEAD=8.0*(FG_BEADS_PER_RES/20.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.9*math.sqrt(20.0/FG_BEADS_PER_RES) # make sure that radius*rest_length scale with sqrt(# beads)
obstacles={}
fgs_to_anchor_coords={}
TUNNEL_RADIUS=375
NUCLEAR_ENVELOPE_WIDTH=250
OBSTACLE_SCALE_FACTOR=1.0 # inflate obstacles a bit to prevent artificial cavities
kap_k=4.75
kap_range=4.5
sigma0_deg=45
sigma1_deg=45
kap_interaction_sites=4

def parse_commandline():
    parser = argparse.ArgumentParser(description='Create a config file loaded from an RMF model of the full NPC')
    parser.add_argument('config_file', metavar='config_file', type=str,
                        default='config.pb',
                        help='output configuration file')
    parser.add_argument('input_rmf_file', metavar='input_rmf_file', type=str,
                        default='config.pb',
                        help='rmf file of NPC model (e.g., InputData/47-35_1spoke.rmf3)')
    parser.add_argument('--diffuser_sizes', metavar='diffuser_sizes', type=int, nargs='*',
                        default=range(14,29,2),
                        help='list of diffuser sizes in A')
    parser.add_argument('--n_diffusers', type=int,
                        default=200,
                        help='number of diffuser molecules for each type of diffuser')
    args = parser.parse_args()
    return args


def get_node_nup_and_range_by_name(node_name, parent_name, grandparent_name):
    if parent_name=='Beads':
        identifier= node_name
    else:
        identifier= parent_name
    result= re.search('^([A-Za-z0-9]*)\.*[0-9]*@*[0-9]*_([0-9]*)-([0-9]*)_(bead|pdb)',
                      identifier)
#    print(identifier, node_name, parent_name)
#    print result.groups()
    if result is not None:
        nup= result.group(1)
        assert(grandparent_name is None or re.match(nup,grandparent_name))
        if parent_name=='Beads':
            res_from= int(result.group(2))
            res_to= int(result.group(3))
        else:
            res_from= int(node_name)
            res_to= int(node_name)
    else:
        if parent_name=='Beads':
            try:
                result= re.search('^(.*)@*[0-9]*$', grandparent_name)
                nup= result.group(1)
                res_from= int(node_name)
                res_to= int(node_name)
            except:
                print("Couldn't parse bead node", node_name, "grandparent", grandparent_name)
                return False
        else:
            return False # not a proper leaf for patsing
    #    assert(len(result.groups(0))==4)
    return (nup, res_from, res_to)

def get_node_nup_and_range(node):
    '''
    Parses node name (for beads) for parent's name (for res1 or res10
    nodes from PDB), and return nup name (first word) and residue
    range (second and third words)

    Returns nup name, first res and last res
    '''
    node_name= node.get_name()
    parent= node.get_parent()
    parent_name= parent.get_name()
    grandparent_name= None
    if(parent_name=='Beads'):
        grandparent_name= parent.get_parent().get_name()
    return get_node_nup_and_range_by_name(node_name, parent_name, grandparent_name)


def is_node_in_fg_domain(node, fg_params):
    try:
        [nup, res_from, res_to]= get_node_nup_and_range(node)
    except:
#        print(node_name, " is not an anchor node")
        return False
    if nup not in fg_params.keys():
        return False
    fg_froms=[fgp.res_from for fgp in fg_params[nup].itervalues()]
    fg_tos=[fgp.res_to for fgp in fg_params[nup].itervalues()]
    fg_from= min(fg_froms)
    fg_to= max(fg_tos)
    return (min(res_to,fg_to)-max(res_from,fg_from)>=0)

def is_anchor_node(node, fgs_regions_to_params):
    try:
#        print("is anchor", node.get_name())
        [nup, res_from, res_to]= get_node_nup_and_range(node)
    except:
        #        print(node_name, " is not an anchor bead")
        raise
        return False
    if nup not in fgs_regions_to_params.keys():
        return False
    anchor_fgp= fgs_regions_to_params[nup]['anchor']
    return  \
        res_from==anchor_fgp.res_from and \
        res_to  ==anchor_fgp.res_to

def get_fgs_regions_to_params(is_remove_gle1_and_nup42):
    # TODO: adjust motif density in different nups, rg?
    fgs_regions_to_params={}
    nan= float('nan')
    default_fgp=FGParamsFactory \
                 ( res_from=nan,
                   res_to= nan,
                   self_k= 1.50,
                   self_range= 6.00,
                   kap_k= 3.58,
                   kap_range= 4.95,
                   nonspec_k= 0.07,
                   nonspec_range= 5.00,
                   backbone_k= 0.0075,
                   backbone_tau= 50 )
    default_disordered_fgp= default_fgp.get_copy()
    default_disordered_fgp.kap_k= None
    default_disordered_fgp.kap_range= None
    default_FSFG= default_fgp.get_copy()
    default_FSFG.self_k= 1.5
    default_FSFG.nonspec_k= 0.01
    default_GLFG= default_fgp.get_copy()
    default_GLFG.self_k= 1.55
    default_GLFG.nonspec_k= 0.13
    fgs_regions_to_params['Nsp1']= {
        'N': default_GLFG.get_copy(1, 180),
        'C': default_FSFG.get_copy(181, 550),
        's': default_disordered_fgp.get_copy(551, 636),
        'anchor': FGParamsFactory(res_from= 637,
                           res_to= 637)
    }
    fgs_regions_to_params['Nup100']= {
        '': default_GLFG.get_copy(1, 550),
        's':  default_disordered_fgp.get_copy(551,800),
        'anchor': FGParamsFactory(res_from= 801,
                           res_to= 815)
    }
    fgs_regions_to_params['Nup116']= {
        '': default_GLFG.get_copy(1, 750),
        's': default_disordered_fgp.get_copy(751, 950),
        'anchor': FGParamsFactory(res_from= 951,
                           res_to= 965)
    }
    fgs_regions_to_params['Nup159']= {
        '': default_FSFG.get_copy(442, 881,
                                        nonspec_k= 0.07),
        's': default_disordered_fgp.get_copy(882, 1116),
        'anchor': FGParamsFactory(res_from= 1117,
                           res_to= 1117)
    }
    if not is_remove_gle1_and_nup42:
        fgs_regions_to_params['Nup42']= {
            '': default_GLFG.get_copy(1, 363),
            'anchor': FGParamsFactory(res_from= 364,
                               res_to= 413)
        }
    fgs_regions_to_params['Nup49']= {
        '': default_GLFG.get_copy(1, 240),
        's': default_disordered_fgp.get_copy(241, 269),
        'anchor': FGParamsFactory(res_from= 270,
                           res_to= 270)
    }
    fgs_regions_to_params['Nup57']= {
        '': default_GLFG.get_copy(1, 200),
        's': default_disordered_fgp.get_copy(201, 286),
        'anchor': FGParamsFactory(res_from= 287,
                           res_to= 287)
    }
    fgs_regions_to_params['Nup145']= {
        '': default_GLFG.get_copy(1, 200),
        's': default_disordered_fgp.get_copy(201,250),
        'anchor': FGParamsFactory(res_from= 251,
                           res_to= 275)
        # Nup145Ns?
    }
    fgs_regions_to_params['Nup1']= {
        'anchor': FGParamsFactory(res_from= 151,
                           res_to= 200),
        's': default_disordered_fgp.get_copy(201, 325),
        'm': default_FSFG.get_copy(326, 815),
        'c': default_GLFG.get_copy(816, 1076)
    }
    fgs_regions_to_params['Nup60']= {
        'anchor': FGParamsFactory(res_from= 251,
                           res_to= 300),
        's': default_disordered_fgp.get_copy(301, 398),
        '': default_FSFG.get_copy(399, 539)
    }
    return fgs_regions_to_params


def test_is_node_in_fg_domain():
    m= IMP.Model()
    p_Nsp1= IMP.Particle(m, 'Nsp1')
    p_Nup60= IMP.Particle(m, 'Nup60')
    p_Nsp1_beads= IMP.Particle(m, 'Beads')
    p_Nup60_beads= IMP.Particle(m, 'Beads')
    p_Nsp1_1_50= IMP.Particle(m, 'Nsp1_1-50_bead_floppy')
    p_Nsp1_550_600= IMP.Particle(m, 'Nsp1_550-600_bead_floppy')
    p_Nsp1_550_636= IMP.Particle(m, 'Nsp1_550-636_bead_floppy')
    p_Nsp1_637_750= IMP.Particle(m, 'Nsp1_637-750_bead_floppy')
    p_Nup60_350_398= IMP.Particle(m, 'Nup60_350-398_bead_floppy')
    p_Nup60_351_397= IMP.Particle(m, 'Nup60_351-397_bead_floppy')
    p_Nup60_1_7= IMP.Particle(m, 'Nup60_1-7_bead_floppy')
    p_Nsp1_1_50_pdb= IMP.Particle(m, 'Nsp1_1-50_pdb')
    p_Nsp1_50= IMP.Particle(m, '50')
    for pi in m.get_particle_indexes():
        IMP.atom.Hierarchy.setup_particle(m, pi)
    IMP.atom.Hierarchy(p_Nsp1).add_child(p_Nsp1_beads)
    IMP.atom.Hierarchy(p_Nup60).add_child(p_Nup60_beads)
    IMP.atom.Hierarchy(p_Nsp1_beads).add_child(p_Nsp1_1_50)
    IMP.atom.Hierarchy(p_Nsp1_beads).add_child(p_Nsp1_550_600)
    IMP.atom.Hierarchy(p_Nsp1_beads).add_child(p_Nsp1_550_636)
    IMP.atom.Hierarchy(p_Nsp1_beads).add_child(p_Nsp1_637_750)
    fgs_regions_to_params= get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42)
    assert(get_node_nup_and_range(IMP.atom.Hierarchy(p_Nsp1_1_50))==('Nsp1',1,50))
    assert(is_node_in_fg_domain(IMP.atom.Hierarchy(p_Nsp1_1_50), fgs_regions_to_params))
    # assert(is_node_in_fg_domain('Nsp1_550-600_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(is_node_in_fg_domain('Nsp1_550-636_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(is_node_in_fg_domain('Nsp1_1-750_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(is_node_in_fg_domain('Nsp1_636-750_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_node_in_fg_domain('Nsp1_637-750_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(is_anchor_node('Nup60_351-398_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60_350-398_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60_351-397_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60_1-7_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(is_anchor_node('Nup60.1_351-398_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60.1_350-398_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60.1_351-397_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(not is_anchor_node('Nup60.1_1-7_bead_floppy', 'Beads', fgs_regions_to_params))
    # assert(get_node_nup_and_range('50', 'Nsp1_1-50_pdb')==('Nsp1',50,50))
    # assert(is_node_in_fg_domain('50', 'Nsp1_1-50_pdb', fgs_regions_to_params))

def get_fg_regions_to_params_ordered(fg_regions_to_params):
    ''' return an OrderedDict object mapping from type to FGParamsFactory,
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
    global fgs_to_anchor_coords
    for child in parent.get_children():
        xyzr=IMP.core.XYZR(child)
        coords=xyzr.get_coordinates()
        radius=round(xyzr.get_radius(),1)
        # Apply 8-fold symmetry
        for i in range(8):
            R=IMP.algebra.get_rotation_about_normalized_axis([0,0,1],
                                                             i*math.pi/4.0)
#            print parent.get_parent().get_name(), parent.get_name(), child_name,
            coords_i=R*coords
            distance=math.sqrt(coords_i[0]**2+coords_i[1]**2)
            if(distance-radius>TUNNEL_RADIUS and
                    abs(coords_i[2])+radius < 0.5*NUCLEAR_ENVELOPE_WIDTH): # Filter
                continue
#            print i, coords_i, radius, distance,
            if is_anchor_node(child, fg_params):
#                print("Anchor bead found", child_name)
                [fg_name,res_from,res_to]= get_node_nup_and_range(child)
                if fg_name in fgs_to_anchor_coords:
                    fgs_to_anchor_coords[fg_name].append(coords_i)
                else:
                    fgs_to_anchor_coords[fg_name]=[coords_i]
            else:
                if is_node_in_fg_domain(child, fg_params):
                    print("# Removing obstacle node that is in fg_domain",
                          child.get_name(), parent.get_name())
                    continue
#                print "Obstacle"
                if radius in obstacles:
                    obstacles[radius].append(coords_i)
                else:
                    obstacles[radius]=[coords_i]

def get_all_anchor_coords():
    global fgs_to_anchor_coords
    ret=[]
    for anchor_coords in fgs_to_anchor_coords.itervalues():
        ret.extend(anchor_coords)
    return ret

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
    '''
    Coarse grain obstacles

    return a dictionary from radii to obstacles coordinates at that radii
    in the coarse grained representation (note the output dictionary is
    of defaultdict(list) type)
    '''
    # bin obstacles by distance from central axis and z-axis, such that
    # coarse-graining of each bin is an decreasing resolution from center outward
    # and increasing resolution for increasing Z values
    anchor_coords= get_all_anchor_coords()
    print(anchor_coords[0:3])
    nn_anchors= IMP.algebra.NearestNeighbor3D(anchor_coords)
    in_obstacles= defaultdict(list)
    for (radius, coords) in obstacles.iteritems():
        for coord in coords:
            nearest_i= nn_anchors.get_nearest_neighbor(coord)
            nearest_anchor= anchor_coords[nearest_i]
#            print(coord, nearest_anchor)
            distance_anchor= IMP.algebra.get_distance(coord, nearest_anchor)
            max_allowed_error_A= 3.0 * (max(distance_anchor-10.0,10.0)/10.0) ** 1.2
            #            dXY= math.sqrt(coord[0]**2+coord[1]**2)
            #            dXY_by_Z= dXY * (500-abs(coord[2]))/500.0
            #            max_allowed_error_A= 1.5 * (dXY_by_Z/150.0)**1.15
            max_allowed_error_A= min(50.0, max_allowed_error_A)
            max_allowed_error_A= math.ceil(max_allowed_error_A/2.0)*2.0
            s= IMP.algebra.Sphere3D(coord,radius)
            in_obstacles[max_allowed_error_A].append(s)
    intermediate_obstacles=[]
    for max_allowed_error_A, spheres in in_obstacles.iteritems():
        intermediate_obstacles=intermediate_obstacles+ \
                                IMP.algebra.get_simplified_from_volume(spheres, max_allowed_error_A)
    out_obstacles=defaultdict(list)
    for sphere in intermediate_obstacles:
        R=max(4.0, round(sphere.get_radius()))
        out_obstacles[R].append(sphere.get_center())
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
        is_gle1_or_nup42= re.match("Gle1", nup_name) or re.match("Nup42", nup_name)
        if is_single_spoke and is_gle1_or_nup42:
            print("Gle1/Nup42 detected in main spoke:", nup_name)
        is_remove= (is_gle1_or_nup42 and IS_REMOVE_GLE1_AND_NUP42)
        if is_single_spoke and not is_remove:
            for nup_node in nup_root.get_children():
                print("Handling", nup_node, "son of", nup_name)
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
    config.backbone_k.lower=0.0075 #kcal/mol/A^2
    config.is_backbone_harmonic=1
    config.backbone_tau_ns.lower=50.0
    config.time_step_factor.lower=2
    config.excluded_volume_k.lower=5
    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
    config.slack.lower = 10
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
args=parse_commandline()
print(args)
print(args.config_file, args.input_rmf_file)
test_is_node_in_fg_domain() # just to test it works correctly
config= get_basic_config()
# Add FGs:
fgs_regions_to_params= get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42)
add_obstacles_from_rmf(config, args.input_rmf_file, fgs_regions_to_params)
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
rrange= args.diffuser_sizes
for radius in rrange:
    inert_name="R%d" % radius
    nonspecifics[radius]= IMP.npctransport.add_float_type(config,
                                                  number=args.n_diffusers,
                                                  radius=radius,
                                                  type_name=inert_name,
                                                  interactions=0)
    kap_name="kap%d" % radius
    kaps[radius]= IMP.npctransport.add_float_type(config,
                                             number=args.n_diffusers,
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
                                             range_sigma1_deg= sigma1_deg)

# Add obstacles


# dump to file
f=open(args.config_file, "wb")
f.write(config.SerializeToString())
print(config)
