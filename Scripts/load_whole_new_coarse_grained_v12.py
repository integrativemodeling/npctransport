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
from FGParamsFactory import *
import my_util
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


IS_TOROID= True
IS_REMOVE_GLE1_AND_NUP42=True
#outfile = sys.argv[1]
#input_rmf = sys.argv[2] #'wholeNPC_0.rmf3')
COARSE_GRAINED_OBSTACLES= True
Z_TRANSFORM=0 #-75.5
FG_RES_PER_BEAD=20 #10
FG_RADIUS_PER_BEAD=8.0*(FG_RES_PER_BEAD/20.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.9*math.sqrt(20.0/FG_RES_PER_BEAD) # make sure that radius*rest_length scale with sqrt(# beads)
obstacles={}
fgs_to_anchor_coords={}
if IS_TOROID:
    TUNNEL_RADIUS=540 # R=540 for toroid from Nature 2018 (not 390 as stated, 390 is R-r)
    NUCLEAR_ENVELOPE_WIDTH= 300 # r=150 A, H=rx2, for toroid from Nature 2018
else:
    TUNNEL_RADIUS= 375
    NUCLEAR_ENVELOPE_WIDTH= 250
OBSTACLE_SCALE_FACTOR=1.0 # inflate obstacles a bit to prevent artificial cavities
kap_k=4.75
kap_range=4.5
sigma0_deg=45
sigma1_deg=45
kap_interaction_sites=4

def parse_commandline():
    parser = argparse.ArgumentParser \
             ( description='Create a config file loaded from an RMF model of the full NPC',
               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config_file', metavar='config_file', type=str,
                        default='config.pb',
                        help='output configuration file')
    parser.add_argument('input_rmf_file', metavar='input_rmf_file', type=str,
                        default='config.pb',
                        help='rmf file of NPC model (e.g., InputData/47-35_1spoke.rmf3)')
    parser.add_argument('--diffusers_radii', metavar='diffuser_radius', type=int, nargs='*',
                        default=range(14,29,2),
                        help='list of diffuser radii in A')
    parser.add_argument('--n_diffusers', type=int,
                        default=200,
                        help='number of diffuser molecules for each type of diffuser')
    parser.add_argument('--box_size', type=float,
                        default=2000,
                        help='size of simulation box edge in ansgtroms')
    parser.add_argument('--only_nup', type=str,
                        help='if specified, creates a box with only a single copy of given nup')
    args = parser.parse_args()
    return args


##
def get_surface_distance_from_axis_aligned_ellipsoid(xyz, sphere_radius,
                                                     origin,
                                                     rv, rh):
    '''
    return the distance between the surface of specified sphere and
    the surface of an axis aligned ellipsoid with vertical semi-axis rv and
    horizontal semi-axis rh, centered at origin.

    Returns:
    If the sphere penetrates the ellipsoid by some distance d, returns -d.

    xyz - center of sphere to be compared to ellipsoid
    sphere_radius - radius of sphere to be compared to ellipsoid
    origin - ellipsoid 3D origin
    rv - vertical ellipsoid radius (z axis)
    rh - horizontal ellipsoid radius (x and y axes)
    '''
    EPS=0.000001
    v= xyz-origin
    dXY2= v[0]**2+v[1]**2
    dZ2 = v[2]**2
    #    theta = atan(dXY/(dZ+EPS))
    dv2= dXY2 + dZ2 + EPS
    sinTheta2=dXY2/dv2
    cosTheta2=dZ2/dv2
    cur_r = math.sqrt(rv**2*cosTheta2 + rh**2*sinTheta2)
    dv=math.sqrt(dv2)
    return dv-cur_r-sphere_radius

#
def is_outside_slab_with_toroidal_pore(xyz, sphere_radius,
                                       R, r,
                                       allowed_overlap_A=0.0):
    '''
    verify that sphere (xyz, sphere_radius) is out of a slab (think
    membrane) with thickness 2xr, and toroidal pore with major radius
    R and minor radius r, whose central axis is along (0,0,z)

    xyz, spehre_radius - sphere coords
    R - major toroid radius
    allowed_overlap_A - penetration allowed that would still result in True. Use negatqive value for a conservative answer that requires additional slack between spheres and slab surface to be considered outside.
    '''
    EPS=0.000001
    if xyz[2] - sphere_radius + allowed_overlap_A > r: # above slab
        return True
    if xyz[2] + sphere_radius - allowed_overlap_A < -r: # below slab
        return True
    # In slab vertically (but still might be in pore = outside slab):
    xy0= IMP.algebra.Vector3D(xyz[0], xyz[1], 0.0)
    rxy0= xy0.get_magnitude()
    if rxy0 + sphere_radius > R:
        return False # in slab outside external radius
    if(rxy0>EPS):
        xy0_major= xy0 * ( R / rxy0 ) # xy0_major is a vector of magnitude R in the same direction as xy0
    else:
        xy0_major=IMP.algebra.Vector3D(0, R, 0.0)
    distance = get_surface_distance_from_axis_aligned_ellipsoid \
               (xyz, sphere_radius, xy0_major, r, r)
    return distance + allowed_overlap_A > 0

def is_outside_slab_with_cylindrical_pore(xyz, sphere_radius,
                                          R, r, allowed_overlap_A):
    '''
    Similar to is_outside_slab_with_toroidal_pore but for a cylinder of radius R
    in a slab of height r*2, whose central axis is along (0,0,z)

    xyz, spehre_radius - sphere coords
    R - cylinder radius
    r - half of slab thickness (along Z)
    allowed_overlap_A - penetration allowed that would still result in True. Use negat
    '''
    distance_xy=math.sqrt(xyz[0]**2 + xyz[1]**2)
    in_pore_interior= distance_xy - sphere_radius - allowed_overlap_A <= R
    outside_slab= abs(xyz[2]) + radius + allowed_overlap_A > r
    return outside_slab or in_pore_interior


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
#    SCALE_SELF_K= 0.88 # not for GLFGs for now
#    SCALE_SELF_K_GLFG= 0.96
    default_fgp=FGParamsFactory \
                 ( res_from=nan,
                   res_to= nan,
                   self_k= 1.28,
                   self_range= 6.00,
                   kap_k= 3.58,
                   kap_range= 4.95,
                   nonspec_k= 0.01,
                   nonspec_range= 5.00,
                   backbone_k= 0.0075,
                   backbone_tau= 50 )
    default_disordered_fgp= default_fgp.get_copy()
    default_disordered_fgp.kap_k= None
    default_disordered_fgp.kap_range= None
    default_FSFG= default_fgp.get_copy()
    default_FSFG.self_k= 1.28
    default_FSFG.nonspec_k= 0.01
    default_GLFG= default_fgp.get_copy()
    default_GLFG.self_k= 1.41
    default_GLFG.nonspec_k= 0.08
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
    fgs_regions_to_params['FSFG_generic']= {
        '': default_FSFG.get_copy(1, 120),
    }
    fgs_regions_to_params['GLFG_generic']= {
        '': default_GLFG.get_copy(1, 120),
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
    IS_REVERSE= False
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
            nbeads = int(math.ceil((region_nres_to_add+0.0) / FG_RES_PER_BEAD))
            for i in range(nbeads):
                bead_suffixes.append(region)
            total_nres= total_nres + region_nres
    return bead_suffixes, is_anchor

def add_fgs(config, type_name, fg_params,
            apply_anchor= True,
            anchor_coordinates= None):
    ''' add fgs of type "type_name" to config, with length and subtypes names
        of each bead set according to fg_params. Use global variable
        FG_RES_PER_BEAD to decide number of residues per bead

        config - protobuf config object to which fgs will be added
        type_name - label of FGs type name
        fg_params - a dictionary from type_name to an internal dictionary
                    of subdomain suffixes (e.g. "C" for "Nsp1C" to FGParam objects,
                    such that each FGParam object describes the residue range of the subdomains.
                    For instance fg_params["Nsp1"]["C"] describes the subdomain "Nsp1C",
                    between fg_params["Nsp1"]["C"].res_from and fg_params["Nsp1"]["C"].res_to
        apply_anchor - if true, anchor the i'th FG chain to anchor_coordinates[i], using the
                    bead fg_params[type_name]["anchor"], which must be at either chain terminal
        anchor_coordinates - a list of anchor coordinates for each FG chain (or None)

    '''
    global FG_RES_PER_BEAD
    bead_suffixes, is_anchor \
        = get_bead_suffixes_and_is_anchor(fg_params, FG_RES_PER_BEAD)
    for fg_param in fg_params.itervalues():
        print(fg_param.backbone_k, config.backbone_k.lower)
        assert(fg_param.backbone_k in [None,config.backbone_k.lower])
        assert(fg_param.backbone_tau in [None, config.backbone_tau_ns.lower])
    nbeads = len(bead_suffixes)
    if apply_anchor:
        number_fgs= len(anchor_coordinates)
    else:
        number_fgs= 1
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= nbeads,
                                      number=number_fgs,
                                      radius=FG_RADIUS_PER_BEAD,
                                      interactions= FG_INTERACTIONS_PER_BEAD,
                                      rest_length_factor = FG_REST_LENGTH_FACTOR)
    if 'anchor' in bead_suffixes and apply_anchor:
        assert(bead_suffixes[0]=='anchor' and 'anchor' not in bead_suffixes[1:])
        for coordinates in anchor_coordinates:
#            print(type_name, "Coords:", coordinates)
            pos =fgs.anchor_coordinates.add()
            pos.x= coordinates[0]
            pos.y= coordinates[1]
            pos.z= coordinates[2] + Z_TRANSFORM
    fgs.type_suffix_list.extend(bead_suffixes)
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
            R= TUNNEL_RADIUS
            r= 0.5*NUCLEAR_ENVELOPE_WIDTH
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
                # print "Obstacle"
                # Filter out obstacles that overlaps nuclear envelope
                if (IS_TOROID and not is_outside_slab_with_toroidal_pore(coords_i, radius, R, r)) or \
                   (not IS_TOROID and not is_outside_slab_with_cylindrical_pore(coords_i, radius, R, r)):
#                    print("# Removing obstacle node that is in nuclear envelope",
#                          child.get_name(), parent.get_name())
                    continue
#                print("# Adding obstacle node",
#                      child.get_name(), parent.get_name())
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
#        print("Res:1", nup_name, ";", node_name)
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

def get_basic_config(cmdline_args):
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
    config.box_side.lower= cmdline_args.box_size
    if IS_TOROID:
        config.slab_is_on.lower=2
    else:
        config.slab_is_on.lower=1 # cylinder
    config.slab_thickness.lower=NUCLEAR_ENVELOPE_WIDTH
    config.tunnel_radius.lower=TUNNEL_RADIUS
    config.is_xyz_hist_stats=1
    return config

def add_kaps_and_inerts(config,
                        n_diffusers,
                        diffusers_radii,
                        kap_interaction_sites,
                        fg_regions_to_params):
    '''
    For each radius in the list diffusers_radii, add
    n_diffusers x kaps and n_diffusers x inerts to config,
    with kap_interaction_sites per kap molecule, and using
    fg_regions_to_params to add interactions with fgs
    '''
    nonspecifics={}
    kaps={}
    SPECIAL_HACK= True
    for radius in diffusers_radii:
        inert_name="inert%d" % radius
        nonspecifics[radius]= IMP.npctransport.add_float_type(config,
                                                              number=n_diffusers,
                                                              radius=radius,
                                                              type_name=inert_name,
                                                              interactions=0)
        # TODO: for now, this is very ad hoc, and should be generalized
        #       also, same binding sites for FGs and kaps, need to think how to handle
        if radius>=80 and SPECIAL_HACK:
            assert(40 in kaps)
            n_interactions= int(math.ceil(4*radius/40.0))
            nonspecifics[radius].interactions.lower= n_interactions
            kaps[40].number.lower= args.n_diffusers * n_interactions
            kaps[40].interactions.lower= kaps[40].interactions.lower + 1
            interaction= IMP.npctransport.add_interaction \
                          ( config,
                            name0= "kap40",
                            name1= inert_name,
                            interaction_k= 10.0,
                            interaction_range= 30.0)
            interaction.active_sites0.append(0)
            continue

        kap_name="kap%d" % radius
        kaps[radius]= IMP.npctransport.add_float_type(config,
                                                      number=args.n_diffusers,
                                                      radius=radius,
                                                      type_name=kap_name,
                                                      interactions=kap_interaction_sites)
        # Add interactions:
        for fg_name, regions_to_params in fgs_regions_to_params.iteritems():
            for region, params in regions_to_params.iteritems():
                # IMP.npctransport.add_interaction(config, name0=fg_name
                #                                  + region,
                #                                  name1=inert_name,
                #                                  interaction_k=0,
                #                                  interaction_range=0)
                interaction= IMP.npctransport.add_interaction \
                             ( config,
                               name0=fg_name + region,
                               name1=kap_name,
                               interaction_k= params.kap_k,
                               interaction_range= params.kap_range,
                               range_sigma0_deg= sigma0_deg,
                               range_sigma1_deg= sigma1_deg)
                if(params.nonspec_k is not None):
                    interaction.nonspecific_k.lower= params.nonspec_k
                if kap_name=="kap40" and SPECIAL_HACK:
                    for site_id in range_inclusive(1, kap_interaction_sites):
                        interaction.active_sites1.append(site_id)



# *************************
# ********* MAIN: *********
# *************************
IS_SCAFFOLD_AND_PORE= True
args=parse_commandline()
print(args)
print(args.config_file, args.input_rmf_file)
test_is_node_in_fg_domain() # just to test it works correctly
config= get_basic_config(args)
# Add FGs:
fgs_regions_to_params= get_fgs_regions_to_params(IS_REMOVE_GLE1_AND_NUP42)
if(args.only_nup is not None):
    fgs_regions_to_params= {args.only_nup : fgs_regions_to_params[args.only_nup]}# discard all other nups
    fgs_to_anchor_coords= {args.only_nup : None }
    config.slab_is_on.lower= 0
    IS_SCAFFOLD_AND_PORE=False
if IS_SCAFFOLD_AND_PORE:
    add_obstacles_from_rmf(config, args.input_rmf_file, fgs_regions_to_params)
#print("FGs to anchor keys", fgs_to_anchor_coords.keys())
for (name,coords) in fgs_regions_to_params.iteritems(): # TODO: get coords right
    apply_anchor= IS_SCAFFOLD_AND_PORE
    if apply_anchor:
        if name not in fgs_to_anchor_coords:
            continue
        anchor_coordinates= fgs_to_anchor_coords[name]
    else:
        anchor_coordinates= None
    add_fgs(config,
            name,
            fgs_regions_to_params[name],
            apply_anchor= apply_anchor,
            anchor_coordinates= anchor_coordinates)
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
                if(params0.nonspec_k is None or params1.nonspec_k is None):
                    nonspec_k=0.01
                else:
#                    nonspec_k= 0.5*(params0.nonspec_k+params1.nonspec_k)
                    nonspec_k= math.sqrt(params0.nonspec_k*params1.nonspec_k) # geometric mean
                interactionFG_FG.nonspecific_k.lower= nonspec_k
# Add kaps and inerts:
if not args.only_nup:
    add_kaps_and_inerts(config,
                        args.n_diffusers,
                        args.diffusers_radii,
                        kap_interaction_sites,
                        fgs_regions_to_params)


# dump to file
f=open(args.config_file, "wb")
f.write(config.SerializeToString())
print(config)
