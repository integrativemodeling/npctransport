#!/usr/bin/bash /Users/barak/imp/fast/setup_environment.sh /anaconda/bin/python

#########################################
# Example: imppy make_cfg.py config.pb > config.txt
#
# Author: Barak Raveh barak@salilab.org
# Data written: May 5, 2017
# Date last udpated: May 5, 2017 or later
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
import argparse

def usage_and_exit(exit_code):
      print sys.argv[0], "<outfile> <n_on> <is_slide> [-b/--boxsize <boxsize>]"
      sys.exit(exit_code)

def is_true(s):
      '''
      returns true if s equals "true", "t" or a non-zero numeric string.
      Case insensitive.
      '''
      s=s.lower()
      try:
            s=int(s)
            return s<>0
      except:
            pass
      return s=="true" or s=="t"

# fetch params
# Usage: <cmd> <outfile>
parser = argparse.ArgumentParser(description="make configuration file for npctransport.")
parser.add_argument("outfile", help="config output file")
parser.add_argument("non", type=int, help="number of FG sites active on chain")
#slide_group = parser.add_mutually_exclusive_group()
parser.add_argument("slide", type=is_true,
                    help="use a sliding interaction potential if equals True, T or a non-zero int (case insensitive)",
                    default=False)
parser.add_argument("-b","--boxscale", type=float,
                   help="scaling factor for box volume", default=1.0)
parser.add_argument("-t", "--temperature_scan", action="store_true",
                    help="create a range of temeprature values")
parser.add_argument("-f", "--time_step_factor", type=float,
                    help="simulation time step factor", default=1.0)
parser.add_argument("-n", "--nbeads", type=int,
                    help="number of beads in FG chains", default=6)
args = parser.parse_args()
#print "[%s]" % outfile
SLIDE=args.slide # False
print "Slide", SLIDE
if(SLIDE):
    kap_k=3.25
    kap_range=4.5
    range_sigma0_deg=45
    range_sigma1_deg=45
else:
    kap_k=2.75
    kap_range=4
    range_sigma0_deg=0
    range_sigma1_deg=0
kap_interaction_sites=4
# NOTE: site potential energy is kap_k*kap_range kCal/mol
FG_BEADS_PER_RES=20 #10
FG_RADIUS_PER_BEAD=6.5*(FG_BEADS_PER_RES/15.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=2.0*math.sqrt(15.0/FG_BEADS_PER_RES) # make sure that radius*rest_length scale with sqrt(# beads), unless we assume strong cohesiveness
#FG_REST_LENGTH_FACTOR=1.75
fgs={}

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=10
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=0.05 # kcal/mol/A^2
    config.is_backbone_harmonic=1
    config.backbone_tau_ns.lower=50.0
#    create_range(config.backbone_tau_ns, 50, 250, steps=3, base=1)
    config.time_step_factor.lower= args.time_step_factor
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=10 # NO-DEBUG
    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
    create_range(config.nonspecific_k, 0.01, 0.4, steps=5, base=1)
    config.slack.lower = 30
    config.number_of_trials= 1
    config.angular_D_factor.lower= 0.3 #increased dynamic viscosity relative to water
    config.statistics_interval_ns= 0.25
    config.output_statistics_interval_ns= 30.0
    ###
    #simulation bounding volumes:
    config.box_is_on.lower= 1
#config.dump_interval_ns=0.01
#config.simulation_time_ns=5
    config.dump_interval_ns= 1000
    config.simulation_time_ns= 1000
    config.box_is_on.lower= 1
    config.box_side.lower=400 * args.boxscale**(1/3.0)
    config.slab_is_on.lower= 0
    return config

def add_fgs(config, type_name, nres, anchor_coordinates=[None]):
    nbeads = int(math.ceil(nres / FG_BEADS_PER_RES))
    if anchor_coordinates[0]<>None:
       nbead= nbead+1 # for anchor
    fgs= IMP.npctransport.add_fg_type(config,
                                      type_name= type_name,
                                      number_of_beads= nbeads,
                                      number=len(anchor_coordinates),
                                      radius=FG_RADIUS_PER_BEAD,
                                      interactions= FG_INTERACTIONS_PER_BEAD,
                                      rest_length_factor = FG_REST_LENGTH_FACTOR)
    for coordinates in anchor_coordinates:
        if coordinates is None:
           continue
        pos =fgs.anchor_coordinates.add()
        pos.x= coordinates[0]
        pos.y= coordinates[1]
        pos.z= coordinates[2] + Z_TRANSFORM
    return fgs





# ********* MAIN: *********
config= get_basic_config()
if(args.temperature_scan):
      IMP.npctransport.create_range(config.temperature_k, lb=275, ub=335, steps=21, base=1)
fgs=[]
fgs.append( add_fgs(config, "Nsp1", args.nbeads * FG_BEADS_PER_RES) )
for i in range(0, args.non):
    fgs[-1].type_suffix_list.append("") # non-interacting species
for i in range(0, args.nbeads-args.non):
    fgs[-1].type_suffix_list.append("_NI") # non-interacting species

# Add kaps / inerts
nonspecifics={}
kaps={}
rrange=[20]
for i in rrange:
    kap_name="kap%d" % i
    kaps[i]= IMP.npctransport.add_float_type(config,
                                             number=1,
                                             radius=i,
                                             type_name=kap_name,
                                             interactions=kap_interaction_sites)
#   create_range(kaps[i].number, lb=1, ub=600, steps=40, base=1.1)
    interaction=IMP.npctransport.add_interaction(config,
        name0="fg0",
        name1=kap_name,
        interaction_k=kap_k,
        interaction_range=kap_range,
        is_on=1,
        range_sigma0_deg=range_sigma0_deg,
        range_sigma1_deg=range_sigma1_deg
    )
    create_range(interaction.interaction_k, 0.8*kap_k, 1.2*kap_k, steps=10)
    create_range(interaction.interaction_range, 0.5*kap_range, 1.5*kap_range, steps=10)


# Add fg-fg interactions:
for fg0 in fgs:
    for fg1 in fgs:
        interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                           name0= fg0.type,
                                                           name1= fg1.type,
                                                           interaction_k= 0.01,
                                                           interaction_range= 6)

        create_range(interactionFG_FG.interaction_k, 0.01, 3.0, steps=5, base=1)
# dump to file
f=open(args.outfile, "wb")
f.write(config.SerializeToString())
print config
