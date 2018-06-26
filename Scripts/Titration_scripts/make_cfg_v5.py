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
parser.add_argument("-f", "--max_fg_concentration", type=float,
                    help="maximal concentration of fg repeat chains in M", default=0.001)
parser.add_argument("-k","--max_kap_concentration", type=float,
                    help="maximal concentration of kaps in M", default=0.005)
parser.add_argument("-r","--ratio_kaps_to_fg", type=float,
                    help="approximate molar ratio of kap to fgs", default=1.0)
parser.add_argument("-t", "--temperature_scan", action="store_true",
                    help="create a range of temeprature values")
parser.add_argument("-s", "--time_step_factor", type=float,
                    help="simulation time step factor", default=1.0)
parser.add_argument("-n", "--nbeads", type=int,
                    help="number of beads in FG chains", default=6)
args= parser.parse_args()
AVOGADRO= 6.022140857E23
ratio= args.ratio_kaps_to_fg
n_fg_chains= 60
ratio_actual= math.ceil(ratio*n_fg_chains)/n_fg_chains if ratio>0.0 else 0.0
n_kaps= int(n_fg_chains*ratio_actual)
for modulus in [2,3,5,2,3,5,2,3,5]:
      while(n_kaps%modulus==0
            and n_fg_chains%modulus==0
            and max(n_kaps/modulus,n_fg_chains/modulus)>=5
            and (min(n_kaps/modulus,n_fg_chains/modulus)>=1 or (ratio_actual==0.0 and n_fg_chains/modulus>=1)) ):
         print "Dividing",n_kaps,n_fg_chains,"by",modulus
         n_kaps= n_kaps/modulus
         n_fg_chains= n_fg_chains/modulus
print "n_fg_chains", n_fg_chains, "n_kaps", n_kaps
n_max= max(n_fg_chains,n_kaps)
min_box_volume_L_fgs= n_fg_chains/AVOGADRO/args.max_fg_concentration
min_box_volume_L_kaps= n_kaps/AVOGADRO/args.max_kap_concentration
box_volume_L= max(min_box_volume_L_fgs, min_box_volume_L_kaps)
box_volume_A3= box_volume_L*1E27
BOX_SIDE_A= box_volume_A3**(1/3.0)
FG_CONCENTRATION= n_fg_chains/AVOGADRO/box_volume_L
KAP_CONCENTRATION= n_kaps/AVOGADRO/box_volume_L
print("FG nups concaentration %.3e [M]" % FG_CONCENTRATION)
print("Kap concaentration %.3e [M]" % KAP_CONCENTRATION)
print("Molar ratio Kap:FG chains: %.3f" % ratio_actual)
#print "[%s]" % outfile
SLIDE=args.slide # False
print "Slide", SLIDE
if(SLIDE):
    kap_k=3.25
    kap_range=4.5
    range_sigma0_deg=45
    range_sigma1_deg=45
else:
    kap_k=2.95
    kap_range=4
    range_sigma0_deg=0
    range_sigma1_deg=0
kap_interaction_sites=4
# NOTE: site potential energy is kap_k*kap_range kCal/mol

FG_BEADS_PER_RES=20 #10
FG_RADIUS_PER_BEAD=8*(FG_BEADS_PER_RES/20.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.9*math.sqrt(20.0/FG_BEADS_PER_RES) # make sure that radius*rest_length scale with sqrt(# beads), unless we assume strong cohesiveness
#FG_REST_LENGTH_FACTOR=1.75
fgs={}

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=10
    config.interaction_range.lower=10
    config.backbone_k.lower=0.05 # kcal/mol/A^2
    config.is_backbone_harmonic=1
    config.backbone_tau_ns.lower=50.0
    config.time_step_factor.lower= args.time_step_factor
    config.excluded_volume_k.lower=10 # NO-DEBUG
    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
    create_range(config.nonspecific_k, 0.01, 0.4, steps=15, base=1)
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
    config.box_side.lower= BOX_SIDE_A
    config.slab_is_on.lower= 0
    return config



# ********* MAIN: *********
config= get_basic_config()
if(args.temperature_scan):
      IMP.npctransport.create_range(config.temperature_k, lb=275, ub=335, steps=21, base=1)
fg= IMP.npctransport.add_fg_type(config,
                                 type_name="fg0",
                                 number_of_beads=args.nbeads,
                                 number=n_fg_chains,
                                 radius=FG_RADIUS_PER_BEAD,
                                 interactions=FG_INTERACTIONS_PER_BEAD,
                                 rest_length_factor = FG_REST_LENGTH_FACTOR)
for i in range(0, args.non):
    fg.type_suffix_list.append("") # non-interacting species
for i in range(0, args.nbeads-args.non):
    fg.type_suffix_list.append("_NI") # non-interacting species

# Add kaps / inerts
nonspecifics={}
kaps={}
rrange=[20]
for i in rrange:
    kap_name="kap%d" % i
    if n_kaps>0:
        kaps[i]= IMP.npctransport.add_float_type(config,
            number=n_kaps,
            radius=i,
            type_name=kap_name,
            interactions=kap_interaction_sites)
        interaction=IMP.npctransport.add_interaction(config,
            name0="fg0",
            name1=kap_name,
            interaction_k=kap_k,
            interaction_range=kap_range,
            is_on=1,
            range_sigma0_deg=range_sigma0_deg,
            range_sigma1_deg=range_sigma1_deg
        )
#        create_range(interaction.interaction_k, 0.9*kap_k, 1.2*kap_k, steps=5, base=1)
#        create_range(interaction.interaction_range, 0.9*kap_range, 1.1*kap_range, steps=5, base=1)


# Add fg-fg interactions:
interactionFG_FG= IMP.npctransport.add_interaction(config,
    name0= fg.type,
    name1= fg.type,
    interaction_k= 0.01,
    interaction_range= 6)
create_range(interactionFG_FG.interaction_k, 0.01, 3.0, steps=15, base=1)

# dump to file
f=open(args.outfile, "wb")
f.write(config.SerializeToString())
print config
