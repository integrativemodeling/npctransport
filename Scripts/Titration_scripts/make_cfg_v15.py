#!/usr/bin/python
# #bash /Users/barak/imp/fast/setup_environment.sh /anaconda/bin/python
#v10 - 2 interaction sites instead of 4
#v11 - custom insteraction site geometry specified in cfg file
#v12 - specify precise kap and FG concentrations + back to 4 interaction sites, n_kaps set to 16
#v13 - same as v12 but fixed number of n_fgs instead of n_kaps

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

import myutils

AVOGADRO= 6.022140857E23

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


if __name__ != "__main__":
    sys.exit(0)
# fetch params
# Usage: <cmd> <outfile>
parser = argparse.ArgumentParser(description="make configuration file for npctransport.")
parser.add_argument("outfile", help="config output file")
parser.add_argument("repeat_string", type=str, help="a string with a character per bead, with S for off sites and F for on sites. The FG chain has len(repeat_string) beads")
parser.add_argument("slide", type=is_true,
                    help="use a sliding interaction potential if equals True, T or a non-zero int (case insensitive)",
                    default=False)
parser.add_argument("-f", "--fg_concentration", type=float,
                    help="concentration of fg repeat chains in M", default=0.001)
parser.add_argument("-k","--kap_concentration", type=float,
                    help="concentration of kaps in M", default=0.005)
parser.add_argument("-t", "--temperature_scan", action="store_true",
                    help="create a range of temeprature values")
parser.add_argument("-s", "--time_step_factor", type=float,
                    help="simulation time step factor", default=1.0)
parser.add_argument("--kap_valency", type=int,
                    help="number of kap interaction sites", default=4)
parser.add_argument("--kap_k", type=float,
                    help="Kap interaction force coefficient in kCal/mol/A^2 (for slide) or kCal/mol/A (otherwise",
                    default=3.96)
parser.add_argument("--kap_range", type=float,
                    help="Kap interaction range in Angstroms",
                    default=5.5)
parser.add_argument("--min_fg_number", type=int,
                    help="Minimal number of FGs",
                    default=2)
parser.add_argument("--min_kap_number", type=int,
                    help="Minimal number of Kaps",
                    default=1)
args= parser.parse_args()
print(args)
[n_kaps, n_fg_chains, box_side_A]= \
    myutils.get_nkaps_nfgs_boxsideA_at_desired_C_M(
        args.kap_concentration,
        args.fg_concentration,
        args.min_kap_number,
        args.min_fg_number,
        verbose= True)
#print "[%s]" % outfile
SLIDE=args.slide # False
print "#Slide", SLIDE
if(SLIDE):
  kap_k=3.96
  kap_range=5.5
  range_sigma0_deg=45
  range_sigma1_deg=45
else:
  kap_k=2.95
  kap_range=4
  range_sigma0_deg=0
  range_sigma1_deg=0
kap_interaction_sites= args.kap_valency
# NOTE: site potential energy is kap_k*kap_range kCal/mol

FG_RES_PER_BEAD=20 #10
FG_RADIUS_PER_BEAD=8.0*(FG_RES_PER_BEAD/20.0) #6
FG_INTERACTIONS_PER_BEAD=1
FG_REST_LENGTH_FACTOR=1.9*math.sqrt(20.0/FG_RES_PER_BEAD) # make sure that radius*rest_length scale with sqrt(# beads)
FG_EXCLUDED_VOLUME_K= 10
fgs={}

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower=1.0
    #config.dump_interval=1
    config.interaction_k.lower=1e-9
    config.interaction_range.lower=10
    config.backbone_k.lower=0.0075 #kcal/mol/A^2
    config.is_backbone_harmonic=1
    config.backbone_tau_ns.lower=50.0
    config.time_step_factor.lower= args.time_step_factor
    config.excluded_volume_k.lower=10 # Original 5
    # non-specific attraction
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
#    create_range(config.nonspecific_k, 0.01, 0.25, steps=10, base=1)
    config.slack.lower = 10
    config.number_of_trials=1
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water
    config.statistics_interval_ns= 0.25
    config.output_statistics_interval_ns= 30.0
    ###
    #simulation bounding volumes:
    config.box_is_on.lower=1
#config.dump_interval_ns=0.01
#config.simulation_time_ns=5
    config.dump_interval_ns=1000
    config.simulation_time_ns=1000
    config.box_is_on.lower=1
    config.box_side.lower= box_side_A
    config.slab_is_on.lower= 0
    config.is_xyz_hist_stats=1
    return config



# ********* MAIN: *********
config= get_basic_config()
if(args.temperature_scan):
      IMP.npctransport.create_range(config.temperature_k, lb=275, ub=335, steps=21, base=1)
fg= IMP.npctransport.add_fg_type(config,
                                 type_name="fg0",
                                 number_of_beads=len(args.repeat_string),
                                 number=n_fg_chains,
                                 radius=FG_RADIUS_PER_BEAD,
                                 interactions=FG_INTERACTIONS_PER_BEAD,
                                 rest_length_factor = FG_REST_LENGTH_FACTOR)
for c in args.repeat_string:
  if c=='F':
    fg.type_suffix_list.append("") # non-interacting species
  elif c=='S':
    fg.type_suffix_list.append("_NI") # non-interacting species
  else:
      raise ValueError("invalid repeat string", args.repeat_string, "- should only contain S or F")

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
        interaction.nonspecific_k.lower= 0.01
#        create_range(interaction.interaction_k, 0.75*kap_k, 1.25*kap_k, steps=50, base=1)
#        create_range(interaction.interaction_range, 0.75*kap_range, 1.25*kap_range, steps=50, base=1)


# Add fg-fg interactions:
interactionFG_FG= IMP.npctransport.add_interaction(config,
    name0= fg.type,
    name1= fg.type,
    interaction_k= 1.32,
    interaction_range= 6.00 )
interactionFG_FG.excluded_volume_k.lower= FG_EXCLUDED_VOLUME_K
interactionFG_FG.nonspecific_k.lower= 0.01
#create_range(interactionFG_FG.interaction_k, 0.01, 2.0, steps=10, base=1)

# dump to file
f=open(args.outfile, "wb")
f.write(config.SerializeToString())
print config
