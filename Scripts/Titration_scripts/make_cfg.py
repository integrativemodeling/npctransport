#!/usr/bin/python
from IMP.npctransport import *
import sys
import math
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
                   help='scaling factor for box folume', default=1.0)
parser.add_argument("-t", "--temperature_scan", action="store_true",
                    help="create a range of temeprature values")
parser.add_argument("-n", "--nbeads", type=int,
                    help="number of beads in FG chains", default=6)
args = parser.parse_args()
#print "[%s]" % outfile
SLIDE=args.slide # False
print "Slide", SLIDE
BOX_SCALE=args.boxscale
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

def get_basic_config():
    config = Configuration()
    IMP.npctransport.set_default_configuration(config)
    # Box and slab:
    config.box_is_on.lower=1
    config.box_side.lower=400 * BOX_SCALE**(1/3.0)
    config.slab_is_on.lower=0
    # config.slab_thickness.lower=300
    # config.tunnel_radius.lower=190
    config.statistics_fraction.lower=1.0
    config.interaction_k.lower=0.001
    config.interaction_range.lower=5
    # create_range(config.backbone_k, .2, 1, 10
    config.backbone_k.lower=0.2
    config.time_step_factor.lower=4.0
#    create_range(config.time_step_factor, .5, 8, 5, 2)
    #create_range(config.rest_length_factor, .5, 1, 10)
    config.excluded_volume_k.lower=5
    # Non-specific attraction:
    config.nonspecific_range.lower= 5.0
    config.nonspecific_k.lower= 0.01
    config.slack.lower = 8
    config.number_of_trials=1
    config.angular_D_factor.lower=0.3 #increased dynamic viscosity relative to water
    ###
    # Simulation stat and output rates and times:
    config.dump_interval_ns=1000.0
    config.output_statistics_interval_ns=1000.0
    config.statistics_interval_ns=0.1 # note events at lower time resolution will be missed
    config.simulation_time_ns=10000
    return config



# ********* MAIN: *********
config= get_basic_config()
if(args.temperature_scan):
      IMP.npctransport.create_range(config.temperature_k, lb=275, ub=335, steps=21, base=1)
#config.simulation_time_ns=5

# fg_cyto= IMP.npctransport.add_fg_type(config,
#                                  number_of_beads=12,
#                                  number=8,
#                                  radius=6,
#                                  interactions=1,
#                                  rest_length_factor = 1.5)
nres_per_bead= 20
fg_nres= args.nbeads*nres_per_bead
fg_radius= 9.5*math.sqrt(nres_per_bead/20.0)
fg_middle= IMP.npctransport.add_fg_type(config,
                                 type_name="fg0",
                                 number_of_beads=args.nbeads,
                                 number=200,
                                 radius=fg_radius,
                                 interactions=1,
                                 rest_length_factor = 1.5)
for i in range(0, args.non):
    fg_middle.type_suffix_list.append("") # non-interacting species
for i in range(0, args.nbeads-args.non):
    fg_middle.type_suffix_list.append("_NI") # non-interacting species
#create_range(fg_middle.number_of_beads,lb=1,ub=6,steps=6,base=1)
# fg_nuclear= IMP.npctransport.add_fg_type(config,
#                                  number_of_beads=12,
#                                  number=8,
#                                  radius=6,
#                                  interactions=1,
#                                  rest_length_factor = 1.5)
#surface_area_ratio_to_R25 = math.pow( floaters_R / 25.0, 2 )
#kap_n_interactions = int( math.ceil ( 12 * surface_area_ratio_to_R25 ) )
#create_range(kaps.radius, lb = 10, ub = 30, steps = 5, base = 1)
nonspecifics={}
kaps={}
rrange=range(20,21,2)
for i in rrange:
    kap_name="kap%d" % i
    kaps[i]= IMP.npctransport.add_float_type(config,
                                             number=8,
                                             radius=i,
                                             type_name=kap_name,
                                             interactions=kap_interaction_sites)
    create_range(kaps[i].number, lb=1, ub=600, steps=40, base=1.1)
    IMP.npctransport.add_interaction(config,
                                     name0="fg0",
                                     name1=kap_name,
                                     interaction_k=kap_k,
                                     interaction_range=kap_range,
                                     is_on=1,
                                     range_sigma0_deg=range_sigma0_deg,
                                     range_sigma1_deg=range_sigma1_deg
                                     )


interactionFG_FG= IMP.npctransport.add_interaction(config,
                                                   name0= "fg0",
                                                   name1= "fg0",
                                                   interaction_k=0.01,
                                                   interaction_range= 6)
#create_range(interactionFG_FG.interaction_k, lb = 0.01, ub = 10, steps = 7, base = 2)
## internal FG-FG
#for i in range(3):
#    for j in range(i,3):
#        interactionFG_FG= IMP.npctransport.add_interaction(config,
#                                                           name0= "fg%d" % i,
#                                                           name1= "fg%d" % j,
#                                                           interaction_k= float(fgs_k),
#                                                           interaction_range= 2)

# dump to file
f=open(args.outfile, "wb")
f.write(config.SerializeToString())
print config
