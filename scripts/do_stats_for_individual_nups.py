#!/bin/python
from __future__ import print_function
#!/usr/bin/env python
from IMP.npctransport import *
import math
import stats_util
import my_util
import subprocess
import sys
import numpy as np
import os

imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp python"
npc_show_stats="$imppy $HOME/imp_git/repository/modules/npctransport/utility/show_statistics.py"

def handle_fg(fg, skip_n_rows, Ns=0):
  config_fname="config_{fg}.txt".format(fg=fg)
  if(not os.path.exists(config_fname)):
    raise(ValueError("{} does not exist - skipping {}".format(config_fname, fg)))
  cmd_n_fgs= "cat {} | ".format(config_fname) + \
        " awk 'BEGIN{is=0}{if(is==1){print $NF; is=2; }" + \
        " else { if($1 ~ /number_of_b/){is=1}}}'"
  n_fgs= int( subprocess.check_output(cmd_n_fgs, shell=True) )
#  Rbead=8.0
  output= Output()
  fname= "output_{fg}.pb".format(fg=fg+Ns*'N')
  with open(fname, "rb") as f:
    output.ParseFromString(f.read())
    Rbead=output.assignment.fgs[0].radius.value
    stats= output.statistics
    Rgs= []
    Rgs_normalized= []
    Rg2s= []
    e2es= []
    e2es_normalized= []
    e2e2s= []
    bond_ds= []
    bond_d2s= []
    assert(len(stats.fgs)==1)
    # Below, correcting Rg from coarse grain to compare directly to SAXS data,
    # using Rg_real^2 = Rg_CG^2 + 0.8*Rbead^2, where 0.8Rbead^2 is the formula
    # for a sphere's Rg
    for op in stats.fgs[0].order_params[skip_n_rows:]:
      Rg= op.mean_radius_of_gyration
      Rg_normalized= math.sqrt(Rg**2 + 0.8*Rbead**2)
      Rg2= op.mean_square_radius_of_gyration
      Rgs.append(Rg)
      Rgs_normalized.append(Rg_normalized)
      Rg2s.append(Rg2)
      e2e= op.mean_end_to_end_distance
      e2e_normalized= e2e + Rbead
      e2e2= op.mean_square_end_to_end_distance
      e2es.append(e2e)
      e2es_normalized.append(e2e_normalized)
      e2e2s.append(e2e2)
      bond_d= op.mean_bond_distance
      bond_d2= op.mean_square_bond_distance
      bond_ds.append(bond_d)
      bond_d2s.append(bond_d2)
    Rg_mean, Rg_sem= stats_util.get_time_series_mean_and_sem(Rgs_normalized)
    Rg_stddev= math.sqrt(np.mean(Rg2s) - np.mean(Rgs)**2)
    e2e_mean, e2e_sem= stats_util.get_time_series_mean_and_sem(e2es_normalized)
    e2e_stddev= math.sqrt(np.mean(e2e2s) - np.mean(e2es)**2)
    bond_d_mean, bond_d_sem= stats_util.get_time_series_mean_and_sem(bond_ds)
    bond_d_stddev= math.sqrt(np.mean(bond_d2s) - bond_d_mean**2)
    # Rg:
    print("{} n_fgs={:d} n_res={:d} ".format(fg,n_fgs,n_fgs*20), end="\t")
    print("Rg  {mean:.2f} +- {conf95:.2f}\tstd-dev {stddev:.2f}\tn= {n}"\
          .format(mean=Rg_mean, conf95=1.96*Rg_sem, stddev=Rg_stddev, n=len(Rgs)),
          end="\t")
    print("For nres=120 at Flory=0.5: {:.2f}".format(Rg_mean*(120/20.0/n_fgs)**0.5),
          end="\t")
    print("Flory=0.4: {:.2f}".format(Rg_mean*(120/20.0/n_fgs)**0.4),
          end="\t")
    print()
    # End to end:
    print("{} n_fgs={:d} n_res={:d} ".format(fg,n_fgs,n_fgs*20), end="\t")
    print("e2e {mean:.2f} +- {conf95:.2f}\tstd-dev {stddev:.2f}\tn= {n}"\
          .format(mean=e2e_mean, conf95=1.96*e2e_sem, stddev=e2e_stddev, n=len(e2es)),
          end="\t")
    print("For nres=120 at Flory=0.5: {:.2f}".format(e2e_mean*(120/20.0/n_fgs)**0.5),
          end="\t")
    print("Flory=0.4: {:.2f}".format(e2e_mean*(120/20.0/n_fgs)**0.4),
          end="\t")
    print()
    # Bond length:
    print("{} n_fgs={:d} n_res={:d} ".format(fg,n_fgs,n_fgs*20), end="\t")
    print("bond-length  {mean:.2f} +- {conf95:.2f}\tstd-dev {stddev:.2f}\tn= {n}"\
          .format(mean=bond_d_mean, conf95=1.96*bond_d_sem, stddev=bond_d_stddev, n=len(bond_ds)),
          end="\t")
    print()
    # Space:
    print()
    f.close()


if(len(sys.argv)>1):
  skip_n_rows= int(sys.argv[1])
else:
  skip_n_rows=100
if(len(sys.argv)>2):
  FGs= sys.argv[2:]
else:
  FGs= ["FSFG_generic", "GLFG_generic", "Nup1", "Nup100", "Nup49", "Nsp1"]
for fg in FGs:
  for Ns in my_util.range_inclusive(3,0,-1):
    try:
      handle_fg(fg, skip_n_rows, Ns=Ns)
      print ("Handled {}{} successfully".format(fg, Ns*'N'))
      break
    except:
      print("Skipping {}{}".format(fg, Ns*'N'))
  print()

  #   # II. End-to-end
  #   echo -n "$fg n_fgs=$n_fgs n_res=$((n_fgs*20)) "
  #   $npc_show_stats output_$fg.pb | grep mean_end_to_end | awk '{if(NR>'"$SKIP_ROWS"'){Rbead='"$Rbead"'; e2e=$2+Rbead; s=s+e2e; s2=s2+e2e*e2e; n=n+1}}END{e=s/n; e2=s2/n; var=e2-e*e; sem=sqrt(var/(n-1)); print e, "+-", 1.96*sem, "[95% conf]  std-dev", sqrt(var), "  n =", n, " For nres=120 at Flory=0.5: ", e*sqrt(120/20.0/'"$n_fgs"'), " at Flory=0.4: ", e*(120/20.0/'"$n_fgs"')^0.4}';
  #   echo
