#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import pickle
import numpy as np

CACHE_FNAME='TMP.d_stats.{:s}.{:s}.p'.format(sys.argv[1].replace("/","_").replace(".",""), 
                                             sys.argv[2])
# confidence interval required
CONF_INTERVAL=0.95
print(sys.argv)
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds
STATS_FROM_NS=STATS_FROM_SEC*1e9

def accumulate(dictionary, key, value, weight=1.0):
    entry=np.array([weight,value*weight,value**2*weight])
    if key in dictionary:
        dictionary[key] = dictionary[key] + entry
    else:
        dictionary[key] = entry

def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def print_stats(table, total_sim_time_sec):
    # print transports per particle per second
    print("Total simulation time: {:.6f} [sec]".format(total_sim_time_sec))
    keys=sorted(table.keys())
    for key in keys:
        #    mean_per_sec=table[key]/n/sim_time_sec;
        n= table[key][0]
        mean= table[key][1]/n
        mean2= table[key][2]/n
        var=mean2-mean**2
        try:
            stddev=math.sqrt(var)
            conf=get_sems_for_conf(CONF_INTERVAL)*stddev/math.sqrt(n)
        except:
            print("# Warning: Can't compute std-dev, probably always same value")
            conf=0.0
        print('{:} {:8.3e} +- {:8.3e} sec'.format(key, mean, conf))

############# Main ############
fnames=set(glob.glob(sys.argv[1]+"*.pb"))
try:
    with open(CACHE_FNAME,'rb') as f:
        [ table,
          n,
          total_sim_time_sec,
          processed_fnames,
          cached_stats_from_sec
          ] = pickle.load(f)
    assert(fnames.issuperset(processed_fnames))
    assert(cached_stats_from_sec == STATS_FROM_SEC)
    print("Cache: processed {:} files".format(len(processed_fnames)))
    print_stats(table, total_sim_time_sec)
except:
    print("NOT USING CACHE")
    table={}
    n=0
    total_sim_time_sec=0.0
    processed_fnames=set()
#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print("Processing", len(fnames), " files that are not present in cache")
for fname in fnames:
    try:
        i = get_output_num(fname)
        #    if (i % 7 <> modulus):
        #        continue
        output= Output()
        with open(fname, "rb") as f:
            output.ParseFromString(f.read())
        start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
        end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
        sim_time_sec = end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
    except ValueError as e:
        print("EXCEPTION ", fname, e)
    except:
        print("EXCEPTION ", fname)
        continue
    processed_fnames.add(fname)
    if sim_time_sec<0:
        print("SKIP", fname, " ", sim_time_sec)
        continue
    total_sim_time_sec= total_sim_time_sec+sim_time_sec
    for floater in output.statistics.floaters:
        prev_time_ns = start_sim_time_sec*1e9
        for op in floater.order_params[1:]:
            if prev_time_ns < STATS_FROM_NS:
                prev_time_ns = op.time_ns
                continue
            dt = op.time_ns - prev_time_ns
            D=op.diffusion_coefficient
            I=op.interacting_fraction
            accumulate(table, "D_" + floater.type, D, dt)
            accumulate(table, "I_" + floater.type, I, dt)
            prev_time_ns = op.time_ns
    n=n+1
    if(n%10>0 or n==0):
        continue
    with open(CACHE_FNAME,'wb') as f:
        pickle.dump([ table,
                      n,
                      total_sim_time_sec,
                      processed_fnames,
                      STATS_FROM_SEC
                      ]
                    , f)

    print_stats(table, total_sim_time_sec)
print_stats(table, total_sim_time_sec)
