#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle

CACHE_FNAME='TMP.d_stats.%s.%s.p' % (sys.argv[1], sys.argv[2])
# confidence interval required
CONF_INTERVAL=0.95
print sys.argv
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds

def accumulate(dictionary, key, value):
    entry=[1,value,value**2]
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
    print "Total simulation time: %.6f [sec]" % total_sim_time_sec
    keys=sorted(table.keys())
    for key in keys:
        #    mean_per_sec=table[key]/n/sim_time_sec;
        n= table[key][0]
        mean= table[key][1]/n
        mean2= table[key][2]/n
        var=mean2-mean**2
        stddev=math.sqrt(var)
        conf=get_sems_for_conf(CONF_INTERVAL)*stddev/math.sqrt(n)
        print key, '%8.2e +- %8.2e sec' % (mean, conf)

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
    print "Cache: processed", len(processed_fnames), "files"
except:
    print "NOT USING CACHE"
    table={}
    n=0
    total_sim_time_sec=0.0
    processed_fnames=set()
#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print "Processing", len(fnames), " files that are not present in cache"
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
        print "EXCEPTION ", fname, e
    except:
        print "EXCEPTION ", fname
        continue
    processed_fnames.add(fname)
    if sim_time_sec<0:
        print "SKIP", fname, " ", sim_time_sec
        continue
    total_sim_time_sec= total_sim_time_sec+sim_time_sec
    for floater in output.statistics.floaters:
#        accumulate(table,floater.type,floater.avg_n_transports);
        count=0.0
        D=floater.diffusion_coefficient
        accumulate(table, floater.type, D)
    n=n+1
    #    total_sim_time_sec=total_sim_time_sec+output.assignment.simulation_time_ns*(1e-9)
    if(n%40>0 or n==0):
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
