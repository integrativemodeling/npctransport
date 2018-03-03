#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle

CACHE_FNAME='TMP.transport_stats_cache.%s.%s.p' % (sys.argv[1], sys.argv[2])
# confidence interval required
CONF=0.95
print sys.argv
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds

def accumulate(dictionary, key, value):
    if key in dictionary:
        dictionary[key] = dictionary[key] + value
    else:
        dictionary[key] = value

def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def print_stats(N, table, total_sim_time_sec):
    # print transports per particle per second
    print "Total simulation time: %.6f [sec/particle]" % total_sim_time_sec
    keys=sorted(table.keys())
    for key in keys:
        #    mean_per_sec=table[key]/n/sim_time_sec;
        mean=table[key]/total_sim_time_sec
        sem=math.sqrt(mean/total_sim_time_sec/N[key])
        print key, '%8.1f +- %8.2f /sec/particle (confidence %d%%);   time for all particles=%.3f [sec]' % \
            (mean, get_sems_for_conf(CONF)*sem, CONF*100, total_sim_time_sec*N[key])

############# Main ############
fnames=set(glob.glob(sys.argv[1]+"*.pb"))
try:
    with open(CACHE_FNAME,'rb') as f:
        [ N,
          table,
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
    N={}
    table={}
    n=0
    total_sim_time_sec=0.0
    processed_fnames=set()
#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print "Processing", len(fnames), " files that are not present in cache"
for fname in fnames:
    is_loaded=False
    try:
#        i = get_output_num(fname)
        #    if (i % 7 <> modulus):
        #        continue
        output= Output()
        with open(fname, "rb") as f:
            output.ParseFromString(f.read())
        start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
        end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
        sim_time_sec = end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
        is_loaded=True
    except ValueError as e:
        print "EXCEPTION ", fname, e
    except:
        print "EXCEPTION ", fname
    if not is_loaded:
        continue
    processed_fnames.add(fname)
    if sim_time_sec<0:
        print "SKIP", fname, " ", sim_time_sec
        continue
    total_sim_time_sec= total_sim_time_sec+sim_time_sec
    for floater in output.assignment.floaters:
        cur_N=floater.number.value
        if floater.type in N:
            assert(N[floater.type]==cur_N)
        else:
            N[floater.type]=cur_N
    for floater in output.statistics.floaters:
#        accumulate(table,floater.type,floater.avg_n_transports);
        count=0.0
        for t_ns in floater.transport_time_points_ns:
            t_sec=t_ns*1e-9
            if t_sec>STATS_FROM_SEC:
                count=count+1.0
        accumulate(table, floater.type, count/N[floater.type])
    n=n+1
    #    total_sim_time_sec=total_sim_time_sec+output.assignment.simulation_time_ns*(1e-9)
    if(n%10>0 or n==0):
        continue
    with open(CACHE_FNAME,'wb') as f:
        pickle.dump([ N,
                      table,
                      n,
                      total_sim_time_sec,
                      processed_fnames,
                      STATS_FROM_SEC
                      ]
                    , f)

print_stats(N, table, total_sim_time_sec)
