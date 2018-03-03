#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle

CACHE_FNAME='TMP.rg_stats_cache.p.%s' % sys.argv[2]
# confidence interval required
CONF=0.95
print sys.argv
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds
STATS_FROM_NS=STATS_FROM_SEC*(1e+9)

def accumulate(dictionary, key, value, weight):
    if weight==0.0:
        return
    if key in dictionary:
        old_entry=  dictionary[key]
        old_value= old_entry[0]
        old_weight= old_entry[1]
        new_weight= weight+old_weight+0.00001
        dictionary[key]= [(old_value*old_weight+value*weight)/new_weight, new_weight]
    else:
        dictionary[key]= [value, weight]

def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def do_stats(RGs, RG2s):
    # print transports per particle per second
    keys=sorted(RGs.keys())
    for key in keys:
        assert(key in RG2s)
        #    mean_per_sec=table[key]/n/sim_time_sec;
        mean= RGs[key][0]
        time_sec= RGs[key][1]
        mean2= RG2s[key][0]
        std= math.sqrt(mean2-mean**2)
#        sem=std/math.sqrt(n)
        print key, 'mean-rg: %8.1f [A] (std-dev: %8.1f [A])  time=%8.1f [sec]' % \
            (mean, std, time_sec)

############# Main ############
fnames=set(glob.glob(sys.argv[1]+"*.pb"))
try:
    with open(CACHE_FNAME,'rb') as f:
        [ RGs,
          RG2s,
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
    RGs={}
    RG2s={}
    n=0
    total_sim_time_sec=0.0
    processed_fnames=set()
#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print "Processing", len(fnames), " files that are not present in cache"
for fname in fnames:
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
    except ValueError as e:
        print "EXCEPTION ", fname, e
        continue
    except:
        print "EXCEPTION ", fname
        continue
    processed_fnames.add(fname)
    if sim_time_sec<0:
        print "SKIP", fname, " ", sim_time_sec
        continue
    total_sim_time_sec= total_sim_time_sec+sim_time_sec
    for fg in output.statistics.fgs:
        prev_time_ns= STATS_FROM_NS
        for fgop in fg.order_params:
            if fgop.time_ns<STATS_FROM_NS:
                continue
            time_interval_ns= fgop.time_ns-prev_time_ns
            rg=fgop.mean_radius_of_gyration
            rg2=fgop.mean_square_radius_of_gyration
            accumulate(RGs, fg.type, rg, time_interval_ns)
#            print "%.0f %.0f %.2f %.0f" % (fgop.time_ns, prev_time_ns, rg, time_interval_ns)
            accumulate(RG2s, fg.type, rg2, time_interval_ns);
            n=n+1
            prev_time_ns= fgop.time_ns

    #    total_sim_time_sec=total_sim_time_sec+output.assignment.simulation_time_ns*(1e-9)
    if(n%50>0 or n==0):
        continue
    with open(CACHE_FNAME,'wb') as f:
        pickle.dump([ RGs,
                      RG2s,
                      n,
                      total_sim_time_sec,
                      processed_fnames,
                      STATS_FROM_SEC
                      ]
                    , f)

do_stats(RGs, RG2s)
