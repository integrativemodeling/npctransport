#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle
import multiprocessing
import os.path

CACHE_FNAME='TMP.transport_stats_cache.{}.{}.p' .format(os.path.basename(sys.argv[1]),
                                                       sys.argv[2])
# confidence interval required
CONF=0.95
print sys.argv
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds

def pickle_globals(cache_fname):
    global N
    global table
    global n
    global total_sim_time_sec
    global processed_fnames
    global STATS_FROM_SEC
    print "UPDATING CACHE"
    with open(cache_fname,'wb') as f:
        pickle.dump([ N,
                      table,
                      n,
                      total_sim_time_sec,
                      set(processed_fnames),
                      STATS_FROM_SEC ]
                    , f)

def unpickle_or_initialize_globals(cache_fname):
    global N
    global table
    global n
    global total_sim_time_sec
    global processed_fnames
    global cached_stats_from_sec
    try:
        with open(cache_fname,'rb') as f:
            [ N,
              table,
              n,
              total_sim_time_sec,
              processed_fnames,
              cached_stats_from_sec
            ] = pickle.load(f)
            assert(fnames.issuperset(processed_fnames))
            assert(cached_stats_from_sec == STATS_FROM_SEC)
            print("Cache: processed {} files".format(len(processed_fnames)))
    except:
        print "NOT USING CACHE"
        N={}
        table={}
        n=0
        total_sim_time_sec=0.0
        processed_fnames=set()

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


def open_file(fname, processed_fnames):
    global STATS_FROM_SEC

    print("Opening {}".format(fname))
    is_loaded=False
    try:
#        i = get_output_num(fname)
        #    if (i % 7 <> modulus):
        #        continue
        output= Output()
        with open(fname, "rb") as f:
            output= Output()
            fstring= f.read()
            output.ParseFromString(fstring)
            print("Parsed {}".format(fname))
        start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
        end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
        sim_time_sec = end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
    except ValueError as e:
        print("EXCEPTION {} {}".format(fname, e))
        raise
    except:
        print("EXCEPTION {}".format(fname))
        raise
    if sim_time_sec<0:
        print("SKIP negative normalized sim-time {} stats from sec {}" \
              .fomrat(sim_time_sec, STATS_FROM_SEC))
        return None
    Ns_dict={}
    counts_dict={}
    for floater_a in output.assignment.floaters:
        Ns_dict[floater_a.type]= floater_a.number.value
    for floater_s in output.statistics.floaters:
        count=0.0
        for t_ns in floater_s.transport_time_points_ns:
            t_sec=t_ns*1e-9
            if t_sec>STATS_FROM_SEC:
                count=count+1.0
        counts_dict[floater_s.type]= count
    file_summary={"fname":fname,
                  "sim_time_sec":sim_time_sec,
                  "Ns_dict": Ns_dict,
                  "counts_dict":counts_dict}
    print("Opened {}".format(fname))
    processed_fnames.append(fname)
    return file_summary




def _sum_output_stats_exception(file_summary):
    global N
    global table
    global n
    global processed_fnames
    global CACHE_FNAME
    global STATS_FROM_SEC
    global total_sim_time_sec
    cache_frequency=500

    if file_summary is None:
        print("Empty output skipped")
    total_sim_time_sec= total_sim_time_sec+file_summary["sim_time_sec"]
    for floater_type, floater_N in file_summary["Ns_dict"].iteritems():
        if floater_type in N:
            assert(N[floater_type]==floater_N)
        else:
            N[floater_type]=floater_N
    for floater_type, count in file_summary["counts_dict"].iteritems():
        accumulate(table, floater_type, count/N[floater_type])
    n=n+1
    if(len(processed_fnames) % 10 == 0):
        print("Number of files processed {0:d}".format(len(processed_fnames)))
    if(len(processed_fnames) % cache_frequency == 0):
        if CACHE_FNAME is not None:
            pickle_globals(CACHE_FNAME)

def sum_output_stats(file_summary):
    try:
        _sum_output_stats_exception(file_summary)
    except:
        if file_summary is not None:
            print("Exception summarizing stats of {0}".format(file_summary["fname"]))
        else:
            print("Unknown exception in sum_output_stats()")
        raise



############# Main ###########
fnames=set(glob.glob(sys.argv[1]+"*.pb"))
unpickle_or_initialize_globals(CACHE_FNAME)
if len(sys.argv)>3 and sys.argv[3]=='report_cache_only':
    print_stats(N, table, total_sim_time_sec)
    sys.exit()

#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print "Processing", len(fnames), " files that are not present in cache"
pool= multiprocessing.Pool(processes=8)
manager= multiprocessing.Manager()
processed_fnames= manager.list(processed_fnames)
print("Starting pool")
R=[]
fnames=list(fnames)
for fname in fnames:
    try:
#        sum_output_stats(open_file(fname, processed_fnames))
        r= pool.apply_async(open_file,
                            args=(fname, processed_fnames),
                            callback=sum_output_stats)
        R.append(r)
    except:
        print("Failed on {0}".format(fname))
pool.close()
pool.join()
print("Pool closed and joined")
print("Verifying status of individual processes...")
for i,r in enumerate(R):
    try:
        r.get()
    except:
        print("Process {:d} failed".format(i))
print("Pickling and reporting stats:")
pickle_globals(CACHE_FNAME)
print_stats(N, table, total_sim_time_sec)
