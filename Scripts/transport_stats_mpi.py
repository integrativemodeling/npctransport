#!/usr/bin/env python
from __future__ import print_function
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import numpy as np
import time
if sys.version_info.major == 3:
    import pickle
else:
    import cPickle as pickle
import multiprocessing
import os.path
try:
    from schwimmbad.mpi import MPIPool, MPI
    MPI=True
except:
    MPI=False
import pandas as pd

print(sys.argv)
OUTPUT_PREFIX= sys.argv[1]
STATS_FROM_SEC_STR= sys.argv[2]# start stats after specified seconds
STATS_FROM_SEC= float(STATS_FROM_SEC_STR)
IS_REPORT_CACHE= len(sys.argv)>3 and sys.argv[3]=='report_cache_only'
CACHE_FNAME= 'TMP.transport_stats_cache.{}.{}.p' .format\
    (OUTPUT_PREFIX.replace("/","_").replace(".",""), \
     STATS_FROM_SEC_STR)
OUTPUT_FILENAME= 'STATS_{}_from_{}_seconds.csv'.format\
    (OUTPUT_PREFIX.replace("/","_").replace(".",""), \
     STATS_FROM_SEC_STR)
# confidence interval required
CONF= 0.95

def is_picklable(obj):
  ''' Return true if obj is picklable '''
  try:
    pickle.dumps(obj)

  except pickle.PicklingError:
    return False
  return True

def pickle_globals(cache_fname):
    global N
    global table
    global n
    global total_sim_time_sec
    global processed_fnames
    global STATS_FROM_SEC
    print("UPDATING CACHE with {:} processed files".format(len(processed_fnames)))
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
             ]= pickle.load(f)
            assert(fnames.issuperset(processed_fnames))
            assert(cached_stats_from_sec == STATS_FROM_SEC)
            print("Cache: processed {} files".format(len(processed_fnames)))
            #   except KeyboardInterrupt as e:
            #           raise
    except:
        print("NOT USING CACHE")
        N= {}
        table= {}
        n= 0
        total_sim_time_sec= 0.0
        processed_fnames= set()

def accumulate(dictionary, key, value):
    if key in dictionary:
        dictionary[key]= dictionary[key] + value
    else:
        dictionary[key]= value

def get_output_num(fname):
    m= re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None:
        raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a= 0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def print_stats(N,
                table,
                total_sim_time_sec,
                output_filename):
    '''
    print transports per particle per second to stdout and to file output_filename in csv format
    '''
    global STATS_FROM_SEC_STR
    print("Total simulation time: {:.6f} [sec/particle]".format( total_sim_time_sec))
    rows= []
    keys= sorted(table.keys())
    for key in keys:
        #    mean_per_sec=table[key]/n/sim_time_sec;
        mean= table[key]/total_sim_time_sec
        sem= math.sqrt(mean/total_sim_time_sec/N[key])
        print(key +
              '{:8.1f} +- {:8.2f} /sec/particle (confidence {:.0f}%);   time for all particles={:.3f} [sec]'.format\
              (mean,
               get_sems_for_conf(CONF)*sem,
               CONF*100,
               total_sim_time_sec*N[key]))
        rows.append({'type':key,
                     'time_sec_mean':mean,
                     'time_sec_sem':sem,
                     'total_sim_time_sec':total_sim_time_sec*N[key]
                     })
    df=pd.DataFrame(rows)
    df.to_csv(output_filename,
              index_label= "From " + STATS_FROM_SEC_STR)


def _open_file_exception(fname):
    global STATS_FROM_SEC

    print("Opening {}".format(fname))
    is_loaded= False
    is_ok= False
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
        sim_time_sec= end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
    if sim_time_sec<0:
        print("SKIP negative normalized sim-time {} stats from sec {}" \
              .format(sim_time_sec, STATS_FROM_SEC))
        ret_value= {"fname":fname,"status":-1}
        assert(is_picklable(ret_value))
        return ret_value
    Ns_dict= {} # a dictionary that maps the number of molecules for each molecule in the simulation
    counts_dict= {}
    for floater_a in output.assignment.floaters:
        Ns_dict[floater_a.type]= floater_a.number.value
    for floater_s in output.statistics.floaters:
        count= 0.0
        for t_ns in floater_s.transport_time_points_ns:
            t_sec= t_ns*1e-9
            if t_sec>STATS_FROM_SEC:
                count= count+1.0
                counts_dict[floater_s.type]= count
    file_summary= {"fname":fname,
                   "sim_time_sec":sim_time_sec,
                   "Ns_dict": Ns_dict,
                   "counts_dict":counts_dict,
                   "status":0}
    print("Opened {}".format(fname))
    assert(is_picklable(file_summary))
#    np.random.RandomState([ord(c) for c in fname])
#    if np.random.rand()<0.5:
#        raise ValueError("DEBUG exception " + fname)
    return file_summary

def open_file(fname):
    try:
        return _open_file_exception(fname)
    except Exception as e:
        print("Exception in file {} - {}".format(fname, e))
        return {"fname":fname, "status":-1}


def _sum_output_stats_exception(file_summary):
    global N
    global table
    global n
    global processed_fnames
    global CACHE_FNAME
    global STATS_FROM_SEC
    global total_sim_time_sec
    cache_frequency= 500

    if file_summary["status"]<0:
        print("Empty or defect output file skipped {}".format(file_summary["fname"]))
        return
    processed_fnames.add(file_summary["fname"])
    total_sim_time_sec= total_sim_time_sec+file_summary["sim_time_sec"]
    for floater_type, floater_N in file_summary["Ns_dict"].items():
        if floater_type in N:
            assert(N[floater_type]==floater_N)
        else:
            N[floater_type]= floater_N
    for floater_type, count in file_summary["counts_dict"].items():
        accumulate(table, floater_type, count/N[floater_type])
        n= n+1
    if(len(processed_fnames) % 10 == 0):
        print("Number of files processed {0:d}".format(len(processed_fnames)))
    if(len(processed_fnames) % cache_frequency == 0):
        if CACHE_FNAME is not None:
            pickle_globals(CACHE_FNAME)

def sum_output_stats(file_summaries):
    if isinstance(file_summaries, dict):
        file_summaries= [file_summaries]
    for file_summary in file_summaries:
        try:
            _sum_output_stats_exception(file_summary)
        except KeyboardInterrupt:
            print("Keyboard interruption caught")
            raise
        except:
            if file_summary is not None:
                print("Exception summarizing stats of {0}".format(file_summary["fname"]))
#                raise
            else:
                print("Unknown exception in sum_output_stats()")
                #        raise

############# Main ###########
fnames= set(glob.glob(OUTPUT_PREFIX+"*.pb"))
unpickle_or_initialize_globals(CACHE_FNAME)
if IS_REPORT_CACHE:
    print_stats(N, table, total_sim_time_sec,
                OUTPUT_FILENAME)
    sys.exit()
fnames= fnames.difference(processed_fnames)
print("Processing {:d} files that are not present in cache".format(len(fnames)))
fnames= list(fnames)
print("Starting pool")
if MPI:
    pool= MPIPool()
    pool.wait(lambda: sys.exit(0))
    results= pool.map(open_file,
                      fnames,
                      callback=sum_output_stats)
    pool.close()
else:
    pool= multiprocessing.Pool(processes=40)
    #manager= multiprocessing.Manager()
    #processed_fnames= manager.list(processed_fnames)
    results= pool.map_async(open_file,
                            fnames,
                            callback=sum_output_stats)
    results.wait()
    pool.close()
    pool.join()
    print("Pool joined")
print("Pool closed")
pickle_globals(CACHE_FNAME)
print_stats(N, table, total_sim_time_sec,
            OUTPUT_FILENAME)
