#!/usr/bin/env python
#from __future__ import print_function
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import numpy as np
import sys
import pickle
import multiprocessing

try:
    from schwimmbad.mpi import MPIPool, MPI
    _tmp_pool = MPIPool()
    MPI = True
except:
    MPI = False
import pandas as pd

def parse_args():
    if sys.argv[1] in ['-h', '--help']:
        print("Usage: <IMP> python transport_stats.py OUTPUT_PREFIX STATS_FROM_SEC [N_PROCESSORS or report_cache_only]")
        print()
        print("       OUTPUT_PREFIX: prefix of output .pb files to process")
        print("       STATS_FROM_SEC: start stats after specified seconds of simulation")
        print("       N_PROCESSORS: number of processors to use, default=1")
        print("       report_cache_only: report stats from cache only")
        print()
        print("Example: transport_stats.py output/N3/output 5e-6 4")
        print("      Process all output files in output/N3/output* from 5 microseconds using 4 CPUs")
        print()
        print("Example: transport_stats.py output/N3/output 5e-6 report_cache_only")
        print("      Report stats from cache for output/N3/output* from 5 microseconds")
        sys.exit()
    print(sys.argv)
    OUTPUT_PREFIX= sys.argv[1]
    STATS_FROM_SEC_STR= sys.argv[2]# start stats after specified seconds
    N_PROCESSORS= int(sys.argv[3]) if len(sys.argv)>3 else 1
    IS_REPORT_CACHE= len(sys.argv)>3 and sys.argv[3]=='report_cache_only'
    return OUTPUT_PREFIX, STATS_FROM_SEC_STR, N_PROCESSORS, IS_REPORT_CACHE

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
    for floater_type in N.keys():
        if floater_type not in table:
            table[floater_type]= 0.0 # make sure all table entires are initialized, in case they weren't in old caches
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
    '''
    add value to dictionary[key] if key existsm otherwise add it and
    set it to value
    '''
    if key in dictionary:
        dictionary[key]= dictionary[key] + value
    else:
        dictionary[key]= value

def get_output_num(fname):
    m= re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None:
        raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

def get_poisson_interval(k, alpha=0.05):
    """
    uses chisquared info to get the 1-alpha poisson confidence
    interval for the rate of an event per time period given k observed
    events over that time period

    :param k: rate (time^-1) for which interval is computed
    :param alpha: one minus the confidence interval
    :return: a numpy array with low and high interval boundaries
    """
    a = alpha
    low, high = (stats.chi2.ppf(a/2, 2*k) / 2,
                 stats.chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0:
        low = 0.0
    return np.array([low, high])

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
        # Compute transport events per particle per second:
        n_per_particle= table[key]
        n_per_particle_per_sec= n_per_particle / total_sim_time_sec # per particle per second
        # Compute confidence intervals for transport events per particle per second:
        n_particles= N[key]
        n= n_per_particle * n_particles # total number of events observed
        n_per_particle_per_sec_interval= get_poisson_interval(n, 1.0 - CONF) \
            / total_sim_time_sec / n_particles
        sem= math.sqrt(n_per_particle_per_sec/total_sim_time_sec/n_particles) # for backward comp.
        print(
              '{:15s} {:7.1f} -> {:8.1f} <- {:<7.1f} /sec/particle (confidence {:.0f}%);   time for all particles={:.3f} [sec]'.format\
              (key,
               n_per_particle_per_sec_interval[0],
               n_per_particle_per_sec,
               n_per_particle_per_sec_interval[1],
               CONF*100,
               total_sim_time_sec*n_particles))
        rows.append({'type':key,
                     'time_sec_mean': n_per_particle_per_sec,
                     'time_sec_sem':sem,
                     'total_sim_time_sec':total_sim_time_sec * n_particles,
                     'n_per_particle_per_sec': n_per_particle_per_sec,
                     'n_per_particle_per_sec_lbound': n_per_particle_per_sec_interval[0],
                     'n_per_particle_per_sec_ubound': n_per_particle_per_sec_interval[1],
                     'confidence_level': CONF,
                     })
    df=pd.DataFrame(rows)
    df.to_csv(output_filename,
              index_label= "From " + STATS_FROM_SEC_STR)


def _open_file(fname):
    ''' 
    handle an individual file
    :return file_summary: a dictionary with file summary information
        with fields:
        "fname" - the file name being handled
        "sim_time_sec" - how much simulation time was accumulated in the file
        "Ns_dict" - a dictionary that maps the number of molecules of each type in the simulation
        "counts_dict" - a dictionary that maps the number of transport events for each
                        molecule type during the simulateion time
        "status" - 0 if successful, -1 if not
        
    '''
    global STATS_FROM_SEC

    print("Opening {}".format(fname))
    output= Output()
    with open(fname, "rb") as f:
        output= Output()
        fstring= f.read()
        output.ParseFromString(fstring)
 #       print("Parsed {}".format(fname))
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
    print("Statistics computed for {}".format(fname))
    assert(is_picklable(file_summary))
    return file_summary

def open_file_no_exception(fname):
    try:
        return _open_file(fname)
    except Exception as e:
        print("Exception in file {} - {}".format(fname, e))
        return {"fname":fname, "status":-1}

def _sum_output_stats(file_summary, cache_frequency=250):
    '''
    Accumulate statistics from the run results stored in file_summary
    :param cache_frequency: how often to cache results
    '''
    global N
    global table
    global n
    global processed_fnames
    global CACHE_FNAME
    global STATS_FROM_SEC
    global total_sim_time_sec

    if file_summary["status"]<0:
        print("Empty or defect output file skipped {}".format(file_summary["fname"]))
        return
    processed_fnames.add(file_summary["fname"])
    total_sim_time_sec += file_summary["sim_time_sec"]
    for floater_type, floater_N in file_summary["Ns_dict"].items():
        if floater_type in N:
            assert(N[floater_type]==floater_N)
        else:
            N[floater_type]= floater_N
    for floater_type, count in file_summary["counts_dict"].items():
        # add number of transport events per particle of type floater_type
        accumulate(table, floater_type, count/N[floater_type])
        n += 1
    if(len(processed_fnames) % 10 == 0):
        print("Number of files processed {0:d}".format(len(processed_fnames)))
    if(len(processed_fnames) % cache_frequency == 0):
        if CACHE_FNAME is not None:
            pickle_globals(CACHE_FNAME)

def sum_output_stats_no_exception(file_summaries):
    ''' 
    accumulate statistics from individual file summaries
    :param file_summaries: a list of file summary outputs from _open_file()
    '''
    print(f"Callback handling {len(file_summaries)} file summaries")
    if isinstance(file_summaries, dict):
        file_summaries= [file_summaries]
    for file_summary in file_summaries:
        try:
            _sum_output_stats(file_summary)
        except KeyboardInterrupt:
            print("Keyboard interruption caught")
            raise
        except Exception as e:
            if file_summary is not None:
                print("Exception summarizing stats of {0} '{1}'".format(file_summary["fname"], e))
            else:
                print("Unknown exception in sum_output_stats() '{0}'", e)
                
                
def process_files(fnames, is_mpi, n_processors=1):
    '''
    process files in fnames in parallel using n_processors if is_mpi is False,
    or the MPI protocol if is_mpi is True.
    '''
    print("Processing {:d} files that are not present in cache".format(len(fnames)))
    fnames= list(fnames)
    if is_mpi:
        print("Starting an MPI pool")
        pool= MPIPool()
        pool.wait(lambda: sys.exit(0))
        results= pool.map(open_file_no_exception,
                        fnames,
                        callback=sum_output_stats_no_exception)
        pool.close()
    else:
        print(f"Starting a pool of {N_PROCESSORS} processors")
        pool= multiprocessing.Pool(processes=N_PROCESSORS)
        results= pool.map_async(open_file_no_exception,
                                fnames,
                                callback=sum_output_stats_no_exception,
                                error_callback=lambda e: print("Error in process_files: {}".format(e)),
                                chunksize=2)
        results.wait()
        pool.close()
        pool.join()
    print("Pool closed")

############# Main ###########
if __name__ == '__main__':
    OUTPUT_PREFIX, STATS_FROM_SEC_STR, N_PROCESSORS, IS_REPORT_CACHE= parse_args()
    STATS_FROM_SEC= float(STATS_FROM_SEC_STR)
    CACHE_FNAME= 'TMP.transport_stats_cache.{}.{}.p'.format\
        (OUTPUT_PREFIX.replace("/","_").replace(".",""), \
        STATS_FROM_SEC_STR)
    OUTPUT_FILENAME= 'STATS_{}_from_{}_seconds.csv'.format\
        (OUTPUT_PREFIX.replace("/","_").replace(".",""), \
        STATS_FROM_SEC_STR)
    CONF= 0.95 # confidence interval required
    # Get going:
    fnames= set(glob.glob(OUTPUT_PREFIX+"*.pb"))
    unpickle_or_initialize_globals(CACHE_FNAME)
    if IS_REPORT_CACHE:
        print_stats(N, table, total_sim_time_sec,
                    OUTPUT_FILENAME)
        sys.exit()
    fnames= fnames.difference(processed_fnames)
    process_files(fnames, is_mpi=MPI, n_processors=N_PROCESSORS)
    pickle_globals(CACHE_FNAME)
    print_stats(N, table, total_sim_time_sec,
                OUTPUT_FILENAME)
