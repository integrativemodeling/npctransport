#TODO: STILL NEED TO BE TWRITTEN
#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle

CACHE_FNAME='TMP.rg_stats_cache.%s_%s.pickle' % (sys.argv[1], sys.argv[2])
# confidence interval required
CONF=0.95
print sys.argv
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds
STATS_FROM_NS=STATS_FROM_SEC*(1e+9)

def accumulate(dictionary, key, value, time_ns):
    new_value=[value*time_ns, time_ns, value**2*time_ns, 1] # mean / time / mean-square / n
    if key in dictionary:
        dictionary[key] = [x+y for x,y in zip (dictionary[key],new_value)]
    else:
        dictionary[key] = new_value

def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def do_stats(table, units_label, is_fraction=False):
    ''' return dictionary with mean and standard-error of mean for each key in table'''
    # print transports per particle per second
    keys=sorted(table.keys())
    ret={}
    for key in keys:
        t_ns= table[key][1]
        mean= table[key][0]/(t_ns+0.0001)
        mean2= table[key][2]/(t_ns+0.0001)
        n= table[key][3]
        stdvar= math.sqrt(mean2-mean**2)
        stderr= stdvar/math.sqrt(n+0.0001) # dividing by n is probably ok when constant stat intervals are used
        conf=get_sems_for_conf(CONF)
        if(is_fraction):
            print key, 'mean-k: %10.2f%% +- %10.2f%% [%s] ; t=%.1f [ns]' % \
                (mean*100, stderr*100*conf, units_label, t_ns)
        else:
            print key, 'mean-k: %10.8f +- %10.8f [%s] ; t=%.1f [ns]' % \
                (mean, stderr*conf, units_label, t_ns)
        ret[key]=(mean, stderr, t_ns)
    return ret

############# Main ############
fnames=set(glob.glob(sys.argv[1]+"*.pb"))
try:
    with open(CACHE_FNAME,'rb') as f:
        [ k_ons_i,
          k_on2s_i,
          k_offs_i,
          k_off2s_i,
          k_ons_ii,
          k_on2s_ii,
          k_offs_ii,
          k_off2s_ii,
          fbounds,
          fbounds2,
          bound_times,
          unbound_times,
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
    k_ons_i={}
    k_on2s_i={}
    k_offs_i={}
    k_off2s_i={}
    k_ons_ii={}
    k_on2s_ii={}
    k_offs_ii={}
    k_off2s_ii={}
    fbounds={}
    fbounds2={}
    bound_times={}
    unbound_times={}
    n=0
    total_sim_time_sec=0.0
    processed_fnames=set()
#modulus = int(sys.argv[2])
fnames=fnames.difference(processed_fnames)
print "Processing", len(fnames), " files that are not present in cache"
for fname in fnames:
    try:
        output= Output()
        with open(fname, "rb") as f:
            output.ParseFromString(f.read())
        start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
        end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
        sim_time_sec = end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
#        print "TIME", sim_time_sec
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
    # for floater in output.statistics.floaters:
    #     ftype = floater.type
    #     was_bound=floater.order_params[0].interacting_fraction>0.0
    #     n_on=0
    #     n_off=0
    #     time_state=0.0
    #     for fop in floater.order_params:
    #         if fop.time_ns<STATS_FROM_NS:
    #             continue
    #         fbound = fop.interacting_fraction
    #         accumulate(fbounds, ftype, fbound, 1)
    #         is_bound=fbound>0.0
    #         if(is_bound != was_bound):
    #             if is_bound:
    #                 accumulate(unbound_times, ftype, time_state, 1) # finished unbound
    #                 n_on=n_on+1
    #             else:
    #                 accumulate(bound_times, ftype, time_state, 1) # finished bound
    #                 n_off=n_off+1
    #             was_bound=is_bound
    #             time_state=0.0
    #         else:
    #             time_state=time_state+1
    ii=0;
    for interaction in output.statistics.interactions:
        iname = interaction.type0 + "-" + interaction.type1
        for iop in interaction.order_params:
            if iop.time_ns<STATS_FROM_NS:
                continue
            k_on_i=  iop.avg_on_per_unbound_i_per_ns
            t_on_i_ns= iop.on_i_stats_period_ns
            k_off_i= iop.avg_off_per_bound_i_per_ns
            t_off_i_ns=iop.off_i_stats_period_ns
            k_on_ii=  iop.avg_on_per_unbound_ii_per_ns
            t_on_ii_ns= iop.on_ii_stats_period_ns
            k_off_ii= iop.avg_off_per_bound_ii_per_ns
            t_off_ii_ns=iop.off_ii_stats_period_ns
            if(k_on_i != float('inf')):
                accumulate(k_ons_i, iname, k_on_i, t_on_i_ns)
            if(k_off_i != float('inf')):
                accumulate(k_offs_i, iname, k_off_i, t_off_i_ns)
            if(k_on_ii != float('inf')):
                accumulate(k_ons_ii, iname, k_on_ii, t_on_ii_ns)
            if(k_off_ii != float('inf')):
                accumulate(k_offs_ii, iname, k_off_ii, t_off_ii_ns)
            fbound2=iop.avg_fraction_bound_particles_ii
            if(fbound2 != float('inf')):
                accumulate(fbounds2, iname, fbound2, 1)
            ii=ii+1
            if(ii%10000==0):
                fbound2_stats=do_stats(fbounds2, "of floats (from interaction order params)", is_fraction=True)
    n=n+1
    if(n%5>0 or n==0):
        continue
    with open(CACHE_FNAME,'wb') as f:
        pickle.dump(
            [ k_ons_i,
              k_on2s_i,
              k_offs_i,
              k_off2s_i,
              k_ons_ii,
              k_on2s_ii,
              k_offs_ii,
              k_off2s_ii,
              fbounds,
              fbounds2,
              bound_times,
              unbound_times,
              n,
              total_sim_time_sec,
              processed_fnames,
              STATS_FROM_SEC
            ]
            , f)

print "K_on_i:"
on_stats_i=do_stats(k_ons_i, "per ns per unbound fg motif")
print "K_off_i:"
off_stats_i=do_stats(k_offs_i, "per ns per bound fg motif")
print "K_on_ii:"
on_stats_ii=do_stats(k_ons_ii, "per ns per unbound float")
print "K_off_ii:"
off_stats_ii=do_stats(k_offs_ii, "per ns per bound float")
print "fraction bounds:"
#do_stats(fbounds, "of floats", is_fraction=True)
fbound2_stats=do_stats(fbounds2, "of floats (from interaction order params)", is_fraction=True)
for key in on_stats_ii.keys():
    [kon,_,t_kon]= on_stats_ii[key]
    [koff,_,t_koff]= off_stats_ii[key]
    if t_kon*t_koff==0.0: continue
    kd=koff/kon
    n_on= kon* t_kon
    n_off=koff*t_koff
    stderr_n_on= math.sqrt(n_on)
    stderr_n_off= math.sqrt(n_off)
    conf= get_sems_for_conf(CONF)
    n_on_high=  n_on +  conf*stderr_n_on
    n_on_low=   n_on -  conf*stderr_n_on
    n_off_high= n_off + conf*stderr_n_off
    n_off_low=  n_off - conf*stderr_n_off
    kon_high=  n_on_high /t_kon
    kon_low=   n_on_low  /t_kon
    koff_high= n_off_high/t_koff
    koff_low=  n_off_low /t_koff
    kd_high= koff_high/kon_low
    kd_low=  koff_low /kon_high
    fb=kon/(kon+koff)
    fb_high=kon_high/(kon_high+koff_low)
    fb_low= kon_low /(kon_low +koff_high)
    print key, "Kd %.4f range (%.4f .. %.4f)" % ( kd, kd_low, kd_high)
    print key, "%% bound estimate from kon/koff times %.1f%% range (%.1f%% .. %.1f%%)" % (100*fb, 100*fb_low, 100*fb_high)
#[A]<-->[A*] # ignoring [B] in these statistics
#kon[A]=koff[A*]
#[A]=kd[A*]
#[A*]/(kd+1)[A*]= 1/(kd+1) = kon/(koff+kon)

# print "bound times:"
# bstat= do_stats(bound_times, "x100 ns")
# print "unbound times (x100 ns):"
# ustat= do_stats(unbound_times, "x100 ns")
# for key in bstat.keys():
#     [mb,sb,_]= bstat[key]
#     [mu,su,_]= ustat[key]
#     conf= get_sems_for_conf(CONF)
#     kd= mu/mb
#     kd_high= (mu+su*conf)/(mb-sb*conf)
#     kd_low=  (mu-su*conf)/(mb+sb*conf)
#     fb= mb/(mu+mb)
#     fb_high= (mb+sb*conf)/(mb+sb*conf+mu-su*conf) # since fb<1.0, take mb+sb*conf rather than mb-sb*conf in both nominator and denominator
#     fb_low=  (mb-sb*conf)/(mb-sb*conf+mu+su*conf) # since fb<1.0, take mb+sb*conf rather than mb-sb*conf in both nominator and denominator
#     print key, "Kd", kd, "range", kd_low, kd_high
#     print key, "fracton bound estimate from bound/unbound times", fb, "range", fb_low, fb_high
