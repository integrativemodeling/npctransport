#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats

# confidence interval required
CONF=0.95
STATS_FROM_SEC=float(sys.argv[2]) # start stats after specified seconds
STATS_FROM_NS=STATS_FROM_SEC*(1E+9)

def update_transport_stats(dictionary, key, passage_time_ns):
    passage_time_ns=float(passage_time_ns)
    if key in dictionary:
        item=dictionary[key]
        item['sum']= item['sum'] +passage_time_ns
        item['sum2']=item['sum2']+passage_time_ns**2
        item['n']=   item['n']+1
        item['list'].append(passage_time_ns)
        print >>item['f'], passage_time_ns
    else:
        f=open('transport_times_ns.%s.txt' % key, 'w')
        dictionary[key] = {'sum':passage_time_ns,
                           'sum2':passage_time_ns**2,
                           'n':1,
                           'list':[passage_time_ns],
                           'f':f}
        print >>f, passage_time_ns


def test_update_transport_stats():
    A={}
    for i in range(1,11):
        update_transport_stats(A,'tmp1',i)
        update_transport_stats(A,'tmp2',i*10)
    assert(A['tmp1']['sum']==11*5)
    assert(A['tmp2']['sum']==110*5)
    assert(A['tmp1']['sum2']==385)
    assert(A['tmp2']['sum2']==38500)
    assert(A['tmp1']['n']==10)
    assert(A['tmp2']['n']==10)


def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

#################### MAIN ###########
test_update_transport_stats()
transport_stats={}
N={}
#modulus = int(sys.argv[2])
n=0
total_sim_time_sec=0.0
excess_ns={}
for fname in glob.glob(sys.argv[1]+"*.pb"):
    try:
        i = get_output_num(fname)
        #    if (i % 7 <> modulus):
        #        continue
        f=open(fname, "rb")
        #    print "File", fname
        output= Output()
        output.ParseFromString(f.read())
        f.close()
    except ValueError as e:
        print "EXCEPTION ", fname, e
    except:
        print "EXCEPTION ", fname
        continue
    start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
    end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
    sim_time_sec = end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
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
    # save sum, sum**2 and n for passage times (time between transport events)
    # for particles by their radius; account for excess time at end of each simulation
    for floater in output.statistics.floaters:
        t0_ns=STATS_FROM_NS
        if not floater.type in excess_ns:
            excess_ns[floater.type]=0.0
        for t1_ns in floater.transport_time_points_ns:
            if(t1_ns<STATS_FROM_NS):
                continue
            passage_time_ns=t1_ns-t0_ns
            if t0_ns==STATS_FROM_NS:
                if excess_ns[floater.type]>=1000.0:
                    passage_time_ns=passage_time_ns+excess_ns[floater.type]
                    update_transport_stats(transport_stats, floater.type,
                                           passage_time_ns)
                excess_ns[floater.type]=0.0
            else:
                update_transport_stats(transport_stats, floater.type,
                                       passage_time_ns)
            t0_ns=t1_ns
        cur_excess_ns=end_sim_time_sec*(1E+9)-t0_ns
        excess_ns[floater.type]=excess_ns[floater.type]+cur_excess_ns
    f.close()
    #    total_sim_time_sec=total_sim_time_sec+output.assignment.simulation_time_ns*(1e-9)

# print transports per particle per second
print "Total simulation time: %.6f [sec/particle]" % total_sim_time_sec
keys=sorted(transport_stats.keys())
for key in keys:
#   E_X is mean transport time in ns over N[key] particles
    E_X_ns  = transport_stats[key]['sum'] /transport_stats[key]['n']
    E_X2_ns2=transport_stats[key]['sum2']/transport_stats[key]['n']
    try:
        stddev_ns=math.sqrt(E_X2_ns2 - E_X_ns**2)
        sem_ns=stddev_ns/math.sqrt(transport_stats[key]['n'])
    except:
        sem_ns=-1.0
    # normalize by number of particles of that type:
    E_X_ns=E_X_ns*N[key]
    stddev_ns=stddev_ns*N[key]
    sem_ns=sem_ns*N[key]
    # convert from ns to sec
    E_X_sec=E_X_ns*1E-9
    stddev_sec=stddev_ns*1E-9
    sem_sec=sem_ns*1E-9
    print key, '%8.6f +- %8.6f [sec] (confidence %d%%);   time for all particles=%.3f [sec]' % \
        (E_X_sec, get_sems_for_conf(CONF)*sem_sec, CONF*100, total_sim_time_sec*N[key])
    transport_stats[key]['f'].close()
