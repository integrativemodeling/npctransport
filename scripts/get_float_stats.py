#!/usr/bin/env python
from IMP.npctransport import *
#import seaborn as sbn
import numpy as np
#import pandas as pd
import sys
import random
import cPickle as pickle

# TODO: add stats time to pickle and to do_stats(), though will invalidate old caches

def do_stats(N,S):
#    print "TOTAL STATS TIME [sec]: ", stats_time_ns*1E-9
    for i in range(N):
        for f_type,H in S[i].iteritems():
            fname="S%d.%s.txt" % (i, f_type)
            F=open(fname,'w')
            for row in H:
                for x in row:
                    print >>F, "%2d " % x,
                print >>F, ""
            F.close()

#import cPickle as p; D=p.load(open('TMP.get_float_stats_cache.p','rb')); print len(D[2]);

############# Main ############
CACHE_FNAME='TMP.get_float_stats_cache.p'
fnames=set(sys.argv[1:])
N=50
try:
    with open(CACHE_FNAME,'rb') as f:
        [ n, S, processed_fnames ] = pickle.load(f)
    assert(fnames.issuperset(processed_fnames))
    assert(n==N)
except:
    print "NOT USING CACHE"
    S=[{} for i in range(N)]
    processed_fnames=set()
stats_time_ns=0.0
fnames=fnames.difference(processed_fnames)
for ifile, fname in enumerate(fnames):
    try:
        F=open(fname, "rb")
        o= Output()
        o.ParseFromString(F.read())
    except:
        print "Skipping ", fname, " due to error"
        continue
    s=o.statistics
    print "Processing", fname
    if(s.bd_simulation_time_ns<100.0 or not s.floaters[0].HasField("zr_hist")):
        continue;
    stats_time_ns=stats_time_ns+s.bd_simulation_time_ns
    for i in range(N):
        for f in (set(s.floaters) | set(s.fgs)):
            if not f.HasField('zr_hist'):
                continue
            H=[]
            for row in f.zr_hist.ints_list:
                row_list=[x for x in row.ints]
                H.append(row_list)
            H=np.array(H)
            if(random.random()<0.99):
                if f.type in S[i]:
                    S[i][f.type]=S[i][f.type]+H
                else:
                    S[i][f.type]=H
    F.close()
    processed_fnames.add(fname)
    print ifile
    if(ifile % 50 == 0 and ifile>0):
        print "UPDATING CACHE"
        with open(CACHE_FNAME,'wb') as f:
            pickle.dump([ N, S, processed_fnames], f)
        do_stats(N,S)

print "UPDATING CACHE"
with open(CACHE_FNAME,'wb') as f:
    pickle.dump([ N, S, processed_fnames], f)
do_stats(N,S)
