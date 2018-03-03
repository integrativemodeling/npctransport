
from IMP.npctransport import *
#import seaborn as sbn
import numpy as np
#import pandas as pd
import sys
import random
import cPickle as pickle

N=1
DISABLE_RANDOM=True or (N==1)

# TODO: add stats time to pickle and to do_stats(), though will invalidate old caches

def do_stats(N,S):
#    print "TOTAL STATS TIME [sec]: ", stats_time_ns*1E-9
    for i in range(N):
        print "hi",i
        for f_type,XYZ in S[i].iteritems():
            print "hello",f_type
            fname="S%d.%s.txt" % (i, f_type)
            F=open(fname,'w')
            for YZ in XYZ:
                for Z in YZ:
                    for n_xyz in Z:
                        print >>F, "%2d " % n_xyz,
                    print >>F, ""
            # print separator line of -1s
            #            if Z is not None:
            #                for n_xyz in Z:
            #                    print >>F, -1,
            #                print >>F, ""
            F.close()

#import cPickle as p; D=p.load(open('TMP.get_float_stats_cache.p','rb')); print len(D[2]);

############# Main ############
CACHE_FNAME='TMP.get_float_stats_cache.p'
fnames=set(sys.argv[1:])
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
    print "Processing", fname
    s=o.statistics
    if s.bd_simulation_time_ns<900.0: # or not o.HasField('xyz_hist_stats'):
           print "Skipping too-short", s.bd_simulation_time_ns, " simulation file",  fname
           continue
#    if not o.is_xyz_hist_stats or not s.floaters[0].HasField("xyz_hist"):
#        continue;
    stats_time_ns=stats_time_ns+s.bd_simulation_time_ns
    for i in range(N):
        for f in s.fgs:#(set(s.floaters) | set(s.fgs)):
            print f.type
            if not f.HasField('xyz_hist'):
                print "Skipping f.type"
                continue
            xyz=f.xyz_hist.ints_lists
            if(len(xyz)==0):
                continue
            yz0=xyz[0].ints_list
            if(len(yz0)==0):
                continue
            z0=yz0[0].ints
            XYZ=np.zeros([len(xyz),len(yz0),len(z0)])
            print np.shape(XYZ)
            for ii,yz_entry in enumerate(xyz):
                for jj,z_row in enumerate(yz_entry.ints_list):
                    for kk,n in enumerate(z_row.ints):
                        XYZ[ii][jj][kk]=n
            #            print "XYZ done", fname
            #            print "XYZ\n", XYZ, "\n=="
            if(random.random()<0.9 or DISABLE_RANDOM):
                #                f_type="FGs"
                f_type=f.type
#                if(f_type<>"Nup116" and f_type<> "Nup100" and f_type <> "Nsp1"):
#                    continue
                if f_type in S[i]:
                    try:
                        S[i][f_type]=S[i][f_type]+XYZ
                    except:
                        print "Couldn't add XYZ of shape", XYZ.shape,
                        print " to S[%d][%s] of shape" % (i, f_type),
                        print S[i][f_type].shape
                        continue
                else:
                    S[i][f_type]=XYZ

    F.close()
    processed_fnames.add(fname)
    print ifile
    if(ifile % 10 == 0 and ifile>0):
        print "UPDATING CACHE"
        with open(CACHE_FNAME,'wb') as f:
            pickle.dump([ N, S, processed_fnames], f)
        do_stats(N,S)

print "UPDATING CACHE"
with open(CACHE_FNAME,'wb') as f:
    pickle.dump([ N, S, processed_fnames], f)
do_stats(N,S)
