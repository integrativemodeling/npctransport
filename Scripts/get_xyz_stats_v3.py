
from IMP.npctransport import *
#import seaborn as sbn
import numpy as np
#import pandas as pd
import sys
import random
import cPickle as pickle
import RMF
import gzip

N=1
DISABLE_RANDOM=True or (N==1)

# TODO: add stats time to pickle and to do_stats(), though will invalidate old caches

def do_stats(N,S):
#    print "TOTAL STATS TIME [sec]: ", stats_time_ns*1E-9
    print "Writing stat files"
    for i in range(N):
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
#    with gzip.open(CACHE_FNAME+".gzip",'rb') as f:
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
        F= RMF.HDF5.open_file(fname)
        G_fgs= F.get_child_group("fg_xyz_hist")
        G_floaters= F.get_child_group("floater_xyz_hist")
    except:
        print "Skipping ", fname, " due to error"
        continue
    print "Processing", fname
    datasets={}
    n_fgs= G_fgs.get_number_of_children()
    for j in range(n_fgs):
        fg_type= G_fgs.get_child_name(j)
        datasets[fg_type]= G_fgs.get_child_int_data_set_3d(fg_type)
    n_floaters= G_floaters.get_number_of_children()
    for j in range(n_floaters):
        floater_type= G_floaters.get_child_name(j)
        datasets[floater_type]= G_floaters.get_child_int_data_set_3d(floater_type)
    for i in range(N):
        for type_name,dataset in datasets.iteritems():
            XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
            print XYZ_as_vector.shape
            dim= int(round(XYZ_as_vector.shape[0]**(1/3.0)))
            print dim
            assert(dim**3==XYZ_as_vector.shape[0])
            XYZ=np.reshape(XYZ_as_vector,[dim,dim,dim])
            print XYZ.shape
            if(random.random()<0.9 or DISABLE_RANDOM):
                if type_name in S[i]:
                    try:
                        S[i][type_name]=S[i][type_name]+XYZ
                    except:
                        print "Couldn't add XYZ of shape", XYZ.shape,
                        print " to S[%d][%s] of shape" % (i, type_name),
                        print S[i][type_name].shape
                        continue
                else:
                    S[i][type_name]=XYZ
    processed_fnames.add(fname)
    print ifile
    if(ifile % 20 == 0 and ifile>0):
        do_stats(N,S)
        print "UPDATING CACHE ifile=%d" % ifile
        with open(CACHE_FNAME,'wb') as f:
            pickle.dump([ N, S, processed_fnames], f)

do_stats(N,S)
print "UPDATING CACHE"
#with gzip.GzipFile(CACHE_FNAME+'.gzip', 'wb') as f:
with open(CACHE_FNAME,'wb') as f:
    pickle.dump([ N, S, processed_fnames], f)
