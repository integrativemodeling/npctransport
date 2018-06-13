#!/bin/python
from IMP.npctransport import *
#import seaborn as sbn
import numpy as np
#import pandas as pd
import sys
import random
import cPickle as pickle
import RMF
import gzip
import multiprocessing

#IS_SKIP_FGS= False
IS_SKIP_FGS= False
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

def dataset_to_xyz(dataset):
    '''
    Convert a dataset from an npctransport HDF5
    output file into a 3-D np.array object '''
    XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
    dim= int(round(XYZ_as_vector.shape[0]**(1/3.0)))
    assert(dim**3==XYZ_as_vector.shape[0])
    XYZ=np.reshape(XYZ_as_vector,[dim,dim,dim])
    return XYZ


def handle_file(fname, processed_fnames):
    ''' Reads an npctransport HDF5 output file and load all the 3D
        matrices of it in a dictionary from type to a 3-D np.array
        object

        fname - file to be handles
        processed_fnames - a shared managed list of all filenames that have been processed
    '''
    global IS_SKIP_FGS
    print("Handling {0}".format(fname))
    try:
        F= RMF.HDF5.open_file(fname)
        if not IS_SKIP_FGS:
            G_fgs= F.get_child_group("fg_xyz_hist")
        G_floaters= F.get_child_group("floater_xyz_hist")
    except:
        print "Skipping ", fname, " due to error"
        return None
    xyzs={}
    if not IS_SKIP_FGS:
        n_fgs= G_fgs.get_number_of_children()
        for j in range(n_fgs):
            fg_type= G_fgs.get_child_name(j)
            dataset= G_fgs.get_child_int_data_set_3d(fg_type)
            xyzs[fg_type]= dataset_to_xyz(dataset)
    n_floaters= G_floaters.get_number_of_children()
    for j in range(n_floaters):
        floater_type= G_floaters.get_child_name(j)
        dataset= G_floaters.get_child_int_data_set_3d(floater_type)
        xyzs[floater_type]= dataset_to_xyz(dataset)
    print("Done handling {0}".format(fname))
    processed_fnames.append(fname)
    return xyzs

def _sum_xyzs_exception(xyzs):
    print("Processing datasets returned by handle_file()")
    global N
    global S
    global DISABLE_RANDOM
    global processed_fnames
    global CACHE_FNAME
    if xyzs is None:
        print("Empty xyzs skipped")
        return
    cache_frequency=500
    for i in range(N):
        for type_name,xyz in xyzs.iteritems():
            if(random.random()<0.9 or DISABLE_RANDOM):
                if type_name in S[i]:
                    try:
                        S[i][type_name]=S[i][type_name]+xyz
                    except:
                        print "Couldn't add XYZ of shape", xyz.shape,
                        print " to S[%d][%s] of shape" % (i, type_name),
                        print S[i][type_name].shape
                        continue
                else:
                    S[i][type_name]= xyz
    print("Number of files processed {0:d}".format(len(processed_fnames)))
    if(len(processed_fnames) % cache_frequency == 0):
        do_stats(N,S)
        if CACHE_FNAME is not None:
            print "UPDATING CACHE"
            with open(CACHE_FNAME,'wb') as f:
                pickle.dump([ N, S, list(processed_fnames)], f)

def sum_xyzs(xyzs):
    try:
        _sum_xyzs_exception(xyzs)
    except:
        print("Unknown exception in sum_xyzs()")
        raise


############# Main ############
if __name__ == '__main__':
    #import cPickle as p; D=p.load(open('TMP.get_float_stats_cache.p','rb')); print len(D[2]);
    CACHE_FNAME='TMP.get_float_stats_cache.p'
#    CACHE_FNAME= None
    fnames=set(sys.argv[1:])
    try:
        #    with gzip.open(CACHE_FNAME+".gzip",'rb') as f:
        if CACHE_FNAME is None:
            assert(False)
        with open(CACHE_FNAME,'rb') as f:
            [ n, S, processed_fnames ] = pickle.load(f)
            assert(fnames.issuperset(processed_fnames))
            assert(n==N)
            print("Using cache", CACHE_FNAME)
    except:
        print "NOT USING CACHE"
        S=[{} for i in range(N)]
        processed_fnames=set()
    stats_time_ns=0.0
    fnames=fnames.difference(processed_fnames)
    pool= multiprocessing.Pool(processes=8)
    manager= multiprocessing.Manager()
    processed_fnames= manager.list(processed_fnames)
    print("Starting pool")
    for fname in fnames:
        try:
            print(fname)
            pool.apply_async(handle_file,
                             args=(fname, processed_fnames),
                             callback=sum_xyzs)
        except:
            print("Failed on {0}.format(fname)")
    # for r in R:
    #     r.wait()
    #     r.get()
#    pool.wait()
    pool.close()
    pool.join()
    do_stats(N,S)
    if CACHE_FNAME is not None:
        print "UPDATING CACHE"
        with open(CACHE_FNAME,'wb') as f:
            pickle.dump([ N, S, list(processed_fnames)], f)
