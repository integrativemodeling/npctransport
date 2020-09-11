#!/bin/python
# Synposis: get statistics on spatial distribution of various molecules
# in the output HDF5 file from an fg_simulation simulation

from IMP.npctransport import *
import RMF
import argparse
import pickle
import glob
import gzip
import multiprocessing
import numpy as np
import os
import os.path
import random
import sys
try:
    import schwimmbad
    SCHWIMMBAD_OK=True
except ImportError:
    SCHWIMMBAD_OK=False
MPI_OK=False
if SCHWIMMBAD_OK:
    try:
        import schwimmbad.mpi
        schwimmbad.mpi.MPIPool()
        MPI_OK=True
        print("MPI is on")
    except ImportError:
        print("MPI is off (no MPI module available)")
    except ValueError as e:
        print("MPI is off (not activated as mpi)")


#IS_SKIP_FGS= False
IS_SKIP_FGS= True
N=1
DISABLE_RANDOM=True or (N==1)

# TODO: add stats time to pickle and to do_stats(), though will invalidate old caches

def do_stats(N,S):
    #    print("TOTAL STATS TIME [sec]: ", stats_time_ns*1E-9)
    if not os.path.isdir("Output"):
        assert(not os.path.exists("Output"))
        os.mkdir("Output")
    print("Writing stat files to Output/")
    for i in range(N):
        for f_type,XYZ in S[i].items():
            fname="Output/S{:d}.{:s}.txt".format(i, f_type)
#            print ("Fname {}".format(fname))
            with open(fname,'w') as F:
                for YZ in XYZ:
                    for Z in YZ:
                        for n_xyz in Z:
                            print("{:2d} ".format(n_xyz), end='', file=F)
                        print("", file=F)

def dataset_to_xyz(dataset):
    '''
    Convert a dataset from an npctransport HDF5
    output file into a 3-D np.array object '''
    XYZ_as_vector=np.array(dataset.get_block([0,0,0],dataset.get_size()))
    dim= int(round(XYZ_as_vector.shape[0]**(1/3.0)))
    assert(dim**3==XYZ_as_vector.shape[0])
    XYZ=np.reshape(XYZ_as_vector,[dim,dim,dim])
    return XYZ


#def handle_file(fname):
#    return True

def handle_file(fname):
    ''' Reads an npctransport HDF5 output file and load all the 3D
        matrices of it in a dictionary from type to a 3-D np.array
        object

        fname - file to be handles
    '''
    global IS_SKIP_FGS
    print("Handling {0}".format(fname))
    try:
        F= RMF.HDF5.open_file(fname)
        if not IS_SKIP_FGS:
            G_fgs= F.get_child_group("fg_xyz_hist")
        G_floaters= F.get_child_group("floater_xyz_hist")
    except:
        print("Skipping " + fname + " due to error")
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
    return (xyzs, fname)

def _sum_xyzs_exception(xyzs_and_fname):
#    print("Processing datasets returned by handle_file()")
    global N
    global S
    global DISABLE_RANDOM
    global processed_fnames # a set
    global CACHE_FNAME
    global MPI_ON
    print("_sum_xyzs_exception")
    if xyzs_and_fname is None:
        print("Empty xyzs skipped")
        return
    xyzs, fname= xyzs_and_fname
    cache_frequency=250
#    print("Checkpoint 1 in {}".format(fname))
    for i in range(N):
        for type_name,xyz in xyzs.items():
            if(random.random()<0.9 or DISABLE_RANDOM):
                if type_name in S[i]:
                    try:
                        S[i][type_name]=S[i][type_name]+xyz
                    except:
                        print("Couldn't add XYZ of shape" + xyz.shape)
                        print(" to S{:d}{:s} of shape".format(i, type_name))
                        print(S[i][type_name].shape)
                        continue
                else:
                    S[i][type_name]= xyz
#    print("Checkpoint 2 in {}".format(fname))
    processed_fnames.add(fname)
#    print("Checkpoint 3 in {}".format(fname))
#    print("Number of files processed {0:d}".format(len(processed_fnames)))
    if(len(processed_fnames) % cache_frequency == 0):
        do_stats(N,S)
        if CACHE_FNAME is not None:
            print("UPDATING CACHE")
            with open(CACHE_FNAME,'wb') as f:
                pickle.dump([ N, S, list(processed_fnames)], f)
    xyzs.clear() # Clear memory since no more use for xyzs - otherwise pool.map is trying to return it
    print("Done summarizing {} statistics".format(fname))

def sum_xyzs(xyzs_and_file):
    # Do not add any exceptions to this function - to prevent hung processes
    try:
        _sum_xyzs_exception(xyzs_and_file)
    except:
        if xyzs_and_file is not None:
            print("Exception summarizing stats of {0}".format(xyzs_and_file[1]))
        else:
            print("Unknown exception in sum_xyzs(None)")

def get_cmdline_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('input_prefixes', metavar='input file prefixes', type=str, nargs='+',
                        help='List of prefixes of input hdf5 file from NPC simulations')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=MPI_OK,
                       action="store_true", help="Run with MPI (default is false, unless using mpirun/exec with schimmbad module installed, which forces it to be true")
    args = parser.parse_args()
    return args


def run_schwimmbad_pool(filenames, is_mpi, n_processes):
    print("Running with Scwimmbad pool")
    pool= schwimmbad.choose_pool(mpi =is_mpi,
                                 processes= n_processes)
    #    pool= schwimmbad.mpi.MPIPool()
    #    if is_mpi and not pool.is_master():
    #       pool.wait(lambda: sys.exit(0))
    results= pool.map(handle_file,
                      filenames,
                      callback=sum_xyzs)
    pool.close()
    return results

def run_multiprocessing_pool(filenames, n_processes):
    print("Running with Multiprocessing pool")
    with multiprocessing.Pool(processes= n_processes) as pool:
        R= []
        #    manager= multiprocessing.Manager()
        for fname in fnames:
            try:
                print("Processing", fname)
                res= pool.apply_async(handle_file,
                                      args=(fname,), #, processed_fnames_managed),
                                 callback=sum_xyzs)
                R.append(res)
            except:
                print("Failed on {0}".format(fname))
            pool.close()
            pool.join()
            for i,r in enumerate(R):
                print(i)
                print(r.ready())
                print(r.successful())

def get_input_fnames_set_from_prefixes(input_prefixes):
    ''' return value is set '''
    input_fnames= []
    for input_prefix in input_prefixes:
        input_fnames.extend(glob.glob(input_prefix + "*.hdf5"))
    return set(input_fnames)


############# Main ############
if __name__ == '__main__':
    if IS_SKIP_FGS:
        print("Skipping FGs")
    
    if SCHWIMMBAD_OK:
        args= get_cmdline_args()
        if(args.mpi):
            assert(MPI_OK)
            #    CACHE_FNAME= None
        input_prefixes= args.input_prefixes
    else:
        input_prefixes= sys.argv[1:]
    CACHE_FNAME='Output/CACHE.get_float_stats_{}.p'.format('__'.join(input_prefixes).replace('/','_'))
    fnames= get_input_fnames_set_from_prefixes(input_prefixes)
    print("Processing {} files with prefixes: {}".format(len(fnames),
                                                         ", ".join(input_prefixes)))
    try:
        #    with gzip.open(CACHE_FNAME+".gzip",'rb') as f:
        if CACHE_FNAME is None:
            assert(False)
        with open(CACHE_FNAME,'rb') as f:
            [ n, S, processed_fnames ] = pickle.load(f)
            processed_fnames= set(processed_fnames)
            assert(fnames.issuperset(processed_fnames))
            assert(n==N)
            print("Using cache file {} with {} pre-processed files".format(CACHE_FNAME,
                                                                           len(processed_fnames)))
    except AssertionError as e:
        print(e)
        raise
    except:
        print("NOT USING CACHE")
        S=[{} for i in range(N)]
        processed_fnames= set()
    stats_time_ns=0.0
    fnames=fnames.difference(processed_fnames)
    print("Fnames")
    print(fnames)
#    for fname in fnames:
#        handle_file(fname)
    if SCHWIMMBAD_OK:
        print("Schwimmbad:")
        results= run_schwimmbad_pool(fnames,
                                     is_mpi= args.mpi,
                                     n_processes= args.n_cores)
    else:
        print("Multiprocessing:")
        results= run_multiprocessing_pool(fnames,
                                          n_processes= 8)
    do_stats(N,S)
    if CACHE_FNAME is not None:
        print("UPDATING CACHE")
        with open(CACHE_FNAME,'wb') as f:
            pickle.dump([ N, S, list(processed_fnames)], f)
