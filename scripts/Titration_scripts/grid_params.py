#!/usr/bin/python
# Configuration file defined in a python script
import numpy as np

class GridParams:
    def __init__(self):
        # Valencies
        self.FG_valencies=[1, 2, 3, 4, 6] # added 6 to Ryo's
        self.NTF2_valencies=[1, 2, 4]
        # For Guest, Ryo used 10000 points between 1nM ~ 10M in
        # logspace (numpy.logspace(-9 1,  10000))
        self.Guest_M= np.logspace(-9,1,500) # using less sample points to make runs feasible
        #
        self.Host_M=np.array([1, 10, 100, 1000])*1.0E-6 # from uM to M
        # KD monos = per-site affinities
        self.KD_mono_M=np.array([40,20,10,5,2.5,1.25])*1E-3 # from mM to M
        # Directions (forward - host=FG and guest=NTF2, reverse - vice versa)
        self.Is_forward=[True, False]
        # Quadratic model for range,k to kD, where other parameters are fixed
        # as in make_cfg_v14/15.py such that for Q being a shortcut for QuadraticLogKDParam
        # log10(kD) = Q[0] + Q[1]*range + Q[2]*k + Q[3]*(range*k) + Q[4]*range^2 + Q[5]*k^2
        self.QuadraticLogKDParams = [ -1.32828572e+00,
                                     5.55793607e-01,
                                     3.36939985e-01,
                                     -1.42452838e-01,
                                     -6.42247161e-02,
                                     -5.76784272e-04]


def compute_k_from_r_KD_and_quadratic_model(r, KD_M, Q, verbose=True):
    '''
    compute k (force coefficient) to obtain desired KD for a given value of r (the range),
    based on quadratic model Q in which:
    log10(KD) = Q[0] + Q[1]*r + Q[2]*k + Q[3]*(rk) + Q[4]*r^2 + Q[5]*k^2
    '''
    a= Q[5]
    b= Q[2] + Q[3]*r
    c= Q[0] + Q[1]*r + Q[4]*(r**2) - np.log10(KD_M) # -log10(KD) to ensure solution equals desired KD
    if verbose:
        print("Solving: {:.5f}k^2 + {:.5f}k + {:.5f} = log10({:.5f})" \
              .format(a, b, c, KD_M))
    # calculate the discriminant
    d = (b**2) - (4*a*c)
    if d<0:
        return None
    # find two solutions
    sol1 = (-b-math.sqrt(d))/(2*a)
    sol2 = (-b+math.sqrt(d))/(2*a)
    if verbose:
        print('The solution are {0} and {1}'.format(sol1,sol2))
    solutions=[]
    if(sol1>0):
        solutions.append(sol1)
    if(sol2>0):
        solutions.append(sol2)
    return solutions

def compute_KD_from_k_r_and_quadratic_model(k, r, Q):
    '''
    compute KD for a given value of  k (force coefficient) and r (the range),
    based on quadratic model Q in which:
    log10(KD) = Q[0] + Q[1]*r + Q[2]*k + Q[3]*(rk) + Q[4]*r^2 + Q[5]*k^2
    '''
    logKD= Q[0] + Q[1]*r + Q[2]*k + Q[3]*(r*k) + Q[4]*r**2 + Q[5]*k**2
    return 10**logKD
