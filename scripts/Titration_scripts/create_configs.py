#!/usr/bin/python
# Configuration file defined in a python script
import numpy as np
import grid_params
import subprocess
import math
import os.path

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

def make_cfg(KD_M, fg_string, kap_valency, C_kap_M, C_fg_M,
             QuadraticLogKDParams):
    range_A=5.5
    kap_ks= compute_k_from_r_KD_and_quadratic_model(range_A, KD_M,
                                                QuadraticLogKDParams,
                                                verbose=False)
    if(len(kap_ks)==0):
        print("Error: no k was found to obtain KD of {} [M] with range {} [A]".format(KD_M, range_A))
        return
    kap_k=Ks[0]
    if len(kap_ks)>1:
        print("Warning: more than one k results in KD of {} [M] with range {} - using first of {}" \
              .format(KD_M, range_A, kap_ks))
    label= "KapRange_{:.4f}__KapK{:.4f}__FG_{}__KapValency_{:d}__CKapM_{:.4e}__C_FGM{:.4e}"\
              .format(range_A,
                      kap_k,
                      fg_string,
                      kap_valency,
                      C_kap_M,
                      C_fg_M)
    cfg_pb_file= "Configs/config_{}.pb".format(label)
    cfg_txt_file= "Configs/config_{}.txt".format(label)
    print("Making {}".format(cfg_pb_file))
    cmd_template= "python MyScripts/make_cfg_v15.py -k {C_kap_M} -f {C_fg_M}" +\
        " --kap_valency {kap_valency} -kap_k {kap_k} {cfg_file} {fg_string} T" +\
        " --time_step_factor {timestep}"
    cmd=cmd_template.format(C_kap_M= C_kap_M,
                            C_fg_M= C_fg_M,
                            kap_valency= kap_valency,
                            kap_k= kap_k,
                            cfg_file= cfg_pb_file,
                            fg_string= fg_string,
                            timestep= 1.0)
    if os.path.exists(cfg_pb_file):
        print("Skipping existing file {}".format(cfg_pb_file))
        return
    with open(cfg_txt_file,"w") as f_cfg_txt_file:
        try:
            subprocess.check_call(cmd, shell=True, stdout=f_cfg_txt_file)
        except:
            print("ERROR: Couldn't create {}".format(cfg_pb_file))

def main():
    params=grid_params.GridParams()
    for KD_M in params.KD_mono_M:
        for fg_valency in params.FG_valencies:
            fg_string= 'F'*fg_valency+'S'*(6-fg_valency)
            for kap_valency in params.NTF2_valencies:
                for C_guest_M in params.Guest_M:
                    for C_host_M in params.Host_M:
                        for is_forward in params.Is_forward:
                            if is_forward:
                                C_kap_M= C_guest_M
                                C_fg_M= C_host_M
                            else:
                                C_fg_M= C_guest_M
                                C_kap_M= C_host_M
                            print("{} {}".format(C_kap_M, C_fg_M))
                            make_cfg(KD_M, fg_string, kap_valency, C_kap_M, C_fg_M,
                                     params.QuadraticiLogKDParams)

if __name__ == "__main__":
    main()
