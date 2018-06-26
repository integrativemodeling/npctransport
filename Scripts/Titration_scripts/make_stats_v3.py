#!/usr/bin/python
from __future__ import print_function
import sys
import subprocess
import tempfile
import random
import re
import glob
from IMP.npctransport import *

N_FGMs=int(sys.argv[1])
Slide=int(sys.argv[2])
AVOGADRO=6.0221409E+23
A3_IN_LITERS=1E-27
FROM_TIME_NS=10000
global FG_CHAINS_MOLARITY
FG_CHAINS_MOLARITY= -1
RANDOM_FRACTION=1.0

def accumulate(dictionary, key, value, weight, custom_data):
    if key in dictionary:
        old_value= dictionary[key][0]
        old_weight= dictionary[key][1]
        new_weight= old_weight + weight
        new_value= (dictionary[key][0]*old_weight + value*weight) / new_weight
        dictionary[key][0]= new_value
        dictionary[key][1]= new_weight
        dictionary[key][2]= custom_data
    else:
        dictionary[key] = [value, weight, custom_data]

def get_indexed_fields(message, message_name="", indexed_fields={}):
    ''' return all fields in a message with an indexed value (index field),
        as a dictionary from nexted message name to value of indexed field '''
    for fd,value in (message.ListFields()):
        if fd.name=="index":
            indexed_fields[message_name]= message.value
        if fd.message_type is not None:
            if(message_name!="" and message_name[-1]!="."):
                message_name=message_name+"."
            if(fd.label == fd.LABEL_REPEATED):
                for i,single_value in enumerate(value):
                    indexed_fields= get_indexed_fields(single_value,
                                                       message_name+fd.name+str(i),
                                                       indexed_fields)
            else:
                indexed_fields= get_indexed_fields(value,
                                                   message_name+fd.name,
                                                   indexed_fields)
    return indexed_fields



def get_assignment_key(assignment):
    return tuple(get_indexed_fields(assignment).values())

def get_mean_energy_from_output(fname, from_time_ns):
    global FG_CHAINS_MOLARITY
    output= Output()
    with open(fname, "rb") as f:
        output.ParseFromString(f.read())
    indexed_fields= get_indexed_fields(output.assignment)
    work_unit= output.assignment.work_unit
    statistics= output.statistics
    s=0.0
    weight=0.0
    prev_time_ns=from_time_ns
    for gop in statistics.global_order_params:
        if gop.time_ns<from_time_ns:
            continue
        elapsed_time_ns= gop.time_ns-prev_time_ns
        s= s+gop.energy*elapsed_time_ns
        weight= weight+elapsed_time_ns
        prev_time_ns= gop.time_ns
    fgs_assignment= output.assignment.fgs[0]
    kaps_assignment= output.assignment.floaters[0]
    n_kaps= float(kaps_assignment.number.value)
    n_sites_kaps= n_kaps * kaps_assignment.interactions.value
    box_volume_L= output.assignment.box_side.value**3 * A3_IN_LITERS
    kaps_molarity= n_kaps / AVOGADRO / box_volume_L
    kap_sites_molarity= n_sites_kaps / AVOGADRO / box_volume_L
    n_fg_chains= float(fgs_assignment.number.value)
    fg_chains_molarity= n_fg_chains / AVOGADRO / box_volume_L
    if(FG_CHAINS_MOLARITY==-1):
        FG_CHAINS_MOLARITY= fg_chains_molarity
#    else:
#        assert(fg_chains_molarity==FG_CHAINS_MOLARITY) # assumed equal in all instances, at least for now
    n_sites_fgs = n_fg_chains * N_FGMs * fgs_assignment.interactions.value
    sites_ratio = n_sites_kaps / n_sites_fgs
    molar_ratio = n_kaps / n_fg_chains
    if weight>0.0:
        return  [s/weight, weight, [fg_chains_molarity, n_fg_chains], [kaps_molarity, n_kaps], indexed_fields]
    else:
        return [0.0, 0.0, [fg_chains_molarity, n_fg_chains], [kaps_molarity, n_kaps], indexed_fields]


def print_M2E_dictionary(M2E_dictionary, indexed_fields):
    row=0
    print("FG_CHAINS_MOLARITY %.3e" % FG_CHAINS_MOLARITY, "(base)")
    print("#", end='')
    for field_name in indexed_fields.keys():
        print(field_name, end=' ')
    print("sites_ratio E time_us M dE N Mt Mt_factor")
    for sites_ratio in sorted(M2E_dictionary.keys()):
        row= row+1
        E= M2E_dictionary[sites_ratio][0] # energy
        [M,N,Mt,Nt]= M2E_dictionary[sites_ratio][2] # molarity of NTF sites; nsites of NTRs
        # print sites_ratio, M,N,Mt,Nt
        Mt_factor=Mt/FG_CHAINS_MOLARITY
        time_us= M2E_dictionary[sites_ratio][1] / 1000.0 # time of stats
        if row==1:
            #        M0=M
            dE=float("nan")
            print("#", end='')
        else:
            #        print E,E1,N,N1
            Nt=float(Nt)
            Nt1=float(Nt1)
            dE = (E/Nt-E1/Nt1) / (N/Nt-N1/Nt1) #(M0/(M-M1)) # emergy per kap site in kcal/mol
        for field_value in indexed_fields.values():
            print("%.3f" % field_value, end=' ')
        print("%6.2e %8.1f %8.1f %6.3f %8.1f %d %6.2e %8.1f" % \
            (M, E, time_us, sites_ratio, dE, N, Mt, Mt_factor))
        E1= E
        M1= M
        N1= N
        Mt1= Mt
        Nt1= Nt
        Mt_factor1 = Mt_factor




####### MAIN ##########
M2E_dictionary={}
#for i in range(1,1921,1):
    #    if((round(i/6)+1)%2 <> Slide) or ((i-1)%6+1 <> N_FGMs)):
    #        continue
if(Slide==1):
    slide_string="True"
else:
    slide_string="False"
fname_match= glob.glob("config%d_Slide%s_ratio*/output*.pb" % (N_FGMs,slide_string))
for fname in fname_match:
    if(random.uniform(0,1)>RANDOM_FRACTION):
        continue
    print("#" + fname)
    try:
        [E_total, time_ns, [fg_chains_molarity, n_fg_chains], [kaps_molarity,n_kaps], indexed_fields]= \
            get_mean_energy_from_output(fname, FROM_TIME_NS)
        molar_ratio= kaps_molarity/fg_chains_molarity
        assignment_key= tuple(indexed_fields.values())
        if assignment_key not in M2E_dictionary.keys():
            M2E_dictionary[assignment_key]=[{}, indexed_fields]
        accumulate( M2E_dictionary[assignment_key][0],
                    key= molar_ratio,
                    value= E_total,
                    weight= time_ns,
                    custom_data= [kaps_molarity, n_kaps, fg_chains_molarity, n_fg_chains]
                    )
    except:
        print("#EXCEPTION", fname)
        raise
        pass

# Print it all:
for assignment_key, cur_M2E_dictionary in M2E_dictionary.iteritems():
    print_M2E_dictionary(cur_M2E_dictionary[0], cur_M2E_dictionary[1])
