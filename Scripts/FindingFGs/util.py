#!/usr/bin/python
import re
import os
import Bio.SeqUtils.ProtParamData
import numpy as np

kd= Bio.SeqUtils.ProtParamData.kd # Kyte-Doolitle hydrophobicity score
min_kd= min(kd.values())
max_kd= max(kd.values())
kd_normalized= { k:(v-min_kd)/(max_kd-min_kd) for k,v in kd.items() }
kd_normalized['X']= 0.0 # rare cases...

def is_fg(seq):
    ''' return true if sequence of amino acid is considerd to be an FG repeat '''
    return re.search("FG.{0,40}?FG.{0,40}?FG.{0,40}?FG.{0,40}?FG.{0,40}?FG.{0,40}?FG", seq)
#    return re.search("DE.*DE.*DE.*DE.*DE.*DE.*DE.*DE", seq)

def remove_whitespace(seq):
    return re.sub("[" + os.linesep + "\r\n]", "", seq)

def get_charge_at_pH_7(aa):
    if aa in ['K','R']:
        return +1.00
    if aa in ['D','E']:
        return -1.00
    if aa=='H':
        return +0.25
    return 0.00

def get_proper_moving_average(vector, w):
    ''' 
    Get moving average of vector using only properly averaged entries for
    window size w.

    Returns a length len(vector)+1-w vector
    '''
    cumsum, moving_avgs = [0], []
    for i, x in enumerate(vector, 1):
        cumsum.append(cumsum[i-1] + x)
        if i>=w:
            moving_avg = (cumsum[i] - cumsum[i-w])/w
            #can do stuff with moving_ave here
            moving_avgs.append(moving_avg)
    return moving_avgs
 

def get_fold_index(seq, w=40):
    '''
    Returns the FoldIndex value for amino acid sequence seq for window size w 
    See https://doi.org/10.1093/bioinformatics/bti537 Priluski et al 2005
    '''
    global kd_normalized
    h= [kd_normalized[aa] for aa in seq] # hydrophobicity
    c= [get_charge_at_pH_7(aa) for aa in seq] # charge
#    print(h)
#    print(c)
    hma= np.array(get_proper_moving_average(h, w))
    cma= np.array(get_proper_moving_average(c, w))
#    print(hma)
#    print(cma)
    net_cma= np.abs(cma)
    fi= 2.785*hma - net_cma - 1.151
    return fi
    
