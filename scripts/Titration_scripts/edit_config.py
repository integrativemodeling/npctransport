#!/usr/bin/env python
from IMP.npctransport import *
import sys

def do_edit(in_fname, out_fname,
            kap_k=None, kap_range=None):
    '''
    in_fname - input filename
    out_fname - output filename
    kap_k - new value for kap-fg interaction_k coefficient (if None, ignore)
    kap_range - new value for kap-fg interaction_range coefficient (if None, ignore)
   '''
    with open(in_fname, "rb") as fin:
        config= Configuration()
        config.ParseFromString(fin.read())
    config.time_step_factor.lower=2.0
    config.output_statistics_interval_ns=250.0
    for interaction in config.interactions:
        if (interaction.type0=='kap20' and interaction.type1=='fg0' or
            interaction.type0=='fg0' and interaction.type1=='kap20'):
            if kap_k is not None:
                interaction.interaction_k.lower= kap_k
            if kap_range is not None:
                interaction.interaction_range.lower= kap_range
            break
    with open(out_fname, "wb") as fout:
        fout.write(config.SerializeToString())
        fout.close()

if __name__ == "__main__":
    do_edit(*sys.argv)
