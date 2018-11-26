#!/usr/bin/env python
from IMP.npctransport import *
import sys
import glob
import re
import math
import scipy.stats as stats
import cPickle as pickle
import numpy as np
import os
import traceback

####### Global params ###
CONF= 0.95 # confidence interval required
AVOGADRO=6.0221409E+23 # molecules/mole
L_per_A3=1E-27 # Liters per A^3
kB_kcal_per_mol_K=0.0019872041
TEMP_MIN_DS=273
TEMP_MAX_DS=335
EPSILON=1e-9
IS_NEW_SITE_STATS= False



def pretty_molarity(molarity):
    ''' Converts molarity in units of [M] to a pretty string, in either M, mM, uM, nM or picoM '''
    if molarity > 5:
        return "{:.1f} M".format(molarity)
    if molarity > 0.5:
        return "{:.2f} M".format(molarity)
    if molarity > 5E-3:
        return "{:.1f} mM".format(molarity*1E+3)
    if molarity > 0.5E-3:
        return "{:.2f} mM".format(molarity*1E+3)
    if molarity > 5E-6:
        return "{:.1f} uM".format(molarity*1E+6)
    if molarity > 0.5E-6:
        return "{:.2f} uM".format(molarity*1E+6)
    if molarity > 5E-9:
        return "{:.1f} nM".format(molarity*1E+9)
    if molarity > 0.5E-9:
        return "{:.2f} nM".format(molarity*1E+9)
    if molarity > 5E-12:
        return "{:.1f} picoM".format(molarity*1E+12)
    if molarity > 0.5E-12:
        return "{:.2f} picoM".format(molarity*1E+12)
    return "{:.2e} M".format(molarity)


def get_indexed_fields(message, message_name="", indexed_fields={}):
    ''' return all fields in a protobuf message with an indexed value (index field),
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

def accumulate(dictionary, key, value, time_ns):
    new_value=[value*time_ns, time_ns, value**2*time_ns, 1] # mean / time / mean-square / n
    if key in dictionary:
        dictionary[key] = [x+y for x,y in zip (dictionary[key],new_value)]
    else:
        dictionary[key] = new_value

def accumulate_check(dictionary, key, value, time_ns):
    ''' accumulate after checking value is not inf or nan,
        return true if added
    '''
    check_ok = (value != float('inf') and not math.isnan(value))
    if check_ok:
        accumulate(dictionary, key, value, time_ns)
    return check_ok


def get_output_num(fname):
    m = re.search('output[_a-zA-Z]*([0-9]*)\.pb', fname)
    if m is None: raise ValueError('cannot parse %s' % fname)
    return int(m.groups(0)[0])

# returns number of sems for symmetric confidence intravel conf_interval
def get_sems_for_conf(conf_interval):
    a=0.5+0.5*conf_interval
    return stats.norm.ppf(a)

def do_stats(table, units_label, is_fraction=False, prefix="", verbose=True):
    ''' return dictionary with mean and standard-error of mean for each key in table'''
    # print transports per particle per second
    keys=sorted(table.keys())
    ret={}
    for key in keys:
        t_ns= table[key][1]
        mean= table[key][0]/(t_ns+0.0001)
        mean2= table[key][2]/(t_ns+0.0001)
        n= table[key][3]
        stddev= math.sqrt(mean2-mean**2)
        stderr= stddev/math.sqrt(n+0.0001) # dividing by n is probably ok when constant stat intervals are used
        conf=get_sems_for_conf(CONF)
        if verbose:
            if(prefix <> ""):
                print prefix,
            if(is_fraction):
                print key, 'mean: %10.2f%% +- %10.2f%% [%s] ; t=%.1f [ns]' % \
                    (mean*100, stderr*100*conf, units_label, t_ns)
            else:
                print key, 'mean: %10.8f +- %10.8f [%s] ; t=%.1f [ns]' % \
                    (mean, stderr*conf, units_label, t_ns)
        ret[key]=(mean, stderr, t_ns)
    return ret

def do_stats_stddev(table, table_squared,
                    units_label,
                    is_fraction=False,
                    prefix="",
                    verbose= True):
    '''
    return standard deviation around mean for each key in table,
    based on table of measurements and table of square measurements
    using var(X)=E(X^2)-E(X)^2,
    table - dictionary of keys and values
    table_square - dictionary of keys and squared values
    is_fraction - are values a fraction between 0 and 1 or not
    prefix - prefix to print before output if not ""
    '''
    # print transports per particle per second
    keys=sorted(table.keys())
    ret={}
    for key in keys:
        t_ns= table[key][1]
        mean= table[key][0]/(t_ns+0.0001)
        mean2= table_squared[key][0]/(t_ns+0.0001)
        stddev= math.sqrt(mean2-mean**2)
        if verbose:
            if(prefix <> ""):
                print prefix,
            if(is_fraction):
                print key, 'std-dev: %10.2f%% [%s] ; t=%.1f [ns]' % \
                    (stddev*100, units_label, t_ns)
            else:
                print key, 'std-dev: %10.8f [%s] ; t=%.1f [ns]' % \
                    (stddev, units_label, t_ns)
        ret[key]=(stddev, t_ns)
    return ret


def get_k_low_high(k_stats):
    [k,_,t]= k_stats
#    print k,t
    if t==0:
        return [-1,-1,-1]
    n= k * t
    stderr_n= math.sqrt(n)
    conf= get_sems_for_conf(CONF)
    n_high=  n + conf*stderr_n
    n_low=   n - conf*stderr_n
    k_high=  n_high / t
    k_low=   n_low / t
    return [k, k_low, k_high]

def get_temperature_from_key(key, key_suffix=""):
    ''' if key_suffix is supplied - return match only if it is included after temperature '''
    m=re.match("(?P<temperature_K>[0-9\.]+) " + key_suffix, key)
    if not m:
        return None
    temperature_K=float(m.group('temperature_K'))
    return temperature_K

def get_dH_and_dS(KDs_dict, key_suffix):
    x=[]
    y=[]
    for key,KDs in KDs_dict.iteritems():
        temperature_K= get_temperature_from_key(key, key_suffix)
        if not temperature_K:
            continue
        if temperature_K<TEMP_MIN_DS or temperature_K>TEMP_MAX_DS:
            continue
        x.append(-1.0/temperature_K/kB_kcal_per_mol_K)
        [KD, KD_low, KD_high] = KDs
        y.append(-math.log(KD)) # natural log of kA
    # ln(kA) = dH*(-1.0/(kB*T)) + dS/kB <==> dH-dS*T = -kB*T*ln(kA) <==> dG = -kB*T*ln(kA)
    if(len(x)<=1):
        return [float('nan'),float('nan')]
    [dH,dS_per_kB]=np.polyfit(x,y,1)
    dS=dS_per_kB*kB_kcal_per_mol_K
    return [dH,dS]

def get_dH_and_dS2(KDs_dict, key_suffix):
    x=[]
    y=[]
    for key,KDs in KDs_dict.iteritems():
        m=re.match("(?P<temperature_K>[0-9\.]+) " + key_suffix, key)
        if not m:
            continue
        temperature_K=float(m.group('temperature_K'))
        if temperature_K<TEMP_MIN_DS or temperature_K>TEMP_MAX_DS:
            continue
        x.append(temperature_K)
        [KD, KD_low, KD_high] = KDs
        y.append(math.log(KD)*kB_kcal_per_mol_K*temperature_K)
    # ln(kA) = dH*(-1.0/(kB*T)) + dS/kB <==> dH-dS*T = -kB*T*ln(kA) <==> dG = -kB*T*ln(kA)
    if len(x)<=1:
        return [float('nan'),float('nan')]
    [minus_dS,dH]=np.polyfit(x,y,1)
    dS= -minus_dS
    return [dH,dS]




############# Main ############

def do_all_stats(fnames, STATS_FROM_SEC, verbose=True, return_outputs=None):
    #TEMPERATURE_K_INDEX=int(sys.argv[3])
    N_KAPS=None # Will be updated from assignment (should be consistent among inputs)
    N_SITES_PER_KAP=None # Will be updated from assignment (should be consistent among inputs)
    N_SITES_PER_FG=None # Will be updated from assignment (should be consistent among inputs)
    N_FGS=None # Will be updated from assignment (should be consistent among inputs)

    STATS_FROM_NS= STATS_FROM_SEC*(1e+9)
    k_ons_i= {}
    k_offs_i= {}
    k_ons_ii= {}
    k_offs_ii= {}
    k_ons_ss= {} # mc = missing contact
    k_offs_ss= {}
    KDs_chains= {}
    fbounds1_floats= {}
    fbounds2_floats= {}
    fbounds2= {}
    fbounds_sites1_old= {}
    fbounds_sites2_old= {}
    fbounds_sites1_new= {}
    fbounds_sites2_new= {}
    bound_times= {}
    unbound_times= {}
    mean_rg= {}
    mean_rg2= {}
    mean_dmax= {}
    mean_dmax2= {}
    mean_bond_rest_length= {}
    mean_bond_rest_length2= {}
    energy= {}
    n= 0
    total_sim_time_sec=0.0
    #modulus = int(sys.argv[2])
    indexed_fields=None # indexed fields that need to be the same for all output files
    if return_outputs is not None and return_outputs==True:
        return_outputs= []
    else:
        return_outputs= None
    for fname in fnames:
        if(not os.path.isfile(fname)):
            print("Skipping {} - not exists or not a regular file".format(fname))
            continue
        if(os.stat(fname).st_size == 0):
            print("Warning: skipping {} - empty file".format(fname))
            continue
        output= IMP.npctransport.Output()
        try:
            f=open(fname, "rb")
            output.ParseFromString(f.read())
        except ValueError as e:
            print "Warning: EXCEPTION ", fname, e
            continue
        except:
            print "Warning: EXCEPTION ", fname
            continue
        if return_outputs is not None:
            return_outputs.append(output)
        if(len(output.statistics.global_order_params)<2):
            print("Warning: skipping file {} - does not contain proper statistics".format(fname))
            continue
        cur_indexed_fields= get_indexed_fields(output.assignment)
        cur_indexed_fields.pop("temperature", None) # make sure temperature is not included in indexed fields
        if indexed_fields is None:
            indexed_fields= cur_indexed_fields
        else:
            assert(indexed_fields.values()==cur_indexed_fields.values()) # TODO: perhaps except temperature
        temperature_k= output.assignment.temperature_k.value
    #        if(output.assignment.temperature_k.index <> TEMPERATURE_K_INDEX):
    #            continue
        start_sim_time_sec= output.statistics.global_order_params[0].time_ns*(1e-9)
        end_sim_time_sec= output.statistics.global_order_params[-1].time_ns*(1e-9)
        sim_time_sec= end_sim_time_sec-max(start_sim_time_sec,STATS_FROM_SEC)
        box_volume_L= output.assignment.box_side.value**3 * L_per_A3
    #        print "TIME", sim_time_sec
        if sim_time_sec<0:
            print "SKIP", fname, " ", sim_time_sec
            continue
        total_sim_time_sec= total_sim_time_sec+sim_time_sec
        # for floater in output.statistics.floaters:
        #     ftype = floater.type
        #     was_bound=floater.order_params[0].interacting_fraction>0.0
        #     n_on=0
        #     n_off=0
        #     time_state=0.0
        #     for fop in floater.order_params:
        #         if fop.time_ns<STATS_FROM_NS:
        #             continue
        #         fbound = fop.interacting_fraction
        #         accumulate(fbounds, ftype, fbound, 1)
        #         is_bound=fbound>0.0
        #         if(is_bound != was_bound):
        #             if is_bound:
        #                 accumulate(unbound_times, ftype, time_state, 1) # finished unbound
        #                 n_on=n_on+1
        #             else:
        #                 accumulate(bound_times, ftype, time_state, 1) # finished bound
        #                 n_off=n_off+1
        #             was_bound=is_bound
        #             time_state=0.0
        #         else:
       #             time_state=time_state+1
        ii=0
        assert(len(output.assignment.fgs)==1)
        assert(len(output.assignment.floaters)==1)
        for fg in output.assignment.fgs:
            if N_FGS is None:
                N_FGS= fg.number.value
                N_SITES_PER_FG= fg.interactions.value
            else:
                assert(N_KAPS==floater.number)
                assert(N_SITES_PER_KAP==floater.interactions)
        for floater in output.assignment.floaters:
            if N_KAPS is None:
                N_KAPS= floater.number.value
                N_SITES_PER_KAP= floater.interactions.value
            else:
                assert(N_KAPS==floater.number)
                assert(N_SITES_PER_KAP==floater.interactions)
        for floater in output.statistics.floaters:
            iname = "%.2f fg0 - %s" % (temperature_k, floater.type)
            prev_time_ns=floater.order_params[0].time_ns
            for fop in floater.order_params:
                elapsed_time_ns= fop.time_ns - prev_time_ns
                prev_time_ns= fop.time_ns
                if fop.time_ns<STATS_FROM_NS:
                    continue
                fb2= fop.interacting_fraction
                if(math.isnan(fb2)):
                    continue
                nb2 = fb2 * N_KAPS
                if nb2>0:
                    fb1= nb2 * fop.chains_per_interacting_floater / (N_FGS+0.0)
                else:
                    fb1=0.0
                accumulate_check(fbounds1_floats, iname, fb1, elapsed_time_ns)
                accumulate_check(fbounds2_floats, iname, fb2, elapsed_time_ns)
    #           print "KD@%.0f ns is %.1e" % (fop.time_ns, KD)

        for interaction in output.statistics.interactions:
            iname = "%.2f %s - %s" % (temperature_k , interaction.type0, interaction.type1 )
            for iop in interaction.order_params:
                if iop.time_ns<STATS_FROM_NS:
                    continue
                k_on_i=  iop.avg_on_per_unbound_i_per_ns
                t_on_i_ns= iop.on_i_stats_period_ns
                k_off_i= iop.avg_off_per_bound_i_per_ns
                t_off_i_ns=iop.off_i_stats_period_ns
                k_on_ii=  iop.avg_on_per_unbound_ii_per_ns
                t_on_ii_ns= iop.on_ii_stats_period_ns
                k_off_ii= iop.avg_off_per_bound_ii_per_ns
                t_off_ii_ns=iop.off_ii_stats_period_ns
                k_on_ss= iop.avg_on_per_missing_contact_per_ns
                t_on_ss_ns= iop.on_stats_period_ns
                k_off_ss= iop.avg_off_per_contact_per_ns
                t_off_ss_ns= iop.off_stats_period_ns
                accumulate_check(k_ons_i, iname, k_on_i, t_on_i_ns)
                accumulate_check(k_offs_i, iname, k_off_i, t_off_i_ns)
                accumulate_check(k_ons_ii, iname, k_on_ii, t_on_ii_ns)
                accumulate_check(k_offs_ii, iname, k_off_ii, t_off_ii_ns)
                accumulate_check(k_ons_ss, iname, k_on_ss, t_on_ss_ns)
                accumulate_check(k_offs_ss, iname, k_off_ss, t_off_ss_ns)
                fbound2= iop.avg_fraction_bound_particles_ii
                accumulate_check(fbounds2, iname, fbound2, 1)
                fbound_sites1_old= iop.avg_contacts_per_particle_i/N_SITES_PER_FG
                accumulate_check(fbounds_sites1_old, iname, fbound_sites1_old, 1)
                fbound_sites2_old= iop.avg_contacts_per_particle_ii/N_SITES_PER_KAP
                accumulate_check(fbounds_sites2_old, iname, fbound_sites2_old, 1)
                if iop.HasField('avg_fraction_bound_particle_sites_i'):
                    IS_NEW_SITE_STATS= True
                    fbound_sites1_new= iop.avg_fraction_bound_particle_sites_i
                    accumulate_check(fbounds_sites1_new, iname, fbound_sites1_new, 1)
                    fbound_sites2_new= iop.avg_fraction_bound_particle_sites_ii
                    accumulate_check(fbounds_sites2_new, iname, fbound_sites2_new, 1)
#                    print("DEBUG: {:} {:.3f} {:.3f}".format(iname,
#                                                            fbound_sites1_new,
#                                                            fbound_sites2_new))

                ii=ii+1
        for fg in output.statistics.fgs:
            fg_name= "%.2f %s" % (temperature_k, fg.type)
            prev_time_ns=fg.order_params[0].time_ns
            for fgop in fg.order_params:
                elapsed_time_ns= fgop.time_ns - prev_time_ns
                prev_time_ns= fgop.time_ns
                if fgop.time_ns<STATS_FROM_NS:
                    continue
                cur_mean_rg= fgop.mean_radius_of_gyration
                cur_mean_rg2= fgop.mean_square_radius_of_gyration
                accumulate_check(mean_rg, fg_name, cur_mean_rg, elapsed_time_ns)
                accumulate_check(mean_rg2, fg_name, cur_mean_rg2, elapsed_time_ns)
                cur_mean_dmax= fgop.mean_end_to_end_distance
                cur_mean_dmax2= fgop.mean_square_end_to_end_distance
                accumulate_check(mean_dmax, fg_name, cur_mean_dmax, elapsed_time_ns)
                accumulate_check(mean_dmax2, fg_name, cur_mean_dmax2, elapsed_time_ns)
                cur_mean_bond_rest_length= fgop.mean_bond_distance
                cur_mean_bond_rest_length2= fgop.mean_square_bond_distance
                accumulate_check(mean_bond_rest_length, fg_name,
                                 cur_mean_bond_rest_length, elapsed_time_ns)
                accumulate_check(mean_bond_rest_length2, fg_name,
                                 cur_mean_bond_rest_length2, elapsed_time_ns)
        for op in output.statistics.global_order_params:
            if op.time_ns<STATS_FROM_NS:
                continue
            if(op.energy < 1E+100 and not math.isnan(op.energy)):
                accumulate_check(energy,"energy", op.energy, 1)
        n=n+1

    if n==0:
        raise RuntimeError \
            ("Error - could not process " +
             "single file {}".format(sys.argv[1]) if len(fnames)==1 \
                 else "{:d} files".format(len(fnames)))
    if verbose:
        print "Indexed Fields:", indexed_fields
        print "K_on_i:"
    on_stats_i=do_stats(k_ons_i, "per ns per unbound fg motif", verbose=verbose)
    if(verbose):
        print "K_off_i:"
    off_stats_i=do_stats(k_offs_i, "per ns per bound fg motif", verbose=verbose)
    if(verbose):
        print "K_on_ii:"
    on_stats_ii=do_stats(k_ons_ii, "per ns per unbound float", verbose=verbose)
    if(verbose):
        print "K_off_ii:"
    off_stats_ii=do_stats(k_offs_ii, "per ns per bound float", verbose=verbose)
    if(verbose):
        print "K_on_ss:"
    on_stats_ss=do_stats(k_ons_ss, "per ns per missing NTF-FG motif contact", verbose=verbose)
    if(verbose):
        print "K_off_ss:"
    off_stats_ss=do_stats(k_offs_ss, "per ns per NTF-FG motif contact", verbose=verbose)
    if(verbose):
        print "fraction bounds:"
    #do_stats(fbounds, "of floats", is_fraction=True)
    fbound2_stats=do_stats(fbounds2, "of floats (from interaction order params)",
                           is_fraction=True, verbose=verbose)
    fbounds_sites1_stats_old=do_stats(fbounds_sites1_old, "of FG-sites (from old interaction order params)",
                                      is_fraction=True, verbose=verbose)
    fbounds_sites2_stats_old=do_stats(fbounds_sites2_old, "of float-sites (from old interaction order params)",
                                      is_fraction=True, verbose=verbose)
    KDs_dict_from_fbounds_old={}
    ret_value_new= {}
    for key in fbounds_sites1_stats_old.iterkeys():
        ret_value_new[key]= {}
    for key in fbounds_sites1_stats_old.iterkeys():
        if not re.search("fg0 - kap20", key):
            continue
        # Assume A is FG sites and B is kap sites
        EPS=1E-9
        B0 = N_KAPS * N_SITES_PER_KAP / AVOGADRO / box_volume_L # [B0] = total kap sites concentrtion
        fbA= max(min(fbounds_sites1_stats_old[key][0],1.0), EPSILON) # [AB]/[A0]
        fbB= max(min(fbounds_sites2_stats_old[key][0],1.0), EPSILON) # [AB]/[B0]
        # Confidence intervals:
        conf95_factor= get_sems_for_conf(0.95)
        fbA_conf95= fbounds_sites1_stats_old[key][1] * conf95_factor
        fbB_conf95= fbounds_sites2_stats_old[key][1] * conf95_factor
        fbA_low=  max(min(fbA-fbA_conf95, 1-1.5*EPSILON), 0.5*EPSILON)
        fbA_high= max(min(fbA+fbA_conf95, 1-0.5*EPSILON), 1.5*EPSILON)
        fbB_low=  max(min(fbB-fbB_conf95, 1-1.5*EPSILON), 0.5*EPSILON)
        fbB_high= max(min(fbB+fbB_conf95, 1-0.5*EPSILON), 1.5*EPSILON)
        # KD computation:
        A_per_AB= (1-fbA)/(fbA) # [A]/[AB]
        B= (1-fbB)*B0 # [B]
        KD = A_per_AB * B # [A]*[B]/[AB]
        KD_low=  B0 * (1-fbA_high) * (1-fbB_high) / fbA_high
        KD_high= B0 * (1-fbA_low) * (1-fbB_low) / fbA_low
        if verbose:
            print key, "KD calculation OLD", "[kap-sites]", pretty_molarity(B0),\
            "%bound A", fbA*100, "%bound B", fbB*100
            print key, "KD from fraction bound {} OLD".format(pretty_molarity(KD))
        ret_value_new[key]['fbound_sitesA_OLD']= fbA
        ret_value_new[key]['fbound_sitesA_lbound_OLD']= min(fbA_low, 1.0)
        ret_value_new[key]['fbound_sitesA_ubound_OLD']= min(fbA_high, 1.0)
        ret_value_new[key]['fbound_sitesB_OLD']= min(fbB, 1.0)
        ret_value_new[key]['fbound_sitesB_lbound_OLD']= min(fbB_low, 1.0)
        ret_value_new[key]['fbound_sitesB_ubound_OLD']= min(fbB_high, 1.0)
        ret_value_new[key]['KD_sites_OLD']= max(KD,1e-12)
        ret_value_new[key]['KD_sites_lbound_OLD']= max(KD_low,1e-12)
        ret_value_new[key]['KD_sites_ubound_OLD']= max(KD_high, 1e-12)
        if(KD<=0.0):
            continue
        KDs_dict_from_fbounds_old[key]=[KD,-1.0,-1.0]
    [dH,dS]= get_dH_and_dS2(KDs_dict_from_fbounds_old, "fg0 - kap20")
    if not math.isnan(dH) and not math.isnan(dS):
        T= 297.15
        if verbose:
            print "fg0-kap20 site-site-from-fbounds-old@%.2fK dH %.2e dS %.2e dS*T %.2e dG %.2e [kcal/mol]" % (T, dH, dS, dS*T, dH-dS*T)

    if IS_NEW_SITE_STATS:
        fbounds_sites1_stats_new=do_stats(fbounds_sites1_new, "of FG-sites (from new interaction order params)",
                                          is_fraction=True, verbose=verbose)
        fbounds_sites2_stats_new=do_stats(fbounds_sites2_new, "of float-sites (from new interaction order params)",
                                          is_fraction=True, verbose=verbose)
        KDs_dict_from_fbounds_new={}
        for key in fbounds_sites1_stats_new.iterkeys():
            if not re.search("fg0 - kap20", key):
                continue
            # A is FG sites and B is kap sites
            EPS=1E-9
            B0 = N_KAPS * N_SITES_PER_KAP / AVOGADRO / box_volume_L # [B0] = total kap sites concentrtion
            fbA= max(min(fbounds_sites1_stats_new[key][0],1.0), EPSILON) # [AB]/[A0]
            fbB= max(min(fbounds_sites2_stats_new[key][0],1.0), EPSILON) # [AB]/[B0]
            conf95_factor= get_sems_for_conf(0.95)
            fbA_conf95= fbounds_sites1_stats_new[key][1] * conf95_factor
            fbB_conf95= fbounds_sites2_stats_new[key][1] * conf95_factor
            fbA_low=  max(min(fbA-fbA_conf95, 1-1.5*EPSILON), 0.5*EPSILON)
            fbA_high= max(min(fbA+fbA_conf95, 1-0.5*EPSILON), 1.5*EPSILON)
            fbB_low=  max(min(fbB-fbB_conf95, 1-1.5*EPSILON), 0.5*EPSILON)
            fbB_high= max(min(fbB+fbB_conf95, 1-0.5*EPSILON), 1.5*EPSILON)
#            A_per_AB= (1-fbA)/(fbA) # [A]/[AB]
#            B= (1-fbB)*B0 # [B]
#            KD =    A_per_AB * B # [A]*[B]/[AB]
            KD=     B0 * (1-fbA) * (1-fbB) / fbA
            KD_low=  B0 * (1-fbA_high) * (1-fbB_high) / fbA_high
            KD_high= B0 * (1-fbA_low) * (1-fbB_low) / fbA_low
            if verbose:
                print key, "KD calculation NEW", "kap C", pretty_molarity(B0),\
                "%bound A", fbA*100, "%bound B", fbB*100
                print key, "KD from fraction bound {} NEW".format(pretty_molarity(KD))
            if(KD<=0.0):
                continue
            KDs_dict_from_fbounds_new[key]=[KD,-1.0,-1.0]
            ret_value_new[key]['fbound_sitesA_NEW']= fbA
            ret_value_new[key]['fbound_sitesA_lbound_NEW']= fbA_low
            ret_value_new[key]['fbound_sitesA_ubound_NEW']= fbA_high
            ret_value_new[key]['fbound_sitesB_NEW']= fbB
            ret_value_new[key]['fbound_sitesB_lbound_NEW']= fbB_low
            ret_value_new[key]['fbound_sitesB_ubound_NEW']= fbB_high
            ret_value_new[key]['KD_sites_NEW']= KD
            ret_value_new[key]['KD_sites_lbound_NEW']= KD_low
            ret_value_new[key]['KD_sites_ubound_NEW']= KD_high
        [dH,dS]= get_dH_and_dS2(KDs_dict_from_fbounds_new, "fg0 - kap20")
        if not math.isnan(dH) and not math.isnan(dS):
            T= 297.15
            print "fg0-kap20 site-site-from-fbounds-new@%.2fK dH %.2e dS %.2e dS*T %.2e dG %.2e [kcal/mol]" % (T, dH, dS, dS*T, dH-dS*T)



    # [A][B]/[AB] = [A][B]/([B0]-[B]) = [A][B]/([A0]-[A])

    #print on_stats_ii
    KDs_kaps_FGMs_dict={}
    for key in on_stats_i.keys():
        try:
            if not key in on_stats_i or not key in off_stats_i:
                continue
            [k_on_i, k_on_i_low, k_on_i_high]=  get_k_low_high(on_stats_i[key])
            [k_off_i, k_off_i_low, k_off_i_high]=  get_k_low_high(off_stats_i[key])
            if k_on_i<=0.0 or k_off_i<0.0:
                continue
            kd_i=k_off_i/k_on_i
            kd_i_high= k_off_i_high/k_on_i_low
            kd_i_low=  k_off_i_low /k_on_i_high
            fb_i=k_on_i/(k_on_i+k_off_i)
            fb_i_high=k_on_i_high/(k_on_i_high+k_off_i_low)
            fb_i_low= k_on_i_low /(k_on_i_low +k_off_i_high)
            KDs= [ x for x in [kd_i, kd_i_low, kd_i_high]]
            if verbose:
                print key, "KD_i %.2e [no-units] (%.2e .. %.2e)" % tuple(KDs)
                print key, "%% bound FG individal beads estimate from k_on/k_off times %.1f%% range (%.1f%% .. %.1f%%)" \
                    % (100*fb_i, 100*fb_i_low, 100*fb_i_high)
            if not key in on_stats_ii or not key in off_stats_ii:
                continue
            [k_on_ii, k_on_ii_low, k_on_ii_high]=  get_k_low_high(on_stats_ii[key])
            [k_off_ii, k_off_ii_low, k_off_ii_high]=  get_k_low_high(off_stats_ii[key])
            if k_on_ii<=0.0 or k_off_ii<0.0:
                continue
            kd_ii= k_off_ii/k_on_ii
            kd_ii_high= k_off_ii_high/k_on_ii_low
            kd_ii_low=  k_off_ii_low /k_on_ii_high
            fb_ii=k_on_ii/(k_on_ii+k_off_ii)
            fb_ii_high=k_on_ii_high/(k_on_ii_high+k_off_ii_low)
            fb_ii_low= k_on_ii_low /(k_on_ii_low +k_off_ii_high)
            KDs= [ x for x in [kd_ii, kd_ii_low, kd_ii_high]]
            if verbose:
                print key, "KD_ii %.2e [no-units] (%.2e .. %.2e)" % tuple(KDs)
                print key, "%% bound floats estimate from k_on/k_off times %.1f%% range (%.1f%% .. %.1f%%)" % \
                    (100*fb_ii, 100*fb_ii_low, 100*fb_ii_high)
            kaps_molar=N_KAPS/AVOGADRO/box_volume_L
            bound_kaps_molar=kaps_molar/(kd_ii+1.0)
            if verbose:
                print "Bound kap-FG motif molar:", bound_kaps_molar
            KDs = [bound_kaps_molar*t for t in \
                       [kd_i*kd_ii, kd_i_low*kd_ii_low, kd_i_high*kd_ii_high]]
            if KDs[0]<=0: # non-positive kd is necessarily an estimation error
                print  "KD_site_site_contacts INVALID-ESTIMATE [M]"
                continue
            pretty_KDs= [pretty_molarity(KD) for KD in KDs]
            if verbose:
                print key, "KD_kaps_FGMs {0} ({1} .. {2})".format(*pretty_KDs)
                print KDs
            KDs_kaps_FGMs_dict[key]= KDs
        except:
            print "EXCEPTION in NTF-FGM stats of key", key
            raise
        pass
    [dH,dS]= get_dH_and_dS(KDs_kaps_FGMs_dict, "fg0 - kap20")
    if not math.isnan(dH) and not math.isnan(dS) and verbose:
        T=297.15
        print "fg0-kap20 kap-FG@%.2fK dH %.2e dS %.2e dS*T %.2e dG %.2e [kcal/mol]" % (T, dH, dS, dS*T, dH-dS*T)

    KDs_dict={}
    ret_value= KDs_dict
    for key in on_stats_ss.keys():
        try:
            if not key in on_stats_ss or not key in off_stats_ss:
                continue
            [k_on, k_on_low, k_on_high]=  get_k_low_high(on_stats_ss[key]) # per ns per missing contact
            [k_off, k_off_low, k_off_high]=  get_k_low_high(off_stats_ss[key]) # per ns per contact
            if k_on<=0.0 or k_off<0.0:
                continue #(k_on can't be zero, negative means invalid)
            KD=k_off/k_on
            KD_high= k_off_high/k_on_low
            KD_low=  k_off_low /k_on_high
            fb=k_on/(k_on+k_off)
            fb_high=k_on_high/(k_on_high+k_off_low)
            fb_low= k_on_low /(k_on_low +k_off_high)
    #        if box_volume_L is None:
    #            print "box_volume_L not set for some reasom"
    #            box_volume_L=814.38
            if KD<=0: # non-positive KD is necessarily an estimation error
                print  "KD_site_site_contacts INVALID-ESTIMATE [M]"
                continue
    #        bound_sites_molar=[4*kaps_molar/(t+1.0) for t in KDs] # Need to correct for FG concentration
            KDs_dict[key]= [ x/AVOGADRO/box_volume_L for x in [KD, KD_low, KD_high]]
            pretty_KDs= [pretty_molarity(KD) for KD in KDs_dict[key]]
            if verbose:
                print key, "KD_site_site_contacts {} ({} .. {})".format(*pretty_KDs)
#            ret_value= KDs_dict
            if re.search("fg0 - kap20", key):
                ret_value_new[key]['k_on_per_ns_per_missing_ss_contact']= k_on
                ret_value_new[key]['k_on_per_ns_per_missing_ss_contact_lbound']= k_on_low
                ret_value_new[key]['k_on_per_ns_per_missing_ss_contact_ubound']= k_on_high
                ret_value_new[key]['k_off_per_ns_per_ss_contact']= k_off
                ret_value_new[key]['k_off_per_ns_per_ss_contact_lbound']= k_off_low
                ret_value_new[key]['k_off_per_ns_per_ss_contact_ubound']= k_off_high
            T=get_temperature_from_key(key)
            dG= math.log(KDs_dict[key][0])*kB_kcal_per_mol_K*T
            if verbose:
                print "dG %.2f" % dG
    #        print key, "%% bound site-site contacts estimate from k_on/k_off times %.1f%% range (%.1f%% .. %.1f%%)" % (100*fb, 100*fb_low, 100*fb_high)
        except:
            print "EXCEPTION in site-site stats of key", key
            raise
        if verbose:
            print KD,T,dG
        pass

    energy_stats= do_stats(energy, "kcal/mol", is_fraction=False, verbose=verbose)
    if ret_value_new is not None:
        ret_value_new["energy"]= energy_stats["energy"]

    [dH,dS]= get_dH_and_dS(KDs_dict, "fg0 - kap20")
    if not math.isnan(dH) and not math.isnan(dS) and verbose:
        T= 297.15
        print "fg0-kap20 site-site@%.2fK dH %.2f dS %.2e dS*T %.2f dG %.2f [kcal/mol]" % (T, dH, dS, dS*T, dH-dS*T)

    chain_KDs_from_float={}
    fbounds1_floats_stats= do_stats(fbounds1_floats, "%FGs bound (from float stats)",
                                    is_fraction=True, verbose=verbose)
    fbounds2_floats_stats= do_stats(fbounds2_floats, "%floats bound (from float stats)",
                                    is_fraction=True, verbose=verbose)
    kaps_molar=N_KAPS/AVOGADRO/box_volume_L
    fg_chains_molar=N_FGS/AVOGADRO/box_volume_L
    if verbose:
        print("FG molecular concentration {}" \
              .format(pretty_molarity(fg_chains_molar)))
        print("NTF ('kap') molecular concentration {}" \
              .format(pretty_molarity(kaps_molar)))
    for key in fbounds2_floats_stats.keys():
        fb1_list=fbounds1_floats_stats[key]
        fb2_list=fbounds2_floats_stats[key]
        fb1= max(min(fb1_list[0], 1-EPSILON), EPSILON)
        fb2= max(min(fb2_list[0], 1-EPSILON), EPSILON)
        conf95_factor= get_sems_for_conf(0.95)
        fb1conf95= fb1_list[1]*conf95_factor
        fb2conf95= fb2_list[1]*conf95_factor
        fb1low=  max(min(fb1-fb1conf95, 1-1.5*EPSILON), 0.5*EPSILON)
        fb1high= max(min(fb1+fb1conf95, 1-0.5*EPSILON), 1.5*EPSILON)
        fb2low=  max(min(fb2-fb2conf95, 1-1.5*EPSILON), 0.5*EPSILON)
        fb2high= max(min(fb2+fb2conf95, 1-0.5*EPSILON), 1.5*EPSILON)
#        print(fb1,fb1low, fb1high)
#        print(fb2,fb2low, fb2high)
        KD =    fg_chains_molar * (1-fb1)     * (1.0-fb2)   / fb2    # [Afree]*[Bfree]/[AB]
        KDlow=  fg_chains_molar * (1-fb1high) * (1-fb2high) / fb2high
        KDhigh= fg_chains_molar * (1-fb1low) * (1-fb2low) / fb2low
        chain_KDs_from_float[key]=[KD,KDlow, KDhigh]
        if verbose:
            print key, "KD_chain_interactions {} ({}..{})".format\
                (*[pretty_molarity(xx) for xx in (KD, KDlow, KDhigh)])
        if key in ret_value_new:
            ret_value_new[key]["fbound_chainsA"]= fb1
            ret_value_new[key]["fbound_chainsB"]= fb2
            ret_value_new[key]["fbound_chainsA_lbound"]= fb1low
            ret_value_new[key]["fbound_chainsB_lbound"]= fb2low
            ret_value_new[key]["fbound_chainsA_ubound"]= fb1high
            ret_value_new[key]["fbound_chainsB_ubound"]= fb2high
            ret_value_new[key]["KD_chains"]= KD
            ret_value_new[key]["KD_chains_lbound"]= KDlow
            ret_value_new[key]["KD_chains_ubound"]= KDhigh
    [dH,dS]= get_dH_and_dS(chain_KDs_from_float, "fg0 - kap20")
    if not math.isnan(dH) and not math.isnan(dS) and verbose:
        print "fg0-kap20 chains@%.2fK dH %.2e dS %.2e dS*T %.2e dG %.2e [kcal/mol]" % (T, dH, dS, dS*T, dH-dS*T)

    mean_rg_stats= do_stats(mean_rg, "A",
             is_fraction=False, prefix="Rg", verbose=verbose)
    if ret_value_new is not None:
        for key in mean_rg_stats:
            ret_value_new["Rg_{}".format(key)]= mean_rg_stats[key][0]
    do_stats_stddev(mean_rg, mean_rg2, "A",
                    is_fraction=False, prefix="Rg", verbose=verbose)
    mean_dmax_stats= do_stats(mean_dmax, "A",
             is_fraction=False, prefix="Dmax", verbose=verbose)
    if ret_value_new is not None:
        for key in mean_dmax_stats:
            ret_value_new["dmax_{}".format(key)]= mean_dmax_stats[key][0]
    do_stats_stddev(mean_dmax, mean_dmax2, "A",
                    is_fraction=False, prefix="Dmax", verbose=verbose)
    do_stats(mean_bond_rest_length, "A",
             is_fraction=False, prefix="bond-rest-length", verbose=verbose)
    do_stats_stddev(mean_bond_rest_length,
                    mean_bond_rest_length2,
                    "A",
                    is_fraction=False,
                    prefix="bond-rest-length",
                    verbose=verbose)
    try:
        if return_outputs is None:
            return ret_value, ret_value_new
        else:
            return ret_value, ret_value_new, return_outputs
    except:
        raise
        pass

def main():
    print sys.argv
    fnames= [sys.argv[1]] # set(glob.glob(sys.argv[1]))#+"*.pb"))
    STATS_FROM_SEC= float(sys.argv[2]) # start stats after specified seconds
    do_all_stats(fnames, STATS_FROM_SEC)

if __name__ == "__main__":
    # execute only if run as a script
    try:
        main()
    except:
        print('%s: %s' % (sys.argv[0], traceback.format_exc()))
        sys.exit(-1)
