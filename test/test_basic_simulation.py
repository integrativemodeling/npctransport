#!/usr/bin/python
# Synopsis:    Unit test for a basic FG-NTR simulation
#
# Overview:    Test that a basic simulation runs without errors
#              Note this simulation is not realistic, just a test using a 
#              simplified NPC model
#
# Author:       Barak Raveh barak.raveh@mail.huji.ac.il 
#               (based on tutorial by Roi Eliasian roi.eliasian@mail.huji.ac.il)
#
# Created:      September 2025
# Last updated: September 2025
#
# Expected running time with default parameters: ~5-60 minutes on a basic laptop
# (Reduce simulation time and initialization time parameters below for faster tests))

import IMP
import IMP.npctransport 
import RMF
import itertools
import numpy as np
import scipy.stats
import os
import pathlib
import subprocess
import sys
import tempfile


def is_conda_python() -> bool:
    # Works for base and named envs; also check base_prefix in case you're in a venv
    prefixes = {sys.prefix, getattr(sys, "base_prefix", ""), getattr(sys, "real_prefix", "")}
    return any((pathlib.Path(p) / "conda-meta").is_dir() for p in prefixes if p)

######################
##### Parameters #####
######################
# FG simulation folder -
# if using IMP source installation, this would be <IMP_BUILD_FOLDER>/modules/npctransport/bin/;
# if using conda installation, install the imp module, and then this would be <CONDA_ENV_FOLDER>/bin/ for the correct environment
PYTHON_PREFIX = sys.prefix
FG_SIMULATION_FOLDER = os.path.join(sys.prefix, "bin") if is_conda_python() else "<PATH/TO/FG/SIMULATION/BINARY>/bin"
CONFIG_FILE_NAME = "config.pb"
OUTPUT_FILE_NAME = "output.pb"
MOVIE_FILE_NAME = "movie.rmf"
FINAL_FILE_NAME = "final.rmf"
SIMULATION_TIME_NS = 200 # total simulation time will greatly affect runtime (bestr not below 10 for this test)
SIMULATION_OUTPUT_STATISTICS_INTERVAL_NS = 1 # interval at which output statistics are recorded
SIMULATION_SHORT_INIT_FACTOR = 0.25 # fraction of the time to use for a short initialization run (best not below 0.25 for this test)
BOUNDING_BOX_SIDE_A = 1000 # side of the cubic bounding box
NUCLEAR_ENVELOPE_THICKNESS_A = 150 # thickness of the slab representing the nuclear envelope
N_NTRS = 400
NTR_RADII_AND_SITES = [(20,4)] # Radius 20, 4 Interaction sites distributed uniformally
FG_CHAIN_LENGTH_AA = 320
FG_AA_PER_BEAD = 20
FG_CHAIN_LENGTH_BEADS = FG_CHAIN_LENGTH_AA // FG_AA_PER_BEAD
# demonstrate two types of FG beads: 3/4 of the beads are type C, 1/4 are type N
FG_N_BEADS_C = FG_CHAIN_LENGTH_BEADS * 3 // 4
FG_N_BEADS_N = FG_CHAIN_LENGTH_BEADS - FG_N_BEADS_C
FG_CHAIN_SUFFIX_LIST = ["_C"] * FG_N_BEADS_C + ["_N"] * FG_N_BEADS_N
FG_BEAD_RADIUS_A = 8.0
FG_REST_LENGTH_FACTOR = 1.9
# Polymer chain parameters for "C"- and "N"-type FG beads (akin to e.g. GLFG and FxFG)
# These parameters are illustrative but they in fact resemble real values fitted
# for GLFG and FxFG type repeats in Raveh et al. PNAS 2025 (see Supplementary Appendix Table S13)
FG_SELF_K = {"_N" : 1.47, "_C" : 1.32}
FG_SELF_RANGE = 6.00
FG_SELF_NONSPEC_K = {"_N" : 0.01, "_C" : 0.08}
FG_NTR_K = 2.64
FG_NTR_RANGE = 5.5
FG_NTR_SIGMA0_DEG = 45.0
FG_NTR_SIGMA1_DEG = 45.0

def add_fgs_to_config(config):
    ''' 
    Add FGs and their internal interactions to the configuration 
    @param config: Configuration object to modify
    @return:       List of FG chains added
    ''' 
    # Anchor coordinates for 32 FG chains in a simplified NPC model (4 layers of 8 FGs each)
    anchor_coordinates = [[185.0, 0.0, -75.0],
                          [130.8147545195113, 130.8147545195113, -75.0], 
                          [1.1327982892113017e-14, 185.0, -75.0], 
                          [-130.81475451951127, 130.8147545195113, -75.0], 
                          [-185.0, 2.2655965784226034e-14, -75.0], 
                          [-130.81475451951133, -130.81475451951127, -75.0], 
                          [-3.398394867633905e-14, -185.0, -75.0], 
                          [130.81475451951127, -130.81475451951133, -75.0], 
                          [120.0480947161671, 0.0, -37.5], 
                          [84.88682184232671, 84.88682184232671, -37.5], 
                          [7.350825746894615e-15, 120.0480947161671, -37.5], 
                          [-84.8868218423267, 84.88682184232671, -37.5],
                          [-120.0480947161671, 1.470165149378923e-14, -37.5], 
                          [-84.88682184232673, -84.8868218423267, -37.5],
                          [-2.2052477240683848e-14, -120.0480947161671, -37.5], 
                          [84.88682184232668, -84.88682184232673, -37.5]]
    # Add another 16 FGs in the same positions but at +37.5 and +75.0 z coordinates
    anchor_coordinates_slice = anchor_coordinates[:]
    for anchor in anchor_coordinates_slice:
        anchor_coordinates.append([anchor[0], anchor[1], -anchor[2]])
    fgs = []
    for frame_i, coordinates in enumerate(anchor_coordinates):
        cur_fg = IMP.npctransport.add_fg_type(config,
                                            type_name=f"fg{frame_i}",
                                            number_of_beads=FG_CHAIN_LENGTH_BEADS,
                                            number=1,
                                            radius=FG_BEAD_RADIUS_A,
                                            interactions=1,
                                            rest_length_factor=FG_REST_LENGTH_FACTOR,
                                            d_factor=1.0,
                                            interaction_k_factor=1.0,
                                            interaction_range_factor=1.0)
        pos = cur_fg.anchor_coordinates.add()
        pos.x = coordinates[0]
        pos.y = coordinates[1]
        pos.z = coordinates[2]
        fgs.append(cur_fg)
        cur_fg.type_suffix_list.extend(FG_CHAIN_SUFFIX_LIST)
    # Add internal FG interactions: 
    unique_suffixes = set(FG_CHAIN_SUFFIX_LIST)
    for frame_i in itertools.combinations_with_replacement(range(len(fgs)), 2):
        for suffix0 in unique_suffixes:
            for suffix1 in unique_suffixes:
                interactionFG_FG = IMP.npctransport.add_interaction(config,
                                                            name0=f"fg{frame_i[0]}{suffix0}",
                                                            name1=f"fg{frame_i[1]}{suffix1}",
                                                            interaction_k=0.5*FG_SELF_K[suffix0] + 0.5*FG_SELF_K[suffix1],
                                                            interaction_range=FG_SELF_RANGE)
                interactionFG_FG.nonspecific_k.lower = np.sqrt(FG_SELF_NONSPEC_K[suffix0] * FG_SELF_NONSPEC_K[suffix1])
    return fgs

def add_NTRs_to_config(config, fgs):
    ''' 
    Add NTRs and their interactions with fgs to the configuration 
    @param config: Configuration object to modify
    @param fgs:    List of FG chains added
    @return:       List of NTRs added
    ''' 
    ntr_vals = NTR_RADII_AND_SITES
    ntrs = []
    for vals in ntr_vals:
        cur_ntr  = IMP.npctransport.add_float_type(config,
                                    number=N_NTRS,
                                    radius=vals[0],
                                    type_name=f"NTR{vals[0]}",
                                    interactions=vals[1],
                                    d_factor=1.0,
                                    interaction_k_factor=1.0,
                                    interaction_range_factor=1.0)
        ntrs.append(cur_ntr)

    ############################
    # Add NTR-FG interactions: #
    ############################
    for frame_i in range(len(fgs)):
        for vals in ntr_vals:
            for suffix in set(FG_CHAIN_SUFFIX_LIST):
                IMP.npctransport.add_interaction(config,
                                        name0=f"fg{frame_i}{suffix}",
                                        name1=f"NTR{vals[0]}",
                                        interaction_k=FG_NTR_K,
                                        interaction_range=FG_NTR_RANGE,
                                        range_sigma0_deg=FG_NTR_SIGMA0_DEG,
                                        range_sigma1_deg=FG_NTR_SIGMA1_DEG)
    return ntrs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

def make_config_file(config_file_name):
    ''' 
    Create a configuration file for a basic FG-NTR simulation in a temporary folder
    @param config_file_name: Name of the configuration file to create
    ''' 
    # See User's manual for interpretation of parameters
    config = IMP.npctransport.Configuration()
    IMP.npctransport.set_default_configuration(config)
    config.statistics_fraction.lower = 1.0
    config.interaction_k.lower = 1e-09
    config.interaction_range.lower = 10
    config.backbone_k.lower = 0.0075
    config.time_step_factor.lower = 2.0
    config.slack.lower = 10
    config.number_of_trials = 1
    config.statistics_interval_ns = 0.1
    config.dump_interval_ns = 1
    config.output_statistics_interval_ns = SIMULATION_OUTPUT_STATISTICS_INTERVAL_NS
    config.simulation_time_ns = SIMULATION_TIME_NS
    config.box_is_on.lower = 1
    config.box_side.lower = BOUNDING_BOX_SIDE_A
    config.slab_is_on.lower = 2
    config.slab_thickness.lower = NUCLEAR_ENVELOPE_THICKNESS_A
    config.tunnel_radius.lower = 185
    config.is_xyz_hist_stats = 1
    config.backbone_tau_ns.lower = 50.0
    config.is_backbone_harmonic = 1
    config.xyz_stats_crop_factor = 1
    config.xyz_stats_voxel_size_a = 10
    config.xyz_stats_max_box_size_a = 700
    config.nonspecific_k.lower = 0.01
    config.excluded_volume_k.lower = 10.0
    config.angular_D_factor.lower=0.3
    config.is_multiple_hdf5s = False
    config.full_output_statistics_interval_factor = 1
    # Add FGs and NTRs:
    fgs = add_fgs_to_config(config)
    ntrs = add_NTRs_to_config(config, fgs)
    # Dump to a file:
    f = open(config_file_name, "wb")
    f.write(config.SerializeToString())
    

def load_time_energy(pb_path):
    ''' load time and energy from a simulation output file '''
    with open(pb_path, "rb") as f:
        output = IMP.npctransport.Output()
        fstring = f.read()
        output.ParseFromString(fstring)
        time_ns = []
        energy = []
        for i in range(len(output.statistics.global_order_params)):
            time_ns.append(output.statistics.global_order_params[i].time_ns)
            energy.append(output.statistics.global_order_params[i].energy)
    return time_ns, energy

def _rmf_has_depth_with_site(root, i):
    """ returns true if node subtree thru first child is at least i
        levels, including the root node itself, and the lead is a site """
#  print root, i, len(root.get_children())
    if (i==1) and root.get_name()=="site":
        return True
    c = root.get_children()
    if len(c) == 0:
        return False
    return _rmf_has_depth_with_site(c[0], i-1)

def _rmf_add_nodes(node, tf, type_prefixes, depth=0):
    '''
    node - rmf node to scan
    tf - typed factory
    type_prefixes - list of full type prefixes (e.g. "Nup1" for "Nup1N")

    adds only nodes whose type name begins with any of the specified type prefixes
    
    @return list of lists of RMF nodes, each list contains all nodes of a given type
            found at the same level (e.g. all NTRs of a given type)
    '''
    children = node.get_children()
    ret = []
    if len(children)==0:
        return ret
    if _rmf_has_depth_with_site(node, 3) and tf.get_is(children[0]):
        child_type = tf.get(children[0]).get_type_name()
        if any([child_type.startswith(tp) for tp in type_prefixes]):
            ret.append(children)
    for c in children:
        ret += _rmf_add_nodes(c, tf,  type_prefixes, depth+1)
    return ret

def _rmf_add_nodes_exact(node, tf, types, depth=0):
    '''
    node - rmf node to scan
    tf - typed factory
    types - list of full types (e.g. "Nup1N")

    adds only nodes whose type name begins with any of the specified type prefixes
    '''
    children = node.get_children()
    ret = []
    if len(children)==0:
        return ret
    if _rmf_has_depth_with_site(node, 3) and tf.get_is(children[0]):
        child_type = tf.get(children[0]).get_type_name()
        if any([child_type == tp for tp in types]):
            ret.append(children)
    for c in children:
        ret += _rmf_add_nodes_exact(c, tf,  types, depth+1)
    return ret

def load_NTR_data_from_RMF(rmf_path, n_NTRs, max_frames_per_file):
    ''' 
    Load NTR trajectories from an RMF file for the first NTR type 
    @param rmf_path:      Path to the RMF file
    @param n_NTRs:        Number of NTRs in the simulation
    @param max_frames_per_file: Maximal number of frames to read from the RMF file 
                            (from the end of the file)
    @return:              Numpy array of shape (n_NTRs, 3, frames_per_file)
                          containing the trajectories of the NTRs in angstrom units
                          that is trajectories[NTR_index, x/y/z, frame_index]
    '''
    print("Max frames per file:", max_frames_per_file)
    in_fh = RMF.open_rmf_file_read_only(rmf_path)
    rff = RMF.ReferenceFrameFactory(in_fh)
    tf = RMF.TypedFactory(in_fh)
    radii = [ntr_radius_and_site[0] for ntr_radius_and_site in NTR_RADII_AND_SITES]
    NTR_types = [f"NTR{radius}" for radius in radii] # name defined in config.pb
    # load data
    type2chains={}
    for NTR_type in NTR_types:
        type2chains[NTR_type] = _rmf_add_nodes(in_fh.get_root_node(), tf, [NTR_type])
    NTR_type = NTR_types[0] # for this test we use only one type
    frames = list(in_fh.get_frames())
    frame_i = 0
    if len(frames) > max_frames_per_file:
        frame_i = len(frames) - max_frames_per_file
        frames = frames[-max_frames_per_file:]
    trajectories = np.zeros(shape=[n_NTRs, 3, len(frames)], dtype=np.float32)
    for fi, f in enumerate(frames):
        in_fh.set_current_frame(f)
        # read data
        for NTR_i in range(n_NTRs):
            coord = rff.get(type2chains[NTR_type][0][NTR_i]).get_translation() 
            trajectories[NTR_i, 0, fi] = coord[0] 
            trajectories[NTR_i, 1, fi] = coord[1] 
            trajectories[NTR_i, 2, fi] = coord[2] 
        frame_i += 1     
    print("Read", len(frames), "frames from", rmf_path)
    return trajectories

if __name__ == "__main__":
    # Make tmp_output folder if it doesn't exist
    tmp_dir = "tmp_output/"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    config_file_path = os.path.join(tmp_dir, CONFIG_FILE_NAME)
    output_file_path = os.path.join(tmp_dir, OUTPUT_FILE_NAME)
    movie_file_path = os.path.join(tmp_dir, MOVIE_FILE_NAME)
    final_file_path = os.path.join(tmp_dir, FINAL_FILE_NAME)
    make_config_file(config_file_path)
    print(f"Configuration file created: {config_file_path}")
    print("Running simulation...")
    try:
        subprocess.run([f"{FG_SIMULATION_FOLDER}/fg_simulation", 
                        "--configuration", config_file_path,
                        "--output", output_file_path,
                        "--conformations", movie_file_path,
                        "--final_conformations", final_file_path,
                        "--short_init_factor", "0.25",
                        "--short_sim_factor", "1.0",
                        "--log_level", "SILENT"])
    except FileNotFoundError:
        print(f"Error: fg_simulation binary not found in {FG_SIMULATION_FOLDER}.")
        print(f"Please check the FG_SIMULATION_FOLDER variable in the script.")
        print(f"If using CONDA, please install the IMP module and make sure you")
        print(f"are using python from the correct environment.")
        sys.exit(1)
    print("Simulation complete")
    time_ns, energy = load_time_energy(output_file_path)
    # crop the first 50% of the data as equilibration
    crop_index = len(energy) // 2
    time_ns = time_ns[crop_index:]
    energy = energy[crop_index:]
    mean_energy = np.mean(energy)
    stddev_energy = scipy.stats.tstd(energy) if len(energy) > 1 else 0.0
    stderr_energy = scipy.stats.sem(energy) if len(energy) > 1 else 0.0
    print(f"Mean energy: {mean_energy} kcal/mol over {len(energy)} frames ({time_ns[-1]} ns simulated)")
    print(f"Std-dev energy: {stddev_energy} kcal/mol")
    print(f"Std-error energy: {stderr_energy} kcal/mol")
    if(SIMULATION_TIME_NS >= 9.99 and SIMULATION_SHORT_INIT_FACTOR >= 0.2499):
        assert(mean_energy < 1000)
        assert(stddev_energy < 200)
        assert(stderr_energy < 50)
    else:
        print("Note: simulation time is short, so energy statistics may not be meaningful")
    trajectories = load_NTR_data_from_RMF(movie_file_path, N_NTRS, max_frames_per_file=50)
    Z = trajectories[:, 2, :].flatten() # all NTRs and all frames
    # Compute mean Z over all NTRs and all frames
    mean_z = np.mean(Z)
    stddev_z = scipy.stats.tstd(Z) if len(trajectories) > 1 else 0.0
    stderr_z = scipy.stats.sem(Z) if len(trajectories) > 1 else 0.0
    print(f"Mean Z of NTRs: {mean_z} A")
    print(f"Std-dev Z of NTRs: {stddev_z} A")
    print(f"Std-error Z of NTRs: {stderr_z} A")
    # Count how many NTRs are within the pore region 
    ntr_in_pore = np.sum(np.abs(Z) < NUCLEAR_ENVELOPE_THICKNESS_A/2) / trajectories.shape[2] # average per frame
    fraction_in_pore = ntr_in_pore / (N_NTRS * 1)
    print(f"Percentage of NTRs in pore (|Z|<NE/2): {100*fraction_in_pore:.2f}%; {ntr_in_pore} out of {N_NTRS} NTRs")
    R = np.sqrt(trajectories[:, 0, :]**2 + trajectories[:, 1, :]**2).flatten()
    R_in_pore = R[np.abs(Z) < NUCLEAR_ENVELOPE_THICKNESS_A/2]
    if len(R_in_pore) == 0:
        print("No NTRs found in pore region for R statistics")
    else:
        mean_R_in_pore = np.mean(R_in_pore) 
        stddev_R_in_pore = scipy.stats.tstd(R_in_pore) if len(R_in_pore) > 1 else 0.0
        stderr_R_in_pore = scipy.stats.sem(R_in_pore) if len(R_in_pore) > 1 else 0.0
        print(f"Mean R of NTRs in pore (|Z|<NE/2): {mean_R_in_pore} A")   
        print(f"Std-dev R of NTRs: {stddev_R_in_pore} A")
        print(f"Std-error R of NTRs: {stderr_R_in_pore} A")
