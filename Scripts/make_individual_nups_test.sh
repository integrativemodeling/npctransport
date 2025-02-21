#!/bin/bash
# NOTE: create soft link to Scripts folder (with load_whole_new_*.py) before running
set -e
TIME_STEP_FACTOR=${1:-2}
imp_folder=$HOME/imp_git/fast
imp=$imp_folder/setup_environment.sh
imppy="$imp /usr/bin/python"
npc_fg_simulation="$imp $imp_folder/module_bin/npctransport/fg_simulation"
for i in Nsp1 Nup1 Nup100 Nup49 FSFG_generic GLFG_generic; do
    $imppy Scripts/load_whole_new_coarse_grained_v12.py --only_nup $i config_$i.pb ../InputData/47-35_1spoke.rmf3  --time_step_factor $TIME_STEP_FACTOR  >& config_$i.txt
    $npc_fg_simulation --configuration config_$i.pb --output output_$i.pb --short_sim_factor 1000.0 >& LOG.$i &
done
