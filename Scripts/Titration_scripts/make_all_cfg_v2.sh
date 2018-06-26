#!/bin/bash
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp /usr/bin/python"
for i in `seq 1 6`; do
    for k in False True; do
        $imppy MyScripts/make_cfg_v2.py --boxscale 1.0 -n 6 config${i}_Slide$k.pb $i $k > config${i}_Slide$k.txt;
    done;
 done
for k in False True; do
    $imppy MyScripts/make_cfg_v2.py --boxscale 1.0 -n 12 config${i}_Slide$k.pb 12 $k > config12_Slide$k.txt;
done
