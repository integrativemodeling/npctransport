#!/bin/bash
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp /usr/bin/python"
VERSION=8
TIME_STEP_FACTOR=$1

for i in `seq 1 6`; do
    for k in True; do
        for r in 0.0 0.6 1.6 3.6; do
            prefix=config_v${VERSION}_n${i}_Slide${k}_ratio${r}
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            $imppy MyScripts/make_cfg_v${VERSION}.py -n 6 -r $r --time_step_factor $TIME_STEP_FACTOR config.pb $i $k > config.txt;
            cd ..;
        done;
    done;
 done
i=12
for k in True; do
    for r in 0.0 0.6 1.6  3.6; do
        prefix=config_v${VERSION}_n${i}_Slide${k}_ratio${r}
        echo $prefix
        if [ ! -e $prefix ]; then mkdir $prefix; fi;
        cd $prefix
        if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
        if [ ! -e Logs/ ]; then mkdir Logs; fi;
        $imppy MyScripts/make_cfg_v${VERSION}.py -n 12 -r $r --time_step_factor $TIME_STEP_FACTOR config.pb $i $k > config.txt
        cd ..;
    done;
done
