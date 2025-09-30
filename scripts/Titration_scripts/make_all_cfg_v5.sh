#!/bin/bash
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp /usr/bin/python"

for i in `seq 1 6`; do
    for k in True; do
        for r in 0.0 0.2 0.6 1.0 1.6 2.4 3.6; do
            prefix=config_v5_n${i}_Slide${k}_ratio${r}
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            $imppy MyScripts/make_cfg_v5.py -n 6 -r $r config.pb $i $k > config.txt;
            cd ..;
        done;
    done;
 done
i=12
for k in True; do
    for r in  0.2 0.6 1.0 1.6 2.4 3.6; do
        prefix=config_v5_n${i}_Slide${k}_ratio${r}
        echo $prefix
        if [ ! -e $prefix ]; then mkdir $prefix; fi;
        cd $prefix
        if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
        if [ ! -e Logs/ ]; then mkdir Logs; fi;
        $imppy MyScripts/make_cfg_v5.py -n 12 -r $r config.pb $i $k > config.txt
        cd ..;
    done;
done
