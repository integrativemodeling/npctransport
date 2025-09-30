#!/bin/bash
# Equivalent to Figure 3 (reverse titration) in Sam and Ryo's paper
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp /usr/bin/python"
VERSION=12
TIME_STEP_FACTOR=$1


# kap_C in Mol:
kap_C=0.000250
for bead_string in FSSSSS; do
    for k in True; do
        for fg_titration in `seq -4 4`; do
            fg2kap=`python -c "print 2**$fg_titration;"`
            fg_C=`python -c "print $kap_C * $fg2kap;"`
            echo $fg_C
            prefix=config_v${VERSION}_${bead_string}_Slide${k}_kapC${kap_C}M_fgC${fg_C}M
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            if [ ! -e config.pb ]; then
                $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
            fi
#            if [ ! -e configMovie.pb ]; then
#                $imppy MyScripts/make_cfg_v${VERSION}movie.py -r $r --time_step_factor $TIME_STEP_FACTOR configMovie.pb $bead_string $k > configMovie.txt;
#            fi
            cd ..;
        done;
    done;
 done


# kap_C in Mol:
kap_C=0.000250
for bead_string in FFFSSS; do
    for k in True; do
        for fg_titration in `seq -4 3`; do
            fg2kap=`python -c "print 2**$fg_titration;"`
            fg_C=`python -c "print $kap_C * $fg2kap;"`
            echo $fg_C
            prefix=config_v${VERSION}_${bead_string}_Slide${k}_kapC${kap_C}M_fgC${fg_C}M
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            if [ ! -e config.pb ]; then
                $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
            fi
#            if [ ! -e configMovie.pb ]; then
#                $imppy MyScripts/make_cfg_v${VERSION}movie.py -r $r --time_step_factor $TIME_STEP_FACTOR configMovie.pb $bead_string $k > configMovie.txt;
#            fi
            cd ..;
        done;
    done;
 done

# kap_C in Mol:
kap_C=0.000500
for bead_string in FFFFFF; do
    for k in True; do
        for fg_titration in `seq -4 1`; do
            fg2kap=`python -c "print 2**$fg_titration;"`
            fg_C=`python -c "print $kap_C * $fg2kap;"`
            echo $fg_C
            prefix=config_v${VERSION}_${bead_string}_Slide${k}_kapC${kap_C}M_fgC${fg_C}M
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            if [ ! -e config.pb ]; then
                $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
            fi
#            if [ ! -e configMovie.pb ]; then
#                $imppy MyScripts/make_cfg_v${VERSION}movie.py -r $r --time_step_factor $TIME_STEP_FACTOR configMovie.pb $bead_string $k > configMovie.txt;
#            fi
            cd ..;
        done;
    done;
 done
