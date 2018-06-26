#!/bin/bash
# Equivalent to Figure 1-2 (forward titration) in Sam and Ryo's paper
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp /usr/bin/python"
VERSION=13
TIME_STEP_FACTOR=$1

for n in `seq 1 6`; do
    # fg_C in Mol:
    fg_C=`python -c "print 0.000120/$n"`
    bead_string=`python -c "print 'F'*$n+'S'*(6-$n)"`
    for k in True; do
        for kap_titration in `seq -4 7`; do
            kap2fg=`python -c "print 2**$kap_titration;"`
            kap_C=`python -c "print $fg_C * $kap2fg;"`
            echo $fg_C
            prefix=config_v${VERSION}_${bead_string}_Slide${k}_fgC${fg_C}M_kapC${kap_C}M
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            if [ ! -e config.pb ]; then
                $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
            fi
            if [ ! -e configMovie.pb ]; then
                sed 's/^output_statistics_interval_ns: 30.0/output_statistics_interval_ns: 0.25/g' config.txt | grep -v '^#' > configMovie.txt
                $imppy ~/npctransport/utility/txt2config.py configMovie.txt configMovie.pb
            fi
            cd ..;
        done;
    done;
 done


# fg_C in Mol:
fg_C=0.000060
for n in 1 2 3 4; do
    bead_string=`python -c "print 'F'+'S'*$n+'F'+'S'*(4-$n)"`
    for k in True; do
        for kap_titration in `seq -4 7`; do
            kap2fg=`python -c "print 2**$kap_titration;"`
            kap_C=`python -c "print $fg_C * $kap2fg;"`
            echo $fg_C
            prefix=config_v${VERSION}_${bead_string}_Slide${k}_fgC${fg_C}M_kapC${kap_C}M
            echo $prefix
            if [ ! -e $prefix ]; then mkdir $prefix; fi;
            cd $prefix
            if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
            if [ ! -e Logs/ ]; then mkdir Logs; fi;
            if [ ! -e config.pb ]; then
                $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
            fi
            if [ ! -e configMovie.pb ]; then
                sed 's/^output_statistics_interval_ns: 30.0/output_statistics_interval_ns: 0.25/g' config.txt | grep -v '^#' > configMovie.txt
                $imppy ~/npctransport/utility/txt2config.py configMovie.txt configMovie.pb
            fi
            cd ..;
        done;
    done;
 done

# fg_C in Mol:
fg_C=0.000025
bead_string=`python -c "print 'F'*12"`
for k in True; do
    for kap_titration in `seq -1 7`; do
        kap2fg=`python -c "print 2**$kap_titration;"`
        kap_C=`python -c "print $fg_C * $kap2fg;"`
        echo $fg_C
        prefix=config_v${VERSION}_${bead_string}_Slide${k}_fgC${fg_C}M_kapC${kap_C}M
        echo $prefix
        if [ ! -e $prefix ]; then mkdir $prefix; fi;
        cd $prefix
        if [ ! -e MyScripts/ ]; then ln -s ../MyScripts/; fi;
        if [ ! -e Logs/ ]; then mkdir Logs; fi;
        if [ ! -e config.pb ]; then
            $imppy MyScripts/make_cfg_v${VERSION}.py -f $fg_C -k $kap_C --time_step_factor $TIME_STEP_FACTOR config.pb $bead_string $k > config.txt;
        fi
        if [ ! -e configMovie.pb ]; then
            sed 's/^output_statistics_interval_ns: 30.0/output_statistics_interval_ns: 0.25/g' config.txt | grep -v '^#' > configMovie.txt
            $imppy ~/npctransport/utility/txt2config.py configMovie.txt configMovie.pb
        fi
        cd ..;
    done;
done;
