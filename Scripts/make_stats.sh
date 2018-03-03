#!/bin/bash
imp=$HOME/imp_git/fast/setup_environment.sh
imppy="$imp python"

for i in cyl_* TEST_*; do
    cd $i;
    echo $i
    $imppy Scripts/get_float_stats.py output*.pb >& LOG.float &
    $imppy Scripts/transport_stats.py output 10E-6 |sed 's/^R//g' >& STATS &
    cd ..;
done;
