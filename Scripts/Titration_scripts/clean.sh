#!/bin/bash
from=$1
to=$2
for j in config_v13_*; do
    echo $j from $from to $to;
    cd $j;
    for i in `seq $from $to`; do
        for f in `find ./Logs/ -maxdepth 1 -type f -name 'scan_Rg.o*.'"$i"''`; do
            if [ `egrep 'terminate|error while loading|Aborted' $f|wc -l` -gt 0 ]; then
                echo Deleting $f
                rm  output$i.* $f;
                MOVIE=movie$i.rmf
                if [ -e $MOVIE ]; then
                    rm $MOVIE
                fi
                continue # no need for more log files for $i if deleted...
            fi
         done;
     done;
     echo;
     cd ..;
done;

echo FINISHED CLEANING
