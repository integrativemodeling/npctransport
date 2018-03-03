#!/bin/bash

for i in cyl_* TEST*; do
    S=$i/STATS
    if [ -e $S ]; then
        echo $i
        cat $i/STATS
        echo ==
    fi;
done;
