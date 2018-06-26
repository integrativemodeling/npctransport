#!/bin/bash
qstat | awk '{if(NR>2){print NR, $0, $NF-x} x=$NF;}'| awk  -v n=$1 -E MyScripts/moving_avg.awk
