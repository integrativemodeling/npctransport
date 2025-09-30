#!/usr/bin/bash
#TAG=2
#RESTART_TAG='restart1\.pbs'
TAG=$1 # the id after N, e.g., 3 for N3output.pb
RESTART_TAG=$2 # script used for creating N${TAG}output.pb, e.g. restart2.pbs for TAG=3

(j=`mktemp`; (for i in `qstat | awk '{print $1}' | sort -u| grep '[0-9]'`; do echo -n "$i "; qstat -j $i | grep  pbs|awk '{print $2}'; done) | sort -k 2 | grep $RESTART_TAG | awk '{print $1}' > $j; for i in `cat $j`; do qstat -j $i | grep usage|awk '{print $2'}|sed 's/://g'; done;) | sort -n > running_N${TAG}s.txt
for i in N${TAG}*.pb; do if [ ! -e $i.hdf5 ]; then j=`echo $i|sed 's/^.*[^0-9]\([0-9]*\)\.pb$/\1/g'`; if [ `grep "^$j\$" running_N${TAG}s.txt | wc -l` -eq 0 ]; then echo rm N${TAG}output$j.pb; fi; fi; done
