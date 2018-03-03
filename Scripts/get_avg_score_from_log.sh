#!/bin/bash
TMP=`mktemp`
touch $TMP;
for i in `seq 1 1000`; do
    cat Logs/cyl*.$i | grep Score   | grep -v BD | cat -n >> $TMP;
done
sort -n $TMP | awk 'BEGIN{prev=1; n=0};{if($1==prev){s=s+$NF;s2=s2+$NF**2;n=n+1;}else{nn=n+0.000001; print prev,n, s/nn,"+-", sqrt(s2/nn-(s/nn)**2)/sqrt(nn); prev=$1; s=0; s2=0; n=0;}}'
