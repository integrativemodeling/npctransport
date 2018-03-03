#!/bin/bash
K=${1:-10}
TMP=${2:-`mktemp`}
FORCE_DO_TEMP=${3:-1}
#MAX=${4:-500}
if [ $FORCE_DO_TEMP -eq 1 ] || [ ! -e $TMP ]; then
#    (for i in `seq 1 $MAX`; do N=`cat Logs/cyl.job*.i | grep -i 'optimizing' | grep -i 'frame' | grep -v ' 0 frames' |  wc -l`; echo $N; done;) | sort -n | uniq -c > $TMP
    (for i in Logs/cyl.job*; do N=`cat $i | grep -i 'optimizing' | grep -i 'frame' | wc -l`; echo $N; done;) | sort -n | uniq -c > $TMP
fi
awk 'function print_stars(num){for(ii=0;ii<5*num;ii++){printf "*"; if(ii%5==4){printf " "}} printf("\n");} function print_index(all,K,idx,s){printf("%5d %5d ",all[(idx-K/2)%K], idx-K/2); print_stars(s/K);} BEGIN{K='"$K"'; nprev=0; for(i=0;i<K;i++){all[i]=0;}} {n=$2; for(i=nprev+1; i<n; i++){all[i%K]=0; print_index(all,K,i,s); s=s-all[(i+1)%K]} nprev=n; all[n%K]=$1; s=s+$1; print_index(all,K,n,s); s=s-all[(n+1)%K];} END{for(i=nprev+1; i<nprev+K; i++) {all[i%K]=0; print_index(all,K,i,s); s=s-all[(i+1)%K]}}' $TMP