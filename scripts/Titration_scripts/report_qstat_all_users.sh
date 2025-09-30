#!/bin/bsash
echo Running
T=`mktemp`;
qstat -u '*' > $T;
(for i in `cat $T | awk '$5~/r/{print $4}' | sort -u`; do
    echo $i `cat $T | awk '$5~/r/ && $0~/'"$i"'/{n=n+$9} END{print n}'` ;
done ) | sort -rnk2
echo TOTAL: `cat $T | awk '$5~/r/ {n=n+$9} END{print n}'`
rm $T

echo

echo In Queue:
T=`mktemp`;
qstat -u '*' > $T;
(for i in `cat $T | awk '$5=="qw"{print $4}' | sort -u`; do
    echo $i `cat $T | awk '$5=="qw" && $0~/'"$i"'/{if($NF==8){n=n+$8}else{split($9,s,"-",sep); from=int(s[1]); split(s[2],ss,":",sep); to=int(ss[1]); n=n+to-from+1; } } END{print n}'` ;
done ) | sort -rnk2
echo TOTAL: `cat $T | awk '$5=="qw"{if($NF==8){n=n+$8}else{split($9,s,"-",sep); from=int(s[1]); split(s[2],ss,":",sep); to=int(ss[1]); n=n+to-from+1; } } END{print n}'` ;
rm $T
