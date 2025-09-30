k=$1;
s=$2;
(for i in config_v13_${s}_*; do
    for j in `seq $k 900 7200`; do
        echo -n "$i $j ";
        grep -n "t$j.pb" $i/STATS | awk '{print $7, $(NF-3), $(NF-6), $NF, $1, $2, $3, $4}';
        echo;
    done;
 done;) | grep pb | awk '{printf "%-60s %5d %8.3f [mM] %5.2f%% %s", $1, $2, $3*1000.0, $4*100, $5*100, $6; printf "\n"}' | sed 's/^config.*M_kapC\([0-9e\.-]*\)M/\1/g' | sort -gk1| awk '{print "[" $1, ",", $5/100., ",", $6/100., ",", $3, "]," }' | grep -v '\-1'
