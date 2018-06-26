#!/bin/bash
if [[ $# -eq 0 ]]; then
    qsub MyScripts/submit.pbs
elif [[ $# -eq 1 ]]; then
    qsub MyScripts/submit.pbs $1
elif [[ $# -lt 3 ]]; then
    echo "Usage: $0 [cfg-file] [from_id to_id]"
    exit -1
else
    rsync -arv MyScripts/submit.pbs MyScripts_tmp/
    sed -i 's/#$ -t.*/#$ -t '"$2-$3/g" MyScripts_tmp/submit.pbs
    qsub MyScripts_tmp/submit.pbs $1
fi
