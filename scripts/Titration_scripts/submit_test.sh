#!/bin/bash
rsync -arv MyScripts/ MyScripts_tmp/
if [ $# -gt 0 ] ; then
    qsub MyScripts_tmp/_submit_test.pbs "$@"
else
    qsub MyScripts_tmp/_submit_test.pbs "."
fi
