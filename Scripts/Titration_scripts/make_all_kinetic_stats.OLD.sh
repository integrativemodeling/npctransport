#!/bin/bash
UPDATE=1 # If true, update an exisiting STATS file in current folder, skipping processed output files
imp="$HOME/imp_git/fast/setup_environment.sh"
imppy="$imp python"

ignore_regexp='DONT_IGNORE_FIRST_LINE' # Add header only if not updating
if [ $UPDATE -eq 1 ] ; then
    STATS_FILE=STATS;
    if [ -e $STATS_FILE ] && [ `cat $STATS_FILE|wc -l` -gt 0 ]; then
        ignore_regexp='#'
    else
        touch $STATS_FILE
    fi
else
    STATS_FILE=/dev/stdout
fi
touch $STATS_FILE
for output_file in "$@"; do
    if [ $UPDATE -eq 1 ]; then
        is_found=`grep $output_file STATS| wc -l`;
        if [ $is_found -gt 0 ]; then
            continue # skip if already processed
        fi
    fi
    $imppy MyScripts/kinetics_stats_by_temperature.py "$output_file" 0 | sed 's/[,\{\}:'"'"'\(\)]//g' | awk -f MyScripts/process_kinetic_stats.awk | grep -v $ignore_regexp | awk '{if ($0 ~ /#/){print $0, "output_file"} else {print $0, "'"$output_file"'"}}' | sed 's/^.*#//g' | sed 's/^ *//g' | sed 's/  */ /g' | sed 's/ *$//g' | grep -v BLAS | awk '$1~/[0-9]/ && $1>0 || $0~/nteraction/' >> $STATS_FILE
   ignore_regexp='#'
done
