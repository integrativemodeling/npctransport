#!/usr/bin/python
from multiprocessing import Pool, TimeoutError
from glob import glob
import os
from os.path import expanduser
from collections import Counter
import subprocess
import sys

def count_lines_matching(lines, pattern):
  lines = Counter([line for line in lines if pattern in line])
  return len(lines)

def process_folder(folder):

  print("Processing folder {}".format(folder))
  cwd = os.getcwd()
  os.chdir(folder)
  IS_UPDATE= True # If true, update an exisiting STATS file in current folder, skipping processed output files
  home= expanduser("~")
  imp= "{home}/imp_git/fast/setup_environment.sh".format(home= home)
  imppy= "{imp} python".format(imp= imp)
  ignore_regexp= 'DONT_IGNORE_FIRST_LINE' # Add header only if not updating
  if IS_UPDATE:
    STATS_FILE='STATS';
    old_stats_lines=[]
    if os.path.exists(STATS_FILE) and os.path.getsize(STATS_FILE)>0:
        ignore_regexp= '#'
        fid= open(STATS_FILE, 'r')
        for line in fid:
          old_stats_lines.append(line)
        fid.close()
    else:
        fid= open(STATS_FILE, 'w')
        fid.close()
  else:
    STATS_FILE= '/dev/stdout/'
  for output_file in glob("output*.pb"):
    if IS_UPDATE:
        if count_lines_matching(old_stats_lines, output_file) > 0:
            continue # skip if already processed
    cmd1= "{imppy} MyScripts/kinetics_stats_by_temperature.py {output_file} 0 > STATS.{output_file}.txt" \
       .format(imppy=imppy, output_file=output_file)
    cmd2= "cat STATS.{}.txt".format(output_file)  \
       + " | sed \"s/[,:'\{\}\(\)]//g\"" \
       + " | awk -f MyScripts/process_kinetic_stats.awk" \
       + " | grep -v '{}'".format(ignore_regexp) \
       + " | awk '{{ if ($0 ~ /#/) {{print $0, \"output_file\"}} else {{print $0, \"{output_file}\"}} }}'" \
       .format(output_file=output_file) \
       + " | sed 's/^.*#//g' | sed 's/^ *//g' | sed 's/  */ /g' | sed 's/ *$//g'" \
       + " | grep -v BLAS | awk '$1~/[0-9]/ && $1>0 || $0~/nteraction/' | awk 'NF>5'" \
       + " >> {}".format(STATS_FILE)
    try:
      subprocess.check_call(cmd1 + " && " + cmd2, shell= True)
    except:
      print("Error processing {}/{}".format(folder, output_file))
      continue
#    subprocess.check_call(cmd2, shell= True)
    ignore_regexp= '#'
  os.chdir(cwd)
  return 0

#folders=glob(sys.argv[1] + "*/")
folders=sys.argv[1:]
print(folders)
pool= Pool(processes=12)
a=pool.map(process_folder, folders)
