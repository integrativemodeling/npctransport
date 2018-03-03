#!/usr/bin/env python
from IMP.npctransport import *
import sys
fin=open(sys.argv[1], "rb")
output= Output()
output.ParseFromString(fin.read())
fin.close()
output.assignment.output_statistics_interval_frames=500000
for obstacle in output.assignment.obstacles:
    obstacle.radius.value= obstacle.radius.value*1.36
fout=open(sys.argv[2], "wb")
fout.write(output.SerializeToString())
fout.close()
