#!/usr/bin/bash /Users/barak/imp/fast/setup_environment.sh /anaconda/bin/python

#########################################
# Example: imppy Scripts/analyze_anchor_distrance.py InputData/wholeNPC_0.rmf3
#
# Author: Barak Raveh barak@salilab.org
# Data written: Feb 14, 2019
# Date last udpated: Feb 14, 2019
#########################################
from IMP.npctransport import *
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.algebra
import RMF
import math
import re
import sys
import numpy as np
import scipy.spatial.distance
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import pandas as pd

input_rmf = sys.argv[1] #'wholeNPC_0.rmf3')
if len(sys.argv)>2:
   n_spokes= int(sys.argv[2])
else:
  n_spokes= 1
if len(sys.argv)>3:
  csv_output_file= sys.argv[3]
fgs= {}
fgs_anchors_regexp={}
fgs_anchors_regexp['Nsp1.*_601-636_']='Nsp1'
fgs_anchors_regexp['Nup1.*_301-350_']='Nup1'
fgs_anchors_regexp['Nup42.*_364-413_']='Nup42'
fgs_anchors_regexp['Nup49.*_201-269_']='Nup49'
fgs_anchors_regexp['Nup57.*_201-286_']='Nup57'
fgs_anchors_regexp['Nup60.*_351-398_']='Nup60'
fgs_anchors_regexp['Nup100.*_551-575_']='Nup100'
fgs_anchors_regexp['Nup116.*_751-775_']='Nup116'
fgs_anchors_regexp['Nup145.*_201-225_']='Nup145'
fgs_anchors_regexp['Nup159.*_1082-1116_']='Nup159'


def add_fg_anchors_to_coords(fg_anchors_to_coords_map, type_name, anchor_coordinates):
    '''
    Update fg_anchors_to_coords_map with a mapping
    from <type_name.i> to anchor_coordinates[i]
    '''
    for i, coordinates in enumerate(anchor_coordinates):
        anchor_tag= "{}.{:d}".format(type_name, i)
        print("{} {:.2f} {:.2f} {:.2f}".format(anchor_tag, *coordinates))
        fg_anchors_to_coords_map[anchor_tag] = coordinates

def handle_xyz_children(parent):
    for p in parent.get_children():
        for anchor_name,fg_name in fgs_anchors_regexp.iteritems():
           if not re.search(anchor_name, p.get_name()):
               continue
           xyzr= IMP.core.XYZR(p)
           coords= xyzr.get_coordinates()
           for i in range(n_spokes):
           # print "FG-Anchor_"+fgs_anchors_regexp[anchor_name]
              R=IMP.algebra.get_rotation_about_normalized_axis([0,0,1],
                  i*math.pi/4.0)
              coords_i= R*coords
              fg_name_spoke= "S{:d}.{}".format(i, fg_name)
              if fg_name_spoke in fgs:
                  fgs[fg_name_spoke].append(coords_i)
              else:
                  fgs[fg_name_spoke]=[coords_i]
           break

def handle_representation(r):
#    print r.get_name()
    if r.get_name()=="Beads":
        handle_xyz_children(r)
    if re.search("Res:10$", r.get_name()):
        for rr in r.get_children():
            handle_xyz_children(rr)

def get_distance_matrix(name_to_coords_map):
  ordered_keys= sorted(name_to_coords_map.keys())
  v_list= []
  for key in ordered_keys:
    v_list.append(name_to_coords_map[key])
  v_ndarray= np.vstack( v_list )
  print(v_ndarray)
  D= scipy.spatial.distance.pdist(v_ndarray)
  D= scipy.spatial.distance.squareform(D)
  return D, ordered_keys

def get_positions_and_labels(ordered_keys, include_copies= True):
  positions= []
  labels= []
  for i,key in enumerate(ordered_keys):
    #    r= re.search("S[0-9]*\.(.*)\.0$", key)
    r= re.search("S[0-9]*\.(.*\.[0-9]*)$", key)
    if r is not None:
      positions.append(i-0.5)
      labels.append(r.group(1))
  return positions, labels

# ********* MAIN: *********
# Load information from RMF:
f=RMF.open_rmf_file_read_only(input_rmf)
m=IMP.Model()
h=IMP.rmf.create_hierarchies(f,m)
IMP.rmf.load_frame(f,0)
for nup in h[0].get_children():
    if not re.search("@", nup.get_name()) or re.search("@11$", nup.get_name()):
        for r in nup.get_children():
            handle_representation(r)
fg_anchors_to_coords_map= {}
for (name,coords) in fgs.iteritems():
    add_fg_anchors_to_coords(fg_anchors_to_coords_map, name, coords)
[D_angstroms, ordered_keys]= get_distance_matrix(fg_anchors_to_coords_map)
D_nm= D_angstroms / 10.0

# Output CSV:
try:
  df = pd.DataFrame(D_nm,
                  columns= ordered_keys)
  df.to_csv(csv_output_file,
            index=False)
except NameError as e:
  pass

# Plot:
fig, ax = plt.subplots(1,1)
im= ax.imshow(D_nm,
             cmap=plt.cm.afmhot_r,
             vmin= 0.0,
             vmax= 20.0)
positions, labels= get_positions_and_labels(ordered_keys)
plt.xticks(positions, labels)
plt.yticks(positions, labels)
ax.tick_params(axis='x', labelrotation=90)
#ax.tick_params(axis='y', labelrotation=-45)
cb_ax= plt.colorbar(im, ax=ax)
cb_ax.set_label('d [nm]')


D= D_nm
condensedD= scipy.spatial.distance.squareform(D, force='tovector')
# Compute and plot first dendrogram.
fig = plt.figure(figsize=(12, 12))
ax1 = fig.add_axes([0.09,0.1,0.2,0.5])
Y = sch.linkage(condensedD, method='centroid')
Z1 = sch.dendrogram(Y, orientation='left')
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.4,0.71,0.5,0.2])
Y = sch.linkage(condensedD, method='single')
Z2 = sch.dendrogram(Y)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.4,0.1,0.5,0.5])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1,:]
D = D[:,idx2]
im = axmatrix.matshow( D,
                       aspect='auto',
                       origin='lower',
                       cmap=plt.cm.afmhot_r, #YlGnBu,
                       vmin= 0.0,
                       vmax=20.0)
ticks= np.linspace( 0,
                    len(ordered_keys),
                    num= len(ordered_keys),
                    endpoint= False )
print(ticks)
print( pd.Series(ordered_keys)[idx1] )
axmatrix.set_xticks( ticks )
axmatrix.set_xticklabels( pd.Series(ordered_keys)[idx1] )
axmatrix.set_yticks( ticks )
axmatrix.set_yticklabels( pd.Series(ordered_keys)[idx2] )
axmatrix.tick_params( axis='x', labelrotation=90 )

# Plot colorbar.
axcolor = fig.add_axes( [0.91,0.1,0.02,0.6] )
cb= plt.colorbar( im, cax=axcolor )
cb.set_label('d [nm]')
fig.show()
fig.savefig('dendrogram.png')
plt.show()
