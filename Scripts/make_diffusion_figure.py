#!/usr/bin/env python
from __future__ import print_function
from IMP.npctransport import *
import numpy as np
import re
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import sys
import matplotlib
matplotlib.use('PDF')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from multiprocessing import cpu_count, Pool

def get_diffuser_data_from_file(fname):
    f=open(sys.argv[1], "rb")
    config= Output()
    config.ParseFromString(f.read())
    s= config.statistics
    diffuser_data= {}
    for f in s.floaters:
        is_first= True
        x= []
        y= []
        w= []
        for od in f.order_params:
            if is_first:
                prev_time_ns= od.time_ns
                is_first= False
                continue
            #        print(f.type, od.time_ns, od.diffusion_coefficient, od.interacting_fraction)
            dT= od.time_ns - prev_time_ns
            prev_time_ns= od.time_ns
            D_A_squared_per_fs= od.diffusion_coefficient
            D_10_to_minus7_cm2_per_sec= D_A_squared_per_fs * 1E+6  # 1E+7 * 1E-16 cm^2/A^2 * 1E+15 fs/sec = 1E+6
            if(np.isnan(od.interacting_fraction) or \
               np.isnan(D_10_to_minus7_cm2_per_sec) or
               D_10_to_minus7_cm2_per_sec < 1E-15):
                continue
            x.append(od.interacting_fraction)
            y.append(od.diffusion_coefficient*1E+6)
            w.append(dT)
        diffuser_data[f.type]={'x':x, 'y':y, 'w':w}
    print("Finished proceccing", fname)
    return diffuser_data

def get_diffuser_data_from_files(fnames):
    dds=[]
    ncores= cpu_count()
    pool= Pool(ncores)
    dds= pool.map(get_diffuser_data_from_file, fnames)
    dd_all= {}
    for dd in dds:
        for key in dd.iterkeys():
            if key in dd_all:
                dd_all[key]['x'].extend(dd[key]['x'])
                dd_all[key]['y'].extend(dd[key]['y'])
                dd_all[key]['w'].extend(dd[key]['w'])
            else:
                dd_all[key]= dd[key]
    return dd_all

def get_diffusion_coefficient_estimates(diffuser_type, diffuser_data, is_plot=True):
    '''
    get diffusion coefficient from diffuser_data, a dictionary with keys
    'x'=fraction bound,
    'y'=diffusion coefficient,
    'w'=weight probably based on statistics time
    '''
    x= diffuser_data['x']
    y= diffuser_data['y']
    w= diffuser_data['w']
    n= len(x)
    assert(len(y)==n and len(w)==n)
    x= np.array(x).reshape(n,1)
    y= np.array(y).reshape(n,1)
    w= np.array(w).reshape(n,1)
    print("Y shape", y.shape)
    if(len(np.unique(x))==1):
        print("Can't regress on a single unique x value", np.unique(x),
              "in type", diffuser_type)
        return None
    regr= linear_model.LinearRegression()
    print(len(x), len(y), len(w))
    regr.fit(x, y)#,w)
    y_pred = regr.predict(x)
    print("Y pred shape", y_pred.shape)
    print('Coefficients: \n', regr.coef_)
    print('Intercept: \n', regr.intercept_)
    print("Mean squared error: %.2f"
          % mean_squared_error(y, y_pred))
    print('Variance score: %.2f' % r2_score(y, y_pred))
    if is_plot:
        plt.figure()
        plt.scatter(x, y, color='black')
        plt.plot(x, y_pred, color='blue', linewidth=2)
        #    plt.xticks(())
        #    plt.yticks(())
        plt.xlabel('fraction bound')
        plt.ylabel('$D\ [10^{-7} cm^2/sec]$')
        plt.title(diffuser_type)
        plt.savefig(pp, format='pdf')
    # get bound and unbound estimates
    xx= np.array([0.0, 1.0]).reshape(2,1)
    yy= regr.predict(xx)
    return yy
#   plt.show()


############# main ###########3
pp= PdfPages('output.pdf')
if len(sys.argv)<=1:
    print("Usage", sys.argv[0]," output_file_1 [output_file2] ...")
    sys.exit(-1)
dd_all= get_diffuser_data_from_files(sys.argv[1:])
unbound_d={}
bound_d={}
for diffuser_type, diffuser_data in dd_all.iteritems():
    yy= get_diffusion_coefficient_estimates(diffuser_type,
                                            diffuser_data,
                                            is_plot=True)
    if yy is None:
        continue
    unbound_d[diffuser_type]= yy[0]
    bound_d[diffuser_type]= yy[1]
print(unbound_d)
print(bound_d)
x=[]
y0=[]
y1=[]
for kap_type in unbound_d.iterkeys():
    m= re.match('kap([0-9]*)$', kap_type)
    if not m:
        print(kap_type, " is not a kap")
        continue
    radius_A= int(m.group(1))
    MW_kDa= (radius_A/6.6)**3
    x.append(MW_kDa) # in nanometers
    y0.append(unbound_d[kap_type])
    y1.append(bound_d[kap_type])
fig=plt.figure()
plt.scatter(x,y0, color='r', label='outside NPC (cyto/nuc)')
plt.scatter(x,y1, color='b', label='inside NPC')
plt.xticks()
plt.yticks()
#plt.title('diffusion coefficient inside and outside the pore')
plt.xlabel('MW [kDa]')
plt.ylabel('$D\ [10^{-7} cm^2/sec]$')
plt.xlim(xmin=0.0)
plt.ylim(ymin=0.0)
plt.legend()
ax2= plt.gca().twinx()
barwidth=2.5
ax2.bar(x,
        [yy1/yy0 for yy0,yy1 in zip(y0,y1)],
        barwidth,
        label='in/out ratio',
        alpha=0.4,
        color='g')
ax2.set_ylabel('in/out diffusion ratio',
               color='g')
ax2.set_ylim(ymax=1.0)
plt.savefig(pp, format='pdf')
pp.close()
