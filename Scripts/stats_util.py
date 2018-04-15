from __future__ import print_function
import numpy as np
import scipy.stats
import math

def _check_autocovariance_input(x):
        if len(x) < 2:
            raise ValueError('Need at least two elements to calculate autocovariance')

def get_autocovariance(x):
    ''' see https://stackoverflow.com/questions/20110590/how-to-calculate-auto-covariance-in-python?answertab=active#tab-top '''
    _check_autocovariance_input(x)
    x_centered = x - np.mean(x)
    return np.correlate(x_centered, x_centered, mode='full')[len(x) - 1:] / len(x)

def get_time_series_mean_and_sem(y):
    '''
    estimate mean and sem of a time series using
    https://stats.stackexchange.com/questions/274635/calculating-error-of-mean-of-time-series

    params:
    y - array-like time series values, assuming constant time intertvals
    '''
    n= float(len(y))
    sqrt_n= int(math.ceil(math.sqrt(n)))
    n2= n*n
    autocov= get_autocovariance(y)
    e= min(sqrt_n, len(autocov)+1)
    varmean=autocov[0]/n
    for k in range(1,e):
        varmean += (2*autocov[k]*(n-k))/n2
    return np.mean(y), math.sqrt(varmean)

if __name__ == "__main__":
        a= [9,3,0,-1,-3,-3, -1, 0, 1, 3,7, 9, 11, 9, 10, 4, 2, 2, 0, -3 , -5 ,-9, -11, -7, -9, -12, -14, -16, -14, -11, -8, -9, -5, 0, 3, 4,7, 4, 9, 11, 9 , 11, 8, 13,8 , 3 , 0, -1, 1, -2, 1, 3, 1, -5, -8, -9, 1, -3, 13, 4, 1, 3, -1]
        mean, sem= get_time_series_mean_and_sem(a)
        print("Autocovariance:")
        for x in get_autocovariance(a):
                print("{:.2f} ".format(x), end='')
        print()
        print("Mean and sem of time-series: {:.3f} {:.3f}".format(mean, sem))
        print("Mean and sem of independent samples: {:.3f} {:.3f}".format(np.mean(a), scipy.stats.sem(a)))
