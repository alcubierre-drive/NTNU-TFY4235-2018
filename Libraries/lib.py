#!/usr/bin/env python3
# coding: utf8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpltc
import scipy.odr as so
import multiprocessing

np.ALLOW_THREADS = multiprocessing.cpu_count()
plt.rcParams['text.latex.preamble']=['\\usepackage{lmodern}\n\\usepackage{amsmath}\n\\newcommand{\\un}[1]{\\,\\mathrm{#1}}\n\\newcommand{\\angstrom}{\\text{\\,\\AA}}']
params={'text.usetex' : True,
       'font.size' : 12,
       'font.family' : 'lmodern',
       'text.latex.unicode' : True
       }
plt.rcParams.update(params)

error_kwargs = {'ls':'', 'marker':'o', 'markersize':2, 'elinewidth':0.5, 'capsize':3}

# functions for the fits
def fit_both(datax,datay,errx,erry,func,beta0,print_info=False):
    dat = so.RealData(datax,datay,sx=errx,sy=erry)
    mod = so.Model(func)
    odr = so.ODR(dat,mod,beta0=beta0)
    res = odr.run()
    #
    if print_info:
        res.pprint()
    #
    sd_beta = np.sqrt(np.diag(res.cov_beta))
    beta = res.beta
    pull_x = res.xplus
    normed_residuals = np.sqrt((res.eps/erry)**2 + (res.delta/errx)**2)
    pull_y = -np.sign(res.eps)*normed_residuals
    return [ [beta,sd_beta], [pull_x,pull_y], [res.res_var,len(datax)-len(beta0)] ]

def fit_y(datax,datay,erry,func,beta0,print_info=False):
    dat = so.RealData(datax,datay,sy=erry)
    mod = so.Model(func)
    odr = so.ODR(dat,mod,beta0=beta0)
    odr.set_job(fit_type=2)
    res = odr.run()
    #
    if print_info:
        res.pprint()
    #
    sd_beta = np.sqrt(np.diag(res.cov_beta))
    beta = res.beta
    res_x = res.xplus
    res_y = res.eps
    reserr = erry
    return [ [beta,sd_beta], [res_x,res_y,reserr], [res.res_var,len(datax)-len(beta0)] ]

# some functions that reduce plotting-related code
def limedit(array,epsilon=0.05,lower=False):
    amin,amax = min(array),max(array)
    delta = epsilon*(amax-amin)
    if lower:
        return (amin,amax+delta)
    else:
        return (amin-delta,amax+delta)

def plotlims(xplot,yplot,lower=False,x_eps=0.05,y_eps=0.05):
    ylims = limedit(yplot,epsilon=y_eps,lower=lower)
    xlims = limedit(xplot,epsilon=x_eps,lower=False)
    return [xlims,ylims]

def help_arr(array,epsilon=0.05,lower=False,size=300):
    lims = limedit(array,epsilon=epsilon,lower=lower)
    new_arr = np.linspace(lims[0],lims[1],size)
    return new_arr

# peakfinder
def peak_ranges(arr,m=None):
    numpyarr=np.asarray(arr)
    if m is None: m=numpyarr.max()/np.sqrt(2)
    elif m<1 and numpyarr.max()>10: m*=numpyarr.max()
    w=np.where(numpyarr>m)[0]
    return np.split(w,1+np.where(w[1:]-w[:-1]!=1)[0])

def find_peaks_thin(x,y,xerr=0,m=None,p=0,ranges=None):
    if m is None: m=0.2*y.max()
    elif m<1: m*=y.max()
    def peaks(ranges):
        for r in ranges:
            r=r[y[r]>y[r].max()*p]
            xr=x[r]
            yr=y[r]
            ys=yr.sum()
            f=(xr*yr).sum()/ys
            sf=np.sqrt(((xerr*yr)**2).sum())/ys
            yield (f,sf)
    return tuple(np.asarray([k for k in peaks(ranges or peak_ranges(y,m))]).T)
