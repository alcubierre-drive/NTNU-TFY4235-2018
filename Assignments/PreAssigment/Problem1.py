#!/usr/bin/env python3
# coding: utf8

#import lib
import matplotlib.pyplot as plt
import numpy as np
import ctypes
import sys
clib = ctypes.cdll.LoadLibrary(sys.argv[0].replace('.py','.dir/')+'calcs.so')
MAXDATA = ctypes.cast(clib.MAXDATA, ctypes.POINTER(ctypes.c_int)).contents.value
MAXTIME = ctypes.cast(clib.MAXTIME, ctypes.POINTER(ctypes.c_int)).contents.value

def get_all_times_array(dist='u'):
    get_times = clib.get_times_array
    get_times.argtypes = [ctypes.c_char]
    get_times.restype = ctypes.POINTER(ctypes.c_int)
    RESULT = np.ctypeslib.as_array(
            get_times(ord(dist)),
            shape=(MAXDATA,)
            )
    return RESULT

def plot(dist='u',savefig=False):
    raw_data = get_all_times_array(dist)
    mod_data = raw_data[raw_data > 1]

    bins = np.logspace(np.log10(2),np.log10(MAXTIME),25)
    widths = bins[1:]-bins[:-1]
    hist = np.histogram(mod_data,bins=bins)
    hist_norm = hist[0]/widths

    if savefig:
        import lib
    fig = plt.figure(figsize=(6,4))

    plt.bar(bins[:-1],hist_norm,widths,alpha=0.8,color='blue',label='data')
    plt.plot(
            bins,
            2.e-3*MAXDATA**1.5*bins**-1.5,
            color='red',
            label='$N(t_\\mathrm{max})^{\\alpha_\\mathrm{th}}$'
    )

    if dist == 'u':
        plt.title(
            'Uniform: $\\alpha_\\mathrm{th} = 1.5$ ($%i$ data points)' % MAXDATA
            )
        savename = 'uniform.pdf'
    else:
        plt.title(
            'Gaussian: $\\alpha_\\mathrm{th} = 1.5$ ($%i$ data points)' % MAXDATA
            )
        savename = 'gauss.pdf'
    plt.loglog()
    plt.grid()
    plt.xlabel('$t_\\mathrm{max}$')
    plt.ylabel('number in a.u.')
    plt.legend()
    plt.tight_layout()
    if savefig:
        fig.savefig(savename)

if __name__ == '__main__':
    plot('u',True)
    plot('g',True)
    #plt.show()
