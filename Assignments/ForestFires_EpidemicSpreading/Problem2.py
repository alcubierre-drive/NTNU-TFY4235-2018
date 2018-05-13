#!/usr/bin/env python3
# coding: utf8

import lib
# this library is only used to make the plots look nice
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

PREFIX = sys.argv[0].replace('.py','.dir/data_plot/')
lambdas = [0.0,0.5,1.0]
qs = [0.0,0.5,1.0]
pvec = np.linspace(0.1,1.01,101)[:-2]

# just some function to generate the filename with the scheme that is used
def gen_fname(prefix,p,q,lam):
    returnstr = prefix+"p%.3f_q%.3f_lambda%.3f.dat"
    return returnstr % (p,q,lam)

# compact function to get the time-average of an array. Uses every 1000th value,
# thus the reshape.
def time_average(all_array):
    return np.mean(all_array[10000:].reshape((1000,90))[0])

# small class that handles each dataset.
class dataset:
    def __init__(self,p,q,lam):
        fname = gen_fname(PREFIX,p,q,lam)
        self.all_data = np.genfromtxt(fname)
        parsestring = os.path.basename(fname)
        self.p = p
        self.q = q
        self.lam = lam
        self.rho_1_eq_avg = time_average(self.all_data.T[1])
    def get_plotdata(self):
        self.t = self.all_data.T[0]
        self.rho_1 = self.all_data.T[1]
        self.rho_2 = self.all_data.T[2]

# doing the plots for task 3.5
def plot_3_5(pvec,q,lam):
    fig = plt.figure(figsize=(8,6))
    dsets = [dataset(p,q,lam) for p in pvec]
    for d in dsets:
        d.get_plotdata()
        plt.plot(d.t,d.rho_1,label='$p = %.3f$' % d.p)
    plt.grid()
    plt.xlabel('$t$')
    plt.ylabel('$\\rho_1$')
    plt.loglog()
    plt.xlim(1e1,1e5)
    plt.ylim(0.5e-2,1.0)
    plt.title('$q = %.3f, \\lambda = %.3f$' % (q,lam))
    plt.legend()
    plt.tight_layout()
    return fig

# doing the plots for task 3.7
def plot_3_7(pvec,qs,lambdas):
    fig = plt.figure(figsize=(6,4))
    for q in qs:
        for lam in lambdas:
            rho_1_eq_avg = [dataset(p,q,lam).rho_1_eq_avg for p in pvec]
            plt.plot(pvec,rho_1_eq_avg,label='$q = %.3f, \\lambda = %.3f$' % (q,lam))
    plt.grid()
    plt.xlabel('$p$')
    plt.ylabel('$\\bar{\\rho_1}$')
    plt.legend()
    plt.tight_layout()
    return fig

# doing the plots for task 3.9
def plot_3_9():
    prefix = sys.argv[0].replace('.py','.dir/')+'data_time/'
    files = os.listdir(prefix)
    dat = np.genfromtxt(prefix+files[0])
    fig_time = plt.figure(figsize=(6,4))
    plt.plot(dat.T[0],dat.T[2])
    plt.ylim(0,45)
    plt.xlabel('$p$')
    plt.ylabel('$\\bar{t}_\\mathrm{max}$')
    plt.grid()
    plt.title('Extinction time')
    plt.tight_layout()
    fig_rho = plt.figure(figsize=(6,4))
    plt.plot(dat.T[0],dat.T[1])
    plt.xlabel('$p$')
    plt.ylabel('$\\bar{\\rho}_2$')
    plt.grid()
    plt.title('Density $\\bar{\\rho}_2$')
    plt.tight_layout()
    return fig_time, fig_rho

# some lazy stuff to generate all the plots needed and save them accordingly.
# to be (un)-commented when needed
out_prefix = sys.argv[0].replace('.py','.dir/plots/')
figs = [None]*12
names = [
        '3_5_q0.0_l0.0.pdf','3_5_q0.5_l0.0.pdf','3_5_q1.0_l0.0.pdf',
        '3_5_q0.0_l0.5.pdf','3_5_q0.5_l0.5.pdf','3_5_q1.0_l0.5.pdf',
        '3_5_q0.0_l1.0.pdf','3_5_q0.5_l1.0.pdf','3_5_q1.0_l1.0.pdf',
        '3_7.pdf','3_9_t.pdf','3_9_rho.pdf'
        ]
figs[0] = plot_3_5(pvec[::10],0.0,0.0)
figs[1] = plot_3_5(pvec[::10],0.5,0.0)
figs[2] = plot_3_5(pvec[::10],1.0,0.0)
figs[3] = plot_3_5(pvec[::10],0.0,0.5)
figs[4] = plot_3_5(pvec[::10],0.5,0.5)
figs[5] = plot_3_5(pvec[::10],1.0,0.5)
figs[6] = plot_3_5(pvec[::10],0.0,1.0)
figs[7] = plot_3_5(pvec[::10],0.5,1.0)
figs[8] = plot_3_5(pvec[::10],1.0,1.0)
#figs[9] = plot_3_7(pvec,[0.5,1.0],[0.5,1.0])
#figs[10], figs[11] = plot_3_9()

for i,n in enumerate(names):
    figs[i].savefig(out_prefix+names[i])
