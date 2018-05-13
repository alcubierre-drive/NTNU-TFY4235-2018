#!/usr/bin/env python3
# coding: utf8

import lib
# this library is only used to make the plots look nice.
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import numpy as np
import scipy.stats as scst
import sys

PLOTS = True
SAVE = True

# read the data (which got stored in this very very precious file
data = np.genfromtxt(
        sys.argv[0].replace('.py','.dir/very_very_precious_data'),
        skip_header=1)

# make the saved data (that is a little messed up due to multicore issues)
# "human accessible"
def generate_human_data(data):
    myrange = list(set(np.asarray(data.T[0],int)))
    myrange.sort()
    myhumandata = []
    for l in myrange:
        selection = np.asarray(data.T[0],int) == l
        myhumandata.append(data.T[1:,selection])
    return np.asarray(myrange), np.asarray(myhumandata)

# actually transform the data
l_range, l_data = generate_human_data(data)

# further data transformations, now it is in a reasonable form for plotting
L, P = np.meshgrid(l_range,l_data[0][0])
RHO = (l_data[:,2]/l_data[:,1]).T
T = l_data[:,3].T

# prepare the data for the fit, and do the fit
def prepare_data_fit():
    if False:
        # just an initial try to do it for some p only... not used (thus "if
        # False")
        used_range = np.logical_and(P<0.75,P>0.55)
        my_size = P[:,0][np.logical_and(P[:,0] < 0.75,P[:,0] > 0.55)].shape
        L_f = L[used_range].reshape((my_size[0],20))
        P_f = P[used_range].reshape((my_size[0],20))
        T_f = T[used_range].reshape((my_size[0],20))
        RHO_f = RHO[used_range].reshape((my_size[0],20))
    else:
        L_f = L
        P_f = P
        T_f = T
        RHO_f = RHO
    PVEC = P_f.T[0]
    r_rho = np.zeros(PVEC.shape)
    r_t = np.zeros(PVEC.shape)
    for i,p in enumerate(P_f):
        my_l = L_f[i]
        my_rho = RHO_f[i]
        my_t = T[i]
        slope, intercept, r_rho[i], p_value, std_err = \
                scst.linregress(np.log(my_l),np.log(my_rho))
        slope, intercept, r_t[i], p_value, std_err = \
                scst.linregress(np.log(my_l),np.log(my_t))
        del slope, intercept, p_value, std_err
    return PVEC, r_rho**2, r_t**2

# get a special grid for plotting
def get_grid(space,ax):
    loc_maj = plticker.MultipleLocator(base=space)
    loc_min = plticker.MultipleLocator(base=space*0.5)
    ax.grid(which='both')
    ax.xaxis.set_major_locator(loc_maj)
    ax.xaxis.set_minor_locator(loc_min)

# make the plots, and save them (according to global vars)
if PLOTS and SAVE:
    fig_density = plt.figure(figsize=(10,5))
    ax_density = fig_density.add_subplot(111)
    d_cf = ax_density.contourf(P,L,RHO,30,cmap=plt.cm.viridis,alpha=0.8)
    ax_density.axvline(0.625,color='red')
    ax_density.set_xlabel("$p$")
    ax_density.set_ylabel("$L$")
    get_grid(0.1,ax_density)
    ax_density.set_title("Density: $p_c \\approx 0.625$")
    fig_density.colorbar(d_cf,ax=ax_density,label='$\\rho$')
    fig_density.tight_layout()

    fig_time = plt.figure(figsize=(10,5))
    ax_time = fig_time.add_subplot(111)
    t_cf = ax_time.contourf(P,L,T,30,cmap=plt.cm.viridis,alpha=0.8)
    ax_time.axvline(0.625,color='red')
    ax_time.set_xlabel("$p$")
    ax_time.set_ylabel("$L$")
    get_grid(0.1,ax_time)
    ax_time.set_title("Time: $p_c \\approx 0.625$")
    fig_time.colorbar(t_cf,ax=ax_time,label='$t$')
    fig_time.tight_layout()

    fig_density.savefig("Problem1.dir/plots/dens.pdf")
    fig_time.savefig("Problem1.dir/plots/time.pdf")

    P,rr,rt = prepare_data_fit()
    fig_fits = plt.figure(figsize=(6,4))
    plt.plot(P,rr,label='$\\rho$')
    plt.plot(P,rt,label='$t$')
    plt.xlabel('$p$')
    plt.ylabel('$r^2$')
    plt.grid()
    plt.legend()
    plt.title('Critical probability: $p_c \\approx 0.635$')
    plt.axvline(0.635,color='red')
    plt.tight_layout()
    fig_fits.savefig("Problem1.dir/plots/fits.pdf")
