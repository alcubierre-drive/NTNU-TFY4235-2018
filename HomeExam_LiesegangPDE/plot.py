#!/usr/bin/env python3
# coding: utf8

import lib
# above library is just used to make the plots look nice and for a peakfinder
# (coded during some labwork in the second year of studies). Finding these peaks
# could have easily be done manually.

################################################################################
## INFORMATION                                                                ##
##  Throughout the program, there are several global statements that all look ##
##  alike:                                                                    ##
##  *                                                                         ##
##  * if True/False:                                                          ##
##  *     # SWITCH                                                            ##
##  *     [do something]                                                      ##
##  *                                                                         ##
##  They are used to control which plots should be generated. If one of these ##
##  "switches" is changed to True, the subsequent block creates and saves     ##
##  some plots.                                                               ##
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc_sig

def get_fname_test(Dfac, bdy, time):
    # function to simply generate the filenames used to save the data to run the
    # test
    fname = 'tests/Dfac0-%.3f_Dfac1-%.3f_a0-%.3f_b0-%.3f_t-%.3e.dat'
    return fname % (*Dfac, *bdy, time)

def plot_test(Dfac, bdy, time):
    # function that takes care of all the plotting needed for the tests
    fname = get_fname_test(Dfac, bdy, time)
    data = np.genfromtxt(fname)
    fig = plt.figure(figsize=(6,4))
    plt.title( ( ("$d_\\alpha = %.0e, d_\\beta = %.0e, a_0 = " + \
            "%.1f\,\mathrm{m^{-1}}, b_0 = %.1f\,\mathrm{m^{-1}}, " + \
            "t = %.0e\,\mathrm{s}$") % (4.e-7, 4.e-7*Dfac[0], *bdy, \
            0.1*time)).replace("e+0","\\cdot 10^").replace( \
            "e-0","\\cdot {10^-}^") )
    plt.plot(data.T[0],data.T[1],label='$\\alpha(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[2],label='$\\beta(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[3],label='$\\gamma(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[4],label='$\\sigma(\\chi,\\tau)$')
    plt.xlabel("$\\chi$")
    plt.ylabel("$f(\\chi,\\tau)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    fig.savefig('plots/tests_dfaca-%.1f_a0-%.1f_b0-%.1f_t-%.0e.pdf' % \
            (Dfac[1], *bdy, time) )

if False:
    # SWITCH
    # change False to True if plots for the tests should be generated

    # configurations used to run the tests. hardcoded.
    DfacTest1 = (0.5,1.0)
    bdyTest1 = (1.0,0.5)
    DfacTest2 = (2.0,1.0)
    bdyTest2 = (1.0,1.0)
    tvec = np.linspace(1.e4,1.e5,10)
    for t in tvec:
        plot_test(DfacTest1,bdyTest1,t)
        plot_test(DfacTest2,bdyTest2,t)

def plot_const(c0,time):
    # function that takes care of the plotting for the "real" problem, i.e. task
    # d)
    fname = 'const_sim/c0-%.3e_t-%.3e.dat' % (c0,time)
    data = np.genfromtxt(fname)
    fig = plt.figure(figsize=(6,4))
    plt.title( ("$c_0 = %.0e\,\mathrm{m^{-1}}, t = %.1e\,\mathrm{s}$" % \
            (c0,0.1*time)).replace("e+0","\\cdot 10^").replace( \
            "e-0","\\cdot {10^-}^") )
    plt.plot(data.T[0],data.T[1],label='$\\alpha(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[2],label='$\\beta(\\chi,\\tau)$')
    plt.plot(data.T[0],100*data.T[3],label='$100\\cdot\\gamma(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[4],label='$\\sigma(\\chi,\\tau)$')
    plt.xlabel("$\\chi$")
    plt.ylabel("$f(\\chi,\\tau)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    fig.savefig('plots/const_sim_c0-%.0e_t-%.1e.pdf' % (c0,time))

def peak_analysis(time):
    # function that does the analysis of the location of the peaks as in task
    # e). makes plots.

    # get data and find peaks
    fname_repl = "const_sim/c0-%.3e_t-%.3e.dat"
    c0s = [1.e-2,3.e-2,6.e-2]
    fnames = [fname_repl % (c0,time) for c0 in c0s]
    datas = [np.genfromtxt(x).T for x in fnames]
    peaks = [lib.find_peaks_thin(datas[i][0],datas[i][4],m=1.0,
        xerr=(datas[i][0][1]-datas[i][0][0])/np.sqrt(12)) for i in range(3)]

    # code for the first plot
    colors = ['red','blue','orange']
    fig = plt.figure(figsize=(6,4))
    for i in range(3):
        plt.plot(datas[i][0],datas[i][4],color=colors[i],
        label='$c_0 = %i\cdot 10^{-2}\,\mathrm{m^{-1}}$' % int(c0s[i]*100))
        for p in peaks[i][0]:
            plt.axvline(p,color=colors[i],ls='--',linewidth=1.0)
    plt.xlim(0.33,0.53)
    plt.legend()
    plt.xlabel('$\\chi$')
    plt.ylabel(('$\\sigma(\\chi,\\tau=%.1e)$' % \
            (0.1*time)).replace('e+0','\\cdot10^'))
    plt.grid()
    plt.title(('Positions of the peaks for different $c_0$ and $\\tau=%.1e$' % \
            (0.1*time)).replace('e+0','\\cdot10^'))
    plt.tight_layout()
    fig.savefig('plots/const_sim_peakpos_t-%.1e.pdf' % time)

    # calculate differences
    delta_chi = [np.diff(p[0]) for p in peaks]
    delta_chi_err = [(p[1]**2+np.roll(p[1],1)**2)[:-1] for p in peaks]
    nvecs = [np.array(range(1,len(dc)+1)) for dc in delta_chi]
    fig = plt.figure(figsize=(6,4))
    for i in range(3):
        plt.errorbar(nvecs[i],delta_chi[i],yerr=1.e5*delta_chi_err[i],color=colors[i],
        label='$c_0 = %i\cdot 10^{-2}\,\mathrm{m^{-1}}$' % int(c0s[i]*100),
        **lib.error_kwargs)
        plt.plot(nvecs[i],delta_chi[i],color=colors[i])
    plt.legend()
    plt.xlabel('$n$')
    plt.ylabel('$\Delta \chi_n$')
    plt.grid()
    plt.title(('Distance of the peaks for $\\tau=%.1e$' % \
            (0.1*time)).replace('e+0','\\cdot10^'))
    plt.tight_layout()
    fig.savefig('plots/const_sim_peakdiff_t-%.1e.pdf' % time)

if False:
    # SWITCH
    # change False to True if plots for the task d) and e) should be done

    # configuration for these runs. hardcoded.
    c0s = [1.e-2,3.e-2,6.e-2]
    tvec = np.linspace(1.e5,2.e6,20)
    for c0 in c0s:
        for t in tvec:
            plot_const(c0,t)
    peak_analysis(2.e6)
    peak_analysis(1.7e6)
    peak_analysis(1.3e6)
    peak_analysis(1.e6)

def plot_var(s0,time,limedit=False):
    # function that takes care of the plotting for the "real" problem, i.e. task
    # d)
    fname = 'var_sim/s0-%.3e_t-%.3e.dat' % (s0,time)
    data = np.genfromtxt(fname)
    fig = plt.figure(figsize=(6,4))
    plt.title( ("$s_0 = %.2f\,\mathrm{m^{-1}}, t = %.1e\,\mathrm{s}$" % \
            (s0,0.1*time)).replace("e+0","\\cdot 10^").replace( \
            "e-0","\\cdot {10^-}^") )
    plt.plot(data.T[0],data.T[1],label='$\\alpha(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[2],label='$\\beta(\\chi,\\tau)$')
    plt.plot(data.T[0],100*data.T[3],label='$100\\cdot\\gamma(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[4],label='$\\sigma(\\chi,\\tau)$')

    if limedit:
        # get the plotting limits (to only plot the relevant part)
        idx_where = np.where(data.T[4]>0.5)[0]
        idx = [data.T[0][min(idx_where)],data.T[0][max(idx_where)]]
        deltax = (idx[1]-idx[0])
        idx[0] -= 1.5*deltax
        idx[1] += 0.5*deltax
        plt.xlim(idx[0],idx[1])

    plt.xlabel("$\\chi$")
    plt.ylabel("$f(\\chi,\\tau)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    fig.savefig('plots/var_sim_s0-%.0e_t-%.1e.pdf' % (s0,time))

if False:
    # SWITCH
    # simple switch to control the plotting for non-constant D. change to true
    # if plots should be made.

    for s0 in [100,10,1,0.1,0.01]:
        for t in [4.e5,8.e5,1.2e6,1.6e6,2.e6]:
            plot_var(s0,t)

def plot_stoichio(a,b,which,time):
    # function that takes the two stoichiometric coefficients and the selection
    # for the function R(a,b) used in the simulation and makes all plots
    fname = 'st_sim/%s_st-a-%i_st-b-%i_t-%.3e.dat' % (which,a,b,time)
    data = np.genfromtxt(fname)
    fig = plt.figure(figsize=(6,4))
    func_table = {
            'simple' : 'R(a,b) = R\\,a(x,t)b(x,t)',
            'powers' : 'R(a,b) = \\bar{R}\\,a^\\mathsf{A}(x,t)b^\\mathsf{B}(x,t)'
            }
    plt.title( ("$%s,\\mathsf{A}=%.1f,\\mathsf{B}=%.1f,t=%.1e\\,\\mathrm{s}$" %\
            (func_table[which],a,b,0.1*time)).replace('e+0','\\cdot10^') )
    plt.plot(data.T[0],data.T[1],label='$\\alpha(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[2],label='$\\beta(\\chi,\\tau)$')
    plt.plot(data.T[0],100*data.T[3],label='$100\\cdot\\gamma(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[4],label='$\\sigma(\\chi,\\tau)$')
    if max(data.T[4][150:350]) < 10:
        plt.ylim(-0.5,10.5)
    plt.xlabel("$\\chi$")
    plt.ylabel("$f(\\chi,\\tau)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    fig.savefig('plots/st_sim_%s_st-a-%i_st-b-%i_t-%.1e.pdf' % (which,a,b,time))

if False:
    # SWITCH
    # simple switch to control the plotting for different stoichiometric models,
    # change to true if plots should be made.

    for st_ab in [(2.,1.),(1.,10.),(10.,1.)]:
        for which in ['powers','simple']:
            for t in np.asarray([0.4,0.8,1.2,1.6,2.0])*1.e6:
                plot_stoichio(*st_ab,which,t)

def plot_parameter_range(fname,save=False):
    # function that plots the simulations made for task f), reads the parameters
    # from the filename (in order to use it with glob()
    def read_params(fname):
        params = []
        for i,s in enumerate(fname):
            if fname[i:i+3] in ['Da-','Db-','Dc-','a0-','b0-']:
                params.append(float(fname[i+3:i+6]))
        return params
    parms = read_params(fname)
    title = ('$f_a = %.1f, f_b = %.1f, f_c = %.1f, a_0 = %.1f\,' + \
            '\mathrm{m^{-1}}, b_0 = %.1f\,\mathrm{m^{-1}}$') % \
            tuple(parms)
    data = np.genfromtxt(fname)
    fig = plt.figure(figsize=(6,4))
    plt.title(title.replace('e+0','\cdot10^'))
    plt.plot(data.T[0],data.T[1],label='$\\alpha(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[2],label='$\\beta(\\chi,\\tau)$')
    plt.plot(data.T[0],100*data.T[3],label='$100\\cdot\\gamma(\\chi,\\tau)$')
    plt.plot(data.T[0],data.T[4],label='$\\sigma(\\chi,\\tau)$')
    plt.xlabel("$\\chi$")
    plt.ylabel("$f(\\chi,\\tau)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    if save:
        fig.savefig(fname.replace('t-2.000e+06.dat','.pdf').replace('D_sim/', \
                'plots/D_sim_'))
    else:
        plt.show()

if True:
    # SWITCH
    # simple switch to control the plotting for a lot of different
    # configurations of the main equation.

    import glob
    import warnings
    import sys
    warnings.filterwarnings("ignore")
    for i,f in enumerate(glob.glob('D_sim/*.dat')):
        plot_parameter_range(f,save=True)
        sys.stdout.write('\rProgress: %.1f %%' % (100*(i+1)/150))
        sys.stdout.flush()
