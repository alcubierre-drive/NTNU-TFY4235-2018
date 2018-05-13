#!/usr/bin/env python3
# coding: utf8

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sc_la
from helpers import initialValues, delta, simulation, make_bc
import sys

'''
the following functions just implement the different linear schemes, i.e.
upwind, downwind, explicit_centered, implicit_centered, lax_friedrichs and
lax_wendroff. The solution is *not* obtained by a linear solver but rather by
matrix inversion because after testing, the linear solver had a smaller range of
initial conditions to work with (due to over/underflow).
'''
def upwind(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF
    t = ival.N_T
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)*(1.0-alpha)) + \
            np.diag(np.ones(ival.N_X-1),-1)*alpha,
            ival.BoundaryType)
    result = np.linalg.matrix_power(matrix,t) @ initial
    return result

def downwind(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF
    t = ival.N_T
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)*(1.0+alpha)) - \
            np.diag(np.ones(ival.N_X-1),1)*alpha,
            ival.BoundaryType)
    result = np.linalg.matrix_power(matrix,t) @ initial
    return result

def explicit_centered(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF/2.0
    t = ival.N_T
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)) - \
            np.diag(np.ones(ival.N_X-1),1)*alpha + \
            np.diag(np.ones(ival.N_X-1),-1)*alpha,
            ival.BoundaryType)
    result = np.linalg.matrix_power(matrix,t) @ initial
    return result

def implicit_centered(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF/2.0
    t = ival.N_T
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)) + \
            np.diag(np.ones(ival.N_X-1)*alpha,1) - \
            np.diag(np.ones(ival.N_X-1)*alpha,-1),
            ival.BoundaryType)
    result = np.linalg.matrix_power(np.linalg.inv(matrix),t) @ initial
    return result

def lax_friedrichs(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF
    t = ival.N_T
    matrix = make_bc(
            0.5*(np.diag(np.ones(ival.N_X-1)*(1-alpha),k=1) + \
            np.diag(np.ones(ival.N_X-1)*(1+alpha),-1)),
            ival.BoundaryType)
    return np.linalg.matrix_power(matrix,t) @ initial

def lax_wendroff(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF
    t = ival.N_T
    matrix = make_bc(
            (1.-alpha**2)*np.diag(np.ones(ival.N_X)) + \
            np.diag(np.ones(ival.N_X-1)*0.5*(alpha**2-alpha),1) + \
            np.diag(np.ones(ival.N_X-1)*0.5*(alpha**2+alpha),-1),
            ival.BoundaryType)
    return np.linalg.matrix_power(matrix,t) @ initial

'''
The following two functions are used to implement the lax_wendroff scheme for
the Hopf Equation, the first function is iterated ´t´ times over in the second
function.
'''
def hopfeq_scheme_once(initial,alpha=1.0):
    return initial - alpha/4*(np.roll(initial,1)**2-np.roll(initial,-1)**2) + \
            alpha**2/8*((np.roll(initial,1)+initial)*(np.roll(initial,1)**2 - \
            initial**2) - (initial+np.roll(initial,-1))*(initial**2 - \
            np.roll(initial,-1)**2))

def lax_wendroff_hopfeq(initial, ival, c=1.0):
    alpha = c*ival.ALPHA_HOPF_EQ
    t = ival.N_T
    out_t = initial
    for ct in range(t):
        out_t = hopfeq_scheme_once(out_t,alpha=alpha)
    print('INFO: Just for one type of BoundaryConditions')
    return out_t

# table to map strings to functions for the simulations.
hopfTable = {
    'up' : upwind,
    'down' : downwind,
    'ex_c' : explicit_centered,
    'im_c' : implicit_centered,
    'lax_f' : lax_friedrichs,
    'lax_w' : lax_wendroff,
    'lax_w_h' : lax_wendroff_hopfeq
    }

if False:
    # dummy if, only used to not execute this part which was used to generate
    # the plots for the hopf equation
    a = 0.1
    dt = 0.001
    maxslope_anti = np.sqrt(np.exp(1))*a
    my_T = (int(maxslope_anti/dt)-1)*dt

    ival = initialValues([None,1.0,301,int(maxslope_anti/dt)-1,dt,'periodic'])
    up = simulation(ival,hopfTable,'lax_w_h')
    u_i = +np.exp(-(ival.xvec-0.5)**2/2/a**2)+1
    up.simulate(u_i,1.0)
    up.plot_init('$t=0$')
    up.plot_final('$t=%.3f$' % my_T,'Breakdown: $a = %.3f$' % a)
    #up.fig.savefig('hopfBreakdown2.pdf')
    plt.show()

# some human-readable output definitely is nice.
schemeNames = {
    'up' : 'upwind scheme',
    'down' : 'downwind scheme',
    'ex_c' : 'explicit centered scheme',
    'im_c' : 'implicit centered scheme',
    'lax_f' : '{\scshape Lax Friedrichs} scheme',
    'lax_w' : '{\scshape Lax Wendroff} scheme'
    }
XNUM = 101
from mpl_toolkits.mplot3d import axes3d

'''
PLOTTING CODE FOR THE ADVECTION EQUATION
'''
# implement the exact solution to the advection equation
def advection_exact_periodic(x,t,c=1.0):
    return np.exp(-(np.mod((x-c*t),XNUM/(XNUM-1))-0.8)**2/2/0.05**2)

def error_plot(dtvec,dxvec,which,tmax):
    '''
    helper function to create the 3d-plot that is made for the error analysis.
    '''
    dt, dx = np.meshgrid(dtvec,dxvec)
    res = np.zeros_like(dt)
    for i,DT in enumerate(dtvec):
        for j,DX in enumerate(dxvec):
            ival = initialValues([None,1.0,int(1/DX),int(tmax/DT),DT,'periodic'])
            obj = simulation(ival,hopfTable,which)
            obj.simulate(10*ival.HopfFunc,1.0)
            res[i,j] = np.sqrt(np.sum((
                10*advection_exact_periodic(ival.xvec,ival.N_T*ival.DELTA_T,1.0)- \
                    obj.result)**2))
    return dt,dx,res

prefix = 'advectionPlots/'
for scheme in ['down','up','down','ex_c','im_c','lax_f','lax_w']:
    for dt in [1/200,1/104,1/96]:
        CFL = dt*(XNUM-1)
        fig = plt.figure(figsize=(8,6))
        # cfl = const, t not const. plot.
        for time in [0.1,0.4,0.9]:
            tnum = int(time/dt)+1
            real_time = tnum*dt
            ival = initialValues([None,1.0,XNUM,tnum,dt,'periodic'])
            sim = simulation(ival,hopfTable,scheme)
            sim.simulate(ival.HopfFunc,1.0)
            plt.plot(ival.xvec,sim.result,'o-',fillstyle='none',
                    label='simulation: $t=%.3f$' % real_time)
            plt.plot(
                ival.xvec,
                advection_exact_periodic(ival.xvec,real_time,c=1.0),
                '-',label='exact: $t=%.3f$' % real_time)
        plt.grid()
        plt.xlabel('$x$')
        plt.ylabel('$u(x,t)/\mathrm{[a.u.]}$')
        plt.title('%s: $\mathsf{CFL} = %.3f$' % (schemeNames[scheme],CFL))
        plt.tight_layout()
        plt.legend()
        fig.savefig('%s%s_CFL%.3f.pdf' % (prefix,scheme,CFL))
    # truncation error plot
    fig_3d_err = plt.figure(figsize=(8,6))
    ax = fig_3d_err.add_subplot(111,projection='3d')
    ax.plot_surface(
            *error_plot(
                np.logspace(-3,-2,30),np.logspace(-1.9,-1.5,30),
                scheme,3e-2),
            cmap=plt.cm.viridis,alpha=0.8
            )
    ax.set_xlabel('$\\Delta t$')
    ax.set_ylabel('$\\Delta x$')
    ax.set_zlabel('$\\mathsf{TE}$')
    ax.set_title('Truncation Error ($\\mathsf{TE}$) for the %s' % schemeNames[scheme])
    plt.tight_layout()
    fig_3d_err.savefig(prefix+'%s_ERRORS.pdf' % scheme)
