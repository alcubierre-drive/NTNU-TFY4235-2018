#!/usr/bin/env python3
# coding: utf8

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sc_la
from helpers import initialValues, delta, simulation, make_bc
import sys
import warnings
import scipy.special as scs

'''
the following functions just implement the different linear schemes, i.e.
euler_explicit, euler_implicit and crank_nicolson. The solution is *not*
obtained by a linear solver but rather by matrix inversion because after
testing, the linear solver had a smaller range of initial conditions to work
with (due to over/underflow).
'''
def euler_explicit(initial, ival, D=1.0):
    alpha = ival.ALPHA_DIFFUSION*D
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)*(1.0-2.0*alpha)) + \
            np.diag(np.ones(ival.N_X-1)*alpha,k=1) + \
            np.diag(np.ones(ival.N_X-1)*alpha,k=-1),
            ival.BoundaryType)
    return np.linalg.matrix_power(matrix,ival.N_T) @ initial

def euler_implicit(initial, ival, D=1.0):
    alpha = ival.ALPHA_DIFFUSION*D
    matrix = make_bc(
            np.diag(np.ones(ival.N_X)*(1.0+2.0*alpha)) - \
            np.diag(np.ones(ival.N_X-1)*alpha,k=1) - \
            np.diag(np.ones(ival.N_X-1)*alpha,k=-1),
            ival.BoundaryType)
    result = np.linalg.matrix_power(np.linalg.inv(matrix),ival.N_T) @ initial
    return result

#CHECK_CRANK_COMMUTE = True
# the commented part in this function was used to check if the two matrices in
# the crank_nicolson scheme commute, they do.
def crank_nicolson(initial, ival, D=1.0):
    alpha = ival.ALPHA_DIFFUSION*D
    Bmat = make_bc(
            np.diag(np.ones(ival.N_X))*2 - \
            np.diag(np.ones(ival.N_X-1),k=1) - \
            np.diag(np.ones(ival.N_X-1),k=-1),
            ival.BoundaryType)
    left_mat = (2*np.diag(np.ones(ival.N_X))+alpha*Bmat)
    right_mat = (2*np.diag(np.ones(ival.N_X))-alpha*Bmat)
    #if CHECK_CRANK_COMMUTE:
    #    closemat = np.isclose(
    #            np.linalg.inv(left_mat) @ right_mat - \
    #                    right_mat @ np.linalg.inv(left_mat),
    #            np.zeros_like(right_mat))
    #    if not closemat.sum() == closemat.shape[0]**2:
    #        # check if matrices commute to only solve the problem once, not
    #        # recursively
    #        warnings.warn('Matrices do not commute. Wrong result.',RuntimeWarning)
    #result = sc_la.solve(np.linalg.matrix_power(left_mat,ival.N_T),
    #        np.dot(np.linalg.matrix_power(right_mat,ival.N_T),initial),
    #        assume_a="sym")
    result = np.linalg.matrix_power(np.linalg.inv(left_mat) @
            right_mat,ival.N_T) @ initial
    return result

'''
This function implements the crank_nicolson scheme for a variable diffusivity
profile (which is stored in the object ival as ival.DiffStep), of course only
the points of the x values are relevant, thus only an array for ival.DiffStep is
needed.
'''
def crank_nicolson_Dfun(initial, ival, D=1.0):
    if ival.BoundaryType != 'periodic':
        warnings.warn('Dfun Crank only supported with PBC.',RuntimeWarning)
    alpha = 0.5*ival.ALPHA_DIFFUSION
    stepfun = ival.DiffStep*D
    upper = 0.5*(3*stepfun-np.roll(stepfun,-1))
    lower = 0.5*(3*stepfun-np.roll(stepfun,1))
    gamma = np.diag(0.5*(np.roll(stepfun,-1)+np.roll(stepfun,1))-3*stepfun) + \
            np.diag(upper[:-1],k=1) + \
            np.diag(lower[1:],k=-1)
    gamma[0][ival.N_X-1] = lower[0]
    gamma[ival.N_X-1][0] = upper[-1]
    left_mat = np.diag(np.ones(ival.N_X))-gamma*alpha
    right_mat = np.diag(np.ones(ival.N_X))+gamma*alpha
    #if CHECK_CRANK_COMMUTE:
    #    closemat = np.isclose(
    #            np.linalg.inv(left_mat) @ right_mat - \
    #                    right_mat @ np.linalg.inv(left_mat),
    #            np.zeros_like(right_mat))
    #    if not closemat.sum() == closemat.shape[0]**2:
    #        # check if matrices commute to only solve the problem once, not
    #        # recursively
    #        warnings.warn('Matrices do not commute. Wrong result.',RuntimeWarning)
    #result = sc_la.solve(np.linalg.matrix_power(left_mat,ival.N_T),
    #        np.dot(np.linalg.matrix_power(right_mat,ival.N_T),initial))
    result = np.linalg.matrix_power(np.linalg.inv(left_mat) @
            right_mat,ival.N_T) @ initial
    return result

'''
The following functions are the exact solutions to the diffusion equation for
different boundary conditions.
'''
def exact_periodic(x,t,D):
    '''
    This is actually not the exact solution. But since there was none given for
    periodical boundaries, this one was used and normalized accordingly. Noted
    similarly in the solutions to the tasks.
    '''
    result = (4.*np.pi*D*t)**(-0.5)*np.exp(-(x-0.5)**2/(4.*D*t))
    normalization = scs.erf(0.5/np.sqrt(4*D*t))
    return result/normalization

def exact_absorbing(x,t,D,NMAX=200,L=1.0):
    result = np.zeros_like(x)
    for n in range(1,NMAX+1):
        result += np.exp(-(np.pi*n/L)**2*D*t)*np.sin(n*np.pi/L*x) * \
                np.sin(n*np.pi*0.5)
    return 2*result/L

def exact_reflecting(x,t,D,NMAX=200,L=1.0):
    result = np.ones_like(x)/L
    for n in range(1,NMAX+1):
        result += 2/L*np.exp(-(np.pi*n/L)**2*D*t)*np.cos(n*np.pi/L*x)* \
                np.cos(n*np.pi*0.5)
    return result

'''
The exact solution for ival.DiffStep just a step function. Taken from the
exercise sheet.
'''
def exact_dfun(x,t,Dfun,x0=0.5):
    fac_left = 2*np.sqrt(Dfun[0]/Dfun[-1])/(1+np.sqrt(Dfun[0]/Dfun[-1]))/ \
            np.sqrt(4*np.pi*Dfun[0]*t)*np.exp(-(x-x0)**2/4/Dfun[0]/t)
    fac_right = 2/(1+np.sqrt(Dfun[0]/Dfun[-1]))/np.sqrt(4*np.pi*Dfun[-1]*t)* \
            np.exp(-(x-x0)**2/4/Dfun[-1]/t)
    result_left = (x<=x0)*fac_left
    result_right = (x>x0)*fac_right
    return result_left+result_right

# table to get the correct exact function for a short key.
diffusionExact = {
        'periodic' : exact_periodic,
        'absorbing' : exact_absorbing,
        'reflecting' : exact_reflecting,
        'Dfun': exact_dfun
        }

# table to get the correct simulation function for a short key.
diffusionTable = {
        'eu_e' : euler_explicit,
        'eu_i' : euler_implicit,
        'cr_n' : crank_nicolson,
        'cr_n_D' : crank_nicolson_Dfun
        }

# table to translate the short keys into human-readable names
translationTitle = {
        'eu_e' : 'explicit \\textsc{Euler}-scheme',
        'eu_i' : 'implicit \\textsc{Euler}-scheme',
        'cr_n' : '\\textsc{Crank-Nicolson}-scheme'
        }

def get_ival(n_t,delta_t,n_x,bdy):
    '''
    small helper function for the initial values used. merely useful.
    '''
    obj = initialValues([None,1.0,51,n_t,delta_t,bdy])
    return obj

def simulate(which,boundary,delta_t):
    '''
    wrapper function for one plot. plots the different times in the so-called
    array for the scheme ´which´ and boundary ´boundary´. The time step has the
    size ´delta_t´.
    '''
    times = np.array([0.01,0.04,0.07])
    how_many_times = np.asarray(times/delta_t,dtype=int)
    xnum = 51
    CFL = delta_t*(xnum+1)**2

    fig = plt.figure(figsize=(8,5))

    for i in range(3):
        ival = get_ival(how_many_times[i],delta_t,xnum,boundary)
        obj = simulation(ival,diffusionTable,which)
        obj.simulate(ival.DiffFunc,1.0)
        plt.plot(ival.xvec,obj.result,'-o',fillstyle='none',label='$t=%.2f$, numerical' % times[i])
        plt.plot(ival.xvec,diffusionExact[boundary](ival.xvec,
            ival.N_T*ival.DELTA_T,1.0),'-',label='$t=%.2f$, exact' % times[i])
    plt.xlabel('$x$')
    plt.ylabel('$u / [\\mathrm{a.u.}]$')
    plt.title('%s, %s boundary, $\\mathsf{CFL} = %.3f$' % \
            (translationTitle[which],boundary,CFL))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    return fig

prefix = 'diffusionPlots/'
if False:
    # save all kinds of different simulations
    for scheme in ['eu_e','eu_i','cr_n']:
        for bdy in ['reflecting','absorbing','periodic']:
            for dt in [1e-4,1e-3,1.818e-4,1.852e-4]:
                fname = prefix + '%s_%s_CFL%.3f.pdf' % (scheme,bdy[:3],52**2*dt)
                simulate(scheme,bdy,dt).savefig(fname)

def error_plot(dtvec,dxvec,which,tmax):
    '''
    helper function to create the 3d-plot that is made for the error analysis.
    Logarithmic scaling of all axes.
    '''
    dt, dx = np.meshgrid(dtvec,dxvec)
    res = np.zeros_like(dt)
    for i,DT in enumerate(dtvec):
        for j,DX in enumerate(dxvec):
            ival = initialValues([None,1.0,int(1/DX),int(tmax/DT),DT,'reflecting'])
            obj = simulation(ival,diffusionTable,which)
            obj.simulate(ival.DiffFunc,1.0)
            res[i,j] = np.sqrt(np.sum((
                diffusionExact['reflecting'](ival.xvec,ival.N_T*ival.DELTA_T,1.0)- \
                    obj.result)**2))
    return np.log10(dt),np.log10(dx),np.log(res)

if False:
    # save the 3d plots for error analysis
    from mpl_toolkits.mplot3d import axes3d
    fig_cr_n = plt.figure(figsize=(8,6))
    ax = fig_cr_n.add_subplot(111,projection='3d')
    ax.plot_surface(
            *error_plot(np.logspace(-4,-3,30),np.logspace(-2,-1.5,30),'cr_n',1e-3),
            cmap=plt.cm.viridis,alpha=0.9
            )
    ax.set_xlabel('$\\mathrm{log_{10}}(\\Delta t)$')
    ax.set_ylabel('$\\mathrm{log_{10}}(\\Delta x)$')
    ax.set_zlabel('$\\mathrm{log}(\\mathsf{TE})$')
    ax.set_title('Truncation Error ($\\mathsf{TE}$) for the {\\scshape Crank-Nicolson}-scheme')
    plt.tight_layout()
    fig_cr_n.savefig(prefix+'cr_n_ERRORS.pdf')

    fig_eu_e = plt.figure(figsize=(8,6))
    ax = fig_eu_e.add_subplot(111,projection='3d')
    ax.plot_surface(
            *error_plot(np.logspace(-4,-3,30),np.logspace(-2,-1.5,30),'eu_e',1e-3),
            cmap=plt.cm.viridis,alpha=0.9
            )
    ax.set_xlabel('$\\mathrm{log_{10}}(\\Delta t)$')
    ax.set_ylabel('$\\mathrm{log_{10}}(\\Delta x)$')
    ax.set_zlabel('$\\mathrm{log}(\\mathsf{TE})$')
    ax.set_title('Truncation Error ($\\mathsf{TE}$) for the explicit {\\scshape Euler}-scheme')
    plt.tight_layout()
    fig_eu_e.savefig(prefix+'eu_e_ERRORS.pdf')

    fig_eu_i = plt.figure(figsize=(8,6))
    ax = fig_eu_i.add_subplot(111,projection='3d')
    ax.plot_surface(
            *error_plot(np.logspace(-4,-3,30),np.logspace(-2,-1.5,30),'eu_i',1e-3),
            cmap=plt.cm.viridis,alpha=0.9
            )
    ax.set_xlabel('$\\mathrm{log_{10}}(\\Delta t)$')
    ax.set_ylabel('$\\mathrm{log_{10}}(\\Delta x)$')
    ax.set_zlabel('$\\mathrm{log}(\\mathsf{TE})$')
    ax.set_title('Truncation Error ($\\mathsf{TE}$) for the implicit {\\scshape Euler}-scheme')
    plt.tight_layout()
    fig_eu_i.savefig(prefix+'eu_i_ERRORS.pdf')

if False:
    # save the CLF=1/6 plot
    dtvec_16 = np.logspace(-4,-3,30)
    dxvec_16 = np.sqrt(6*dtvec_16)
    fig_eu_e_16 = plt.figure(figsize=(8,6))
    ax = fig_eu_e_16.add_subplot(111)
    X,Y,Z = error_plot(dtvec_16,dxvec_16,'eu_e',1e-3)
    ind = np.diag_indices_from(X)
    x = X[ind]
    y = Y[ind]
    z = Z[ind]
    ax.plot(x,z)
    ax.set_xlabel('$\\mathrm{log_{10}}(\\Delta t)$')
    ax.set_ylabel('$\\mathrm{log}(\\mathsf{TE})$')
    ax.set_title('Truncation Error ($\\mathsf{TE}$) for the explicit {\\scshape Euler}-scheme; $\\mathsf{CLF} = 1/6$')
    plt.grid()
    plt.tight_layout()
    fig_eu_e_16.savefig(prefix+'eu_e_ERRORS16.pdf')

T = 1000
eps = 5e-3
DT = 1e-6
ival = initialValues([None,1.0,701,T,DT,'periodic'],eps=eps)
u_init = eps/np.pi/((ival.xvec-0.5)**2+eps**2)
# general check for the cr_n_Dfun-simulation.
if False:
    cnsd = simulation(ival,diffusionTable,'cr_n_D')
    cnsd.simulate(eps/np.pi/((ival.xvec-0.5)**2+eps**2),1.0)
    cnsd.plot_init(plot_default=False)
    plt.plot(ival.xvec,cnsd.result,label='simulated')
    plt.plot(ival.xvec,exact_dfun(ival.xvec,T*DT,ival.DiffStep,x0=0.5),label='exact')
    plt.grid()
    plt.xlabel('$x$')
    plt.ylabel('$u(x,t)/\\mathrm{[a.u.]}$')
    plt.title('Step profile, mathematical comparison')
    plt.legend()
    plt.tight_layout()
    cnsd.fig.savefig(prefix+'cr_n_D_COMPARE.pdf')

from helpers import fourier_rect,many_steps,ran_profile
'''
plot some of the rather obscure profiles. Mostly plotting code.
'''
fig_steps = plt.figure(figsize=(7,4))
ax_steps = fig_steps.add_subplot(111)
fig_results = plt.figure(figsize=(7,4))
ax_results = fig_results.add_subplot(111)
nran = [1,5,10,50]
results = [None,None,None,None]
T = 500
eps = 1e-4
DT = 1e-5
ival = initialValues([None,1.0,101,T,DT,'periodic'],eps=eps)
u_init = \
        eps/np.pi/((ival.xvec-0.0)**2+eps**2) + \
        eps/np.pi/((ival.xvec-0.2)**2+eps**2) + \
        eps/np.pi/((ival.xvec-0.4)**2+eps**2) + \
        eps/np.pi/((ival.xvec-0.6)**2+eps**2) + \
        eps/np.pi/((ival.xvec-0.8)**2+eps**2)
def make_sim(nmax):
    ival.DiffStep = fourier_rect(ival.xvec,h=0.5,N=nmax) + \
            ran_profile(ival.xvec,h=1.0,N=nmax) + 0.5
    nfive = simulation(ival,diffusionTable,'cr_n_D')
    nfive.simulate(u_init,1.0)
    return nfive.result
for i,n in enumerate(nran):
    results[i] = make_sim(n)
    ax_steps.plot(ival.xvec,ival.DiffStep,label='$n=%i$' % n)
    ax_results.plot(ival.xvec,results[i],label='$n=%i$' % n)
ax_steps.legend()
ax_steps.set_xlabel('$x$')
ax_steps.set_ylabel('$D_n(x)$')
ax_steps.grid()
ax_steps.set_title('Step profiles of different smoothness')
fig_steps.tight_layout()
ax_results.legend()
ax_results.set_xlabel('$x$')
ax_results.set_ylabel('$u_n(x,t=%.3f)$' % (T*DT))
ax_results.grid()
ax_results.set_title('Simulations for step profiles of different smoothness')
fig_results.tight_layout()

fig_steps.savefig(prefix+'step_smooth.pdf')
fig_results.savefig(prefix+'step_smooth_res.pdf')
