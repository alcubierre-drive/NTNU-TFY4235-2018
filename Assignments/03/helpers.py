#!/usr/bin/env python3
#coding: utf8

import lib
import numpy as np
import matplotlib.pyplot as plt

def delta(x,x0,N=1):
    '''
    should only be used on-grid. Implementation of the Dirac-delta.
    '''
    return N*np.isclose(x,x0)

class initialValues:
    '''
    Just a storage for initial values. Can be used in both simulations, hopf and
    diffusion.
    '''
    def __init__(self,argv,eps=1e-2):
        self.LENGTH = float(argv[1])
        self.N_X = int(argv[2])
        self.DELTA_X = self.LENGTH/(self.N_X-1)
        self.N_T = int(argv[3])
        self.DELTA_T = float(argv[4])
        self.ALPHA_HOPF = self.DELTA_T/self.DELTA_X
        self.ALPHA_HOPF_EQ = self.DELTA_T/self.DELTA_X
        self.ALPHA_DIFFUSION = self.DELTA_T/self.DELTA_X**2
        self.BoundaryType = argv[5]
        self.xvec = np.linspace(0,self.LENGTH,self.N_X)
        self.HopfFunc = np.exp(-(self.xvec-0.8)**2/2/0.05**2)
        self.DiffFunc = delta(self.xvec,self.xvec[int(self.N_X/2)],N=self.N_X-1)
        #self.DiffStep = 0.5*np.ones_like(self.xvec) + \
        #        0.5*np.heaviside(self.xvec-0.5,1.0)
        self.DiffStep = 0.5*np.ones_like(self.xvec) + \
                0.5/np.pi*(np.pi/2+np.arctan((self.xvec-0.5)/eps))

def fourier_rect(t,h=1.0,N=8):
    '''
    get the fourier sum until term ´N´ for a rectangular function.
    '''
    res = np.zeros_like(t)
    for n in range(1,N+1):
        res += np.sin((2*n-1)*2*np.pi*t)/(2*n-1)
    return res*4*h/np.pi

def many_steps(t,h=1.0,N=8):
    '''
    generate an approximation to a linear function with ´N´ steps.
    '''
    stepwidth = (t.max()-t.min())/N
    stepsize = h/N
    res = np.ones_like(t)*stepsize
    for n in range(1,N):
        res += np.heaviside(t-n*stepwidth,0.5)*stepsize
    return res

def ran_profile(t,h=1.0,N=2):
    '''
    return a median-filtered random profile.
    '''
    import scipy.signal as sc_sig
    filtnum = int(len(t)/N/2)*2+1
    result = sc_sig.medfilt(np.random.rand(len(t)),filtnum)
    return result*h

def make_bc(matrix,which):
    '''
    Applies the boundary conditions of choice to the matrix of the pde. Can
    either treat periodic or reflecting boundary conditions. The scheme works
    like it is presented in the slides for the assignment (and the written
    solution to the assignment).
    '''
    size_x, size_y = matrix.shape
    if which == 'periodic':
        matrix[size_x-1][0] = matrix[0][1]
        matrix[0][size_y-1] = matrix[1][0]
    elif which == 'reflecting':
        matrix[0][1] *= 2
        matrix[size_x-1][size_y-2] *= 2
    elif which == 'absorbing':
        pass
    else:
        print('ERROR: No such option for the type of the boundary condition.')
        # technically not the correct error, but too lazy for custom handlers.
        raise ValueError
    return matrix

class simulation:
    '''
    An almost useless class initially thought to gather all types of
    simulations, but in fact not really used in that way (but just as storage
    for results).
    '''
    def __init__(self, init, table, which):
        self.simfun = table[which]
        self.initialValues = init
    def simulate(self, initial, const=1.0):
        self.initial = initial
        self.const = const
        self.result = self.simfun(initial,self.initialValues,const)
    def plot_init(self,label='initial',plot_default=True):
        self.fig = plt.figure(figsize=(6,4))
        self.ax = self.fig.add_subplot(111)
        if plot_default:
            self.ax.plot(self.initialValues.xvec,self.initial,label=label)
    def plot(self,ax,label):
        ax.plot(self.initialValues.xvec,self.result,label=label)
    def plot_final(self,label,title):
        self.ax.plot(self.initialValues.xvec,self.result,label=label)
        self.ax.set_xlabel('$x$')
        self.ax.set_ylabel('$u(x,t)$')
        self.ax.grid()
        self.ax.set_title(title)
        self.ax.legend()
        self.fig.tight_layout()
