#!/usr/bin/env python3
# coding: utf8

import lib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import path
import matplotlib.mlab as mpl_ml
import numpy as np
import ctypes
import sys
import scipy.linalg as scl
import scipy.sparse as scs
import scipy.sparse.linalg as scs_la
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MaxNLocator
from datetime import datetime

clib = ctypes.cdll.LoadLibrary('Problem1.dir/calcs.so')

LMAX = 2
l = 1
kmax = 4

DELTA = 4**(-LMAX)
EPSILON = 0.1*DELTA

MODE = 'biharmonic'
# can be "five", "nine", "nine-dir", "biharmonic"
BVAL = 0.7
SAVE = True
SAVE_BDY = False
CALCULATE = True
PLOT_MODES = False
PLOT_OMEGA = False
SAVE_FIT = False

def print_time():
    ctime = datetime.now().strftime("%d.%m.%Y - %H:%M:%S")
    print(ctime)

separation = "---------------------------------"

print_time()
print("starting calculations / program")
print(separation)

# function to get the boundary generated in the C-library
def gen_boundary(l):
    c_boundary = clib.create_boundary
    c_boundary.argtypes = [ ctypes.c_int ]
    c_boundary.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
    result = c_boundary(l)
    return np.ctypeslib.as_array(result[0], shape=(4*8**l,)), \
            np.ctypeslib.as_array(result[1], shape=(4*8**l,))

# function that reads values from the C library, is not used after all.
def limiting_values(l):
    c_lim = clib.print_limit_memory
    c_lim.argtypes = [ ctypes.c_int ]
    c_lim.restype = None
    c_lim(l)

# function to tell whether a point is inside or outside a polygon
def winding_check(points, boundary, radius=EPSILON):
    p = path.Path(np.stack(boundary).T)
    return p.contains_points(points,radius=radius)

# get the grid of the simulation
def grid():
    bdy_x, bdy_y = gen_boundary(LMAX)
    xgrid = np.arange(min(bdy_x),max(bdy_x)+4**(-LMAX),4**(-LMAX))
    return np.meshgrid(xgrid,xgrid)

# now some global vaiables. Get the grid, the grid's size, a flattened version
# of these arrays, the boundary and a matrix that tells whether a point lies in-
# or outside the fractal (IN_MATRIX)
GRID_X, GRID_Y = grid()
GRIDSIZE = len(GRID_X)
X_FLAT = GRID_X.flatten()
Y_FLAT = GRID_Y.flatten()
BOUNDARY = gen_boundary(l)
IN_MATRIX = winding_check(np.asarray([X_FLAT,Y_FLAT]).T,BOUNDARY).reshape(GRID_X.shape)

# code to generate a boolean matrix containing the boundary.
xaddidx = 0
yaddidx = 0
if MODE == 'biharmonic':
    XCLOSE = np.isclose(GRID_X,min(BOUNDARY[0]))
    YCLOSE = np.isclose(GRID_Y,min(BOUNDARY[1]))
    class Found(Exception): pass
    try:
        for i in range(GRIDSIZE):
            for j in range(GRIDSIZE):
                if XCLOSE[i,j] and YCLOSE[i,j]:
                    xaddidx = i
                    yaddidx = j
                    raise Found
    except Found:
        pass

# for each point on the boundary, get the corresponding index of the matrix.
def boundary_points_to_idx(boundary):
    x_idx = boundary[0]/DELTA
    y_idx = boundary[1]/DELTA
    x_idx -= np.min(x_idx) - xaddidx
    y_idx -= np.min(y_idx) - yaddidx
    xres = np.asarray(np.around(x_idx),np.int)
    yres = np.asarray(np.around(y_idx),np.int)
    return xres, yres

# just a small helper to sort two indices
def sort_sign(i,j):
    if i < j:
        return i,j
    else:
        return j,i

# final code to get the boundary matrix
def boundary_matrix(boundary):
    x, y = boundary_points_to_idx(boundary)
    mtx = np.zeros((GRIDSIZE,GRIDSIZE),dtype=np.bool)
    for i in range(len(x)):
        mtx[y[i],x[i]] = True
        idx_delta_x = x[i]-x[i-1]
        idx_delta_y = y[i]-y[i-1]
        if not idx_delta_x == 0:
            for newx in range(*sort_sign(x[i-1],x[i])):
                mtx[y[i],newx] = True
        elif not idx_delta_y == 0:
            for newy in range(*sort_sign(y[i-1],y[i])):
                mtx[newy,x[i]] = True
    return mtx

# this is the boolean matrix that contains the boundary
BOUNDARY_MATRIX = boundary_matrix(BOUNDARY)
if MODE == 'biharmonic':
    # in this case we need to know how many next-to-nearest neighbors lie
    # outside the boundary. They are stored in the matrix NOT_INSIDE
    NOT_INSIDE = np.zeros(IN_MATRIX.shape)
    for i in range(GRIDSIZE):
        for j in range(GRIDSIZE):
            if IN_MATRIX[i,j]:
                NOT_INSIDE[i,j] = BOUNDARY_MATRIX.astype(int)[i,j+2]+ \
                        BOUNDARY_MATRIX.astype(int)[i+2,j]+ \
                        BOUNDARY_MATRIX.astype(int)[i-2,j]+ \
                        BOUNDARY_MATRIX.astype(int)[i,j-2]

# just some more definitions, such as how many points lie iniside the boundary,
# and their coordinates
LENGTH_INSIDE = np.sum(1*IN_MATRIX.flatten())
INSIDE_COORDS_X = X_FLAT[IN_MATRIX.flatten()]
INSIDE_COORDS_Y = Y_FLAT[IN_MATRIX.flatten()]
PLOTLIMS = (-0.4,1.4)

# important function mapping the rectangular indices to the ones that are used
# in the computation
def map_grid_to_final_vec(i,j):
    final_idx = np.sum(1*IN_MATRIX.flatten()[:np.ravel_multi_index((i,j),(GRIDSIZE,GRIDSIZE))])
    return final_idx-1

# change the matrix elements of the finite difference matrix
def index_gymnastics(i,j,operator,di,dj,c=1.):
    if IN_MATRIX[i+di,j+dj]:
        nidx_old = map_grid_to_final_vec(i,j)
        nidx_new = map_grid_to_final_vec(i+di,j+dj)
        operator[nidx_new,nidx_old] = c
    return operator

# generate the laplacian / biharmonic finite difference matrix, different
# implementations, explained in the report.
def get_laplace(order='five', b=0.5):
    my_operator = scs.dok_matrix((LENGTH_INSIDE,LENGTH_INSIDE))
    if order == 'five':
        for i in range(GRIDSIZE):
            for j in range(GRIDSIZE):
                flatindex = np.ravel_multi_index((i,j),(GRIDSIZE,GRIDSIZE))
                if IN_MATRIX[i,j]:
                    nidx = map_grid_to_final_vec(i,j)
                    my_operator[nidx,nidx] = -4
                    index_gymnastics(i,j,my_operator,1,0)
                    index_gymnastics(i,j,my_operator,-1,0)
                    index_gymnastics(i,j,my_operator,0,1)
                    index_gymnastics(i,j,my_operator,0,-1)
    elif order == 'nine':
        for i in range(GRIDSIZE):
            for j in range(GRIDSIZE):
                flatindex = np.ravel_multi_index((i,j),(GRIDSIZE,GRIDSIZE))
                if IN_MATRIX[i,j]:
                    nidx = map_grid_to_final_vec(i,j)
                    my_operator[nidx,nidx] = -4*b+(-4/2)*(1-b)
                    index_gymnastics(i,j,my_operator,1,0,b)
                    index_gymnastics(i,j,my_operator,-1,0,b)
                    index_gymnastics(i,j,my_operator,0,1,b)
                    index_gymnastics(i,j,my_operator,0,-1,b)
                    index_gymnastics(i,j,my_operator,1,1,(1-b)/2)
                    index_gymnastics(i,j,my_operator,-1,-1,(1-b)/2)
                    index_gymnastics(i,j,my_operator,-1,1,(1-b)/2)
                    index_gymnastics(i,j,my_operator,1,-1,(1-b)/2)
    elif order == 'nine_dir':
        for i in range(GRIDSIZE):
            for j in range(GRIDSIZE):
                flatindex = np.ravel_multi_index((i,j),(GRIDSIZE,GRIDSIZE))
                if IN_MATRIX[i,j]:
                    nidx = map_grid_to_final_vec(i,j)
                    my_operator[nidx,nidx] = -60/12
                    index_gymnastics(i,j,my_operator,1,0,16/12)
                    index_gymnastics(i,j,my_operator,0,1,16/12)
                    index_gymnastics(i,j,my_operator,0,-1,16/12)
                    index_gymnastics(i,j,my_operator,-1,0,16/12)
                    index_gymnastics(i,j,my_operator,2,0,-1/12)
                    index_gymnastics(i,j,my_operator,0,2,-1/12)
                    index_gymnastics(i,j,my_operator,0,-2,-1/12)
                    index_gymnastics(i,j,my_operator,-2,0,-1/12)
    elif order == 'biharmonic':
        for i in range(GRIDSIZE):
            for j in range(GRIDSIZE):
                flatindex = np.ravel_multi_index((i,j),(GRIDSIZE,GRIDSIZE))
                if IN_MATRIX[i,j]:
                    nidx = map_grid_to_final_vec(i,j)
                    my_operator[nidx,nidx] = 20+NOT_INSIDE[i,j]
                    index_gymnastics(i,j,my_operator,0,1,-8)
                    index_gymnastics(i,j,my_operator,1,0,-8)
                    index_gymnastics(i,j,my_operator,0,-1,-8)
                    index_gymnastics(i,j,my_operator,-1,0,-8)
                    index_gymnastics(i,j,my_operator,0,2,1)
                    index_gymnastics(i,j,my_operator,2,0,1)
                    index_gymnastics(i,j,my_operator,0,-2,1)
                    index_gymnastics(i,j,my_operator,-2,0,1)
                    index_gymnastics(i,j,my_operator,1,1,2)
                    index_gymnastics(i,j,my_operator,1,-1,2)
                    index_gymnastics(i,j,my_operator,-1,1,2)
                    index_gymnastics(i,j,my_operator,-1,-1,2)
    else:
        print("No such option ´%s´ for kwarg ´order´" % order,file=sys.stderr)
        raise ValueError
    if order == 'biharmonic':
        return my_operator/DELTA**4
    else:
        return -my_operator/DELTA**2

# get the solution of the eigenvalue problem.
def solution(operator, kmax):
    v0 = np.zeros(LENGTH_INSIDE)
    v0[int(LENGTH_INSIDE/2)] = 1
    sol = scs_la.eigsh(scs.csc_matrix(operator), k=kmax, which='SM', v0=v0)
    Ufunc = sol[1]
    if MODE == 'biharmonic':
        omega = np.abs(sol[0])**(1/4)
    else:
        omega = np.abs(sol[0])**(1/2)
    return omega, Ufunc

# get an interpolation of the function U, only used for plotting
def modify_arrays(x,y,z,length=300):
    xv = np.linspace(*PLOTLIMS,length)
    X,Y = np.meshgrid(xv,xv)
    Z = mpl_ml.griddata(x,y,z,X,Y,interp='linear').filled(0)
    return X,Y,Z

# save the arrays used for plotting to a file
def save_plot_arrays(x,y,z,savename):
    X,Y,Z = modify_arrays(x,y,z,length=200)
    XYZ = np.stack([X,Y,Z])
    np.savetxt(savename,XYZ[:,:,0].T)
    f = open(savename,"ab")
    for i in range(1,200):
        f.write(b"\n")
        np.savetxt(f,XYZ[:,:,i].T)
    f.close()

# save the points of the boundary to a file
def save_plot_boundary(l,fname):
    bnd = np.stack(gen_boundary(l))
    np.savetxt(fname,bnd.T)
    f = open(fname,"ab")
    np.savetxt(f,[bnd.T[0]])
    f.close()

# generate a plot of a mode
def plot_U_func(X,Y,Z,title):
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111)
    cf = ax.contourf(X,Y,Z,30,cmap=plt.cm.viridis)
    fig.colorbar(cf,ax=ax,label='$U/\\mathrm{[a.u.]}$')
    plt.fill(*BOUNDARY, edgecolor='red', fill=False)
    plt.xlabel('$x/L$')
    plt.ylabel('$y/L$')
    plt.title(title)
    plt.tight_layout()
    return fig

# some logic to use with the parameters set in the beginning of this file ...
if CALCULATE:
    laplace = get_laplace(order=MODE,b=BVAL)
    print_time()
    if MODE == 'biharmonic':
        print("gotten the $\\nabla^4$")
    else:
        print("gotten the laplacian")
    print(separation)
    omega, U = solution(laplace,kmax)
    print_time()
    print("gotten the solution up to mode %i" % kmax)
    print(separation)
    if SAVE:
        np.savetxt('Data/omega_%s_k%i_lm%i_l%i.txt' % (MODE,kmax,LMAX,l),omega)
        np.savetxt('Data/Ufunc_%s_k%i_lm%i_l%i.txt' % (MODE,kmax,LMAX,l),U)
else:
    print_time()
    print("plots...")
    print(separation)
    omega = np.genfromtxt('Data/omega_%s_k%i_lm%i_l%i.txt' % (MODE,kmax,LMAX,l))
    U = np.genfromtxt('Data/Ufunc_%s_k%i_lm%i_l%i.txt' % (MODE,kmax,LMAX,l))

# ... continues ...
if SAVE:
    for j in range(10):
        save_plot_arrays(INSIDE_COORDS_X,INSIDE_COORDS_Y,U[:,j],'Data/pgf_%s_lm%i_l%i_mode%i.dat'
            % (MODE,LMAX,l,j) )
    if SAVE_BDY:
        save_plot_boundary(1,'Data/pgf_boundary_l1.dat')
        save_plot_boundary(2,'Data/pgf_boundary_l2.dat')
        save_plot_boundary(3,'Data/pgf_boundary_l3.dat')
        save_plot_boundary(4,'Data/pgf_boundary_l4.dat')
        save_plot_boundary(5,'Data/pgf_boundary_l5.dat')
        save_plot_boundary(6,'Data/pgf_boundary_l6.dat')
        save_plot_boundary(7,'Data/pgf_boundary_l7.dat')

# ... and on ...
if PLOT_OMEGA or SAVE_FIT:
    delta_N = -np.asarray(range(1,len(omega)+1))+omega**2/4/np.pi
    print_time()
    print("Fitting to $\\Delta N(\\omega)$:")
    fit = lib.fit_y(omega,delta_N,0.1,lambda beta,x: beta[0]*x**beta[1],[1,1])
    omega_help = lib.help_arr(omega,size=100)
    delta_N_help = fit[0][0][0]*omega_help**fit[0][0][1]
    print('d = %.3f' % fit[0][0][1])
    print(separation)
    fig = plt.figure(figsize=(5,4))
    plt.plot(omega,delta_N,label='numerical data')
    plt.plot(omega_help,delta_N_help,label='fit')
    plt.loglog()
    plt.grid(which='both')
    plt.xlabel('$\\omega$')
    plt.ylabel('$\\Delta N(\\omega)$')
    plt.title('$d = %.2f, l_\\mathrm{max} = %i, l = %i, k_\\mathrm{max} =%i$' %
        (fit[0][0][1],LMAX, l, kmax))
    plt.tight_layout()
    if SAVE_FIT:
        print('What?')
        fig.savefig('Data/matplotlib_fit_%s_k%i_lm%i_l%i.pdf' % 
            (MODE,kmax,LMAX,l))

# ... and on ...
if PLOT_MODES:
    for i in range(min([10,kmax])):
        plot_U_func(
            *modify_arrays(INSIDE_COORDS_X,INSIDE_COORDS_Y,U[:,i]),
            "$\\omega = %.2f$, mode $%i$" % (omega[i],i)
        )

# ... and that's the last one.
if PLOT_MODES or PLOT_OMEGA:
    plt.show()
