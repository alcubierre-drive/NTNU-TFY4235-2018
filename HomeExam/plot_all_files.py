#!/usr/bin/env python3
#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import glob
for f in glob.glob('*.dat'):
    d = np.genfromtxt(f).T
    plt.figure()
    plt.title(f.replace('_',' ').replace('-',' '))
    plt.plot(d[0],d[1])
    plt.plot(d[0],d[2])
    plt.plot(d[0],d[3]*300)
    plt.plot(d[0],d[4])
    plt.show()
