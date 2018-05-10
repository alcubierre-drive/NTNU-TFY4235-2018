#!/usr/bin/env python3
# coding: utf8

import numpy as np
import matplotlib.pyplot as plt
from raycast import *

numpoints = 400

# demonstrate that the algorithms for finding points inside a polygon work
Poly = Polygon()
Poly.AddPoint(Point(-1,0))
Poly.AddPoint(Point(-0.5,0))
Poly.AddPoint(Point(0,1))
Poly.AddPoint(Point(0,0.5))
Poly.AddPoint(Point(1,0))
Poly.AddPoint(Point(0.5,0))
Poly.AddPoint(Point(0,-1))
Poly.AddPoint(Point(0,-0.5))

points_ray = -1+2*np.random.rand(int(numpoints/2),2)
points_ang = -1+2*np.random.rand(int(numpoints/2),2)

# use the external file for the ray-casting algorithm
def raycast(polygon, point):
    testline_left = Segment(Point(-999999999,point.y), point)
    testline_right = Segment(point, Point(-999999999,point.y))
    count_left = 0
    count_right = 0
    for e in polygon.GetEdges():
        if EdgesIntersect(testline_left, e):
            count_left += 1
        if EdgesIntersect(testline_right, e):
            count_right += 1
    if count_left % 2 == 0 and count_right % 2 == 0:
        return 0
    else:
        return 1

# just the other algorithm, see report.
def angle_check(Poly,point):
    ang = 0
    for i,x in enumerate(Poly.points):
        tmp = np.arctan2(
                    *(x.arr-point.arr)
                ) - np.arctan2(
                    *(Poly.points[i-1].arr-point.arr)
                )
        if tmp > np.pi:
            ang += tmp-2*np.pi
        elif tmp < -np.pi:
            ang += tmp+2*np.pi
        else:
            ang += tmp
    if np.abs(ang) < np.pi:
        return 0
    else:
        return 1

PLOT=False

# logic to execute this
if PLOT:
    plt.fill(*Poly.arraypoints().T,fill=False,edgecolor='green')
    for p in points_ang:
        if angle_check(Poly,Point(*p)):
            plt.plot(*p,'r.')
        else:
            plt.plot(*p,'b.')
    for p in points_ray:
        if raycast(Poly,Point(*p)):
            plt.plot(*p,'y.')
        else:
            plt.plot(*p,'m.')
    plt.show()
else:
    class savefile:
        bdy = "../Data/pgf_check_inside_bdy.dat"
        ray_in = "../Data/pgf_check_inside_pts_ray_in.dat"
        ray_out = "../Data/pgf_check_inside_pts_ray_out.dat"
        wind_in = "../Data/pgf_check_inside_pts_wind_in.dat"
        wind_out = "../Data/pgf_check_inside_pts_wind_out.dat"
    np.savetxt(savefile.bdy,np.array([*Poly.arraypoints(),Poly.arraypoints()[0]]))
    ray_in = []
    ray_out = []
    wind_in = []
    wind_out = []
    for p in points_ang:
        if angle_check(Poly,Point(*p)):
            wind_in.append(p)
        else:
            wind_out.append(p)
    for p in points_ray:
        if raycast(Poly,Point(*p)):
            ray_in.append(p)
        else:
            ray_out.append(p)
    np.savetxt(savefile.ray_in,np.asarray(ray_in))
    np.savetxt(savefile.ray_out,np.asarray(ray_out))
    np.savetxt(savefile.wind_in,np.asarray(wind_in))
    np.savetxt(savefile.wind_out,np.asarray(wind_out))
