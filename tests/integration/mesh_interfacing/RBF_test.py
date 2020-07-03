#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:28:32 2020

@author: cfse2
"""


import numpy as np
from scipy import interpolate
import pylab as py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def func(x,y):
    return (x+y)*np.exp(-5.0*(x**2 + y**2))
N = 60
x = np.random.uniform(-1.0,1.0, size = N)
y = np.random.uniform(-1.0,1.0, size = N)
z = func(x,y)

newfunc = interpolate.Rbf(x,y,z, function='multiquadric')
xnew, ynew = np.mgrid[-1:1:100j, -1:1:100j]
znew = newfunc(xnew,ynew)
xexa = func(xnew,ynew)

# Generate image plot
py.figure(1)
py.clf()
py.imshow(znew,extent=[-1,1,-1,1], cmap=py.cm.jet)
py.scatter(x,y,30,z,cmap=py.cm.jet)

# Generates a surface plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(xnew,ynew,znew,cmap=cm.coolwarm)
ax.set_xlim((-1,1))
ax.set_ylim((-1,1))
ax.set_zlim((-1,1))
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_title("Airplane")
plt.show()