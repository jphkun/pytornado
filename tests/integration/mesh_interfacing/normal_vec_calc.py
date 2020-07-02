#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:15:57 2020

@author: cfse2
"""


import numpy as np

v = 1/np.sqrt(3)
a = np.array([[v,v,v],[0,0,1],[0,0,1],[0,0,1]])
b = np.array([[v,v,v],[0,0,1],[0,0,1],[0,0,1]])

dot = np.dot(a,b.T)
tensordot = np.tensordot(a, b, axes=((0),(0)))
einsum = np.einsum("ij,ij->i",a,b)
print(dot)
print(tensordot)
print(einsum)