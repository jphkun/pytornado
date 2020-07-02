#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------
# Copyright 2017-2020 Airinnova AB and the PyTornado authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ----------------------------------------------------------------------

# Authors:
# * Aaron Dettmann
# * Jean-Philippe Kuntzer

"""
Reading the aircraft deformation file.

Developed at Airinnova AB, Stockholm, Sweden.
"""

import os
import logging

from commonlibs.logger import truncate_filepath
# from aeroframe.fileio.serialise import load_json_def_fields

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
# from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R
# from numpy.core.umath_tests import inner1d
# import pickle

logger = logging.getLogger(__name__)


def load(aircraft, settings):
    """
    Loads the aircraft deformation file if it exitsts.

    Args:
        :aircraft: (obj) aircraft
        :settings: (obj) settings
    """

    filepath = settings.paths('f_deformation')
    logger.info(f"Reading deformation from file '{truncate_filepath(filepath)}'")

    if not os.path.exists(filepath):
        raise IOError(f"file '{filepath}' not found")
    # File is empty or as good as (this also catches empty JSON file: '{}')
    elif os.stat(filepath).st_size < 10:
        logger.warning(f"Empty deformation file. No deformations are modelled.")
        return

    # def_fields = load_json_def_fields(filepath)
    # for wing_uid, def_field in def_fields.items():
    #     ####
    #     # Convention: Deformation fields starting with '_' will be ignore by CFD
    #     if wing_uid.startswith('_'):
    #         continue
    #     ####

    #     # Convention: Deformation fields ending with '_m' belongs to a 'mirrored' component
    #     if wing_uid.endswith('_m'):
    #         aircraft.wings[wing_uid.replace('_m', '')].def_field_mirror = def_field
    #         continue

    #     aircraft.wings[wing_uid].def_field = def_field

    # TODO: Handle exceptions
    # TODO: Check deformation continuity



def plotting_p(lattice,c):
    """
    Plots the points of the mesh
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lattice[:,0],lattice[:,1],lattice[:,2],color = c)
    var = 10
    ax.set_xlim((-var,var))
    ax.set_ylim((-var,var))
    ax.set_zlim((-var,var))
    ax.set_title("Mesh points")

def plotting_v(lattice,c):
    """
    Plots the points of the mesh
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lattice[:,0],lattice[:,1],lattice[:,2],color = c)
    ax.set_xlim((-5,5))
    ax.set_ylim((-5,5))
    ax.set_zlim((-5,5))
    ax.set_title("HSV points")

def plotting_c(lattice,c):
    """
    Plots the points of the mesh
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lattice[:,0],lattice[:,1],lattice[:,2],color = c)
    ax.set_xlim((-5,5))
    ax.set_ylim((-5,5))
    ax.set_zlim((-5,5))
    ax.set_title("cell centroids")

def plotting_n(lattice_c,lattice_n,c):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(lattice_c[:,0],lattice_c[:,1],lattice_c[:,2],
              lattice_n[:,0],lattice_n[:,1],lattice_n[:,2],color = c)
    ax.set_xlim((-5,5))
    ax.set_ylim((-5,5))
    ax.set_zlim((-5,5))
    ax.set_title("Speed normal vectors")

def plotting_r(lattice_c,rot,c):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(lattice_c[:,0],lattice_c[:,1],lattice_c[:,2],
              rot[:,0],rot[:,1],rot[:,2],color = c)
    ax.set_xlim((-5,5))
    ax.set_ylim((-5,5))
    ax.set_zlim((-5,5))
    ax.set_title("Speed normal vectors")

def shape_2(y):
    """
    Shapte function2. This function computes the slope at each y location for
    a cantilever beam of length "L", with a distributed load of "q". The Young
    modulus is imposed for steel and "I" the second moment of inertia is
    also imposed.
    """
    
    # [N/m] Distributed load
    q = 500
    # [m] Wing span
    L = 5
    logger.debug("L = " + str(L) )
    # [Pa] Elasticity modulus
    E = 210e9
    # [m**4] Second moment of area
    I = 1.330e-6
    # Computes beam deformation
    uz = q*y**2 *(6*L**2 -4*L*y +y**2)/(24*E*I)
    
    return uz

def shape_1(y,settings):
    """
    Shape function 1. This functions respresents line of slope 1. The only 
    purpouse of this shape function is to test if the deformation is done 
    correctly.
    """
    logger.info(settings.settings["aircraft"])
    dea = 0
    if settings.settings["aircraft"] == "1_flat_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0   
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "2_dih_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "3_anh_funcActivated.json":
        if dea == 0: m = 0
        else: m = -0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "4_flat_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "5_dih_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "6_anh_funcActivated.json":
        if dea == 0: m = 0
        else: m = -0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "7_flat_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "8_dih_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "9_anh_funcActivated.json":
        if dea == 0: m = 0
        else: m = -0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "10_flat_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "11_dih_funcActivated.json":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "12_anh_funcActivated.json":
        if dea == 0: m = 0
        else: m = -0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "D150_AGILE_Hangar_funActivated.xml":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    elif settings.settings["aircraft"] == "Boxwing_AGILE_Hangar_funActivated_v3.1.xml":
        if dea == 0: m = 0
        else: m = 0.1
        logger.info("wing slope m = "+str(m))
    else:
        logger.error("Deformation input is wrong!")
        
    h = 0.0
    uz = m*y + h
    
    return uz

def shape_3(z,settings):
    h = 0.0
    m = 0.1
    uy = m*z + h
    return uy

def deformation(lattice,settings):
    """
    Questions
        - Do we need to tilt the HSV?
        - How to tilt the vectors only?
    
    This function deforms the mesh, computes the normal vector and feeds back
    the computed parameters into the lattic class variable.
    
    The stdrun.run function will then continue to compute the rest.
    
    Parameters
    ----------
    lattice : class variable
        Variable of class lattice. This variable contains all the mesh points, 
        collocations points, normal vectors, horseshoe vortex points, surface 
        area of each cell, and some other non-essential information in this 
        part.

    Returns
    -------
    None.
    """
    
    logger.info("=== Starts deformation function ===")
    # reshapes lattice._ in a more user-friendly way
    # Mesh points
    s_p = lattice.p.shape
    # 4 points of the vortex horseshoe
    s_v = lattice.v.shape
    # Collocation point: point of calculation
    s_c = lattice.c.shape
    # Normal vector of each cell
    s_n = lattice.n.shape
    
    """
    Since the program is embeded in classes, it looks like a copy of the 
    variables is necessary in order to modifiy only the variables and not the 
    class variable.
    """
    # Mesh points
    lattice_p = np.copy(lattice.p.reshape((s_p[0]*s_p[1],s_p[2])))
    # HSV
    lattice_v = np.copy(lattice.v.reshape((s_v[0]*s_v[1],s_v[2])))
    # Collocation point
    lattice_c = np.copy(lattice.c)
    # Normal vector
    lattice_n = np.copy(lattice.n)
    # Lattice bml
    lattice_blm = np.copy(lattice.bound_leg_midpoints)

    # Deformation of the wing, kept absolute to keep xz plane symmetry
    y_p = np.abs(lattice_p[:,1])
    y_v = np.abs(lattice_v[:,1])
    y_c = np.abs(lattice_c[:,1])
    y_blm = np.abs(lattice_blm[:,1])
    
    # Deformation of the wing, kept absolute to keep xz plane symmetry
    z_p = np.abs(lattice_p[:,2])
    z_v = np.abs(lattice_v[:,2])
    z_c = np.abs(lattice_c[:,2])
    z_blm = np.abs(lattice_blm[:,2])
    
    # Computes the deformed surface
    shape = 1
    if shape == 1:
        uz_p = shape_1(y_p,settings)
        uz_v = shape_1(y_v,settings)
        uz_c = shape_1(y_c,settings)
        uz_blm = shape_1(y_blm,settings)
    else:
        uz_p = shape_2(y_p)
        uz_v = shape_2(y_v)
        uz_c = shape_2(y_c)
        uz_blm = shape_2(y_blm)
    
    shape3 = True
    if shape3 == True:
        uy_p = shape_3(z_p,settings)
        uy_v = shape_3(z_v,settings)
        uy_c = shape_3(z_c,settings)
        uy_blm = shape_1(z_blm,settings)
    else:
        uy_p = 0
        uy_v = 0
        uy_c = 0
        uy_blm = 0
    # logger.debug("uz_p" + str(uz_p))
    
    # Computes the normal vector by using 3 points
    # BA = lattice.p[:,1]-lattice.p[:,0]
    # BC = lattice.p[:,2]-lattice.p[:,1]
    # n_i = np.cross(BA,BC)
    # norm_n = LA.norm(n_i,axis=1)
    # n_i = n_i/norm_n[:,np.newaxis]
    
    # TODO delete after
    # Computes the initial normal vector
    G = np.concatenate((lattice.c, lattice.c, lattice.c, lattice.c), axis=1)
    mat = lattice.p - G.reshape(s_p[0],s_p[1],s_p[2])
    u, s, vh_i = np.linalg.svd(mat)
    
    # Adds the displacement on all the mesh points
    lattice_p[:,2] = lattice_p[:,2] + uz_p
    lattice_v[:,2] = lattice_v[:,2] + uz_v
    lattice_c[:,2] = lattice_c[:,2] + uz_c
    lattice_blm[:,2] = lattice_blm[:,2] + uz_blm
    lattice_p[:,1] = lattice_p[:,1] + uy_p
    lattice_v[:,1] = lattice_v[:,1] + uy_v
    lattice_c[:,1] = lattice_c[:,1] + uy_c
    lattice_blm[:,1] = lattice_blm[:,1] + uy_blm
    # Reshapes array, this permits pytornado be able to read it again
    lattice.p = lattice_p.reshape((s_p[0],s_p[1],s_p[2]))
    lattice.v = lattice_v.reshape((s_p[0],s_v[1],s_v[2]))
    lattice.c = lattice_c
    lattice.bound_leg_midpoints = lattice_blm
    
    # Computes the normal vector by using 3 points
    # BA = lattice.p[:,1]-lattice.p[:,0]
    # BC = lattice.p[:,2]-lattice.p[:,1]
    # n_f = np.cross(BA,BC)
    # norm_n = LA.norm(n_f,axis=1)
    # n_f = n_f/norm_n[:,np.newaxis]
    
    # TODO delete after
    # Computes the initial normal vector by using SVD
    G = np.concatenate((lattice.c, lattice.c, lattice.c, lattice.c), axis=1)
    mat = lattice.p - G.reshape(s_p[0],s_p[1],s_p[2])
    u, s, vh_f = np.linalg.svd(mat)
    
    # Compares the two methods
    # logger.debug(n)
    # logger.debug(vh_f[:,2])
    # logger.debug(n-vh_f[:,2])
    
    # TODO delete after
    # Computes the roation vector. This vector will be used in the quaternion
    # part and many others
    rot_g = np.cross(vh_f[:,2,:],vh_i[:,2,:])
    rot = rot_g/np.linalg.norm(rot_g,axis=1,)[:,np.newaxis]
    rot[np.isnan(rot)] = 0.0
    
    # Computes the roation vector. This vector will be used in the quaternion
    # part and many others
    # rot_g = np.cross(n_f,n_i)
    # rot = rot_g/LA.norm(rot_g,axis=1)[:,np.newaxis]
    # logger.debug(rot)
    # TODO verify if this one is needed
    # rot[np.isnan(rot)] = 0.0
    
    # TODO get rid of the useless variables when debug is finished
    # Computes the angle between the intial and deformed normal vector
    # ab = lattice.n.dot(new_norm.T)
    ab = np.einsum('ij,ij->i',vh_f[:,2,:],vh_i[:,2,:])
    # ab = np.dot(n_i,n_f)
    # ab = np.tensordot(n_i,n_f,axes=((0),(1)))
    # dot product of vector a and b
    # ab = np.einsum("ij,ij->i",n_i,n_f)
    # logger.debug(ab)
    a = np.linalg.norm(vh_f[:,2,:], axis=1)
    b = np.linalg.norm(vh_i[:,2,:], axis=1)
    # a = np.linalg.norm(n_i, axis=1)
    # b = np.linalg.norm(n_f, axis=1)
    angle = np.empty(len(ab))
    ############################################################
    # selects the correct angle (clockwise or counterlockwise) depending on
    # which side the wing is (y>0 or y<0) and the orientation of the rotation
    # vector (rot) in the "x" direction
    for i in range(len(ab)):
        if   rot[i,0]<=0 and lattice_c[i,1]<0:
            angle[i] = np.arccos(ab[i]/(a[i]*b[i]))  
        elif rot[i,0]>0 and lattice_c[i,1]<0:
            angle[i] = (np.pi + np.arccos(ab[i]/(a[i]*b[i])))
        elif rot[i,0]<=0 and lattice_c[i,1]>=0:
            angle[i] = (np.pi + np.arccos(ab[i]/(a[i]*b[i])))
        else:
            angle[i] = np.arccos(ab[i]/(a[i]*b[i]))
    angle[np.isnan(angle)] = 0.0
    
    # initializes all the quaternion rotations
    quat = np.einsum('i,ij->ij',angle,rot)
    r = R.from_rotvec(quat)
    if lattice.n[:,2].mean() >= 0:
        quat = np.einsum('i,ij->ij',angle,rot)
        r = R.from_rotvec(quat)
        logger.debug("normal vector point UP")
    else:
        quat = np.einsum('i,ij->ij',-angle,rot)
        r = R.from_rotvec(quat)
        logger.debug("normal vector point DOWN")
        
    lattice_n = r.apply(lattice.n)
    logger.debug("  lattice Z direction:\n"+str(lattice_n[0,2]))
    
    if lattice_n[:,2].mean() >= 0:
        lattice.n = np.array((1.0,-1.0,-1.0))*lattice_n
        logger.debug("normal vector point UP")
    else:
        lattice.n = lattice_n
        logger.debug("normal vector point DOWN")
    
    # Corrects the area
    s1 = lattice.p[:,0]-lattice.p[:,1]
    s2 = lattice.p[:,2]-lattice.p[:,3]
    c1 = lattice.p[:,0]-lattice.p[:,3]
    c2 = lattice.p[:,1]-lattice.p[:,2]
    c = 0.5*(np.linalg.norm(c1,axis=1)+np.linalg.norm(c2,axis=1))
    s = 0.5*(np.linalg.norm(s1,axis=1)+np.linalg.norm(s2,axis=1))
    area = s*c
    lattice.a = area