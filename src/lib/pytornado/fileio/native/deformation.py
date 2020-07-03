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
import pandas as pd
from numpy import linalg as LA
# from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R
from scipy import interpolate



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

# def shape_2(y):
#     """
#     Shapte function2. This function computes the slope at each y location for
#     a cantilever beam of length "L", with a distributed load of "q". The Young
#     modulus is imposed for steel and "I" the second moment of inertia is
#     also imposed.
#     """
#     logger.info("Shape function 2 is selected")
#     # [N/m] Distributed load
#     q = 500
#     # [m] Wing span
#     L = 5
#     # logger.debug("L = " + str(L) )
#     # [Pa] Elasticity modulus
#     E = 210e9
#     # [m**4] Second moment of area
#     I = 1.330e-6
#     # Computes beam deformation
#     uz = q*y**2 *(6*L**2 -4*L*y +y**2)/(24*E*I)
    
#     return uz

def shape_1(u,y,settings,lattice):
    """
    Shape function 1. This functions respresents line of slope 1. The only 
    purpouse of this shape function is to test if the deformation is done 
    correctly.
    """
    csv_save = True
    logger.info("Shape function 1 is selected")
    case = settings.settings["aircraft"]
    if case == "1_flat_funcActivated.json": m = 0   
    elif case == "2_dih_funcActivated.json": m = 0.1
    elif case == "3_anh_funcActivated.json": m = -0.1
    elif case == "4_flat_funcActivated.json": m = 0
    elif case == "5_dih_funcActivated.json": m = 0.1
    elif case == "6_anh_funcActivated.json": m = -0.1
    elif case == "7_flat_funcActivated.json": m = 0
    elif case == "8_dih_funcActivated.json": m = 0.1
    elif case == "9_anh_funcActivated.json": m = -0.1
    elif case == "10_flat_funcActivated.json": m = 0
    elif case == "11_dih_funcActivated.json": m = 0.1
    elif case == "12_anh_funcActivated.json": m = -0.1
    elif case == "D150_AGILE_Hangar_funActivated.xml": m = 1e-5
    elif case == "Boxwing_AGILE_Hangar_funActivated_v3.1.xml": m = 1e-5
    elif case == "Optimale_Tornado_SU2_funActivated.xml": m = 0.1
    else: logger.warning("Deformation input UNEXPECTED")
    
    h = 0
    u[0][:,2] = m*y[0] + h
    u[1][:,2] = m*y[1] + h
    u[2][:,2] = m*y[2] + h
    u[3][:,2] = m*y[3] + h
    
    headers = ["x","y","z","dx","dy","dz"]
    points = np.concatenate((lattice.c,u[2]),axis=1)
    
    if csv_save == True: 
        # TODO, change the way the path is selected
        filepath = settings.paths('f_deformation')
        dataset = pd.DataFrame(points,columns=headers)
        absolute_path = "/home/cfse2/Documents/pytornado/tests/integration/mesh_interfacing/wkdir/23_OptiMale/deformation/"
        logger.error(filepath)
        pd.DataFrame(dataset).to_csv(absolute_path+"dataset.csv",
                                     index=False,
                                     float_format='%.18E')
    return u

# def shape_3(z,settings):
#     logger.info("Shape function 3 is selected")
#     h = 0.0
#     m = 0.1
#     uy = m*z + h
#     return uy

def csv_def_load(u,lattice,path):
    """
    Loads a displacement file of format .csv and up

    Returns
    -------
    RBF explained
    https://www.youtube.com/watch?v=OOpfU3CvUkM
    
    None.
    TODO: check if the csv matches the mesh
    TODO: apply deformation to the mesh
    """
    logger.debug("====== load displacement function called ==========")
    logger.debug(path)
    displacements = pd.read_csv(path)
    logger.debug("displacements: \n" + str(displacements))
    x = displacements["x"] + displacements["dx"]
    y = displacements["y"] + displacements["dy"]
    z = displacements["z"] + displacements["dz"]
    newShape = interpolate.Rbf(x,y,z, function='multiquadric')
    # logger.debug("displacement file:\n" + str(displacements.head))
    # logger.debug("actual points: \n" + str(lattice.p))

def csv_builder():
    pass

def mapper(lattice,file):
    pass
    
def deformation(lattice,settings):
    """    
    This function deforms the mesh, computes the normal vector and feeds back
    the computed parameters into the lattice class variable.
    
    The stdrun.run function will then continue to compute the simulation with
    the deformed mesh.
    
    Parameters
    ----------
    lattice : class variable
        Variable of class lattice. This variable contains all the mesh points, 
        collocations points, normal vectors, horseshoe vortex points, surface 
        area of each cell, and some other non-essential information in this 
        part.
    settings : class variable
        Variable of class settings. This variable is used for checking which
        simulation should be done especially during the debug and testing 
        phase. It probably will loose it's importance when the program is 
        fully functionnal.
    Returns
    -------
    None.
    
    TODO: clean the way deformation is made by using an np.array
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
    # Normal vector of each cell
    s_b = lattice.bound_leg_midpoints.shape
    
    # Since the program is embeded in classes, it looks like a copy of the 
    # variables is necessary in order to modifiy only the variables and not 
    # the class variable. just putting var = class_var and then modifying
    # class var will modify var and class_var.
        
    # Mesh points
    lattice_p = np.copy(lattice.p.reshape((s_p[0]*s_p[1],s_p[2])))
    # HSV (Horseshoe vortex points)
    lattice_v = np.copy(lattice.v.reshape((s_v[0]*s_v[1],s_v[2])))
    # Collocation point
    lattice_c = np.copy(lattice.c)
    # Normal vector
    lattice_n = np.copy(lattice.n)
    # Lattice blm (boundleg midpoint)
    lattice_b = np.copy(lattice.bound_leg_midpoints)
    
    # Computes the initial normal vector. SVD has a proprety that in the vh
    # matrix, all the vectors are unitary orthonormal vectors. This permits
    # to have the new reference frame and compute the angles between old
    # and new reference frame.
    G = np.concatenate((lattice.c, lattice.c, lattice.c, lattice.c), axis=1)
    mat = lattice.p - G.reshape(s_p[0],s_p[1],s_p[2])
    u, s, vh_i = np.linalg.svd(mat)
    
    # Deformation of the wing, kept absolute to keep xz plane symmetry
    y_p = np.abs(lattice_p[:,1])
    y_v = np.abs(lattice_v[:,1])
    y_c = np.abs(lattice_c[:,1])
    y_b = np.abs(lattice_b[:,1])
    y = np.array([y_p,y_v,y_c,y_b])
    
    # Computes the deformed surface. u_ will be the displacement that needs to
    # be added to the given point.
    # TODO reduce this 24 lines to 10
    u_p = np.zeros((s_p[0]*s_p[1],s_p[2]))
    u_v = np.zeros((s_v[0]*s_v[1],s_v[2]))
    u_c = np.zeros((s_c[0],s_c[1]))
    u_b = np.zeros((s_b[0],s_b[1]))
    def_u = np.array([u_p,u_v,u_c,u_b])
    
    # user input choice TODO: put this into the settings file
    shape_func = 1
    filename = "/deformation/dataset.csv"
    path = str(settings.project_dir) + filename
    if shape_func == 1:   def_u = shape_1(def_u,y,settings,lattice)
    # elif shape_func == 2: shape_2(u,y,settings)
    # elif shape_func == 3: shape_3(u,y,settings)
    elif shape_func == 4: def_u = csv_def_load(def_u,settings,path)
    else: logger.info("No shape function selected")
    
    # Adds the displacement on all the mesh points
    lattice_p = lattice_p + def_u[0]
    lattice_v = lattice_v + def_u[1]
    lattice_c = lattice_c + def_u[2]
    lattice_b = lattice_b + def_u[3]
    
    # Reshapes array, and feeds the values back to pytornado format for the
    # solver to continue its task
    lattice.p = lattice_p.reshape((s_p[0],s_p[1],s_p[2]))
    lattice.v = lattice_v.reshape((s_p[0],s_v[1],s_v[2]))
    lattice.c = lattice_c
    lattice.bound_leg_midpoints = lattice_b
    
    # Computes the deformed reference frame values by using the same SVD 
    # proprety as before.
    G = np.concatenate((lattice.c, lattice.c, lattice.c, lattice.c), axis=1)
    mat = lattice.p - G.reshape(s_p[0],s_p[1],s_p[2])
    u, s, vh_f = np.linalg.svd(mat)
    
    # Computes the roation pivot vector. (kind of a hinge axis for rotation)
    rot_g = np.cross(vh_f[:,2,:],vh_i[:,2,:])
    rot = rot_g/np.linalg.norm(rot_g,axis=1,)[:,np.newaxis]
    rot[np.isnan(rot)] = 0.0
    
    # Computes the angle between the intial and deformed normal vector
    # dot product of vector a (initial state) and b (deformed state).
    ab = np.einsum('ij,ij->i',vh_f[:,2,:],vh_i[:,2,:])
    a = LA.norm(vh_f[:,2,:], axis=1)
    b = LA.norm(vh_i[:,2,:], axis=1)
    angle = np.empty(len(ab))
    
    # selects the correct angle (clockwise or counterlockwise) depending on
    # which side the wing is (y>0 or y<0) and the spacial orientation of the 
    # rotation vector (rot) in the "x" direction
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
    
    # Corrects symmetry issues. SVD has a tendency to make the normal vector 
    # point upwards which is not what pytornado uses as reference 
    # (pytornado -> z axis pointing down).
    if lattice.n[:,2].mean() >= 0:
        quat = np.einsum('i,ij->ij',angle,rot)
        r = R.from_rotvec(quat)
        logger.debug("normal vector point UP")
    else:
        quat = np.einsum('i,ij->ij',-angle,rot)
        r = R.from_rotvec(quat)
        logger.debug("normal vector point DOWN")
        
    lattice_n = r.apply(lattice.n)
    if lattice_n[:,2].mean() >= 0:
        lattice_n = np.array((1.0,-1.0,-1.0))*lattice_n
    
    # Computes the new surface area by using a first order method. Could be
    # impoved but is consistent with how the C code computes the surface area.
    s1 = lattice.p[:,0]-lattice.p[:,1]
    s2 = lattice.p[:,2]-lattice.p[:,3]
    c1 = lattice.p[:,0]-lattice.p[:,3]
    c2 = lattice.p[:,1]-lattice.p[:,2]
    c = 0.5*(np.linalg.norm(c1,axis=1)+np.linalg.norm(c2,axis=1))
    s = 0.5*(np.linalg.norm(s1,axis=1)+np.linalg.norm(s2,axis=1))
    area = s*c
    # feeds result back to pytornado in the correrct format
    lattice.a = area