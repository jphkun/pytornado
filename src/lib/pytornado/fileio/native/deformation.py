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
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R
from scipy import interpolate
from numpy.core.umath_tests import inner1d


logger = logging.getLogger(__name__)


# def load(aircraft, settings):
#     """
#     Loads the aircraft deformation file if it exitsts.

#     Args:
#         :aircraft: (obj) aircraft
#         :settings: (obj) settings
#     """

#     filepath = settings.paths('f_deformation')
#     logger.info(f"Reading deformation from file '{truncate_filepath(filepath)}'")

#     if not os.path.exists(filepath):
#         raise IOError(f"file '{filepath}' not found")
#     # File is empty or as good as (this also catches empty JSON file: '{}')
#     elif os.stat(filepath).st_size < 10:
#         logger.warning(f"Empty deformation file. No deformations are modelled.")
#         return

#     # def_fields = load_json_def_fields(filepath)
#     # for wing_uid, def_field in def_fields.items():
#     #     ####
#     #     # Convention: Deformation fields starting with '_' will be ignore by CFD
#     #     if wing_uid.startswith('_'):
#     #         continue
#     #     ####

#     #     # Convention: Deformation fields ending with '_m' belongs to a 'mirrored' component
#     #     if wing_uid.endswith('_m'):
#     #         aircraft.wings[wing_uid.replace('_m', '')].def_field_mirror = def_field
#     #         continue

#     #     aircraft.wings[wing_uid].def_field = def_field

#     # TODO: Handle exceptions
#     # TODO: Check deformation continuity

class Mesh_Def:

    def __init__(self, lattice):
        """
        *_p : lattice points
        *r_ p: lattice points reshaped array for ease of use
        *_v : horseshoe vortex points
        *r_v : horseshoe vortex points reshaped for ease of use
        *_c : cell collocation points
        *_n : normal vector directions
        *_b : bound leg midpoints
        f_a : final area of each individual panels
        u_* : defomation of the given type of point. follows the same naming
              pattern as for the previous 5 items
            
        ----------
        lattice : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # stores lattice shapes (only theses two are needed, the others are
        # of identical shape)
        self.s_p = lattice.p.shape
        self.s_v = lattice.v.shape
        self.s_c = lattice.c.shape
        self.s_b = lattice.bound_leg_midpoints.shape
        # stores lattice intitial (i_) data
        self.i_p = np.copy(lattice.p)
        self.ir_p = np.copy(lattice.p.reshape((self.s_p[0]*self.s_p[1],self.s_p[2])))
        self.i_v = np.copy(lattice.v)
        self.ir_v = np.copy(lattice.v.reshape((self.s_v[0]*self.s_v[1],self.s_v[2])))
        self.i_c = np.copy(lattice.c)
        self.i_n = np.copy(lattice.n)
        self.i_b = np.copy(lattice.bound_leg_midpoints)
        # stores lattice final (f_) data
        self.f_p = np.zeros([self.s_p[0],self.s_p[1],self.s_p[2]])
        self.fr_p = np.zeros([self.s_p[0]*self.s_p[1],self.s_p[2]])
        self.f_v = np.zeros([self.s_p[0],self.s_p[1],self.s_p[2]])
        self.fr_v = np.zeros([self.s_p[0]*self.s_p[1],self.s_p[2]])
        self.f_c = np.zeros([self.s_c[0]*self.s_c[1]])
        self.f_n = np.zeros([self.s_c[0]*self.s_c[1]])
        self.f_b = np.zeros([self.s_c[0]*self.s_c[1]])
        self.f_a = np.zeros([self.s_c[0]*self.s_c[1]])
        # Cells absolute y corrdinates (needed for testing and debug)
        # Deformation of the wing, kept absolute to keep xz plane symmetry
        self.y_p = np.abs(self.ir_p[:,1])
        self.y_v = np.abs(self.ir_v[:,1])
        self.y_c = np.abs(self.i_c[:,1])
        self.y_b = np.abs(self.i_b[:,1])
        # Mesh displacement
        self.u_p = np.zeros((self.s_p[0]*self.s_p[1],self.s_p[2]))
        self.u_v = np.zeros((self.s_p[0]*self.s_p[1],self.s_p[2]))
        self.u_c = np.zeros((self.s_c[0],self.s_c[1]))
        self.u_b = np.zeros((self.s_c[0],self.s_c[1]))

    def shape_2(self):
        """
        Shapte function2. This function computes the slope at each y location for
        a cantilever beam of length "L", with a distributed load of "q". The Young
        modulus is imposed for steel and "I" the second moment of inertia is
        also imposed.
        """
        logger.info("Shape function 2 is selected")
        # [N/m] Distributed load
        q = 500
        # [m] Wing span
        L = 5
        # logger.debug("L = " + str(L) )
        # [Pa] Elasticity modulus
        E = 210e9
        # [m**4] Second moment of area
        I = 1.330e-6
        # Computes beam deformation
        uz = q*self.y**2 *(6*L**2 -4*L*self.y +self.y**2)/(24*E*I)
        
        return uz
    
    def shape_1(self,settings):
        """
        Shape function 1. This functions respresents line of slope 1. The only 
        purpouse of this shape function is to test if the deformation is done 
        correctly.
        """
        csv_save = True
        const = 1e-4
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
        elif case == "AircraftMalo-std_funActivated.xml": m = const
        elif case == "B7772VSP_v3.1_funActivated.xml": m = const
        elif case == "BWB_102_VTP1_v3.1_funActivated.xml": m = const
        elif case == "BWB_ACFA_cpacs_v3.1_funActivated.xml": m = const 
        elif case == "Circlewing_Test.v_3.1_funActivated.xml": m = const
        elif case == "D150_AGILE_Hangar_funActivated.xml": m = const
        elif case == "Boxwing_AGILE_Hangar_funActivated_v3.1.xml": m = const
        elif case == "Optimale_Tornado_SU2_funActivated.xml": m = const
        else: logger.warning("Deformation input UNEXPECTED")
        h = 0
        self.u_p[:,2] = m * self.y_p + h
        self.u_v[:,2] = m * self.y_v + h
        self.u_c[:,2] = m * self.y_c + h
        self.u_b[:,2] = m * self.y_b + h

        
        self.fr_p = self.ir_p + self.u_p
        self.f_p = self.i_p + self.u_p.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        
        self.fr_v = self.ir_v + self.u_v
        self.f_v = self.i_v + self.u_v.reshape(self.s_v[0],self.s_v[1],self.s_v[2])
        self.f_c = self.i_c + self.u_c
        self.f_b = self.i_b + self.u_b
        
        # headers = ["x","y","z","dx","dy","dz"]
        # points = np.concatenate((self.i_c,self.u_c),axis=1)
        # if csv_save == True: 
        #     # TODO, change the way the path is selected
        #     filepath = settings.paths('f_deformation')
        #     dataset = pd.DataFrame(points,columns=headers)
        #     absolute_path = "/home/cfse2/Documents/pytornado/tests/integration/mesh_interfacing/wkdir/23_OptiMale/deformation/"
        #     pd.DataFrame(dataset).to_csv(#absolute_path+"dataset.csv",
        #                                   filepath,
        #                                   index=False,
        #                                   float_format='%.18E')
    
    def csv_def_load(self,u,lattice,path):
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
        # # Uploads deformed states into class
        # self.u_p = 
        # self.u_v = 
        # self.u_c =
        # self.u_b = 
        # Addition of initial value and deformation.
        self.fr_p = self.ir_p + self.u_p
        self.f_p = self.fr_p.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        self.fr_v = self.ir_v + self.u_v
        self.f_v = self.fr_v.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        self.f_c = self.i_c + self.u_c
        self.f_b = self.i_b + self.u_b
    
    def csv_builder():
        pass
    
    def mapper(lattice,file):
        pass
        
    def deformation(self,settings):
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
        """
        
        logger.info("=== Starts deformation function ===")
        
        # Computes the initial normal vector. SVD has a proprety that in the vh
        # matrix, all the vectors are unitary orthonormal vectors. This permits
        # to have the new reference frame and compute the angles between old
        # and new reference frame.
        G = np.concatenate((self.i_c, self.i_c, self.i_c, self.i_c), axis=1)
        mat = self.i_p - G.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        u, s, vh_i = LA.svd(mat)
        
        # user input choice 
        # TODO: put this into the settings file
        shape_func = 1
        if shape_func==1:   self.shape_1(settings)
        # elif shape_func == 2: shape_2(u,y,settings)
        else: logger.error("No shape function selected")

        # Computes the deformed reference frame values by using the same SVD 
        # proprety as before.
        G = np.concatenate((self.f_c, self.f_c, self.f_c, self.f_c), axis=1)
        mat = self.f_p - G.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        u, s, vh_f = LA.svd(mat)
        
        # Computes the roation pivot vector. (kind of a hinge axis for rotation)
        rot_g = np.cross(vh_f[:,2,:],vh_i[:,2,:])
        rot = rot_g/np.linalg.norm(rot_g,axis=1,)[:,np.newaxis]
        rot[np.isnan(rot)] = 0.0
        
        # Computes the angle between the intial and deformed normal vector
        # dot product of vector a (initial state) and b (deformed state).
        # ab = np.einsum('ij,ij->i',vh_f[:,2,:],vh_i[:,2,:])
        ab = inner1d(vh_f[:,2,:],vh_i[:,2,:])
        a = LA.norm(vh_f[:,2,:], axis=1)
        b = LA.norm(vh_i[:,2,:], axis=1)
        
        # Computes the new normal vector
        angle = np.arccos(ab/(a*b))
        angle[np.isnan(angle)] = 0.0
        # Some angles might be greater than pi/2 hence a correction is needed
        corrector = np.zeros(angle.shape)
        corrector[angle>np.pi/2]=np.pi
        angle = angle-corrector
        # Rotates the vector using a quaternion
        quat = np.einsum('i,ij->ij',-angle,rot)
        r = R.from_rotvec(quat)
        self.f_n = r.apply(self.i_n)
        
        # Computes the new surface area by using a first order method. Could be
        # impoved but is consistent with how the C code computes the surface area.
        s1 = 0.5*(self.f_p[:,1] - self.f_p[:,0] + self.f_p[:,2] - self.f_p[:,3])
        s = LA.norm(s1,axis=1)
        c1 = 0.5*(self.f_p[:,3] - self.f_p[:,0] + self.f_p[:,2] - self.f_p[:,1])
        c = LA.norm(c1,axis=1)
        # Computes area
        self.f_a = s*c