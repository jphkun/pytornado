#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:54:23 2020

@author: cfse2
"""


import pickle

# path = dir_path+"/wkdir/"+folder
# path_d_dea = path + "/data_defDeactivated.pkl"
# path_l_dea = path + "/data_defActivated.pkl"
d_dea = open("data_defDeactivated.pkl", 'rb')
l_act = open("lattice_defActivated.pkl", 'rb')
l_dea = open("lattice_defDeactivated.pkl", 'rb')
lattice_v_act = pickle.load(l_dea)
print(l_dea)
data_v_dea = pickle.load(d_dea)
lattice_v_dea = pickle.load(l_dea)