#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:17:34 2020

@author: cfse2
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
import os
import pandas as pd 
import numpy as np

def plot_vec_mult(dea_c,act_c,dea_n,act_n):
    """
    Plots the normal vectors
    """
    N = len(act_n)
    Nr = 2
    Nc = N/Nr
    fig = plt.figure(figsize=plt.figaspect(2/N))
    for i in range(N):
        ax = fig.add_subplot(Nr,Nc,i+1, projection='3d')
        ax.scatter(dea_c[i].c[:,0], dea_c[i].c[:,1], dea_c[i].c[:,2],color="blue")
        ax.quiver( dea_c[i].c[:,0], dea_c[i].c[:,1], dea_c[i].c[:,2],
                   dea_n[i].n[:,0], dea_n[i].n[:,1], dea_n[i].n[:,2],color="blue")
        ax.scatter(act_c[i].c[:,0], act_c[i].c[:,1], act_c[i].c[:,2],color="red")
        ax.quiver( act_c[i].c[:,0], act_c[i].c[:,1], act_c[i].c[:,2],
                   act_n[i].n[:,0], act_n[i].n[:,1], act_n[i].n[:,2],color="red")
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_title("case " + str(i+1))
    plt.show()

def plot_vec(dea_c,act_c,dea_n,act_n):
    """
    Plots the normal vectors
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dea_c[:,0],dea_c[:,1],dea_c[:,2],color="blue")
    ax.quiver( dea_c[:,0],dea_c[:,1],dea_c[:,2],
               dea_n[:,0],dea_n[:,1],dea_n[:,2],color="blue")
    ax.scatter(act_c[:,0],act_c[:,1],act_c[:,2],color="red")
    ax.quiver( act_c[:,0],act_c[:,1],act_c[:,2],
               act_n[:,0],act_n[:,1],act_n[:,2],color="red")

    if max(dea_c[:,0])> 7:
        ax.set_xlim((10,50))
        ax.set_ylim((-20,20))
        ax.set_zlim((-20,20))
    else:
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_title("airplane")
    plt.show()

def plot_points_mult(dea,act):
    """
    Plots the centroids
    """
    N = len(dea)
    Nr = 2
    Nc = N/Nr
    fig = plt.figure(figsize=plt.figaspect(2/N))
    for i in range(N):
        ax = fig.add_subplot(Nr,Nc,i+1, projection='3d')
        ax.scatter(dea[i].c[:,0], dea[i].c[:,1], dea[i].c[:,2],color="blue")
        ax.scatter(act[i].c[:,0], act[i].c[:,1], act[i].c[:,2],color="red")
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_title("case " + str(i+1))
        ax.view_init(elev=0., azim=i)
    plt.show()

def plot_points(dea,act):
    """
    Plots the centroids
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dea[:,0], dea[:,1], dea[:,2],color="blue")
    ax.scatter(act[:,0], act[:,1], act[:,2],color="red")
    if max(dea[:,0])> 7:
        ax.set_xlim((10,50))
        ax.set_ylim((-20,20))
        ax.set_zlim((-5,35))
    else:
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_title("Airplane")
    plt.show()
    
def load_displacement_file(filename):
    """
    Loads a displacement file of format .csv and up

    Returns
    -------
    None.

    """
    displacements = pd.read_csv(filename)
    print("====== load displacement function called ==========")
    return displacements

def main():
    arr = os.listdir('wkdir')
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    folders = []
    for folder_name in arr:
        if folder_name[0].isdigit():
            folders.append(folder_name)
    folders = sorted(folders)
    ### TODO delete this once it works
    folders2 = folders[:3]#+folders[9:]
    folders = folders2
    # print(arr)
    ###
    
    # Allocates memory to upload the different results into it
    data_dea = []
    data_act = []
    lattice_dea = []
    lattice_act = []
    folders = [#"1_pytornado_flat",
               #"2_pytornado_dih",
               #"3_pytornado_anh",
               #"4_pytornado_flat_flaps",
               #"5_pytornado_dih_flaps",
               #"6_pytornado_anh_flaps",
               #"7_pytornado_flat_flaps_swept",
               #"8_pytornado_dih_flaps_swept",
               #"9_pytornado_anh_flaps_swept",
               #"10_pytornado_flat_flaps_aft",
               #"11_pytornado_dih_flaps_aft",
               #"12_pytornado_anh_flaps_aft",
               #"21_D150",
               #"22_Boxwing",
               "23_OptiMale"
               ]
    for folder in folders:
        # Builds all the needed paths
        path = dir_path+"/wkdir/"+folder
        path_d_dea = path + "/data_defDeactivated.pkl"
        path_l_dea = path + "/lattice_defDeactivated.pkl"
        path_d_act =path + "/data_defActivated.pkl"
        path_l_act = path + "/lattice_defActivated.pkl"
        # Opens up the .pkl files
        d_dea = open(path_d_dea, 'rb')
        l_dea = open(path_l_dea, 'rb')
        d_act = open(path_d_act, 'rb')
        l_act = open(path_l_act, 'rb')
        # Loads the .pkl values into variables
        data_v_dea = pickle.load(d_dea)
        lattice_v_dea = pickle.load(l_dea)
        dat_v_act = pickle.load(d_act)
        lattice_v_act = pickle.load(l_act)
        # Closes the .pkl files
        d_dea.close()
        l_dea.close()
        d_act.close()
        l_act.close()
        # Concatenate all the .pkl for an easy access
        data_dea.append(data_v_dea)
        lattice_dea.append(lattice_v_dea)
        data_act.append(dat_v_act)
        lattice_act.append(lattice_v_act)
    
    # Plots normal vectors and cell centroids
    # print(lattice_dea)
    for i in range(len(folders)):
        # plot_points(lattice_dea[i][2],lattice_act[i][2])
        plot_vec(lattice_dea[i][2],
                  lattice_act[i][2],
                  lattice_dea[i][3],
                  lattice_act[i][3])
    

    # # Compares results with deformation function activated or not
    items = ["x","y","z","D","C","L","l","m","n"]
    data_dea = [data_dea[i].forces for i in range(len(folders))]
    data_act = [data_act[i].forces for i in range(len(folders))]
    index =[folders[i][12:] for i in range(len(folders))]
    pd.set_option('display.max_columns', None)
    df1 = pd.DataFrame(data_dea,index,columns =items)
    df2 = pd.DataFrame(data_act,index,columns =items)
    error = 100*((df1-df2)/df2)
    
    print(df1)
    print(df2)
    print(error)
    
    
    
    path = '/home/cfse2/Documents/pytornado/tests/integration/mesh_interfacing/wkdir/23_OptiMale/deformation/'
    filename = 'Optimale_disp.csv'
    full_path = path + filename
    dlr = load_displacement_file(full_path)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dlr["x"],dlr["y"],dlr["z"],color="blue")
    ax.scatter(lattice_dea[0][2][:,0],
               lattice_dea[0][2][:,1],
               lattice_dea[0][2][:,2],color="red") 
    
if __name__ == "__main__":
    main()