#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:17:34 2020

@author: cfse2
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
import os
import pandas as pd 

def plot_vec_mult(dea,act):
    """
    Plots the normal vectors
    """
    N = len(dea)
    Nr = 2
    Nc = N/Nr
    fig = plt.figure(figsize=plt.figaspect(2/N))
    for i in range(N):
        ax = fig.add_subplot(Nr,Nc,i+1, projection='3d')
        ax.scatter(dea[i].c[:,0], dea[i].c[:,1], dea[i].c[:,2],color="blue")
        ax.quiver( dea[i].c[:,0], dea[i].c[:,1], dea[i].c[:,2],
                   dea[i].n[:,0], dea[i].n[:,1], dea[i].n[:,2],color="blue")
        ax = fig.add_subplot(Nr,Nc,i+1, projection='3d')
        ax.scatter(act[i].c[:,0], act[i].c[:,1], act[i].c[:,2],color="red")
        ax.quiver( act[i].c[:,0], act[i].c[:,1], act[i].c[:,2],
                   act[i].n[:,0], act[i].n[:,1], act[i].n[:,2],color="red")
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_title("case " + str(i+1))
    plt.show()

def plot_vec(dea,act):
    """
    Plots the normal vectors
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dea.c[:,0],dea.c[:,1], dea.c[:,2],color="blue")
    ax.quiver( dea.c[:,0],dea.c[:,1], dea.c[:,2],
               dea.n[:,0],dea.n[:,1], dea.n[:,2],color="blue")
    ax.scatter(act.c[:,0],act.c[:,1],act.c[:,2],color="red")
    ax.quiver( act.c[:,0],act.c[:,1],act.c[:,2],
                act.n[:,0],act.n[:,1],act.n[:,2],color="red")
    ax.set_xlim((10,50))
    ax.set_ylim((-20,20))
    ax.set_zlim((-15,25))
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
    ax.scatter(dea.c[:,0], dea.c[:,1], dea.c[:,2],color="blue")
    ax.scatter(act.c[:,0], act.c[:,1], act.c[:,2],color="red")
    ax.set_xlim((10,50))
    ax.set_ylim((-20,20))
    ax.set_zlim((-5,35))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_title("Airplaine")
    plt.show()

def main():
    arr = os.listdir('wkdir')
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    folders = []
    for folder_name in arr:
        if folder_name[0].isdigit():
            folders.append(folder_name)
    folders = sorted(folders)
    ### TODO delete this once it works
    # folders2 = folders[:3]+folders[9:]
    # folders = folders2
    print(arr)
    ###
    
    # Allocates memory to upload the different results into it
    data_dea = []
    data_act = []
    lattice_dea = []
    lattice_act = []
    folders = ["22_Boxwing"]
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
    
    # plot_points(lattice_dea[0],lattice_act[0])
    plot_vec(lattice_dea[0],lattice_act[0])
    

    # plot_points_mult(lattice_dea,lattice_act)
    # plot_vec_mult(lattice_dea,lattice_act)
    # Compares results with deformation function activated or not
    items = ["x","y","z","D","C","L","l","m","n"]
    
    # TODO compute the error here!
    data_dea = [data_dea[i].forces for i in range(len(folders))]
    data_act = [data_act[i].forces for i in range(len(folders))]
    index =[folders[i][12:] for i in range(len(folders))]
    
    df1 = pd.DataFrame(data_dea,index,columns =items)
    df2 = pd.DataFrame(data_act,index,columns =items)
    
    df1.to_csv('dea.csv')
    df2.to_csv('act.csv')
    
    # print(folders)
    # print(df1)
    # print(df2)

if __name__ == "__main__":
    main()





























# diffla  = lajson - lafunc
# difflp  = lpjson - lpfunc
# difflv  = lvjson - lvfunc
# difflc  = lcjson - lcfunc
# diffln  = lnjson - lnfunc
# diffDwM = Defjson- Deffunc
# diffrhs = rhsjson - rhsfunc
# diffgam = gammajson - gammafunc
# diffblm = blmjson - blmfunc

# print("\n max: \n")
# print("la ",diffla.max())
# print("lp ",difflp.max())
# print("lv ",difflv.max())
# print("lc ",difflc.max())
# print("ln ",diffln.max())
# print("DwM",diffDwM.max())
# print("rhs",diffrhs.max())
# print("gam",diffgam.max())
# print("blm",diffblm.max())

# print("\n min: \n")
# print("la ",diffla.min())
# print("lp ",difflp.min())
# print("lv ",difflv.min())
# print("lc ",difflc.min())
# print("ln ",diffln.min())
# print("DwM",diffDwM.min())
# print("rhs",diffrhs.min())
# print("gam",diffgam.min())
# print("blm",diffblm.min())

# print("\n mean: \n")
# print("la ",diffla.mean())
# print("lp ",difflp.mean())
# print("lv ",difflv.mean())
# print("lc ",difflc.mean())
# print("ln ",diffln.mean())
# print("DwM",diffDwM.mean())
# print("rhs",diffrhs.mean())
# print("gam",diffgam.mean())
# print("blm",diffblm.mean())

# plt.matshow(diffDwM)
# plt.matshow(difflp)
# plt.matshow(difflv)
# plt.matshow(difflc)
# plt.matshow(diffln)

# # reshapes lattice._ in a more user-friendly way
# # Mesh points
# s_p = lpjson.shape
# # 4 points of the vortex horseshoe
# s_v = lvjson.shape
# # Collocation point: point of calculation
# s_c = lcjson.shape
# # Normal vector of each cell
# s_n = lnjson.shape