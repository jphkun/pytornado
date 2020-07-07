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
        ax.set_xlim((0,40))
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
        ax.set_xlim((0,40))
        ax.set_ylim((-20,20))
        ax.set_zlim((-5,35))
    else:
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_title("Airplane points")
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
    folders = ["1_pytornado_flat",
               "2_pytornado_dih",
               "3_pytornado_anh",
               "4_pytornado_flat_flaps",
               "5_pytornado_dih_flaps",
               "6_pytornado_anh_flaps",
               "7_pytornado_flat_flaps_swept",
               "8_pytornado_dih_flaps_swept",
               "9_pytornado_anh_flaps_swept",
               "10_pytornado_flat_flaps_aft",
               "11_pytornado_dih_flaps_aft",
               "12_pytornado_anh_flaps_aft",
               "21_AircraftMalo-std",
               "22_B7772VSP",
               "23_Boxwing",
               "24_BWB_102_VTP1",
               "25_BWB_ACFA_cpacs",
               "26_Circlewing",
               "28_D150",
               "31_OptiMale"
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
    
    # for i in range(len(folders)):
    #     plot_vec(lattice_dea[i][2],
    #               lattice_act[i][2],
    #               lattice_dea[i][3],
    #               lattice_act[i][3])
    
    # for i in range(len(folders)):
    #     plot_points(lattice_dea[i][0],lattice_act[i][0]) 
    # for i in range(len(folders)):
    #     plot_points(lattice_dea[i][1],lattice_act[i][1]) 
    # for i in range(len(folders)):
    #     plot_points(lattice_dea[i][5],lattice_act[i][5]) 
    
    
    # # Compares results with deformation function activated or not
    items = ["x","y","z","D","C","L","l","m","n"]
    # epsilon_dea = [data_dea[i].epsilon for i in range(len(folders))]
    # epsilon_act = [data_dea[i].epsilon for i in range(len(folders))]
    matrix_downwash_dea = [data_dea[i].matrix_downwash for i in range(len(folders))]
    matrix_downwash_act = [data_act[i].matrix_downwash for i in range(len(folders))]
    array_rhs_dea = [data_dea[i].array_rhs for i in range(len(folders))]
    array_rhs_act = [data_act[i].array_rhs for i in range(len(folders))]
    data_dea = [data_dea[i].forces for i in range(len(folders))]
    data_act = [data_act[i].forces for i in range(len(folders))]
    index =[folders[i][0:] for i in range(len(folders))]
    
    pd.options.display.float_format = '{0:4.2f}'.format
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 500)
    df1 = pd.DataFrame(data_dea,index,columns=items)
    df2 = pd.DataFrame(data_act,index,columns=items)
    error = 100*((df1-df2)/df2)
    
    print(df1.round(1))
    print(df2.round(1))
    print("Error %")
    print(error.round(1))
    
    print("\n")
    print("Matrix downwash error")
    print(np.max(matrix_downwash_dea[0]-matrix_downwash_act[0]))
    print(np.min(matrix_downwash_dea[0]-matrix_downwash_act[0]))
    
    print("array rhs error")
    print(np.max(array_rhs_dea[0]-array_rhs_act[0]))
    print(np.min(array_rhs_dea[0]-array_rhs_act[0]))
    
    # print(epsilon_dea)
    # print(epsilon_act)
    # path = '/home/cfse2/Documents/pytornado/tests/integration/mesh_interfacing/wkdir/23_OptiMale/deformation/'
    # filename = 'Optimale_disp.csv'
    # full_path = path + filename
    # print(np.shape(lattice_act[0]))
    
    # header = ["p","v","c","n","a","blm"]
    # for i in range(len(lattice_act)):
    #     print("\n Avion:" + str(i))
    #     for j in range(len(lattice_act[0])):
    #         # print( ("Avion 0, (pvca) " + "%10.3E"% str(i)))
    #         diff1 = np.max(lattice_act[i][j] - lattice_dea[i][j])
    #         diff2 = np.min(lattice_act[i][j] - lattice_dea[i][j])
    #         print(header[j])
    #         print("%2.6f"% (diff1))
    #         print("%2.6f"% (diff2))
    
    # print(lattice_act[0][4][:5])
    # print(lattice_dea[0][4][:5])

if __name__ == "__main__":
    main()