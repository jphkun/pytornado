#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:17:34 2020

@author: cfse2
"""

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import pickle
import os
import pandas as pd
import numpy as np


def plot_vec_mult(dea_c, act_c, dea_n, act_n):
    """
    Plots the normal vectors
    """
    N = len(act_n)
    Nr = 2
    Nc = N/Nr
    fig = plt.figure(figsize=plt.figaspect(2/N))
    for i in range(N):
        ax = fig.add_subplot(Nr, Nc, i+1, projection='3d')
        ax.scatter(dea_c[i].c[:,0], dea_c[i].c[:,1], dea_c[i].c[:,2],
                   color="blue")
        ax.quiver(dea_c[i].c[:,0], dea_c[i].c[:,1], dea_c[i].c[:,2],
                  dea_n[i].n[:,0], dea_n[i].n[:,1], dea_n[i].n[:,2],
                  color="blue")
        ax.scatter(act_c[i].c[:,0], act_c[i].c[:,1], act_c[i].c[:,2],
                   color="red")
        ax.quiver(act_c[i].c[:,0], act_c[i].c[:,1], act_c[i].c[:,2],
                  act_n[i].n[:,0], act_n[i].n[:,1], act_n[i].n[:,2],
                  color="red")

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
    ax.quiver(dea_c[:,0],dea_c[:,1],dea_c[:,2],
              dea_n[:,0],dea_n[:,1],dea_n[:,2],color="blue")
    ax.scatter(act_c[:,0],act_c[:,1],act_c[:,2],color="red")
    ax.quiver(act_c[:,0],act_c[:,1],act_c[:,2],
              act_n[:,0],act_n[:,1],act_n[:,2],color="red")

    if max(dea_c[:,0]) > 7:
        ax.set_xlim((-10,30))
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


def plot_hsv(dea,act,title):
    """
    Plots the centroids
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dea[:,0], dea[:,1], dea[:,2],color="blue")
    ax.scatter(act[:,0], act[:,1], act[:,2],color="red")

    if max(dea[:,0]) > 7:
        ax.set_xlim((0,20))
        ax.set_ylim((-20,20))
        ax.set_zlim((-5,5))
    else:
        ax.set_xlim((-5,5))
        ax.set_ylim((-5,5))
        ax.set_zlim((-5,5))

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.set_title(title)
    plt.show()


def plot_points(dea,act,title):
    """
    Plots the centroids
    """
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    ax.scatter(dea[:,0], dea[:,1], dea[:,2],color="blue")
    ax.scatter(act[:,0], act[:,1], act[:,2],color="red")

    if max(dea[:,0]) > 7:
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
    ax.set_title(title)
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

def panelwise(i,data_dea,data_act,name,plot):
    # Individual fz panel data
    dea = data_dea[i].panelwise[name]
    act = data_act[i].panelwise[name]
    err = 100 * (act - dea) / dea
    
    if plot:
        # print(rhs_dea.shape)
        plt.matshow(np.concatenate([[dea],[dea]]).T)
        plt.title(name + " vector def. DEACTIVATED")
        plt.colorbar()
        plt.show()
        plt.matshow((np.concatenate([[act],[act]]).T))
        plt.title(name + " vector def. ACTIVATED")
        plt.colorbar()
        plt.show()
        plt.matshow((np.concatenate([[err],[err]]).T))
        plt.title(name + " vector def. ERROR")
        plt.show()
        plt.colorbar()
    return err

def main():
    # When to start and stop or plot
    start = 0
    stop = 20
    plot_p = False
    plot_v = False
    plot_n = True
    plot_blm = False
    downwash = False
    rhs_plot = False
    gamma_plot = False
    vm_plot = False
    vx_plot = False
    vy_plot = False
    vz_plot = False
    fx_plot = False
    fy_plot = False
    fz_plot = False

    # Pandas visal settings
    pd.set_option('expand_frame_repr', False)
    pd.options.display.float_format = '{:,.2e}'.format
    pd.set_option('max_colwidth', 10)

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
    folders = folders[start:stop]

    dir_path = os.path.dirname(os.path.realpath(__file__))
    for folder in folders:
        # Builds all the needed paths
        path = dir_path+"/wkdir/"+folder
        path_d_dea = path + "/data_defDeactivated.pkl"
        path_l_dea = path + "/lattice_defDeactivated.pkl"
        path_d_act = path + "/data_defActivated.pkl"
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

    # plots mesh points
    if plot_p:
        for i in range(len(folders)):
            l_dea = np.array(lattice_dea[i][0])
            s = l_dea.shape
            plot_points(lattice_dea[i][0].reshape((s[0]*s[1],s[2])),
                        lattice_act[i][0].reshape((s[0]*s[1],s[2])),
                        title="lattice.p")

    # Plots Horseshoe points
    if plot_v:
        for i in range(len(folders)):
            l_dea = np.array(lattice_dea[i][1])
            s = l_dea.shape
            plot_hsv(lattice_dea[i][1].reshape(s[0]*s[1],s[2]),
                     lattice_act[i][1].reshape(s[0]*s[1],s[2]),
                     title="HSV")

    # Plots normal vectors and cell centroids
    if plot_n:
        for i in range(len(folders)):
            plot_vec(lattice_dea[i][2],
                     lattice_act[i][2],
                     lattice_dea[i][3],
                     lattice_act[i][3])

    # Plots bounleg midpoints
    if plot_blm:
        for i in range(len(folders)):
            plot_points(lattice_dea[i][5],
                        lattice_act[i][5],
                        title="lattice.blm")

    # TODO: Implement everything with a for loop
    # TODO: have a recap table for forces that has the same size everywhre
    # TODO: have a recap table for p v c n a b for each cases
    # TODO: have the min values, the max values, the % error and the absolute
    #       error for all the p v c n a b cases.

    # Compares results with deformation function activated or not
    items = ["x", "y", "z", "D", "C", "L", "l", "m", "n"]
    index = [folders[i][0:] for i in range(len(folders))]
    forces_dea = pd.DataFrame(columns=items)
    forces_act = pd.DataFrame(columns=items)
    coeffs_dea = pd.DataFrame(columns=items)
    coeffs_act = pd.DataFrame(columns=items)

    # Allocates names to forces and coefficient variables
    items2 = ["Gamma max","Gamma min",
              "DWM max","DWM min",
              "rhs max","rhs min"]
    sol_err = pd.DataFrame(columns=items2)

    # Uploads VLM data into diffenrent variables for ease of use
    for i in range(0,len(folders)):
        # Downwash matrix
        matrix_dea = data_dea[i].matrix_downwash
        matrix_act = data_act[i].matrix_downwash
        matrix_dea_sign = np.empty(matrix_dea.shape)
        matrix_dea_sign[matrix_dea < 0] = 1
        matrix_dea_sign[matrix_dea > 0] = -1
        if downwash:
            plt.matshow(matrix_dea_sign)
            plt.title("Downwash matrix def. DEACTIVATED")
            plt.show()
        
        matrix_act_sign = np.empty(matrix_act.shape)
        matrix_act_sign[matrix_act < 0] = 1
        matrix_act_sign[matrix_act > 0] = -1
        if downwash:
            plt.matshow(matrix_act_sign)
            plt.title("Downwash matrix def. ACTIVATED")
            plt.show()
        if downwash:
            plt.matshow(matrix_act_sign - matrix_dea_sign)
            plt.title("Sign difference")
            plt.colorbar()
            plt.show()
        err_m = 100 * (matrix_dea - matrix_act) / matrix_dea
        if downwash:
            plt.matshow(err_m)
            plt.title("Downwash matrix ERROR")
            plt.colorbar()
            plt.show()
            plt.spy(err_m)
        
        # RHS vector
        rhs_dea = data_dea[i].array_rhs
        rhs_act = data_act[i].array_rhs
        err_rhs = 100 * (rhs_dea - rhs_act) / rhs_dea
        
        if rhs_plot:
            plt.matshow(np.concatenate([[rhs_dea],[rhs_dea]]).T)
            plt.title("RHS matrix def. DEACTIVATED")
            plt.colorbar()
            plt.show()
            plt.matshow((np.concatenate([[rhs_act],[rhs_act]]).T))
            plt.title("RHS matrix def. ACTIVATED")
            plt.colorbar()
            plt.show()
            plt.matshow((np.concatenate([[err_rhs],[err_rhs]]).T))
            plt.title("RHS matrix def. ERROR")
            plt.colorbar()
            plt.show()
        
        # print(data_dea[i].matrix_lu)
        # print(data_dea[i].array_pivots)

        # Individual gamma panel data
        gamma_dea = data_dea[i].panelwise['gamma']
        gamma_act = data_act[i].panelwise['gamma']
        err_g = 100 * (gamma_act - gamma_dea) / gamma_act
        sol_err = sol_err.append({"Gamma max": np.max(err_g),
                                  "Gamma min": np.min(err_g),
                                  "DWM max": np.max(err_m),
                                  "DWM min": np.min(err_m),
                                  "rhs max": np.max(err_rhs),
                                  "rhs min": np.min(err_rhs)
                                  },
                                 ignore_index=True)
        if gamma_plot:
            # print(rhs_dea.shape)
            plt.matshow(np.concatenate([[gamma_dea],[gamma_dea]]).T)
            plt.title("Gamma vector def. DEACTIVATED")
            plt.colorbar()
            plt.show()
            plt.matshow((np.concatenate([[gamma_act],[gamma_act]]).T))
            plt.title("Gamma vector def. ACTIVATED")
            plt.colorbar()
            plt.show()
            plt.matshow((np.concatenate([[err_g],[err_g]]).T))
            plt.title("Gamma vector def. ERROR")
            plt.colorbar()
            plt.show()
        
        panelwise(i,data_dea,data_act,'vmag',vm_plot)
        panelwise(i,data_dea,data_act,'vx',vx_plot)
        panelwise(i,data_dea,data_act,'vy',vy_plot)
        panelwise(i,data_dea,data_act,'vz',vz_plot)
        panelwise(i,data_dea,data_act,'fx',fx_plot)
        panelwise(i,data_dea,data_act,'fy',fy_plot)
        panelwise(i,data_dea,data_act,'fz',fz_plot)
        

        # Whole aircraft data
        forces_dea = forces_dea.append(data_dea[i].forces,ignore_index=True)
        forces_act = forces_act.append(data_act[i].forces,ignore_index=True)
        coeffs_dea = coeffs_dea.append(data_dea[i].coeffs,ignore_index=True)
        coeffs_act = coeffs_act.append(data_act[i].coeffs,ignore_index=True)

    # Computation of the error
    err_forces = 100 * (forces_act - forces_dea)/forces_dea
    err_coeffs = 100 * (coeffs_act - coeffs_dea)/coeffs_dea

    # Forces
    print("Forces when displacement function is DEACTIVATED:")
    print(forces_dea)
    print("Forces when displacement function is ACTIVATED:")
    print(forces_act)
    print("Difference in % :")
    pd.options.display.float_format = '{:,.4f}'.format
    print(err_forces)
    print("\n")

    # # Coefficients, the errors are exactly the same as the others, hence
    # # no need to print them
    # print("Coeficients when displacement function is DEACTIVATED:")
    # print(coeffs_dea)
    # print("Coeficients when displacement function is ACTIVATED:")
    # print(coeffs_act)
    # print("Difference in % :")
    # print(err_coeffs)
    # print("\n")

    # Coefficients
    print("Solution error in %")
    print(sol_err)
    print("\n")

    header = ["p max","p min",
              "v max","v min",
              "c max","c min",
              "n max","n min",
              "a max","a min",
              "b max","b min"]
    mesh_err = pd.DataFrame(columns=header)

    # Uploads lattice values into separate variables for ease of use
    for i in range(len(lattice_act)):
        err_p = lattice_act[i][0] - lattice_dea[i][0]
        err_v = lattice_act[i][1] - lattice_dea[i][1]
        err_c = lattice_act[i][2] - lattice_dea[i][2]
        err_n = lattice_act[i][3] - lattice_dea[i][3]
        err_a = lattice_act[i][4] - lattice_dea[i][4]
        err_b = lattice_act[i][5] - lattice_dea[i][5]
        mesh_err = mesh_err.append({"p max": np.max(err_p),
                                    "p min": np.min(err_p),
                                    "v max": np.max(err_v),
                                    "v min": np.min(err_v),
                                    "c max": np.max(err_c),
                                    "c min": np.min(err_c),
                                    "n max": np.max(err_n),
                                    "n min": np.min(err_n),
                                    "a max": np.max(err_a),
                                    "a min": np.min(err_a),
                                    "b max": np.max(err_b),
                                    "b min": np.min(err_b)},ignore_index=True)

    # pd.options.display.float_format = '{:,.3e}'.format
    # print("Absolute mesh error")
    # print(mesh_err)
    # print("\n")

    for i in range(len(lattice_act)):
        err_p = lattice_act[i][0][0] - lattice_dea[i][0][0]
        err_v = lattice_act[i][1][0] - lattice_dea[i][1][0]
        err_c = lattice_act[i][2][0] - lattice_dea[i][2][0]
        err_n = lattice_act[i][3][0] - lattice_dea[i][3][0]
        err_a = lattice_act[i][4][0] - lattice_dea[i][4][0]
        err_b = lattice_act[i][5][0] - lattice_dea[i][5][0]
        mesh_err = mesh_err.append({"p max": np.max(err_p),
                                    "p min": np.min(err_p),
                                    "v max": np.max(err_v),
                                    "v min": np.min(err_v),
                                    "c max": np.max(err_c),
                                    "c min": np.min(err_c),
                                    "n max": np.max(err_n),
                                    "n min": np.min(err_n),
                                    "a max": np.max(err_a),
                                    "a min": np.min(err_a),
                                    "b max": np.max(err_b),
                                    "b min": np.min(err_b)},ignore_index=True)

    # pd.options.display.float_format = '{:,.3e}'.format
    # print("Absolute mesh error in the x direction")
    # print(mesh_err)
    # print("\n")

    for i in range(len(lattice_act)):
        err_p = lattice_act[i][0][1] - lattice_dea[i][0][1]
        err_v = lattice_act[i][1][1] - lattice_dea[i][1][1]
        err_c = lattice_act[i][2][1] - lattice_dea[i][2][1]
        err_n = lattice_act[i][3][1] - lattice_dea[i][3][1]
        err_a = lattice_act[i][4][1] - lattice_dea[i][4][1]
        err_b = lattice_act[i][5][1] - lattice_dea[i][5][1]
        mesh_err = mesh_err.append({"p max": np.max(err_p),
                                    "p min": np.min(err_p),
                                    "v max": np.max(err_v),
                                    "v min": np.min(err_v),
                                    "c max": np.max(err_c),
                                    "c min": np.min(err_c),
                                    "n max": np.max(err_n),
                                    "n min": np.min(err_n),
                                    "a max": np.max(err_a),
                                    "a min": np.min(err_a),
                                    "b max": np.max(err_b),
                                    "b min": np.min(err_b)},ignore_index=True)

    pd.options.display.float_format = '{:,.3e}'.format
    # print("Absolute mesh error in the y direction")
    # print(mesh_err)
    # print("\n")


if __name__ == "__main__":
    main()
