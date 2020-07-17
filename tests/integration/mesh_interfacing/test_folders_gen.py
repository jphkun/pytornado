#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 15:33:08 2020

@author: cfse2
"""

import os
import errno
import json
import glob
import shutil

def main():
    # Genreates validation folder containing all the cases
    directory = "01_validation"
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    
    # Reads all the CPACS file in the current folder and get the names
    
    aircrafts = glob.glob("*.xml")
    aircrafts = sorted(aircrafts)
    
    # Opens-up the settings templates a modify it accordingly
    settings_file = "Tuto_example/settings/Optimale.json"
    with open(settings_file, "r") as fp:
        settings1 = json.load(fp)
    with open(settings_file, "r") as fp:
        settings2 = json.load(fp)
    
    for i in range(len(aircrafts)):
        # generates the validation case
        path1 = directory + "/0"+ str(i+1) +"_case"
        try:
            os.makedirs(path1)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # Generates aircraft folder
        path2 = path1 + "/aircraft"
        try:
            os.makedirs(path2)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
                # Generates airfoils folder
        path3 = path1 + "/airfoils"
        try:
            os.makedirs(path3)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # Generates aircraft deformation folder
        path4 = path1 + "/deformation"
        try:
            os.makedirs(path4)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        # Generates aircraft settings folder
        path5 = path1 + "/settings"
        try:
            os.makedirs(path5)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise 
        # Generates aircraft state folder
        path6 = path1 + "/state"
        try:
            os.makedirs(path6)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        # Deactivated deformation function
        settings_dea = settings1
        settings_dea['deformation'] = "false"
        settings_dea['deformation_method'] = "false"
        settings_dea["aircraft"] = aircrafts[i][:-4]+"_dea.xml"
        
        # Activated deformation function
        settings_act = settings2
        settings_act['deformation'] = "true"
        # Shape function options are:
        #   shape_1
        #   shape_2
        #   load_from_csv
        settings_act['deformation_method'] = "shape_2"
        settings_act["aircraft"] = aircrafts[i][:-4]+"_act.xml"
        
        # Saves file the correct directory
        with open(path5 + "/settings_dea.json", "w") as fp1:
            json.dump(settings_dea, fp1)    
        with open(path5 + "/settings_act.json", "w") as fp2:
            json.dump(settings_act, fp2)
        
        # Copies state template
        template_path = "Tuto_example/state/template.json"
        shutil.copyfile(template_path, path1 + '/state/template.json')
        
        # Copies CPACS
        # CPACS_path = "Tuto_example/state/template.json"
        shutil.copyfile(aircrafts[0], 
                        path1 + '/aircraft/0' + str(i+1) + aircrafts[0][2:-4] + "_dea.xml")
        shutil.copyfile(aircrafts[i], 
                        path1 + '/aircraft/' + aircrafts[i][:-4] + "_act.xml")


if __name__ == "__main__":
    main()

