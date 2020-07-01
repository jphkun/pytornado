#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
# import shutil
from pathlib import Path
import pytornado.database.tools as dbtools

def test_basic_usage():
    """
    Test all the cases
    """
    # if true : Activates deformation and selects the undefored file in the
    #           aircraft file.
    # if false: Deactivates the deformation function and selects the deformed
    #           file. This mode is selected as reference to test the 
    #           deformation function.
    deform = True
    
    # Working directory path
    absolute_path = "/home/cfse2/Documents/pytornado/tests/integration/mesh_interfacing/wkdir/"
    
    # list of all the tests
    paths = [#"1_pytornado_flat",
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
             "22_Boxwing",
             # "11_OptiMale"
             ]
    
    tempActiv = [#"1_flat_funcActivated.json",
                 #"2_dih_funcActivated.json",
                 #"3_anh_funcActivated.json",
                 #"4_flat_funcActivated.json",
                 #"5_dih_funcActivated.json",
                 #"6_anh_funcActivated.json",
                 #"7_flat_funcActivated.json",
                 #"8_dih_funcActivated.json",
                 #"9_anh_funcActivated.json",
                 #"10_flat_funcActivated.json",
                 #"11_dih_funcActivated.json",
                 #"12_anh_funcActivated.json",
                 #"D150_AGILE_Hangar_funcActivated.json",
                 "Boxwing_AGILE_Hangar_funActivated_v3.1.json",
                ]
    
    tempDeactiv = [#"1_flat_funcDeactivated.json",
                   #"2_dih_funcDeactivated.json",
                   #"3_anh_funcDeactivated.json",
                   #"4_flat_funcDeactivated.json",
                   #"5_dih_funcDeactivated.json",
                   #"6_anh_funcDeactivated.json",
                   #"7_flat_funcDeactivated.json",
                   #"8_dih_funcDeactivated.json",
                   #"9_anh_funcDeactivated.json",
                   #"10_flat_funcDeactivated.json",
                   #"11_dih_funcDeactivated.json",
                   #"12_anh_funcDeactivated.json",
                   #"D150_AGILE_Hangar_funcDeactivated.json",
                   "Boxwing_AGILE_Hangar_funDeactivated_v3.1.json",
                   ]
    
    print(paths)
    for i in range(len(paths)):
        # Paths
        project_dir = Path(absolute_path + paths[i])
        
        if deform:
            # deformation function will deform the undeformed airplane mesh
            settings_file = Path(os.path.join(project_dir, 'settings', tempActiv[i]))
        else:
            # since def. function is not actived airplane mesh needs to be deformed by the file
            settings_file = Path(os.path.join(project_dir, 'settings', tempDeactiv[i]))
        
        print(settings_file)
        
        results_dir = Path(os.path.join(project_dir, '_results'))
        plot_dir = Path(os.path.join(project_dir, '_plots'))
        
        # ------ Make sure it runs -----
        with open(settings_file, "r") as fp:
            settings = json.load(fp)
        settings['plot']['results']['show'] = False
        settings['plot']['results']['save'] = False
        with open(settings_file, "w") as fp:
            json.dump(settings, fp)
        os.system(f"pytornado -d --clean --run {settings_file}")
        assert plot_dir.is_dir()
        assert results_dir.is_dir()

test_basic_usage()