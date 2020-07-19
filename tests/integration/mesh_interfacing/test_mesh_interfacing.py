#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
import json
import time
import glob
from pathlib import Path
logging.basicConfig(level=logging.DEBUG)
__prog_name__ = 'test_mesh_interfacing'
logger = logging.getLogger(__prog_name__+"."+__name__)


def test_basic_usage():
    """
    Test all the cases
    """
    # Bechmark time start
    t0 = time.time()

    # if true : Activates deformation and selects the undefored file in the
    #           aircraft file.
    # if false: Deactivates the deformation function and selects the deformed
    #           file. This mode is selected as reference to test the
    #           deformation function.
    deform = False
    # Choices are:
    #   shape1
    #   shape2
    #   csv_def_load
    deformation_method = "shape_2"
    start = 0
    end = 1
    wkdir = "01_validation"
    # Working directory path
    current_directory = os.path.dirname(os.path.abspath(__file__))
    absolute_path = current_directory + "/" + wkdir + "/"
    os.chdir(absolute_path)

    # Cleans up the unnecessary data
    paths = next(os.walk('.'))[1]
    r = [".","__"]
    new = []
    for path in paths:
        if r[0] not in path and r[1] not in path:
            new.append(path)
    paths = sorted(new)

    paths_act = glob.glob(absolute_path + "**/settings/*act*.json")
    paths_dea = glob.glob(absolute_path + "**/settings/*dea*.json")

    for i in range(len(paths[start:end])):
        # Paths
        project_dir = Path(absolute_path + paths[start+i])

        if deform:
            # deformation function will deform the undeformed airplane mesh.
            settings_file = paths_act[i]
            # Path(os.path.join(project_dir,'settings',tempActiv[start+i]))
        else:
            # since def. function is not actived airplane mesh needs to be
            # deformed by the file.
            settings_file = paths_dea[i]
            # Path(os.path.join(project_dir,'settings',tempDeactiv[start+i]))

        # print(settings_file)

        results_dir = Path(os.path.join(project_dir, '_results'))
        plot_dir = Path(os.path.join(project_dir, '_plots'))

        # ------ Make sure it runs -----
        with open(settings_file, "r") as fp:
            settings = json.load(fp)
        settings['plot']['results']['show'] = False
        settings['plot']['results']['save'] = True
        settings['deformation'] = deform
        settings['deformation_method'] = deformation_method
        with open(settings_file, "w") as fp:
            json.dump(settings, fp)
        #   -v, --verbose; -d, --debug; -q, --quiet
        os.system(f"pytornado -v --clean --run {settings_file}")
        # assert plot_dir.is_dir()
        # assert results_dir.is_dir()

    # End of bechmark
    t1 = time.time()
    total_execution_time = t1 - t0
    logger.info("Execution time: " + str(total_execution_time))


if __name__ == '__main__':
    test_basic_usage()