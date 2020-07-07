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
# * Alessandro Gastaldi
# * Aaron Dettmann
# * Jean-Philippe Kuntzer

"""
PyTornado standard functions

Developed for Airinnova AB, Stockholm, Sweden.
"""

import logging
import commonlibs.logger as hlogger
from pytornado.__version__ import __version__
from pytornado.objects.vlm_struct import VLMData
import pytornado.aero.vlm as vlm
import pytornado.fileio as io
import pytornado.plot.makeplots as makeplots
# import pytornado.fileio.native.deformation as deform

### TODO delete once debug is done
import numpy as np
import pickle
###
logger = logging.getLogger(__name__)
__prog_name__ = 'pytornado'


class StdRunArgs:
    """
    Arguments used in 'standard_run'

    Attributes:
        :run: name of file to be loaded
        :verbose: boolean flag, true for verbose logger setting
        :debug: boolean flag, true for debug logger setting
        :quiet: boolean flag, true for quiet logger setting
    """

    def __init__(self, run=None, verbose=False, debug=False, quiet=False):
        self.run = run
        self.verbose = verbose
        self.debug = debug
        self.quiet = quiet


def get_settings(settings_filepath):
    """
    Read settings and return a settings instance
    """

    logger.info("Getting configuration file...")
    settings = io.native.settings.load(settings_filepath)
    return settings


def clean_project_dir(settings):
    """
    Remove old files in project directory

    Args:
        :settings: settings instance
    """

    logger.info("Removing old files...")
    settings.clean()


def standard_run(args):
    """
    Run a standard analysis. This function is the brain of the program. The
    function does the following:
        1) (Setup) Reads the setting file and Sets the logging level for all
            the functions and classes for the rest of the simulation

        2) (Setup aircraft model and flight state) self explanatory

        3) (Generate lattice) Allocates a variable and stores all the mesh
            data into it. The file has the following structure:
            * :lattice.p: panel corner points
            * :lattice.v: panel vortex filament endpoints
            * :lattice.c: panel collocation point
            * :lattice.n: panel normal vector
            * :lattice.a: panel surface area
        4) (Computes VLM) Generates the LHS and RHS of the VLM. LHS being the
            "downwash" matrix and RHS being "boundary" matrix. Solves the Ax=b
            system and then computes the aircrafts results. The last step is
            just summing up all the panels contribution in oder to have the
            aircraft parameters.

    Args:
        :args: arguments (see StdRunArgs())
    """
    # Sets logging level for the whole simulation by reading the "args"
    # variable and displays program version.
    settings = get_settings(settings_filepath=args.run)
    hlogger.init(settings.paths('f_log'), level=args)
    logger = logging.getLogger(__name__)
    logger.info(hlogger.decorate(f"{__prog_name__} {__version__}"))

    # ===== Setup aircraft model and flight state =====
    # Chooses where (CPACS of JSON) to read the aircraft information and
    # assing the values to the aircraft variable.
    if settings.aircraft_is_cpacs: aircraft = io.cpacs.aircraft.load(settings)
    else:                          aircraft = io.native.aircraft.load(settings)

    # Sets state for the simulation (Mach number, altitude, AoA, etc...)
    if settings.state_is_cpacs: state = io.cpacs.state.load(settings)
    else:                       state = io.native.state.load(settings)

    # TODO: load as part of aircraft definition
    # if settings.settings['deformation']:
    #     io.native.deformation.load(aircraft, settings)

    # ===== Generate lattice =====
    vlmdata = VLMData()
    vlm.set_autopanels(aircraft, settings)

    # ----- Iterate through the flight states -----
    for i, cur_state in enumerate(state.iter_states()):

        settings.paths.counter = i

        # TODO: Temporary workaround!
        settings.paths('d_results', make_dirs=True, is_dir=True)
        settings.paths('d_plots', make_dirs=True, is_dir=True)

        # TODO: Don't set refs here. Find better solution!
        cur_state.refs = aircraft.refs
        ##########################################################

        ##########################################################
        # TODO: Find better solution for pre_panelling() function
        make_new_subareas = True if i == 0 else False
        ##########################################################

        lattice = vlm.gen_lattice(aircraft, cur_state, settings, make_new_subareas)
        
        # ===== Mesh Deformation =====    
        logger.info(settings.settings["aircraft"])
        if  "Activated" in settings.settings["aircraft"]:
            # Deforms the mesh and uploads the deformed one into the code
            logger.info("===== Mesh deformation function activated =====")
            mesh_def = io.native.deformation.Mesh_Def(lattice)
            mesh_def.deformation(settings)
            lattice.p = mesh_def.f_p
            lattice.v = mesh_def.f_v
            lattice.c = mesh_def.f_c
            lattice.bound_leg_midpoints = mesh_def.f_b
            lattice.n = mesh_def.f_n
            lattice.a = mesh_def.f_a
        else:
            logger.info("===== Mesh deformation function deactivated =====")

        # ===== VLM =====
        vlm.calc_downwash(lattice, vlmdata)
        vlm.calc_boundary(lattice, cur_state, vlmdata)  # right-hand side terms
        vlm.solver(vlmdata)
        vlm.calc_results(lattice, cur_state, vlmdata)
        
        # TODO delete once the debugging phase is done
        # Saves the results for comparison with the debugger.py function
        path = str(settings.project_dir)
        if  "Activated" in settings.settings["aircraft"]:
            with open(path + '/lattice_defActivated.pkl', 'wb') as l:
                var = [lattice.p,
                       lattice.v,
                       lattice.c,
                       lattice.n,
                       lattice.a,
                       lattice.bound_leg_midpoints]
                pickle.dump(var, l)
            l.close()
            
            with open(path + '/data_defActivated.pkl', 'wb') as d:
                pickle.dump(vlmdata, d)
            d.close()
        
        else:
            with open(path + '/lattice_defDeactivated.pkl', 'wb') as l:
                var = [lattice.p,
                       lattice.v,
                       lattice.c,
                       lattice.n,
                       lattice.a,
                       lattice.bound_leg_midpoints]
                pickle.dump(var, l)
            l.close()
            logger.debug(lattice)
            
            with open(path + '/data_defDeactivated.pkl', 'wb') as d:
                pickle.dump(vlmdata, d)
            d.close()
        
        # ===== Create plots and result files =====
        io.native.results.save_all(settings, aircraft, cur_state, vlmdata)
        makeplots.make_all(settings, aircraft, cur_state, vlmdata, lattice)
        
        ################################################
        # TODO: Find better solution
        ################################################
        # Save AeroPerformance map results
        state.results['Fx'].append(vlmdata.forces['x'])
        state.results['Fy'].append(vlmdata.forces['y'])
        state.results['Fz'].append(vlmdata.forces['z'])
        state.results['FD'].append(vlmdata.forces['D'])
        state.results['FC'].append(vlmdata.forces['C'])
        state.results['FL'].append(vlmdata.forces['L'])
        state.results['Mx'].append(vlmdata.forces['l'])
        state.results['My'].append(vlmdata.forces['m'])
        state.results['Mz'].append(vlmdata.forces['n'])
        ####
        state.results['Cx'].append(vlmdata.coeffs['x'])
        state.results['Cy'].append(vlmdata.coeffs['y'])
        state.results['Cz'].append(vlmdata.coeffs['z'])
        state.results['CD'].append(vlmdata.coeffs['D'])
        state.results['CC'].append(vlmdata.coeffs['C'])
        state.results['CL'].append(vlmdata.coeffs['L'])
        state.results['Cl'].append(vlmdata.coeffs['l'])
        state.results['Cm'].append(vlmdata.coeffs['m'])
        state.results['Cn'].append(vlmdata.coeffs['n'])
        ################################################

    # ---------- Save aeroperformance map ----------
    if settings.aircraft_is_cpacs and settings.state_is_cpacs:
        io.cpacs.results.save_aeroperformance_map(state, settings)

    if settings.settings['save_results']['aeroperformance']:
        io.native.results.save_aeroperformance_map(state, settings)

    logger.info(f"{__prog_name__} {__version__} terminated")

    # ---------- Return data to caller ----------
    results = {
        "lattice": lattice,
        "vlmdata": vlmdata,
        "state": state,
        "settings": settings,
    }
    return results
