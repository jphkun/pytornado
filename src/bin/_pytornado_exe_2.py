#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------
# Copyright 2020 Airinnova AB and the PyTornado authors
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

# Author: Aaron Dettmann

"""
Command line interface
"""

import argparse
from pathlib import Path
import sys

from commonlibs.logger.logger import truncate_filepath

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/cfse2/Documents/pytornado/src/lib')

import pytornado
from pytornado.__version__ import __version__
import pytornado.database.tools as dbtools
import pytornado.stdfun.run as stdrun
import pytornado.stdfun.setup as project_setup


def main():
    """
    Command line interface
    """
    HELP_RUN = f"Settings file (entry point for {stdrun.__prog_name__})"
    HELP_CPACS2JSON = "Convert a CPACS file to native JSON format"
    HELP_MAKE_EXAMPLE = "Generate a minimal working example"
    HELP_LIST_EXAMPLE_FROM_DB = "List example aircraft in the database"
    HELP_MAKE_EXAMPLE_FROM_DB = "Pick an aircraft from the database and \
                                    generate a project directory"
    
    # Default value for debut purposes. There is a trial folder set up to make
    # a chosen simulation. This permits to run a test the program behaviour.
    def_dir = "/home/cfse2/Documents/1_PDM/FE_test/2_Pytornado_test/\
1_pytornado_flat/settings/template.json"
#     def_dir = "/home/cfse2/Documents/1_PDM/FE_test/2_Pytornado_test/\
# 2_pytornado_dihedral/settings/template.json"
    # def_c2j =
    # def_mex =
    # def_lef =
    # def_mef = 
    
    parser = argparse.ArgumentParser(prog=f'{stdrun.__prog_name__} {__version__}')

    group = parser.add_mutually_exclusive_group()
    # Generates a C++ equivalent to BOOST. The idea is to have a help menu when
    # running in the CLI the program and have a set of choices. Some default
    # choices were added to help debug the program.
    group.add_argument(
        '-r', '--run',
        # TODO suppress when debug phase is finished
        nargs='?',
        # default=def_dir,
        metavar='<Settings file>',
        type=str, 
        help=HELP_RUN)
    group.add_argument(
        '--cpacs2json',
        metavar='<CPACS file>', 
        type=str,help=HELP_CPACS2JSON)
    group.add_argument(
        '--make-example',
        action='store_true',
        help=HELP_MAKE_EXAMPLE)
    group.add_argument(
        '--list-example-from-db', '--ldb', 
        action='store_true',
        help=HELP_LIST_EXAMPLE_FROM_DB)
    group.add_argument(
        '--make-example-from-db', 
        '--mdb',metavar='<Aircraft name>', 
        type=str,help=HELP_MAKE_EXAMPLE_FROM_DB)

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '-v', '--verbose', 
        action='store_true')
    group.add_argument(
        '-d', '--debug', 
        action='store_true')
    group.add_argument(
        '-q', '--quiet', 
        action='store_true')

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-c", "--clean", 
        help="remove old project files",
        action="store_true")
    group.add_argument(
        "--clean-only",
        help="clean and exit",
        action="store_true")
    
    args = parser.parse_args()
    
    ## TODO suppress when debug phase is finished
    # All the entry arguments are assigned here. This permits to have a clear 
    # view of what is happening during the comprhention/debug phase.
    args.clean = True
    args.clean_only = False
    args.cpacs2json = False
    args.debug = True
    args.list_example_from_db = False
    args.make_example = False
    args.make_example_from_db = False
    args.quiet = False
    args.verbose = False
    ##
    
    if args.clean or args.clean_only:
        # TODO: more general way (maybe --run must not be provided)
        if args.run is None:
            raise RuntimeError("Settings file must be specified with '--run'")
        stdrun.clean_project_dir(stdrun.get_settings(args.run))
        if args.clean_only:
            return
    
    if args.run:
        stdrun.standard_run(args)
    elif args.cpacs2json:
        project_setup.cpacs2pytornado(args.cpacs2json)
    elif args.make_example:
        project_setup.setup_wkdir()
    elif args.make_example_from_db is not None:
        aircraft_name = Path(args.make_example_from_db).stem
        print(f"Loading '{aircraft_name}' from database...")

        try:
            filepath_aircraft = dbtools.get_aircraft_file_path(aircraft_name)
        except FileNotFoundError:
            dbtools.print_available_aircraft()
            print(f"\nAircraft '{aircraft_name}' file not found. Exit...", 
                  file=sys.stderr)
            sys.exit(1)

        print(f"Aircraft found at '{truncate_filepath(filepath_aircraft)}'...")
        project_setup.setup_wkdir(filepath_aircraft)
    elif args.list_example_from_db:
        dbtools.print_available_aircraft()
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
