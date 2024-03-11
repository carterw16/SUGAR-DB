from __future__ import division, print_function

import pyximport

pyximport.install(language_level=3)

# Python Library Imports
import argparse
import os
from os import getcwd
from sys import path

path.append(getcwd() + "/classes")
path.append(getcwd() + "/lib")
path.append(getcwd())

import logging
from pathlib import Path

from colorama import init
from termcolor import colored

from classes.Methods.PowerFlow import PowerFlow
from classes.Methods.InfeasibilityAnalysis import InfeasibilityAnalysis
from classes.Methods.Multi import Multiperiod

from lib.argparse_actions import ReadJSON

from lib.parser import parser
from lib.initialize_voltage_bounding import initialize_voltage_bounding

def main(TESTCASE, SETTINGS=None, FEATURES=None, path_to_output='output/'):

    # ==========================================================================================
    #     ----------------------- SET UP LOGGER -----------------------------------
    # ==========================================================================================

    Path('log/').mkdir(parents=True, exist_ok=True)
    casename = None
    if 'gridlabd' in TESTCASE:
        casename = TESTCASE.replace('gridlabd/', '')
    elif 'opendss' in TESTCASE:
        casename = TESTCASE.replace('opendss/', '')
    elif 'cyme' in TESTCASE:
        casename = TESTCASE.replace('cyme/', '')

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m-%d %H:%M',
        filename='log/' + casename + '.log',
        filemode='w')

    # create console handler and set level to debug
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # set console formatter
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)

    # add console handler to the root logger
    logging.getLogger('').addHandler(console)

    # ==========================================================================================
    #     -----------------------SOLVER OPTIONS AND SETTINGS -----------------------------------
    # ==========================================================================================
    path_to_output = os.path.normpath(path_to_output)

    # use Colorama to make Termcolor work on Windows too
    init()
    """----------------------Read and Parse Data from Input file-------------------------"""
    if SETTINGS["Run Type"] == "Power Flow":
        casedata, node_key, node_index_ = parser(TESTCASE, SETTINGS, FEATURES)

    elif SETTINGS["Run Type"] == "Infeas":
        casedata, node_key, node_index_ = parser(TESTCASE, SETTINGS, FEATURES)
        if SETTINGS['voltage bound settings']["enforce voltage bounds"]:
            voltage_bound_settings = SETTINGS['voltage bound settings']
            node_index_ = initialize_voltage_bounding(casedata, node_key, node_index_, voltage_bound_settings)

    elif SETTINGS["Run Type"] == "Multi": #TODO: load forecasted inputs
        casedata, node_key, node_index_ = parser(TESTCASE, SETTINGS, FEATURES)
        if SETTINGS['voltage bound settings']["enforce voltage bounds"]:
            voltage_bound_settings = SETTINGS['voltage bound settings']
            node_index_ = initialize_voltage_bounding(casedata, node_key, node_index_, voltage_bound_settings)
    logging.info(
        colored(
            '======================= Data Parse and Load Complete ==============================',
            'white'))
    breakpoint()
    if SETTINGS["Run Type"] == "Multi":
        multi = Multiperiod(TESTCASE, casedata, FEATURES, SETTINGS, path_to_output, node_key, node_index_)
        multi.runMP()
    elif SETTINGS["Run Type"] == "Power Flow":
        powerflow = PowerFlow(TESTCASE, casedata, FEATURES, SETTINGS, path_to_output, node_key, node_index_)
        powerflow.run_pf()

    elif SETTINGS["Run Type"] == "Infeas":
        infeasibility = InfeasibilityAnalysis(TESTCASE, casedata, FEATURES, SETTINGS, path_to_output, node_key, node_index_)
        infeasibility.run_infeasibility_analysis()

    """---------------------- End -------------------------"""



if __name__ == "__main__":
    features = {
        'IBDGs': {
        },
        'Tap Controls': {
            'Fixed': True
        },
        'Current Meas': False
    }

    settings = {
        'Tolerance': 1E-5,
        'Max Iters': 10000,
        'Save File': True,
        'Voltage Limiting':  False,
        'Diode Limiting': False,
        'Load Factor': 1,
        'Homotopy': False,
        'G_homotopy': 400,
        'B_homotopy': 400,
        'Run Type': 'Power Flow',
        'infeas settings': None,
        'Stamp Dual': False,
        'Warm Start': None
    }

    cli_parser = argparse.ArgumentParser(description='Process SUGAR3 inputs.')
    cli_parser.add_argument('--case',
                            help='the name of the test case',
                            type=str)
    cli_parser.add_argument('--settings',
                            help='the optional settings for SUGAR',
                            default=settings,
                            action=ReadJSON,
                            type=str)
    cli_parser.add_argument('--features',
                            help='the optional features for SUGAR',
                            default=features,
                            action=ReadJSON,
                            type=str)
    cli_parser.add_argument('--out',
                            help='the output path',
                            default='output/',
                            type=str)

    args = cli_parser.parse_args()
    main(args.case, args.settings, args.features, args.out)
