from __future__ import division, print_function
from lib.parser import parser

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

import classes.GlobalVars


# CASE NAME
#case = 'gridlabd/R1-12.47-3_debug_test'
# case = 'gridlabd/R4-12.47-1'
# case = 'gridlabd/13node_ieee_NR_SUGAR'
# case = 'opendss/d23802'
case = 'gridlabd/twobus_case_bal_Y'

# Create settings for parser
infeas_settings_dict = {}
infeas_settings_dict['source type'] = 'current' # Initial values sent in lib/intialize.py
infeas_settings_dict['obj'] = 'L2'
infeas_settings_dict['obj type'] = 'power'
infeas_settings_dict['obj scalar'] = 1e-6

# VOLTAGE UNBALANCE SETTINGS
voltage_unbalance_settings_dict = {}
# REQUIRED - either True or False
voltage_unbalance_settings_dict['enforce voltage unbalance'] = False

# VOLTAGE BOUND SETTINGS
voltage_bound_settings_dict = {}
# REQUIRED - either True or False
voltage_bound_settings_dict['enforce voltage bounds'] = False

# WARM START SETTINGS
# If performing a warm start to a test, use the dictionary to pass settings
warm_start_settings_dict = {}
# REQUIRED - either True or False
warm_start_settings_dict['initialize'] = False

# feature: battery settings
infeas_settings_dict['battery_node_list'] = [
                                        {"ID": "B1", "node":"l4", "P_max":300000, "P_min":0,
                                         "Mch": 0.00000005, "Md": 0.00000009, "type":"P", "Bt_prev":0.1, "C_ch":1, "C_d":-0.5, "single_phase":"A"}
]

SETTINGS = {
    'infeas settings': infeas_settings_dict,
    'Stamp Dual': True,
    'voltage unbalance settings': voltage_unbalance_settings_dict,
    'Warm Start': warm_start_settings_dict,
    'voltage bound settings': voltage_bound_settings_dict,
        }


def main(TESTCASE, SETTINGS):
    casedata, node_key, node_index_ = parser(TESTCASE, SETTINGS)

    # casedata namespace contains the following objects:
    # node: electrical node for each devices (same as a bus) 
    # load: electrical load
    # ohline: overhead lines
    # slack: slack generator 
    breakpoint()

main(case, SETTINGS)
