from __future__ import division, print_function
# Python Library Imports
import argparse
import os
import sys

parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.append(parent_dir + "/classes")
sys.path.append(parent_dir + "/lib")
sys.path.append(parent_dir)

from lib.parser import parser

import pyximport

pyximport.install(language_level=3)


#import GlobalVars


# CASE NAME
#case = 'gridlabd/R1-12.47-3_debug_test'
# case = 'gridlabd/R4-12.47-1'
case = parent_dir + '/testcases/gridlabd/GC-12.47-1'
# case = 'opendss/d23802'
# case = 'gridlabd/twobus_case_bal_Y'

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

FEATURES = {
    'IBDGs': {
    },
    'Tap Controls': {
        'Fixed': True
    },
    'Current Meas': False,
}

def main(TESTCASE, SETTINGS, FEATURES):
    casedata, node_key, node_index_ = parser(TESTCASE, SETTINGS, FEATURES)

    # casedata namespace contains the following objects:
    # node: electrical node for each devices (same as a bus) 
    # load: electrical load - to, from, constant power, nominal voltage, name is {}_wind = wind plant, name is {}_PV = PV plant
    # ohline: overhead lines - to, from
    # ugline: underground liens - to, from
    # slack: slack generator 
    # xmfr: transformer - to, from, primary voltage, secondary voltage, power rating
    print(casedata.xfmr[0].__dict__)

main(case, SETTINGS, FEATURES)
