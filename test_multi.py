import os
import sys
from pathlib import Path

from classes.Methods.Multi import Multiperiod
from lib.parser import parser

import pandas as pd
import numpy as np
import math
import lib.parser_post_json as parser_post_json

class TestMP:
    _solved = []
    _results = {
        "Solved": None,
        "Iters": None,
        "V": None
    }

    def runMP(self, case: str, post_process_file: str, settings: dict, features: dict):
        Path("output/").mkdir(parents=True, exist_ok=True)
        path_to_output = os.path.normpath("output/")
        # each period's output file path is "output/{testcase}/{testcase}_p{period}_{type}.csv"

        # parse and store the case data
        casedata, node_key, node_index_ = parser(case, settings["Stamp Dual"], settings["Obj norm"], features)

        # post-process
        casedata = parser_post_json.post_process_json(casedata, post_process_file)

        # run infeasibility
        mp = Multiperiod(case, casedata, features, settings, path_to_output, node_key, node_index_)

        # Add battery
        # mp.battery = casedata.battery


        # Return the return codes to verify the case was solved at each iteration

        for period in range(len(self._solved)):
            if sol == 0:
                case_name = case.replace("gridlabd/","")

                if period == 0:
                    print("case name", case_name)

                print(f"================Period {period}========================")

                voltage_results_path = path_to_output + '/' + case_name +'/' + case_name + '_p' + str(period) + '_voltages_SUGAR.csv'
                # skip first row because it does not contain the column labels
                test_results = pd.read_csv(voltage_results_path, delimiter = ',', skiprows = 1)
                print("Results:", test_results)
            else:
                raise RuntimeError("Infeasibility failed.")

if __name__ == "__main__":

    warm_start_settings = []
    # Do you want to initialize with a warm start?
    warm_start_settings.append(False)

    settings = {
        'Tolerance': 1E-5,
        'MP Tolerance': 1E-5,
        'Max Iters': 10000,
        'Save File': True,
        'Voltage Limiting': False,
        'Diode Limiting': False,
        'Load Factor': 1,
        'Homotopy': False,
        'G_homotopy': 400,
        'B_homotopy': 400,
        'Run Type': 'Infeas',
        'Stamp Dual': True,
        'Obj norm': 2,
        'Warm Start': warm_start_settings
    }

    features = {
        'IBDGs': {
        },
        'Tap Controls': {
            'Fixed': True
        },
        'Current Meas': False,
    }

    #'gridlabd/ieee_4node_stepdown_unbal_D-D' 'testcases/jsons/ieee4bat.json' 'test_suite/ieee4node/ieee4node_voltages.csv'
    case_name = sys.argv[1]
    post_process_file = sys.argv[2]
    benchmark_path = sys.argv[3]
    test = TestMP()
    test.runMP(case_name, post_process_file, settings, features)
    # Note: I do not think there is an infeasibility source at the Slack node - need to modify to include.
    # Save statement in dual_analysis is not for Slack nodes
