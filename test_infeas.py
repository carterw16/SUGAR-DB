import os
import sys
from pathlib import Path

# from classes.Infeasibility import Infeasibility
from classes.Methods.InfeasibilityAnalysis import InfeasibilityAnalysis
from lib.parser import parser

import pandas as pd
import numpy as np
import math

class TestInfeas:
    _solved = False
    _results = {
        "Solved": None,
        "Iters": None,
        "V": None
    }
    def adjust_angle(self, angle_value):
        if angle_value > 360:
            angle_value = angle_value - 360
        elif angle_value < -5:
            angle_value = angle_value + 360
        elif angle_value >= -5 and angle_value <= 5:
            angle_value = angle_value + 360
        
        return angle_value

    def calc_percent_error(self, old_val, new_val):
        if old_val != 0:
            err = (old_val - new_val)/old_val * 100
        else:
            err = 0
        
        return err
    
    def print_error_violations(self, err_type, comparison_value, benchmark_value, test_value, node_name):
        print('Node Name:', node_name)
        print(err_type, comparison_value)
        print('test value: ', test_value)
        print('benchmark value: ', benchmark_value)

    def run_infeas(self, case: str, settings: dict, features: dict) -> dict:        
        Path("output/").mkdir(parents=True, exist_ok=True)
        path_to_output = os.path.normpath("output/")

        # parse and store the case data
        casedata, node_key, node_index_ = parser(case, settings, features)

        # run infeasibility
        infeasibility = InfeasibilityAnalysis(case, casedata, features, settings, path_to_output, node_key, node_index_)

        # Return the infeas return code to verify the case was solved
        self._solved = infeasibility.run_infeasibility_analysis()

        if self._solved == 0:
            # TODO: update and return the results dict with the solved status, the number of iterations, and the solution
            # vector
            self._results["Solved"] = self._solved
            case_name = case.strip("gridlabd/")
            voltage_results_path = path_to_output + '/' + case_name +'/' + case_name + '_' + str(settings['Load Factor']) + '_voltages_SUGAR.csv'
            # skip first row because it does not contain the column labels
            test_results = pd.read_csv(voltage_results_path, delimiter = ',', skiprows = 1)
            return test_results
        
        else:
            raise RuntimeError("Infeasibility failed.")
            
    def parse_benchmark_results(self, benchmark_path: str, test_results: dict, percent_error_threshold: float, print_error_flag: bool) -> dict:
        # Read in results to pandas dataframe
        benchmark_results = pd.read_csv(benchmark_path, delimiter = ',')
        
        comparison_columns = ['Node','A_mag', 'B_mag', 'C_mag', 'A_ang', 'B_ang', 'C_ang']
        comparison_df = pd.DataFrame(columns = comparison_columns, index = range(0, len(test_results)))
        
        violations = 0
        violation_index = []
        
        for test_index in range(0, len(test_results)):
            benchmark_index = int(np.where(test_results['Bus'][test_index] == benchmark_results['node_name'])[0])
            
            test_mag_ang = [test_results['Magnitude1'][test_index], test_results['Magnitude2'][test_index], 
                            test_results['Magnitude3'][test_index],
                            self.adjust_angle(test_results['Angle1'][test_index]), 
                            self.adjust_angle(test_results['Angle2'][test_index]),
                            self.adjust_angle(test_results['Angle3'][test_index])]

            benchmark_mag_ang = [benchmark_results['voltA_mag'][benchmark_index], benchmark_results['voltB_mag'][benchmark_index],
                            benchmark_results['voltC_mag'][benchmark_index], 
                             self.adjust_angle(math.degrees(benchmark_results['voltA_angle'][benchmark_index])), 
                             self.adjust_angle(math.degrees(benchmark_results['voltB_angle'][benchmark_index])),
                             self.adjust_angle(math.degrees(benchmark_results['voltC_angle'][benchmark_index]))]
                
            # percent error - magnitude
            violation_index_labels = ['_phase_A', '_phase_B', '_phase_C', '_phase_A', '_phase_B', '_phase_C']
            print_error_labels = ['V_A_mag_error: ', 'V_B_mag_error: ', 'V_C_mag_error: ', 
                                  'V_A_ang_error: ', 'V_B_ang_error: ', 'V_C_ang_error: ']
            
            for index in range(0, len(comparison_columns)):
                if index == 0:
                    comparison_df[comparison_columns[index]][test_index] = test_results['Bus'][test_index]
                else:
                    comparison_df[comparison_columns[index]][test_index] = self.calc_percent_error(benchmark_mag_ang[index - 1], test_mag_ang[index - 1])
                    if comparison_df[comparison_columns[index]][test_index] > percent_error_threshold:
                        violations += 1
                        violation_index.append(comparison_df['Node'][test_index] + violation_index_labels[index - 1])
                        if print_error_flag:
                            self.print_error_violations(print_error_labels[index - 1], comparison_df[comparison_columns[index]][test_index], 
                                                benchmark_mag_ang[index - 1], test_mag_ang[index - 1], comparison_df['Node'][test_index])
                            print("\n")
        return violations

    def test_R1_12_47_3(self, case, benchmark_path, settings, features, percent_error_threshold, print_error_flag):
        case_name = case
        #benchmark_path = 'test_suite/R1-12.47-3/R1-12.47-3_voltages.csv'
        
        test_results = self.run_infeas(case_name, settings, features)

        # Parse benchmark results and return benchmark results
        violations = self.parse_benchmark_results(benchmark_path, test_results, percent_error_threshold, print_error_flag)
        return violations

    
    def test_accuracy(self, case_name, benchmark_path, error_threshold):
        print_error_flag = False
        # WARM START SETTINGS
        # If performing a warm start to a test, use the dictionary to pass settings
        warm_start_settings_dict = {}
        # REQUIRED - either True or False
        warm_start_settings_dict['initialize'] = False

        if warm_start_settings_dict['initialize'] == True:
            # REQUIRED: starting from a power flow (0), L1 norm (1), or L2 norm (2) solution?
            warm_start_settings_dict['initial solution type'] = 'L1'
            # REQUIRED: load factor from the previous solution
            warm_start_settings_dict['previous lf'] = 1
            # REQUIRED: objective for new test
            warm_start_settings_dict['new objective'] = 'L1'
            # REQUIRED: file path of the previous solution:
            warm_start_settings_dict['solution file path'] = 'Initial_Conditions/gridlabd/R1-12.47-3'

        # INFEASIBILITY ANALYSIS SETTINGS
        # If performing an infeasibility analysis, use the dictionary to pass analysis settings
        infeas_settings_dict = {}
        # REQUIRED - either True or False
        infeas_settings_dict['run infeas'] = True 
        if infeas_settings_dict['run infeas']:
            # OPTIONAL - options are: 'current', 'PQ', 'GB', 'Q', or 'B'; default is 'current'
            infeas_settings_dict['source type'] = 'PQ' 
            # OPTIONAL - options are: 'L1' or 'L2'; default is 'L2'
            infeas_settings_dict['obj'] = 'L2' 
            # OPTIONAL- options are: 'power' or 'current'; default is 'power'
            infeas_settings_dict['obj type'] = 'power'
            # OPTIONAL - either True or False; this setting determines if infeasibility sources are added at the Slack bus
            infeas_settings_dict['stamp slack bus'] = False
            # OPTIONAL - default value is 1e-6 
            infeas_settings_dict['comp slack tol'] = 1e-6 
            # OPTIONAL - either True or False; this setting determines if infeasibility values are printed to the terminal
            infeas_settings_dict['report infeasibility'] = True 
            # OPTIONAL - either True or False; this settings determines whether or not to use NetworkX to draw a map of infeasibility in the network
            # TODO: create maps for PQ and GB variables
            infeas_settings_dict['draw infeasibility maps'] = False
            # OPTIONAL - either True or False; this settings is used in mapping - it normalizes infeasibility current
            # TODO: create for PQ and GB variables
            infeas_settings_dict['normalize mapped infeas var'] = False
            # OPTIONAL - either True or False; this setting determines whether infeasibility sources are added at triplex nodes; default value is False
            infeas_settings_dict['triplex sources'] = True

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
            'Run Type': 'Infeas', #'Power Flow'
            'infeas settings': infeas_settings_dict,
            'Stamp Dual': True,
            'Warm Start': warm_start_settings_dict
        }
        
        features = {
            'IBDGs': {
            },
            'Tap Controls': {
                'Fixed': True
            },
            'Current Meas': False,
        }

        test_results = self.run_infeas(case_name, settings, features)

        # Parse benchmark results and return benchmark results
        violations = self.parse_benchmark_results(benchmark_path, test_results, error_threshold, print_error_flag)

        if violations != 0:
            print(case_name, ': Test failed. The number of violations was ', violations)
        else: 
            print(case_name,': Test passed')
       

if __name__ == "__main__":
    # Test Suite:
    # cases: 
    # 'gridlabd/R1-12.47-3_debug_test' 
    # 'gridlabd/R2-25.00-1_NR'
    # 'gridlabd/R3-12.47-3_NR_SUGAR'
    # 'gridlabd/R4-12.47-1_NR_SUGAR'
    # 'gridlabd/R5-12.47-3_NR_SUGAR' 'gridlad/13node_ieee_NR_SUGAR';

    # benchmark paths:
    # 'test_suite/R1-12.47-3/R1-12.47-3_voltages.csv'
    # 'test_suite/R2-25.00-1/R2-25.00-1_voltages.csv'
    # 'test_suite/R3-12.47-3/R3-12.47-3_voltages.csv'
    # 'test_suite/R4-12.47-1/R4-12.47-1_voltages.csv'
    # 'test_suite/R5-12.47-3/R5-12.47-3_voltages.csv'
    # 'test_suite/ieee13node/ieee13node_voltages.csv'

    case_name = sys.argv[1] 
    benchmark_path = sys.argv[2] 
    error_threshold = float(sys.argv[3])
    test = TestInfeas()
    test.test_accuracy(case_name, benchmark_path, error_threshold)

    # Note: I do not think there is an infeasibility source at the Slack node - need to modify to include. 
    # Save statement in dual_analysis is not for Slack nodes.
