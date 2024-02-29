from contextlib import redirect_stdout
from itertools import count
from termcolor import colored
from scipy.sparse.linalg import spsolve
from types import MethodType
from scipy.sparse import csr_matrix
import numpy as np
from classes.Methods.InfeasibilityAnalysis import InfeasibilityAnalysis

import time

class Multiperiod():
    def __init__(self, testcase, casedata, features, settings, path_to_output, node_key, node_index_):
        self.features = features
        self.settings = settings
        self.case_name = testcase
        self.casedata = casedata
        self.node_key = node_key
        self.node_index_ = node_index_
        self.output = path_to_output

        # Convergence Settings
        self.MP_tolerance = settings['multi settings']["tolerance"]
        self.verbose = settings['multi settings']["verbose"]


        self.err_max = 1.1 * self.MP_tolerance
        self.periods = settings['multi settings']["periods"] # TODO: add as setting
        self.num_loads = 10 ## TODO: calculate number of loads/generators from gridlab file
        self.load_curve = settings['multi settings']['load curve']
        self.num_batteries = len(self.casedata.battery_sources)


    def runMP(self):
        mp_start_time = time.perf_counter()

        # Initialize Coupling Variables
        B = np.zeros((self.periods, self.num_batteries))
        L = np.zeros((self.periods, self.num_batteries))
        M_low = np.zeros((self.periods, self.num_batteries))
        M_up = np.zeros((self.periods, self.num_batteries))

        # Slack Bus Power - FOR DEBUG PURPOSES
        P = np.zeros(self.periods)

        # Initialize infeasibility object
        infeasibility = InfeasibilityAnalysis(self.case_name, self.casedata, self.features, self.settings, self.output, self.node_key, self.node_index_)
        infeasibility.multi_period = True

        # Outer Loop - Epochs
        #while err_max > MP_tolerance:
        #    iteration_count += 1
            # Inner loop -> Run SUGAR3 at each time step

        for t in range(0, self.periods):

            # Modify infeasibility for this period
            infeasibility.load_factor = self.load_curve[t]

            # Initialize with last iteration's solution
            if t > 0:
                infeasibility.Vinit = Vsol
                infeasibility.use_Vinit = True

            # RUN INFEASIBILITY #
            if self.verbose:
                    infeasibility.run_infeasibility_analysis()
            else:
                with open('infeas_out.txt', 'a') as f:
                    with redirect_stdout(f):
                        infeasibility.run_infeasibility_analysis()

            #P[t] = infeasibility.slack

            # Extract Battery coupling terms from solution vector
            for idx in range(0, self.num_batteries):
                bat = infeasibility.battery_sources[idx]

                # SOC constraints
                B[t][idx] = infeasibility.Vsol[bat.Bt]
                L[t][idx] = infeasibility.Vsol[bat.lambda_Bt]
                M_low[t][idx] = infeasibility.Vsol[bat.node_dual_ineq_Bt]
                M_up[t][idx] = infeasibility.Vsol[bat.node_dual_ineq_Bt_upper]
                print(f'Period: {t}, SOC: {B[t]}')


            # Update battery object SOC in infeasibility for next iteration
            for idx in range(0, self.num_batteries):
                infeasibility.battery_sources[idx].Bt_prev = B[t][idx]

            # extract vsol to use as initialization for next solution
            Vsol = infeasibility.Vsol

        # DEBUG #
        print("=======================\nBattery SOC by Period\n======================")
        for row in range(len(B)):
            print(f'Period: {row}, SOC: {B[row]}')


        #TODO properly extract slack bus power
        '''print("=======================\nSlack Bus Power by Period\n======================")
        for row in range(len(B)):
            print(f'Period: {row}, SOC: {B[row]}')'''


