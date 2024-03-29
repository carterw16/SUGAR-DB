from contextlib import redirect_stdout
from itertools import count
from termcolor import colored
from scipy.sparse.linalg import spsolve
from types import MethodType
from scipy.sparse import csr_matrix
import numpy as np
from classes.Methods.InfeasibilityAnalysis import InfeasibilityAnalysis
import ipdb
import pickle
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
        self.PV_capacity = settings['multi settings']["PV capacity"]
        self.wind_capacity = settings['multi settings']["wind capacity"]
        self.max_epochs = settings['multi settings']["max epochs"]


        self.err_max = 1.1 * self.MP_tolerance
        self.periods = settings['multi settings']["periods"]
        self.num_loads = 10 ## TODO: calculate number of loads/generators from gridlab file
        self.load_curve = settings['multi settings']['load curve']
        self.wind_capacity = settings['multi settings']['wind capacity']
        self.PV_capacity = settings['multi settings']['PV capacity']
        self.num_batteries = len(self.casedata.battery_sources)


    def runMP(self):
        mp_start_time = time.perf_counter()


        # Multi-Epoch Results Dictionary - for a single battery
        self.results_dict = {'B': {}, 'B_res': {}, 'P_ch': {}, 'P_ch_res': {}, 
                'P_d': {}, 'P_d_res': {}, 'L': {}, 'L_res': {}, 
                'M_low': {}, 'M_low_res': {}, 'M_up': {}, 'M_up_res': {},
                'P_g': {'slack': []}, 'S_line': {}, 'V': {}, 'total_load': [],
                'total_gen': []}

        # Initialize infeasibility object
        infeasibility = InfeasibilityAnalysis(self.case_name, self.casedata, self.features, self.settings, self.output, self.node_key, self.node_index_)
        infeasibility.multi_period = True

        # Initialize Coupling Variable array for DDP
        L_DDP = np.zeros((self.periods, self.num_batteries, 3))
        # 
        '''L[1][0] = [0.5,0.5,0.5]
        L[2][0] = [0.5,0.5,0.5]
        L[3][0] = [0.5,0.5,0.5]'''

        # Store the solution for each period, compare with previous epoch
        sizeY = int(repr(infeasibility.node_index_)[6:-1])
        self.Vsol_mp_old = np.zeros((self.periods, sizeY))
        self.Vsol_mp = np.zeros((self.periods, sizeY))

        # =====================================================================
        #                       Outer Loop - Epochs
        # =====================================================================
        epoch = 0
        while self.err_max > self.MP_tolerance and epoch < self.max_epochs:


            # =================================================================
            #              Foward Pass - Run SUGAR3 at each time step
            # =================================================================
            for t in range(0, self.periods):

                # Modify infeasibility for this period
                for load in infeasibility.load:
                    if 'wind' in load.name:
                        load._cP = [self.wind_capacity[t]*load.cP_A.real]*3
                    if 'PV' in load.name:
                        load._cP = [self.PV_capacity[t]*load.cP_A.real]*3

                # Initialize with last iteration's solution
                if t > 0:
                    infeasibility.Vinit = Vsol
                    infeasibility.use_Vinit = True
                
                # =================================================================
                #                    BACKWARD PASS - update Lb_next
                # =================================================================
                for idx in range(0, self.num_batteries):
                    bat = infeasibility.battery_sources[idx]
                    
                    # Update lambda_t+1 in each battery for each phase
                    if t < self.periods - 1:
                        bat.Lb_next = L_DDP[t+1][idx]

                # RUN INFEASIBILITY #
                if self.verbose:
                        infeasibility.run_infeasibility_analysis()
                else:
                    with open('infeas_out.txt', 'a') as f:
                        with redirect_stdout(f):
                            infeasibility.run_infeasibility_analysis()

                # EXTRACT RESULTS
                self.extract_results(infeasibility, epoch, t)
                Vsol = infeasibility.Vsol
                res = infeasibility.res_eqn

                if t < self.periods - 1:
                    print(f'Lb_next: {L_DDP[t+1][0]}')
  
                # FOWARD PASS UPDATE
                for idx in range(0, self.num_batteries):
                    # Update Bt_prev
                    infeasibility.battery_sources[idx].Bt_prev = Vsol[infeasibility.battery_sources[idx].Bt_nodes]

                    # Update Lb_next
                    L = Vsol[[bat.lambda_BtA,bat.lambda_BtB,bat.lambda_BtC]].ravel()
                    L_DDP[t][idx] = L
               
                # extract vsol to use as initialization for next solution
                self.Vsol_mp[t] = Vsol.reshape(-1)

                # extract 1st period vsol to use for next epoch
                # if t == 0:
                #    Vsol_period1 = infeasibility.Vsol

 
            if epoch > 0:
                # Calculate Max Error over all periods
                err = self.Vsol_mp.ravel() - self.Vsol_mp_old.ravel()
                self.err_max = np.amax(abs(err))

                # Figure out what node has the error
                argmax_ind = np.argmax(abs(err)) % sizeY
                err_max_node = [
                    ele.name for ele in infeasibility.node if argmax_ind in ele.node_set
                ]
                if err_max_node:
                    print(
                        colored(
                            'Maximum error from epoch %s is %f at node %s' %
                            (epoch, self.err_max, err_max_node[0]), 'red'))
                else:
                    print(
                        colored('Maximum error from iteration %s is %f' % (epoch, self.err_max),
                                'red'))

            # Increment Epoch
            epoch += 1

            # Store Full Multi-Period Solution
            np.copyto(self.Vsol_mp_old, self.Vsol_mp)


            # Re-Initialize Battery
            for idx in range(0, self.num_batteries):
                Bt_INIT = infeasibility.infeas_settings['battery_node_list'][idx]['Bt_prev']
                infeasibility.battery_sources[idx].Bt_prev = [Bt_INIT, Bt_INIT, Bt_INIT]

            # DEBUG #
            #self.debug_print(B)
            print("="*30)
            print(f"EPOCH {epoch}")
            print("="*30)

        # After all epochs, store results dict as pickle for later plotting
        self.write_results()

    def write_results(self): 
        """
        writes time-variant battery variables and entire Multi object to 2 pickles

        saved battery variables and their residuals:
            - Bt SOC lagrange equation: Lb, mu_up, mu_low
            - Lb lagrange equation: Bt, P_d, P_ch
        """
        file_path = "multi_outputs.pkl"
        with open(file_path, 'wb') as pickle_file:
            pickle.dump(self.results_dict, pickle_file)

    def extract_results(self, infeasibility, epoch, t):
        """ Extracts the following grid state variables from a SUGAR solution
        and stores in the Multi period results dictionary
            - node: electrical bus - V = voltage (complex phasor in volts)
            - load: constant power PQ loads (W) - includes generation of PV and wind
              as negative PQ loads
            - lines: S_line = power flow (VA)
            - slack: slack generator - P_g = generation power (W)
            - battery has 3 phases, following variables in length 3 list [PhaseA, PhaseB, PhaseC]: 
                B = state of charge (% of total)  
                P_ch = charge power, 
                P_d = discharge power (W)
            - total load (W)
            - total generation (W)

        Also extracts saved battery variables and their residuals:
            - Bt SOC lagrange equation: Lb, M_up, M_low
            - Lb lagrange equation: Bt, P_d, P_ch
        """

        Vsol = infeasibility.Vsol
        res = infeasibility.res_eqn

        # INITIALIZE DICTIONARIES
        if epoch == 0:
            # initialize line powers
            for ele in infeasibility.curr_meas:
                if ele.name not in self.results_dict['S_line']:
                    self.results_dict['S_line'][ele.name] = []

            # initialize loads
            for load in infeasibility.load:
                if ('wind' in load.name) or ('PV' in load.name):
                    if load.name not in self.results_dict['P_g']:
                        self.results_dict['P_g'][load.name] = []

            # initialize battery variables
            for var in ['B', 'B_res', 'P_ch', 'P_ch_res',
                'P_d', 'P_d_res', 'L', 'L_res', 'M_up', 'M_up_res',
                'M_low','M_low_res']: 
                for bat in infeasibility.battery_sources:
                    if bat.ID not in self.results_dict[var]:
                        self.results_dict[var][bat.ID] = []

        # INITIALIZE EPOCH ARRAY
        if t == 0:
            for ele in infeasibility.curr_meas: # initialize line powers
                self.results_dict['S_line'][ele.name].append(np.zeros((self.periods,3), dtype=np.complex64))

            for var in ['B', 'B_res', 'P_ch', 'P_ch_res',
                'P_d', 'P_d_res', 'L', 'L_res', 'M_up', 'M_up_res',
                'M_low','M_low_res']: 
                for bat in infeasibility.battery_sources:
                    self.results_dict[var][bat.ID].append(np.zeros((self.periods,3)))

            for load in infeasibility.load:
                if ('wind' in load.name) or ('PV' in load.name):
                    self.results_dict['P_g'][load.name].append(np.zeros((self.periods,3)))

            self.results_dict['total_load'].append(np.zeros((self.periods)))
            self.results_dict['total_gen'].append(np.zeros((self.periods)))
            self.results_dict['P_g']['slack'].append(np.zeros((self.periods, 3), dtype=np.complex64))


        # EXTRACT BATTERY VARIABLES
        for idx in range(0, self.num_batteries):
            bat = infeasibility.battery_sources[idx]

             
            # SOC constraints
            
            B = Vsol[[bat.BtA,bat.BtB,bat.BtC]].ravel()
            B_res = res[[bat.BtA,bat.BtB,bat.BtC]].ravel()
            L = Vsol[[bat.lambda_BtA,bat.lambda_BtB,bat.lambda_BtC]].ravel()
            L_res = res[[bat.lambda_BtA,bat.lambda_BtB,bat.lambda_BtC]].ravel()
            M_low = Vsol[[bat.nodeA_dual_ineq_Bt,bat.nodeB_dual_ineq_Bt,bat.nodeC_dual_ineq_Bt]].ravel()
            M_low_res = res[[bat.nodeA_dual_ineq_Bt,bat.nodeB_dual_ineq_Bt,bat.nodeC_dual_ineq_Bt]].ravel()
            M_up = Vsol[[bat.nodeA_dual_ineq_Bt_upper,bat.nodeB_dual_ineq_Bt_upper,bat.nodeC_dual_ineq_Bt_upper]].ravel()
            M_up_res = res[[bat.nodeA_dual_ineq_Bt_upper,bat.nodeB_dual_ineq_Bt_upper,bat.nodeC_dual_ineq_Bt_upper]].ravel()

            Bt_prev = bat.Bt_prev
            P_ch = Vsol[[bat.nodeA_p_plus,bat.nodeB_p_plus,bat.nodeC_p_plus]].ravel()
            P_ch_res = res[[bat.nodeA_p_plus,bat.nodeB_p_plus,bat.nodeC_p_plus]].ravel()
            P_d = Vsol[[bat.nodeA_p_minus,bat.nodeB_p_minus,bat.nodeC_p_minus]].ravel()
            P_d_res = res[[bat.nodeA_p_minus,bat.nodeB_p_minus,bat.nodeC_p_minus]].ravel()

            # Store battery variables from this period
            self.results_dict['B'][bat.ID][epoch][t][:] = B
            self.results_dict['B_res'][bat.ID][epoch][t][:] = B_res
            self.results_dict['L'][bat.ID][epoch][t][:] = L
            self.results_dict['L_res'][bat.ID][epoch][t][:] = L_res
            self.results_dict['M_up'][bat.ID][epoch][t][:] = M_up
            self.results_dict['M_up_res'][bat.ID][epoch][t][:] = M_up_res
            self.results_dict['M_low'][bat.ID][epoch][t][:] = M_low
            self.results_dict['M_low_res'][bat.ID][epoch][t][:] = M_low_res
            self.results_dict['P_ch'][bat.ID][epoch][t][:] = P_ch
            self.results_dict['P_ch_res'][bat.ID][epoch][t][:] = P_ch_res
            self.results_dict['P_d'][bat.ID][epoch][t][:] = P_d
            self.results_dict['P_d_res'][bat.ID][epoch][t][:] = P_d_res

        # EXTRACT NODE VOLTAGES
        for node in infeasibility.node:
           pass 

        # EXTRACT LINE POWERS
        for ele in infeasibility.curr_meas:
            # store line power in dictionary
            self.results_dict['S_line'][ele.name][epoch][t] = [ele.Sa, ele.Sb, ele.Sc]

        # EXTRACT GENERATION POWERS (SLACK, PV, WIND)
        slack = infeasibility.node[0]
        self.results_dict['P_g']['slack'][epoch][t] = [slack.Sa, slack.Sb, slack.Sc]
        for load in infeasibility.load:
            if 'wind' in load.name or 'PV' in load.name:
                self.results_dict['P_g'][load.name][epoch][t] = np.abs(load._cP)

        # TOTAL LOAD
        P_total = 0
        Q_total = 0
        for load in infeasibility.load:
            if ('wind' not in load.name) and ('PV' not in load.name):
                P = np.sum(np.asarray(load._cP))
                Q = np.sum(np.asarray(load._cQ))
                P_total += P
                Q_total += Q

        self.results_dict['total_load'][epoch][t] = np.abs(complex(P_total, Q_total))
        
        # TOTAL GEN
        # Total Solar Gen
        P_total = 0
        for load in infeasibility.load:
            if ('wind' not in load.name) and ('PV' not in load.name):
                P = np.sum(np.asarray(load._cP))
                P_total += P

        self.results_dict['total_load'][epoch][t] = np.abs(P_total)



    def debug_print(self, B, L, P_ch, P_d, P):
        """
        Prints Battery SOC, Bt Dual Variable, Charge Rates, and Slack bus Power
        """

        if self.num_batteries > 0:
            print("=======================\nBattery Params by Period\n======================")
            for row in range(len(B)):
                print(f"Period: {row}, SOC: {B[row][0][0]}, L: {L[row][0][0]}, Pch: {P_ch[row][0][0]}, Pd: {P_d[row][0][0]}")


        print("=======================\nSlack Bus Power by Period\n======================")
        for row in range(len(P)):
            print(f'Period: {row}, S: {P[row]}')


