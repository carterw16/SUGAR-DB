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
        #self.wind_capacity = settings['multi settings']['wind capacity']
        #self.PV_capacity = settings['multi settings']['PV capacity']
        self.num_batteries = len(self.casedata.battery_sources)


    def runMP(self):
        mp_start_time = time.perf_counter()


        # Multi-Epoch Results Dictionary
        self.results_dict = {'B': [], 'B_res': [], 'P_ch': [], 'P_ch_res': [], 
                'P_d': [], 'P_d_res': [], 'L': [], 'L_res': [], 
                'M_low': [], 'M_low_res': [], 'M_up': [], 'M_up_res': []}

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

            # Extend results dict
            for var in self.results_dict:
                self.results_dict[var].append(np.zeros((self.periods, self.num_batteries, 3)))

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

                Vsol = infeasibility.Vsol
                res = infeasibility.res_eqn
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
                    print(f'P_d_res: {P_d_res[0]}')
                    if t < self.periods - 1:
                        print(f'Lb_next: {L_DDP[t+1][0]}')
                    
                    # Update running Lb_next array for DDP
                    L_DDP[t][idx] = L

                    # Store battery variables from this period
                    self.results_dict['B'][epoch][t][idx][:] = B
                    self.results_dict['B_res'][epoch][t][idx][:] = B_res
                    self.results_dict['L'][epoch][t][idx][:] = L
                    self.results_dict['L_res'][epoch][t][idx][:] = L_res
                    self.results_dict['M_up'][epoch][t][idx][:] = M_up
                    self.results_dict['M_up_res'][epoch][t][idx][:] = M_up_res
                    self.results_dict['M_low'][epoch][t][idx][:] = M_low
                    self.results_dict['M_low_res'][epoch][t][idx][:] = M_low_res
                    self.results_dict['P_ch'][epoch][t][idx][:] = P_ch
                    self.results_dict['P_ch_res'][epoch][t][idx][:] = P_ch_res
                    self.results_dict['P_d'][epoch][t][idx][:] = P_d
                    self.results_dict['P_d_res'][epoch][t][idx][:] = P_d_res


                
                # FOWARD PASS UPDATE
                for idx in range(0, self.num_batteries):
                    infeasibility.battery_sources[idx].Bt_prev = self.results_dict['B'][epoch][t][idx]
                
                # extract vsol to use as initialization for next solution
                self.Vsol_mp[t] = Vsol.reshape(-1)

                # extract 1st period vsol to use for next epoch
                # if t == 0:
                #    Vsol_period1 = infeasibility.Vsol

 
            if epoch > 0:
                # Calculate Max Error
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
            #self.debug_print(B, L, P_ch, P_d, P)
            print("="*30)
            print(f"EPOCH {epoch}")
            print("="*30)

        # After all epochs, store results dict as pickle for later plotting
        self.write_results()

    def write_results(self): 
        """
        writes time-variant battery variables to a pickle

        saved variables and their residuals:
            - Bt SOC lagrange equation: Lb, mu_up, mu_low
            - Lb lagrange equation: Bt, P_d, P_ch
        """
        file_path = "battery_outputs.pkl"
        with open(file_path, 'wb') as pickle_file:
            pickle.dump(self.results_dict, pickle_file)
        

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


