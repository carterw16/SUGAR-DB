#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 14:50:14 2021

@author: emfoster
"""

import numpy as np
import math

from itertools import count
from termcolor import colored
from scipy.sparse.linalg import spsolve

from lib.stamp_linear_homotopy_test import stamp_linear_homotopy_test
from lib.stamp_nonlinear_test import stamp_nonlinear_test
from lib.matrix_conversion import matrix_conversion
from lib.check_zero_rows import check_zero_rows
from lib.apply_voltage_limiting import apply_voltage_limiting
from lib.apply_power_stepping import apply_power_stepping
from lib.apply_Tx_Stepping import apply_Tx_Stepping
from lib.calc_infeas_obj import calc_infeas_obj
from lib.calc_residuals import calc_residuals

from classes.Nodes import Nodes
from classes.Map import Map
from classes.Limiting import Limiting
 
import matplotlib.pyplot as plt
import networkx as nx

import pickle
    
def run_pf(Vinit, Ylin, Ylin_row, Ylin_col, sizeY, idx_Ynlin, Jlin, Jlin_h_row, Jlin_h_val, 
           Jlin_h_col, idx_Ylin_H, overheadline, undergroundline, triplexline, xfmr, 
           node, load, triplexload, ibdg, regulator, switch, maxIter, features, curr_meas, 
           enable_CM, lf, gmin_stepping, tol, settings, stamped_ground, opt_on, infeasibility, casename):

    
    # ==============================================================================
    #                             HELPER FUNCTIONS
    # ==============================================================================
    def almostEqual(d1, d2, epsilon= 10**-7):
            return abs(d2 - d1) < epsilon 


    # Variable definitions - previously in main
    # Raw Variables
    errMax_Store = [1e-1]
    lfStore = []
    num_outer = 0
    _iteration_count = count(0)
    V = np.copy(Vinit)
    
    # Flags
    flag_maxiter = False
    InnerLoopComplete = True
    
    residual_check_count = 0
    
    # Power Stepping
    power_stepping = settings["Power Stepping"]
    pstep__ = 1
    if power_stepping:
        min_pstep = 0.01
        pstepHardLimit = 20
        pstep_mod = 2
        lambdas = 1
        flag_ps_break = False
        flag_ps_continue = False
        not_converged = False
    
    # Homotopy
    homotopy_enabled = settings["Homotopy"]
    h_factor_old = 2
    if homotopy_enabled:
        tol = 1e-4
        turnoffVlim = False
        h_factor = 1
        G_homotopy = 400
        B_homotopy = 400
        h_count = 0
        h_step = 0.1
        h_factorStore = []
        flag_homotopy_failed = False
        flag_GBincrease = False
        V_prev_dict = dict() # Tx_stepping
        large_err_step = 0 # Tx_stepping
    else:
        turnoffVlim = False
        h_factor = -1
        h_step = -1
        G_homotopy = 0
        B_homotopy = 0
        h_count = 0
        h_factorStore = []
        past_h_factor = []
        h_osci_count = 0
        flag_homotopy_failed = False
        
    # Voltage Limiting
    if settings["Voltage Limiting"]:
        voltage_limiting_dict = {'maxstep': 0.1, 'minstep': -0.1, 
                                 'Vr Limit': 2, 'Vi Limit': 2, 'dual': opt_on, 'Vr Norm': None, 'Vi Norm': None}
        voltage_limiting = Limiting('voltage', voltage_limiting_dict)
        voltage_limiting_flag = True
    else:
        voltage_limiting_flag = False
    
    if settings["Diode Limiting"]:
        if_index = Nodes.if_index
        mu_index = Nodes.mu_index
        diode_limiting_dict = {'xmin': np.zeros((np.size(if_index,0),1)), 'x_index': if_index, 'mu_min_index': mu_index}
        diode_limiting = Limiting('diode', diode_limiting_dict)
        diode_limiting_flag = True
    else:
        diode_limiting_flag = False


    obj = []
    
    while True:
        if power_stepping:           
            not_converged = True
        else:
            not_converged = False

        if gmin_stepping:
            pass
        else:
            pass

        if homotopy_enabled:
            homotopyRunning = True
        else:
            homotopyRunning = False
        """ Variable for the While Loop"""
        num_outer += 1
        err_max = 1.1 * tol
        innerLoopCount = 0
        pstep_Store = []
        pstep_count = int(0)
        pstep_osci = int(0)

        # Homotopy Variables
        V_previous = np.copy(V)
        past_h_factor = []
        h_osci_count = 0
        err_dict = []

        while err_max > tol or not_converged or homotopyRunning:
            """ =========IF POWER_STEPPING VOLTAGE DEPENDENT WITH Tx Stepping"""
            # Without homotopy use power_steppingV with power_stepping
            # With homotopy use power_steppingV without power_stepping
            # TODO: Add PowerSteppingV ???

            # Increase iteration count
            innerLoopCount += 1
            iterationCount = next(_iteration_count)

            if turnoffVlim is True:
                voltage_limiting = False

            if iterationCount > maxIter:
                flag_maxiter = True
                break

            # TODO: Add Function to Check Generator Limits to Generators Themselves

            # TODO: Incorporate Linear Stamp Updates
            """ ===========================STAMP Non-Linear ==================="""
            # Stamp Non-Linear Elements
            (Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin,
              idx_Jnlin) = stamp_nonlinear_test(node, load, triplexload, ibdg,
                                          curr_meas, enable_CM, V, idx_Ynlin,
                                          lf, pstep__, 1, opt_on, infeasibility) # removed pstepV__

            # Create Y and J nonlin matrix
            Jnlin_col = np.zeros((len(Jnlin_row)), dtype=np.int64)

            Ynlin, Jnlin = matrix_conversion(Ynlin_val, Ynlin_row, Ynlin_col,
                                              Jnlin_row, Jnlin_col, Jnlin_val,
                                              sizeY)
            """===============================End==============================="""
            """---------------------Homotopy Stepping---------------------------"""
            # Stamp Homotopy Circuits
            if homotopy_enabled:
                if not almostEqual(h_factor, h_factor_old) or flag_GBincrease:
                    flag_GBincrease = False
                    h_factor_old = h_factor
                    (Ylin_h_val, Ylin_h_row, Ylin_h_col, Jlin_h_row, Jlin_h_val,
                     idx_Ylin_H, idx_Jlin_H, stamped_ground) = stamp_linear_homotopy_test(
                         overheadline, undergroundline, triplexline, xfmr,
                         curr_meas, features, node, V, curr_meas, stamped_ground,
                         idx_Ylin_H, h_factor,
                         G_homotopy, B_homotopy, False)
                    Jlin_h_col = np.zeros((len(Jlin_h_row)), dtype=np.int64)
                    Ylin_H, Jlin_H = matrix_conversion(Ylin_h_val, Ylin_h_row,
                                                       Ylin_h_col, Jlin_h_row,
                                                       Jlin_h_col, Jlin_h_val,
                                                       sizeY)
            else:
                if not almostEqual(h_factor, h_factor_old):
                    h_factor_old = h_factor
                    h_factor_const = 0
                    (Ylin_h_val, Ylin_h_row, Ylin_h_col, Jlin_h_row, Jlin_h_val,
                      idx_Ylin_H, idx_Jlin_H, stamped_ground) = stamp_linear_homotopy_test(
                          overheadline, undergroundline, triplexline, xfmr,
                          curr_meas, features, node, V, curr_meas, stamped_ground,
                          idx_Ylin_H, h_factor_const,
                          G_homotopy, B_homotopy, False)
                    Jlin_h_col = np.zeros((len(Jlin_h_row)), dtype=np.int64)
                    Ylin_H, Jlin_H = matrix_conversion(Ylin_h_val, Ylin_h_row,
                                                        Ylin_h_col, Jlin_h_row,
                                                        Jlin_h_col, Jlin_h_val,
                                                        sizeY)
            """====================Combine The linear and Non-linear parts======"""
            # TO DO: Nonlinear components currently out 
            Y = Ylin + Ynlin + Ylin_H 
            J = Jlin + Jnlin + Jlin_H 
            
            # Check Y = Y.T
            # Y_temp = abs(Y.todense().round(decimals = 6)) == abs(Y.todense().round(decimals = 6).T)

            # # Diagonal Dominance Check
            # # D = np.abs(Y.diagonal())
            # # S = np.sum(np.abs(Y),axis=1) - D
            # # if np.all(D > S):
            # # 	print("Diagonally Dominant")
            # # else:
            # # 	print("Not Diagonally Dominant")
            #
            # # Eigenvalue Check
            # w, v = eigs(Y)

            # Reformat the Y matrix to not include rows and cols of zeros
            (Y_red, J_red, sol_index) = check_zero_rows(Y, J)
            
            """ ==================Solve the system and calculate the error======"""
            Vsol_red = spsolve(Y_red, J_red)
            Vsol_red = np.reshape(Vsol_red, (Vsol_red.shape[0], 1))
            # Map the Vsol solution to the original V vector
            Vsol = np.zeros((sizeY, 1), dtype=float)
            Vsol[sol_index] = Vsol_red
            
            # m = np.delete(Nodes.if_index, np.array([184, 200, 216, 1048]))
            # Vsol[m] = 1e-10
            # Calculate the error for only the voltage variables
            Verr = V#[Nodes.voltage_index]
            
            Vsol_err = Vsol#[Nodes.voltage_index]
            err = Verr - Vsol_err
            err_max = np.amax(abs(err))
            if np.argmax(abs(err)) in Nodes.voltage_index:
                #err_argmax = Nodes.voltage_index[np.argmax(abs(err))]
                err_argmax = np.argmax(abs(err))
            elif np.argmax(abs(err)) in Nodes.Lr_index:
                err_argmax = np.argmax(abs(err))
            elif np.argmax(abs(err)) in Nodes.Li_index:
                err_argmax = np.argmax(abs(err))
            else:
                err_argmax = len(V) + 1
 
            errMax_Store.append(err_max)
            lfStore.append(lf)
            # Figure out what node has the error
            err_max_node = [
                ele.name for ele in node if err_argmax in ele.node_set
            ]
            if err_max_node:
                print(
                    colored(
                        'Maximum error from this iteration is %f at node %s' %
                        (err_max, err_max_node[0]), 'red'))
            else:
                print(
                    colored('Maximum error from this iteration is %f' % err_max,
                            'red'))
            if math.isnan(err_max):
                break
            
            err_dict.append(err_max)
            
            
            
            
            # Map.make_dual_map(node, overheadline, undergroundline, xfmr, regulator, switch)
            # Map.map_V(V, Ylin_row, Ylin_col, Ynlin_row, Ynlin_col, Ylin_h_row, Ylin_h_col, node, h_factor, xfmr, regulator, switch)
            """-------------------- Calculating Objective (Debugging) --------"""
            # if opt_on:
            #     obj.append(calc_infeas_obj(Vsol, node, xfmr, regulator, switch))

            """--------------------  Tx Steppping ------------------------------"""
            flagVupdate = False

            if homotopy_enabled:
                # if round(h_factor, 6) == 0.000371:
                #     name_node = "Initial_Conditions/" + casename + "_L1node_Homotopy_" + str(lf) + ".pkl"
                #     with open(name_node, "wb") as f:
                #         pickle.dump(node, f)
                    
                #     name_V =  "Initial_Conditions/" + casename + "_L1_V_Homotopy_" + str(lf) + ".npy"
                #     np.save(name_V, Vsol)
                # increase the GB values if little or no progress
                if h_factor == 1:
                    if innerLoopCount > 145 and err_max > 0.1 or innerLoopCount > 145:
                        G_homotopy += 400
                        B_homotopy += 400
                        innerLoopCount = 0
                        flag_GBincrease = True
                        V = np.copy(Vinit)
                        print(colored("Homotopy G and B increased", 'green'))
                        continue
                # Implementation of Vmax to check h_Factor algorithm
                if err_max < tol:
                    # Check homotopy factor almost equal to zero
                    if almostEqual(h_factor, 0.0):
                        tol = 1e-4
                    #TODO (nturnerb@cmu.edu): Fix max and min voltage calculations
                    # (Vmax_node, Vmax, Vmin_node, Vmin, Vmax_ang_node, Vmax_ang,
                    #  Vmin_ang_node, Vmin_ang) = calcVrange(node, Vsol)
                    # print(
                    #     colored(
                    #         'Max - Vmag: %0.5f Vang:%0.5f Node: %s' %
                    #         (Vmax, Vmax_ang, node[Vmax_node].name), 'blue'))
                    # print(
                    #     colored(
                    #         'Min - Vmag: %0.5f Vang:%0.5f Node: %s' %
                    #         (Vmin, Vmin_ang, node[Vmin_node].name), 'blue'))
                if h_factor == 0 and err_max < tol:
                    homotopyRunning = False
                    break
                else:
                    (h_factor, h_step, 
                     past_h_factor, Vsol, V_previous, h_count,
                     h_osci_count, flagVupdate, h_factorStore, large_err_step,
                     Tx_stepping_failed, V_prev_dict) = \
                     apply_Tx_Stepping(err_max, h_factor, h_count, h_osci_count,
                                       past_h_factor, h_step, V_previous, tol, Vsol,
                                       flagVupdate, h_factorStore, large_err_step,
                                       V_prev_dict)
                    if Tx_stepping_failed:
                        flag_homotopy_failed = True
                        break
                if flagVupdate:
                    V = np.copy(Vsol)
                    continue
            """--------------------  Variable Limiting -------------------------"""
            # TODO: Finish Variable Limiting Implementation
            # if variable_limiting and err_max > tol:
            # Find Vmax, Vmin, etc.
            # (Vmax_node, Vmax, Vmin_node, Vmin, Vmax_ang_node,
            #  Vmax_ang, Vmin_ang_node, Vmin_ang) = calcVrange(node, Vsol)
            #     apply_variable_limiting(Vsol, Vmax_lim, Vmin_lim, sigma, sigma_count, errMax_Store, node, Nodes)
            
            """--------------------  Voltage Limiting -------------------------"""
            if voltage_limiting_flag and err_max > tol:
                V = voltage_limiting.apply_voltage_limiting(Vsol, V)
                # V = apply_voltage_limiting(Vsol, V, Nodes, maxStep__, minStep__,
                #                        Vr_norm, Vi_norm, VrLimit__, ViLimit__, opt_on)

            """--------------------  Voltage Limiting -------------------------"""  
            if diode_limiting_flag and err_max > tol:
                V = diode_limiting.apply_diode_limiting(Vsol, V)
                if np.size(np.where(V[Nodes.mu_index] < 0)[0]) != 0:
                    print('Mus are negative')
            else:
                V = np.copy(Vsol)    
            """--------------------  Residuals -------------------------""" 
            if iterationCount % 10 == 0 and opt_on == True:
                res_eqn = calc_residuals(V, Ylin, Jlin, node, load, triplexload, lf, opt_on, infeasibility.obj, infeasibility.epsilon)
            elif iterationCount % 10 == 0 and opt_on == False:
                res_eqn = calc_residuals(V, Ylin, Jlin, node, load, triplexload, lf, opt_on, None, None)
               # if abs(np.sum(res_eqn)) >= 3000:
               #     V = Vinit
                   # V[Nodes.if_index] = np.rand()
                   # V[Nodes.mu_index] = 1e-6
                   # V[Nodes.L_index] = 1
               # print('Sum of residuals at iteration',iterationCount, 'is', np.sum(res_eqn))
            # if residual_check_count == 300:
            #     res_eqn = calc_residuals(V, Ylin, Jlin, node, load, triplexload, lf, opt_on, infeasibility.obj, infeasibility.epsilon)
            #     residual_check_count = 0
            # else:
            #     residual_check_count += 0
            
            """=========================End====================================="""

            if not power_stepping and not homotopy_enabled and err_max < tol:
                break
            """ ======================Apply Power Stepping======================"""
            if power_stepping:
                (pstep__, flag_ps_continue, flag_ps_break, pstep_mod, V,
                 lambdas, innerLoopCount, not_converged, pstep_count,
                 pstep_osci, pstep_Store) = apply_power_stepping(
                     innerLoopCount, pstepHardLimit, errMax_Store, pstep_count,
                     pstep_osci, err_max, tol, pstep__, min_pstep, Vinit, V,
                     False, lambdas, not_converged, pstep_Store,
                     pstep_mod)
                if flag_ps_continue:
                    continue
                elif flag_ps_break:
                    print('Case did not converge - Power Stepping Failed')
                    break
        """##############OUTER LOOP OF POWERFLOW OPTIONS#################"""
        """Mininum power_step or minimum iter fatal stop"""
        if flag_maxiter:
            error_occurred = True
            print("Simulation has failed in power flow. Maximum iteration count reached.")
            outputs = {'inside iterations': iterationCount, 'node info': node, 
                       'Converged V': Vsol, 'outer iters': num_outer, 'primal objective value': obj}
            return outputs
        elif homotopy_enabled:
            if flag_homotopy_failed:
                error_occurred = True
                print("Simulation has failed in power flow with homotopy enabled.")
                break
        
        if InnerLoopComplete:
            power_stepping = False
            homotopy_enabled = False
            turnoffVlim = False
            InnerLoopComplete = False

        # Finish The Simulation
        if err_max < tol:
            # (Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin,
            #   idx_Jnlin) = stamp_nonlinear_test(node, load, triplexload, ibdg,
            #                               curr_meas, enable_CM, Vsol, idx_Ynlin,
            #                               lf, pstep__, 1, opt_on) 
            # # Create Y and J nonlin matrix
            # Jnlin_col = np.zeros((len(Jnlin_row)), dtype=np.int64)

            # Ynlin, Jnlin = matrix_conversion(Ynlin_val, Ynlin_row, Ynlin_col,
            #                                   Jnlin_row, Jnlin_col, Jnlin_val,
            #                                   sizeY)
            
            # Y = Ylin + Ynlin + Ylin_H
            outputs = {'inside iterations': iterationCount, 'node info': node, 
                       'Converged V': Vsol, 'outer iters': num_outer, 'primal objective value': obj, 'lastY': Y}
        else:
            outputs = []
        
        # if opt_on:
        #     plt.plot(obj)
        #     plt.title('Objective value per iteration')
        #     plt.show()
        error_name = "Initial_Conditions/error" + casename + "_L1node.npy"
        np.save(error_name, err_dict)
        print('Number of iterations:', iterationCount)
        
        Yname_node = "Initial_Conditions/" + casename + "_L1_Y_" + str(lf) + ".pkl"
        with open(Yname_node, "wb") as f:
            pickle.dump(Y, f)
        
        Ylinname_node = "Initial_Conditions/" + casename + "_L1_Ylin_" + str(lf) + ".pkl"
        with open(Ylinname_node, "wb") as f:
            pickle.dump(Ylin, f)
        return outputs
