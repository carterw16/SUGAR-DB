#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 17:04:41 2021

@author: emfoster
"""

"""
Created on Wed May 12 10:42:26 2021

@author: emfoster
"""

from classes.Map import Map
from classes.Nodes import Nodes
from classes.GlobalVars import connected_nodes

from lib.normalize_if import normalize_if
from lib.assign_nodes import assign_nodes
from lib.save_output import save_output
from lib.dual_analysis import dual_analysis
from lib.initialize import initialize
from lib.initialize_warm_start import initialize_warm_start
from lib.stamp_linear import stamp_linear
from lib.stamp_nonlinear import stamp_nonlinear
from lib.matrix_conversion import matrix_conversion
from lib.check_zero_rows import check_zero_rows

from classes.CurrentMeasurements import create_current_meas
from classes.Homotopy import Homotopy
from classes.Limiting import Limiting

import time
import pickle
import numpy as np
import math
import logging

from itertools import count
from termcolor import colored
from scipy.sparse.linalg import spsolve
from types import MethodType
from scipy.sparse import csr_matrix
from classes.Limiting import Limiting

import ipdb

from .InfeasibilityStamps import (stamp_ineq_duals, stamp_ineq_duals_tplx, 
                                  stamp_infeas_current_L2, stamp_infeas_current_L2_tplx,
                                  stamp_infeas_current_L1, stamp_infeas_current_L1_tplx,
                                  stamp_infeas_dual_relationship, stamp_infeas_dual_relationship_tplx,
                                  stamp_linear_infeasibility, stamp_nonlinear_infeasibility)


class Infeasibility():
    def __init__(self, testcase, casedata, features, settings, path_to_output, node_key):
        '''

        Parameters
        ----------
        stamp_dual : boolean variable that flags whether or not to incorporate dual variables
        obj : integer value of either 1 for the L1-norm or 2 for the L2-norm that denotes which sets of equations to stamp

        Returns
        -------
        None.

        '''
        self.features = features
        self.settings = settings
        self.path_to_output = path_to_output
        self.case_name = testcase

        # Infeasibility Settings
        self.stamp_dual = settings["Stamp Dual"]
        self.obj = settings["Obj norm"]
        self.normalize_if = False
        self.report_if = True

        # One Norm Settings
        self.cs_eps = 1e-6
        self.epsilon = self.cs_eps 
        self.maxIter = settings["Max Iters"]

        # Newton-Raphson Settings
        self.NR_tolerance = settings["Tolerance"]

        # Initialization Settings
        self.warm_start = settings["Warm Start"][0]
        self.load_factor = settings["Load Factor"]
        
        # Output Settings
        self.save_output_conditions = False
        self.draw_network_maps = False
        self.draw_infeasibility_maps = False
        self.save_results = settings["Save File"]
        
        """---------------------- Initialize Heuristic Settings -------------------------"""        
        # Homotopy
        if settings["Homotopy"]:
            self.homotopy_enabled = True
            g_homotopy = settings["G_homotopy"]
            b_homotopy = settings["B_homotopy"]
            h_factor_init = 1
            homotopy_tolerance = 1e-3
            h_step_size_init = 0.1   
            self.homotopy = Homotopy(homotopy_tolerance, g_homotopy, b_homotopy, h_factor_init, h_step_size_init) 
        else:
            self.homotopy = None
            self.homotopy_enabled = False
       
        # Voltage Limiting
        if settings["Voltage Limiting"]:
            voltage_limiting_dict = {'maxstep': 0.1, 'minstep': -0.1, 
                                    'Vr Limit': 2, 'Vi Limit': 2, 'dual': True, 'Vr Norm': None, 'Vi Norm': None}
            self.voltage_limiting = Limiting('voltage', voltage_limiting_dict)
            self.voltage_limiting_enabled = True
        else:
            self.voltage_limiting_enabled = False
        
        # Diode Limiting 
        if settings["Diode Limiting"]:
            self.diode_limiting_enabled = True
        else:
            self.diode_limiting_enabled = False

        (node_index_, casedata.node,
         casedata.regulator, casedata.xfmr, casedata.switch,
         casedata.fuse, casedata.reactor, casedata.ibdg) = assign_nodes(node_key, casedata.node, casedata.regulator,
                                                                        casedata.xfmr, casedata.switch,
                                                                        casedata.fuse, casedata.reactor,
                                                                        casedata.ibdg, self.obj)
        
        self.node = casedata.node
        self.load = casedata.load
        self.slack = casedata.slack
        self.triplex_load = casedata.triplex_load
        self.xfmr = casedata.xfmr
        self.line_oh = casedata.ohline
        self.line_ug = casedata.ugline
        self.line_tplx = casedata.triplex_line
        self.capacitor = casedata.capacitor
        self.regulator = casedata.regulator
        self.switch = casedata.switch
        self.fuse = casedata.fuse
        self.reactor = casedata.reactor
        self.ibdg = casedata.ibdg
        self.simulation_stats = casedata.stats

        self.node_index_ = node_index_
        self.node_key = node_key

        enable_CM = self.features['Current Meas']
        self.curr_meas = []
        if enable_CM:
            self.activate_current_measures(enable_CM, self.node_index_, self.node_key)
        
    def almostEqual(d1, d2, epsilon= 10**-7):
            return abs(d2 - d1) < epsilon 
    

    def activate_current_measures(self, enable_CM, node_index_, node_key):
        self.curr_meas = create_current_meas(enable_CM, self.oh_lines, self.ug_lines,
                                             self.triplex_lines)

        for ele in self.curr_meas:
            node_index_ = ele.assign_nodes(node_key, node_index_)

        self.node_index_ = node_index_
    
    def calc_residuals(self, V_out, Ylin, Jlin, nodes, node_key, load, triplex_load, lf, stamp_dual, obj, epsilon):
        '''

        Parameters
        ----------
        V_out : converged solution from Newton Raphson iterations
        Ylin : linear Y matrix
        Jlin : linear J vector
        nodes : list of node elements
        load : list of load elements
        triplex_load : list of triplex load elements
        lf : scalar load factor
        stamp_dual : boolean variable for whether or not to incorporate dual equations
        obj : scalar of value 1 or 2; built out for either L1 norm or L2 norm infeasibility analysis
        epsilon : scalar value (should be small) used as relaxation in complementary slackness

        Returns
        -------
        res_eqn : vector the size of V_out; the result of every original KKT condition 
            evaluated at the optimal values from V; individual elements should equal 0

        '''
    # Initializing the vector which houses the residual equations
        res_eqn = np.zeros_like(V_out)
        
        violation = 0
        # Stationarity Constraints (L2 norm)
        # dL/dif = 0 (located in if node) [all elements for this should be linear and this in Ylin, Jlin]
        # dL/dV = 0 (located in lambda nodes) [some elements in this equation are nonlinear, so these components
            # are directly added to res_eqn for load and triplex load]
        
        # Primal Feasibility (Equality Constraints) [some elements in this equation are nonlinear, so these components
        # are directly added to res_eqn for load and triplex load]
        # eq constr = 0 (located in V nodes) 
        for ele in range(len(load)):
            res_eqn = load[ele].calc_residual(V_out, node_key, nodes, res_eqn, lf)


        for ele in range(len(triplex_load)):
            res_eqn = triplex_load[ele].calc_residual(V_out, node_key, nodes, res_eqn, lf)
        
        dual_feas_flag = np.zeros(len(nodes))
        prim_feas_flag = np.empty(len(nodes))

        if obj == 1 and stamp_dual == True:
            for ele in range(len(nodes)):
                if nodes[ele].bustype != 3:  
                    if nodes[ele].isTriplex:
                        mu_r_lb = [nodes[ele].node1_mu_r_lb, nodes[ele].node2_mu_r_lb, nodes[ele].nodeN_mu_r_lb]
                        mu_i_lb = [nodes[ele].node1_mu_i_lb, nodes[ele].node2_mu_i_lb, nodes[ele].nodeN_mu_i_lb] 
                        mu_r_ub = [nodes[ele].node1_mu_r_ub, nodes[ele].node2_mu_r_ub, nodes[ele].nodeN_mu_r_ub] 
                        mu_i_ub = [nodes[ele].node1_mu_i_ub, nodes[ele].node2_mu_i_ub, nodes[ele].nodeN_mu_i_ub] 
                        
                        if_r_plus =  [nodes[ele].node1_if_r_plus,  nodes[ele].node2_if_r_plus,  nodes[ele].nodeN_if_r_plus]
                        if_r_minus = [nodes[ele].node1_if_r_minus, nodes[ele].node2_if_r_minus, nodes[ele].nodeN_if_r_minus]
                        if_i_plus =  [nodes[ele].node1_if_i_plus,  nodes[ele].node2_if_i_plus,  nodes[ele].nodeN_if_i_plus] 
                        if_i_minus = [nodes[ele].node1_if_i_minus, nodes[ele].node2_if_i_minus, nodes[ele].nodeN_if_i_minus]           
        
                        for i in range(0,3):
                            # Complementary Slackness
                            res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]] * V_out[if_r_minus[i]] - epsilon
                            res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]] * V_out[if_i_minus[i]] - epsilon           
                            res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]] * V_out[if_r_plus[i]] - epsilon
                            res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]] * V_out[if_i_plus[i]] - epsilon
                            
                            if dual_feas_flag[ele] == False:
                                if V_out[mu_r_lb[i]] >= 0 and V_out[mu_i_lb[i]] >= 0 and V_out[mu_r_ub[i]] >= 0 and V_out[mu_i_ub[i]] >= 0:
                                    dual_feas_flag[ele] = False
                                else:
                                    dual_feas_flag[ele] = True
                                                    
                    else:
                        mu_r_lb = [nodes[ele].nodeA_mu_r_lb, nodes[ele].nodeB_mu_r_lb, nodes[ele].nodeC_mu_r_lb, nodes[ele].nodeN_mu_r_lb]
                        mu_i_lb = [nodes[ele].nodeA_mu_i_lb, nodes[ele].nodeB_mu_i_lb, nodes[ele].nodeC_mu_i_lb, nodes[ele].nodeN_mu_i_lb] 
                        mu_r_ub = [nodes[ele].nodeA_mu_r_ub, nodes[ele].nodeB_mu_r_ub, nodes[ele].nodeC_mu_r_ub, nodes[ele].nodeN_mu_r_ub] 
                        mu_i_ub = [nodes[ele].nodeA_mu_i_ub, nodes[ele].nodeB_mu_i_ub, nodes[ele].nodeC_mu_i_ub, nodes[ele].nodeN_mu_i_ub] 
                    
                        if_r_plus = [nodes[ele].nodeA_if_r_plus, nodes[ele].nodeB_if_r_plus, nodes[ele].nodeC_if_r_plus, nodes[ele].nodeN_if_r_plus]        
                        if_r_minus = [nodes[ele].nodeA_if_r_minus, nodes[ele].nodeB_if_r_minus, nodes[ele].nodeC_if_r_minus, nodes[ele].nodeN_if_r_minus]
                        if_i_plus = [nodes[ele].nodeA_if_i_plus, nodes[ele].nodeB_if_i_plus, nodes[ele].nodeC_if_i_plus, nodes[ele].nodeN_if_i_plus]   
                        if_i_minus = [nodes[ele].nodeA_if_i_minus, nodes[ele].nodeB_if_i_minus, nodes[ele].nodeC_if_i_minus, nodes[ele].nodeN_if_i_minus]          
        
                        for i in range(0,4):
                            # Complementary Slackness
                            res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]] * V_out[if_r_minus[i]] - epsilon
                            res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]] * V_out[if_i_minus[i]] - epsilon           
                            res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]] * V_out[if_r_plus[i]] - epsilon
                            res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]] * V_out[if_i_plus[i]] - epsilon
                            
                            if dual_feas_flag[ele] == 0:
                                if V_out[mu_r_lb[i]] >= 0 and V_out[mu_i_lb[i]] >= 0 and V_out[mu_r_ub[i]] >= 0 and V_out[mu_i_ub[i]] >= 0:
                                    dual_feas_flag[ele] = 0
                                else:
                                    dual_feas_flag[ele] = 1
                
                    if V_out[if_r_plus[i]] >= 0 and V_out[if_i_plus[i]] >= 0 and V_out[if_r_minus[i]] >= 0 and V_out[if_i_minus[i]] >= 0:
                        prim_feas_flag[ele] = False
                    else:
                        prim_feas_flag[ele] = True
                else:
                    prim_feas_flag[ele] = False
                    dual_feas_flag[ele] = False
        
        # If res_eqn = 0, then stationarity, primal feasibility (equality), and complementary slackness are met
        res_eqn += (Ylin@V_out - Jlin)

        if any(abs(res_eqn[np.append(Nodes.Vr_index, Nodes.Vi_index)]) > 1e-6):
            print('Primal feasibility (equality) condition violated')
            violation = 1
        
        if stamp_dual == True:
            if obj == 1:
                if any(res_eqn[np.append(np.append(Nodes.if_r_plus_index, Nodes.if_i_plus_index), np.append(Nodes.if_r_minus_index, Nodes.if_i_minus_index))] > 1e-6):
                # if any(res_eqn[np.append(Nodes.if_r_index, Nodes.if_i_index)] > 1e-3):
                    print('Stationarity condition (dL/dif) violated')   
                    violation = 1
                
                if any(res_eqn[Nodes.L_index] > 1e-6):
                    print('Stationarity condition (dL/dVr) violated')
                    violation = 1
                
                if np.any(prim_feas_flag) == True:
                    print('Primal feasibility (inequality) violated at nodes', np.where(prim_feas_flag == True))
                    violation = 1
                if np.any(dual_feas_flag) == 1:
                    print('Dual feasibility (inequality) violated at nodes', np.where(dual_feas_flag == True)) 
                    violation = 1
        if violation == 0:
            print("No violations from residuals")
            
        if violation == 1:
            list_of_nodes = []
            violated_indices = np.where(abs(res_eqn) >= 1e-8)[0]
            
            for i in range(len(nodes)):
                node_set = np.array(nodes[i].node_set)
                if any(node_num in violated_indices for node_num in node_set):
                    list_of_nodes.append(nodes[i].name)
        return res_eqn
               
    def run_infeasibility(self):
        returncode = {'success': 0, 'fail': 1, 'error': 2}
        flag_tx = False
        
        simulation_start_time = time.perf_counter()
        # Size of the solution vector
        sizeY = int(repr(self.node_index_)[6:-1])   
        print("size Y is ", sizeY)
        # Initialize V vector
        Vinit = np.zeros([sizeY, 1], dtype=np.float64)
        """------------------------------End----------------------------------------"""
        """-----------------------Initialize V and Q Vector-------------------------"""

        Vinit = initialize(self.node, self.node_key, self.regulator, self.ibdg, Vinit, self.load_factor, self.stamp_dual, self.obj)
        if self.warm_start:
            warm_start_settings = self.settings["Warm Start"]
            Vinit = initialize_warm_start(Vinit, self.case_name, self.load_factor, self.node, warm_start_settings)  

        # Initializing initial size
        if sizeY < 1000:
            idx_Ylin = 100 * sizeY
            idx_Ynlin = 100 * sizeY
            idx_Ylin_H = 100 * sizeY
        else:
            idx_Ylin = 10 * sizeY
            idx_Ynlin = 10 * sizeY
            idx_Ylin_H = 10 * sizeY
    
        print(
            '======================= Variable Initialization Complete =============================='
        )
        """------------------------------End----------------------------------------"""
    
        """ ===========================STAMP Linear ================================"""
        stamped_ground = set()
    
        # Stamp Linear Elements
        (Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Ylin,
         idx_Jlin, stamped_ground) = stamp_linear(self.node_key, self.slack, self.line_oh, 
                                self.line_ug, self.line_tplx,
                                self.node, self.xfmr, self.regulator, self.switch, self.fuse, self.capacitor,
                                self.reactor, self.curr_meas, self.features, stamped_ground, idx_Ylin, Vinit, 
                                self.stamp_dual)
        
        for ele in range(len(self.node)):
            if self.node[ele].bustype != 3:
                Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Ylin, idx_Jlin = self.stamp_linear_infeasibility(self.node[ele], Ylin_val,
                                                                                                            Ylin_row, Ylin_col, 
                                                                                                            idx_Ylin, Jlin_val, 
                                                                                                            Jlin_row, idx_Jlin)

        Jlin_col = np.zeros((len(Jlin_row)), dtype=np.int64)
        Ylin, Jlin = matrix_conversion(Ylin_val, Ylin_row, Ylin_col, Jlin_row,
                                       Jlin_col, Jlin_val, sizeY)
        
        # Raw Variables
        errMax_Store = [1e-1]
        num_outer = 0
        _iteration_count = count(0)
        
        # Flags
        flag_maxiter = False
        flag_isnan = False
        
        # Map System
        if self.draw_network_maps:
            Map.make_KVA_map(self.node, self.line_oh, self.line_ug, self.line_tplx, self.xfmr, self.regulator, self.fuse, self.switch)
            Map.make_map(self.node, self.line_oh, self.line_ug, self.line_tplx, self.xfmr, self.regulator, self.fuse, self.switch)
            
        if self.diode_limiting_enabled:
            if_index = Nodes.if_index
            mu_index = Nodes.mu_index
            diode_limiting_dict = {'xmin': np.zeros((np.size(if_index,0),1)), 'x_index': if_index, 'mu_min_index': mu_index}
            self.diode_limiting = Limiting('diode', diode_limiting_dict)
        
        if self.homotopy_enabled:
            self.homotopy.Vinit = Vinit
        V = np.copy(Vinit)
        V[self.node[0].sr_index] = 1e-3
        V[self.node[0].si_index] = 1e-3
        # V[self.node[0].if_r_minus_index] += 1e-6
        # V[self.node[0].if_i_minus_index] += 1e-6
        rng = np.random.default_rng(0)

        while True:
            """ Variable for the While Loop"""
            num_outer += 1
            err_max = 1.1 * self.NR_tolerance
            innerLoopCount = 0
            res_eqn = self.calc_residuals(Vinit, Ylin, Jlin, self.node, self.node_key, self.load, self.triplex_load, self.load_factor, True, self.obj, self.cs_eps)
    
            while err_max > self.NR_tolerance:
    
                # Increase iteration count
                innerLoopCount += 1
                iteration_count = next(_iteration_count)
    
                # if turnoffVlim is True:
                #     voltage_limiting = False
    
                if iteration_count > self.maxIter:
                    flag_maxiter = True
                    break
    
                """ ===========================STAMP Non-Linear ==================="""
                # Stamp Non-Linear Elements
                (Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin,
                  idx_Jnlin) = stamp_nonlinear(self.node_key, self.node, self.regulator, self.load, self.triplex_load, self.ibdg,
                                              self.curr_meas, self.features, V, idx_Ynlin,
                                              self.load_factor, self.homotopy_enabled, self.homotopy, True) 
                for ele in range(len(self.node)):
                    if self.node[ele].bustype != 3:
                        Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Ynlin, idx_Jnlin = self.stamp_nonlinear_infeasibility(self.node[ele], V, 
                                                                                                                    Ynlin_val, Ynlin_row, Ynlin_col, 
                                                                                                                    idx_Ynlin, Jnlin_val, Jnlin_row, 
                                                                                                                    idx_Jnlin)

                # Create Y and J nonlin matrix
                Jnlin_col = np.zeros((len(Jnlin_row)), dtype=np.int64)
    
                Ynlin, Jnlin = matrix_conversion(Ynlin_val, Ynlin_row, Ynlin_col,
                                                  Jnlin_row, Jnlin_col, Jnlin_val,
                                                  sizeY)
                """===============================End==============================="""
                """---------------------Homotopy Stepping---------------------------"""
                # Create linear homotopy vectors
                
                Ylin_val_H = np.zeros(2*idx_Ylin_H, dtype = float)
                Ylin_col_H = np.zeros(2*idx_Ylin_H, dtype = int)
                Ylin_row_H = np.zeros(2*idx_Ylin_H, dtype = int)    

                Jlin_val_H = np.zeros(2*idx_Ylin_H, dtype = float)  
                Jlin_row_H = np.zeros(2*idx_Ylin_H, dtype = int) 
                Jlin_col_H = np.zeros(2*idx_Ylin_H, dtype = int)  
                
                # Stamp homotopy
                if self.homotopy_enabled and self.homotopy.h_factor != 0:
                    Ylin_row_H, Ylin_col_H, Ylin_val_H, idx_Ylin_H, stamped_ground = self.homotopy.stamp_tx_homotopy(
                                               self.node, self.node_key, self.line_oh, self.line_ug, 
                                               self.line_tplx, self.xfmr, stamped_ground, Ylin_val_H, 
                                               Ylin_row_H, Ylin_col_H, self.obj, Jlin_val_H, Jlin_row_H, V)
                
                #Ylin_H = csr_matrix((Ylin_val_H, (Ylin_row_H, Ylin_col_H)), shape=(sizeY, sizeY), dtype=np.float64)
                Ylin_H, Jlin_H = matrix_conversion(Ylin_val_H, Ylin_row_H, Ylin_col_H,
                                                  Jlin_row_H, Jlin_col_H, Jlin_val_H,
                                                  sizeY)
               
                """====================Combine The linear and Non-linear parts======"""
                Y = Ylin + Ynlin + Ylin_H 
                J = Jlin + Jnlin + Jlin_H
    
                # Reformat the Y matrix to not include rows and cols of zeros
                (Y_red, J_red, sol_index) = check_zero_rows(Y, J)
                
                """ ==================Solve the system and calculate the error======"""
                Vsol_red = spsolve(Y_red, J_red)
                Vsol_red = np.reshape(Vsol_red, (Vsol_red.shape[0], 1))
                # Map the Vsol solution to the original V vector
                Vsol = np.zeros((sizeY, 1), dtype=float)
                Vsol[sol_index] = Vsol_red
                
                # Calculate the error
                Verr = V
                
                Vsol_err = Vsol
                err = Verr - Vsol_err
                err_max = np.amax(abs(err))
                if np.argmax(abs(err)) in Nodes.voltage_index:
                    err_argmax = np.argmax(abs(err))
                elif np.argmax(abs(err)) in Nodes.Lr_index:
                    err_argmax = np.argmax(abs(err))
                elif np.argmax(abs(err)) in Nodes.Li_index:
                    err_argmax = np.argmax(abs(err))
                else:
                    err_argmax = len(V) + 1
     
                errMax_Store.append(err_max)
                
                # Figure out what node has the error
                err_max_node = [
                    ele.name for ele in self.node if err_argmax in ele.node_set
                ]
                if err_max_node:
                    print(
                        colored(
                            'Maximum error from iteration %s is %f at node %s' %
                            (iteration_count, err_max, err_max_node[0]), 'red'))
                else:
                    print(
                        colored('Maximum error from iteration %s is %f' % (iteration_count, err_max),
                                'red'))
                if math.isnan(err_max):
                    flag_isnan = True
                    break
                
                
                err = V - Vsol
                err_max = np.amax(abs(err))
                
    
                """--------------------  Tx Steppping ------------------------------"""
                # Run tx_stepping
                if self.homotopy_enabled:
                    Vsol = self.homotopy.run_tx_stepping(err_max, Vsol)
                    
                    # If homotopy failed, break the loop
                    if self.homotopy.tx_stepping_failed == True:
                        flag_tx = True
                        break

                    if self.homotopy.perturb_conditions == True:
                        random_ones = rng.random(len(Vsol))
                        random_ones[random_ones > 0.5] = 1
                        random_ones[random_ones <= 0.5] = -1
                        random_ones = random_ones.reshape((len(random_ones),1))

                        if self.obj == 1:
                            # temp_indices = np.where(Vsol[self.node[0].if_index] > 1e-6)
                            # new_Vinit = Vsol + random_ones*0.15*Vsol
                            #new_Vinit[self.node[0].if_index] += 1
                            random_infeas_ones = rng.random(len(Vsol[self.node[0].if_index]))
                            random_infeas_ones[random_infeas_ones > 0.25] = 10
                            random_infeas_ones[random_infeas_ones <= 0.25] = 0
                            new_Vinit = Vsol # + random_ones * 0.15 * Vsol
                            new_Vinit[self.node[0].if_index] += random_infeas_ones.reshape((len(random_infeas_ones),1)) #np.reshape(random_infeas_ones, (len(random_infeas_ones), 1))
                        elif self.obj == 2:
                            new_Vinit = Vsol + random_ones*0.15*Vsol + 1e-3
                        else:
                            new_Vinit = Vsol + random_ones*0.15*Vsol
                        Vsol = new_Vinit

                        self.homotopy.perturb_conditions = False
                
                """--------------------  Voltage Limiting -------------------------"""
                if self.voltage_limiting_enabled and err_max > self.NR_tolerance:
                    print("voltage limiting")
                    V = self.voltage_limiting.apply_voltage_limiting(Vsol, V)
    
                """--------------------  Voltage Limiting -------------------------"""  
                if self.diode_limiting_enabled and err_max > self.NR_tolerance:
                    ipdb.set_trace()
                    print("diode liminting")
                    V = self.diode_limiting.apply_diode_limiting(Vsol, V)
                    if np.size(np.where(V[Nodes.mu_index] < 0)[0]) != 0:
                        print('Mus are negative')
                else:
                    V = np.copy(Vsol)    
                """--------------------  Residuals -------------------------""" 
                # if iteration_count % 10 == 0:
                #	res_eqn = self.calc_residuals(V, Ylin, Jlin, self.node, self.node_key, self.load, self.triplex_load, self.load_factor, True, self.obj, self.cs_eps)
                
                """=========================End====================================="""
    
                if not self.homotopy_enabled and err_max < self.NR_tolerance:
                    break
           
            """##############OUTER LOOP OF POWERFLOW OPTIONS#################"""
            """Mininum power_step or minimum iter fatal stop"""

            if flag_maxiter:
                logging.debug(
                    "Simulation has failed. Maximum iteration count reached.")
                exit()
                break
            
            if flag_tx:
                logging.debug(
                    "Simulation has failed with homotopy.")
                exit()
                break
            
            if flag_isnan:
                logging.debug("Simulation has failed.")
                exit()
                break
            
            if self.homotopy_enabled:
                if err_max < self.NR_tolerance and self.homotopy.tx_stepping_successful == True:
                    break
            else:
                if err_max < self.NR_tolerance:
                    break

        simulation_end_time = time.perf_counter()
        sim_total_time = simulation_end_time - simulation_start_time
        self.simulation_stats.append(sim_total_time)
        logging.info(colored('Simulation time is %f seconds' % (sim_total_time), 'blue'))
    
        if flag_maxiter:
            logging.error('Case did not converge - Maximum Iteration Reached')
            return returncode["fail"]
        elif flag_tx:
            logging.error('Case did not converge - Homotopy failed')
            return returncode["fail"]          
        else:
            iteration_count = next(_iteration_count)
            logging.info(
                colored(
                    'Case converged in %d iterations - SUCCESS!' % iteration_count,
                    'green'))        
        
            if self.save_output_conditions:
                if self.obj == 2:
                    name_node = "Initial_Conditions/" + self.case_name + "_L2node_" + str(self.load_factor) + ".pkl"
                    with open(name_node, "wb") as f:
                        pickle.dump(self.node, f)
                    
                    name_V =  "Initial_Conditions/" + self.case_name + "_L2_V_" + str(self.load_factor) + ".npy"
                    np.save(name_V, Vsol)
                
                if self.obj == 1:
                    name_node = "Initial_Conditions/" + self.case_name + "_L1node_" + str(self.load_factor) + ".pkl"
                    with open(name_node, "wb") as f:
                        pickle.dump(self.node, f)
                    
                    name_V =  "Initial_Conditions/" + self.case_name + "_L1_V_" + str(self.load_factor) + ".npy"
                    np.save(name_V, Vsol)
                
                if self.obj == 0:
                    name_node = "Initial_Conditions/" + self.case_name + "_pf_node_" + str(self.load_factor) + ".pkl"
                    with open(name_node, "wb") as f:
                        pickle.dump(self.node, f)
                    
                    name_V =  "Initial_Conditions/" + self.case_name + "_pf_V_" + str(self.load_factor) + ".npy"
                    np.save(name_V, Vsol)
                    
            res_eqn = self.calc_residuals(Vsol, Ylin, Jlin, self.node, self.node_key, self.load, self.triplex_load, self.load_factor, True, self.obj, self.cs_eps)
            
            if self.normalize_if:
                normalized_if, max_if = normalize_if(self.node, Vsol, self.obj)
            
            if self.report_if and self.stamp_dual:
                dual_info_ABC, dual_info_tplx = dual_analysis(Vsol, self.node, self.obj)
            else:
                dual_info_ABC = None
                dual_info_tplx = None
            
            if self.draw_infeasibility_maps:
                Map.map_infeasibility(self.node, self.xfmr, self.regulator, self.switch, self.load, 
                                        self.line_oh, self.line_ug, self.line_tplx, self, self.fuse, Vsol)
                if self.normalize_if:
                    Map.map_normal_infeasibility(self.node, self.xfmr, self.regulator, self.switch, self.load, 
                                            self.line_oh, self.line_ug, self.line_tplx, self, self.fuse, Vsol, max_if)            
            if self.save_results:
                self.write_results(V, iteration_count, self.simulation_stats, dual_info_ABC, dual_info_tplx)

            return returncode["success"]

    def write_results(self, V, iteration_count, simulation_stats, dual_info_ABC, dual_info_tplx):

        for ele in self.node:
            ele.print_slack = True
            ele.calcMagAng(V, False)

        if self.regulator:
            for ele in self.regulator:
                ele.update_values(V)

        enable_IBDGs = True if self.ibdg else False
        if self.ibdg:
            for ele in self.ibdg:
                ele.calc_ibdg_outputs(self.node, self.node_key, V)

        for ele in self.curr_meas:
            ele.calc_currents(V)

        for ele in self.xfmr:
            ele.calc_currents(V)

        for ele in self.fuse:
            ele.calc_currents(V)

        for ele in self.switch:
            ele.calc_currents(V)

        outputs = [
            self.load, self.ibdg, self.regulator, self.curr_meas, self.xfmr, self.fuse, self.switch,
            self.triplex_load, dual_info_ABC, dual_info_tplx
        ]

        # Write csv file with results
        save_output(self.path_to_output, self.case_name, self.node, outputs, iteration_count,
                    simulation_stats, self.settings, enable_IBDGs, self.load_factor, self.stamp_dual)

    
    
    

    
    
    
