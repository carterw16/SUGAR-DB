#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 21:16:19 2021

@author: emfoster
"""

import numpy as np
from classes.Nodes import Nodes

def calc_residuals(V_out, Ylin, Jlin, nodes, node_key, load, triplex_load, lf, stamp_dual, obj, epsilon):
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
                    
                    # vr_index = [nodes[ele].nodeA1_Vr, nodes[ele].node2_Vr, nodes[ele].nodeN_Vr]
                    # vi_index = [nodes[ele].node1_Vi, nodes[ele].node2_Vi, nodes[ele].nodeN_Vi]
            
                    # sr = [nodes[ele].node1_sr, nodes[ele].node2_sr, nodes[ele].nodeN_sr]
                    # si = [nodes[ele].node1_si, nodes[ele].node2_si, nodes[ele].nodeN_si]           
    
                    for i in range(0,3):
                        # Complementary Slackness
                        res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]] * V_out[if_r_minus[i]] - epsilon
                        res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]] * V_out[if_i_minus[i]] - epsilon           
                        res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]] * V_out[if_r_plus[i]] - epsilon
                        res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]] * V_out[if_i_plus[i]] - epsilon
                        
                        # res_eqn[vr_index[i]] += - V_out[if_r_plus[i]] - V_out[if_r_minus[i]]
                        # res_eqn[vi_index[i]] += - V_out[if_i_plus[i]] - V_out[if_i_minus[i]]
                        
                        # CS code for alternative version of L1 norm (two variable boundaries)
                        # res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]]*(-V_out[sr[i]] - V_out[if_r[i]]) - epsilon
                        # res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]]*(-V_out[si[i]] - V_out[if_i[i]]) - epsilon              
                        # res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]]*(-V_out[sr[i]] + V_out[if_r[i]]) - epsilon
                        # res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]]*(-V_out[si[i]] + V_out[if_i[i]]) - epsilon   
                        
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
                    
                    # vr_index = [nodes[ele].nodeA_Vr, nodes[ele].nodeB_Vr, nodes[ele].nodeC_Vr, nodes[ele].nodeN_Vr]
                    # vi_index = [nodes[ele].nodeA_Vi, nodes[ele].nodeB_Vi, nodes[ele].nodeC_Vi, nodes[ele].nodeN_Vi]   
            
                    # sr = [nodes[ele].nodeA_sr, nodes[ele].nodeB_sr, nodes[ele].nodeC_sr, nodes[ele].nodeN_sr]
                    # si = [nodes[ele].nodeA_si, nodes[ele].nodeB_si, nodes[ele].nodeC_si, nodes[ele].nodeN_si]           
    
                    for i in range(0,4):
                        # Complementary Slackness
                        res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]] * V_out[if_r_minus[i]] - epsilon
                        res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]] * V_out[if_i_minus[i]] - epsilon           
                        res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]] * V_out[if_r_plus[i]] - epsilon
                        res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]] * V_out[if_i_plus[i]] - epsilon
                        
                        # res_eqn[vr_index[i]] += - V_out[if_r_plus[i]] - V_out[if_r_minus[i]]
                        # res_eqn[vi_index[i]] += - V_out[if_i_plus[i]] - V_out[if_i_minus[i]]
                
                        # CS code to alternative version of L1 norm (two variable boundaries)
                        # res_eqn[mu_r_lb[i]] += V_out[mu_r_lb[i]]*(-V_out[sr[i]] - V_out[if_r[i]]) - epsilon
                        # res_eqn[mu_i_lb[i]] += V_out[mu_i_lb[i]]*(-V_out[si[i]] - V_out[if_i[i]]) - epsilon              
                        # res_eqn[mu_r_ub[i]] += V_out[mu_r_ub[i]]*(-V_out[sr[i]] + V_out[if_r[i]]) - epsilon
                        # res_eqn[mu_i_ub[i]] += V_out[mu_i_ub[i]]*(-V_out[si[i]] + V_out[if_i[i]]) - epsilon
                        
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
     
    s = 12
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