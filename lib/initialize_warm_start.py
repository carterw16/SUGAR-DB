#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 11:27:34 2021

@author: emfoster
"""
import pickle
import numpy as np

def initialize_from_amrit_solution(V, node_list):
    with open('Vsol_R4-12.47-1.pkl', 'rb') as f:
        data = pickle.load(f)
    for ele in node_list:
        if ele.name not in list(data.keys()) or ele.bustype == 3:
            continue
        if ele.isTriplex:
            voltage_inds = [ele.node1_Vr, ele.node1_Vi, ele.node2_Vr, ele.node2_Vi]
            lambda_inds = [ele.node1_dual_eq_var_r, ele.node1_dual_eq_var_i, ele.node2_dual_eq_var_r, ele.node2_dual_eq_var_i]
            current_inds = [ele.node1_if_r, ele.node1_if_i, ele.node2_if_r, ele.node2_if_i]
            vmag2_inds = [ele.node1_vmag2_index, ele.node2_vmag2_index]
            Lvmag2_inds = [ele.node1_Lvmag2_index, ele.node2_Lvmag2_index]
            umax_inds = [ele.node1_umax_vmag2_index, ele.node2_umax_vmag2_index]
            umin_inds = [ele.node1_umin_vmag2_index, ele.node2_umin_vmag2_index]
        else:
            voltage_inds = [ele.nodeA_Vr, ele.nodeA_Vi, ele.nodeB_Vr, ele.nodeB_Vi, ele.nodeC_Vr, ele.nodeC_Vi]
            lambda_inds = [ele.nodeA_dual_eq_var_r, ele.nodeA_dual_eq_var_i, ele.nodeB_dual_eq_var_r, ele.nodeB_dual_eq_var_i, ele.nodeC_dual_eq_var_r, ele.nodeC_dual_eq_var_i]
            current_inds = [ele.nodeA_if_r, ele.nodeA_if_i, ele.nodeB_if_r, ele.nodeB_if_i, ele.nodeC_if_r, ele.nodeC_if_i]
            vmag2_inds = [ele.nodeA_vmag2_index, ele.nodeB_vmag2_index, ele.nodeC_vmag2_index]
            Lvmag2_inds = [ele.nodeA_Lvmag2_index, ele.nodeB_Lvmag2_index, ele.nodeC_Lvmag2_index]
            umax_inds = [ele.nodeA_umax_vmag2_index, ele.nodeB_umax_vmag2_index, ele.nodeC_umax_vmag2_index]
            umin_inds = [ele.nodeA_umin_vmag2_index, ele.nodeB_umin_vmag2_index, ele.nodeC_umin_vmag2_index]
        V[voltage_inds,0] = data[ele.name]['voltages']
        V[lambda_inds,0] = data[ele.name]['lambdas']
        V[current_inds,0] = -1*np.array(data[ele.name]['feas_currents'])
        V[vmag2_inds,0] = data[ele.name]['vmag2s']
        V[Lvmag2_inds,0] = data[ele.name]['Lvmag2s']
        V[umax_inds,0] = data[ele.name]['umax']
        V[umin_inds,0] = data[ele.name]['umin']
    return V

def initialize_warm_start(Vinit, casename, lf, node, warm_start_settings):
    obj_old = warm_start_settings['initial solution type']    
    lf_old = warm_start_settings['previous lf']
    obj_new = warm_start_settings['new objective']
    old_filename = warm_start_settings['solution file path']
    
    
    if obj_old == 1:
        name_node = old_filename + "_L1node_" + float(lf_old) + ".pkl"
        name_V = old_filename + "_L1_V_" + float(lf_old) + ".npy"

        with open(name_node, 'rb') as f:
            old_node = pickle.load(f)
        
        old_V = np.load(name_V)
        Vinit = old_V
        Vinit[abs(Vinit) == 0] = 1e-9

    if obj_old == 2:
        name_node = old_filename + "_L2node_" + float(lf_old) + ".pkl"
        name_V = old_filename + "_L2_V_" + float(lf_old) + ".npy"
        
        if obj_new == 2:
            old_V = np.load(name_V)
            Vinit = old_V
        else:
            with open(name_node, 'rb') as f:
                old_node = pickle.load(f)
            
            old_V = np.load(name_V)
            
            for i in range(1,len(node)):
                if node[i].isTriplex == 0:
                    Vinit[node[i].nodeA_Vr] = old_V[old_node[i].nodeA_Vr]
                    Vinit[node[i].nodeB_Vr] = old_V[old_node[i].nodeB_Vr]
                    Vinit[node[i].nodeC_Vr] = old_V[old_node[i].nodeC_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].nodeA_Vi] = old_V[old_node[i].nodeA_Vi]
                    Vinit[node[i].nodeB_Vi] = old_V[old_node[i].nodeB_Vi]
                    Vinit[node[i].nodeC_Vi] = old_V[old_node[i].nodeC_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                        if_r_old = [old_V[old_node[i].nodeA_if_r], old_V[old_node[i].nodeB_if_r], old_V[old_node[i].nodeC_if_r]]
                        if_i_old = [old_V[old_node[i].nodeA_if_i], old_V[old_node[i].nodeB_if_i], old_V[old_node[i].nodeC_if_i]]
                        
                        if_r_plus_node = [node[i].nodeA_if_r_plus, node[i].nodeB_if_r_plus, node[i].nodeC_if_r_plus]
                        if_i_plus_node = [node[i].nodeA_if_i_plus, node[i].nodeB_if_i_plus, node[i].nodeC_if_i_plus]
                        if_r_minus_node = [node[i].nodeA_if_r_minus, node[i].nodeB_if_r_minus, node[i].nodeC_if_r_minus]
                        if_i_minus_node = [node[i].nodeA_if_i_minus, node[i].nodeB_if_i_minus, node[i].nodeC_if_i_minus]
                        
                        mu_r_lb = [node[i].nodeA_mu_r_lb, node[i].nodeB_mu_r_lb, node[i].nodeC_mu_r_lb]
                        mu_r_ub = [node[i].nodeA_mu_r_ub, node[i].nodeB_mu_r_ub, node[i].nodeC_mu_r_ub]
                        mu_i_lb = [node[i].nodeA_mu_i_lb, node[i].nodeB_mu_i_lb, node[i].nodeC_mu_i_lb]
                        mu_i_ub = [node[i].nodeA_mu_i_ub, node[i].nodeB_mu_i_ub, node[i].nodeC_mu_i_ub]
                        
                        Lr = [node[i].nodeA_dual_eq_var_r, node[i].nodeB_dual_eq_var_r, node[i].nodeC_dual_eq_var_r]
                        Li = [node[i].nodeA_dual_eq_var_i, node[i].nodeB_dual_eq_var_i, node[i].nodeC_dual_eq_var_i] 
                        
                        old_Lr = [old_node[i].nodeA_dual_eq_var_r, old_node[i].nodeB_dual_eq_var_r, old_node[i].nodeC_dual_eq_var_r]
                        old_Li = [old_node[i].nodeA_dual_eq_var_i, old_node[i].nodeB_dual_eq_var_i, old_node[i].nodeC_dual_eq_var_i] 
                        for phase in range(0,3):
                            Vinit[Lr[phase]] = old_V[old_Lr[phase]]
                            Vinit[Li[phase]] = old_V[old_Li[phase]]
                            
                            if if_r_old[phase] > 1e-3:
                                Vinit[if_r_plus_node[phase]] = if_r_old[phase]
                                Vinit[if_r_minus_node[phase]] = 1e-10
                                Vinit[mu_r_ub] = 1e-6
                            elif -if_r_old[phase] > 1e-3:
                                Vinit[if_r_plus_node[phase]] = 1e-10
                                Vinit[if_r_minus_node[phase]] = abs(if_r_old[phase])  
                                Vinit[mu_r_lb] = 1e-6
                            else:
                                Vinit[if_r_plus_node[phase]] = 1e-10
                                Vinit[if_r_minus_node[phase]] = 1e-10
                                Vinit[mu_r_lb] = 1e-6      
                            if if_i_old[phase] > 1e-3:
                                Vinit[if_i_plus_node[phase]] = if_i_old[phase]
                                Vinit[if_i_minus_node[phase]] = 1e-10
                                Vinit[mu_i_ub] = 1e-6
                            elif -if_i_old[phase] > 1e-3 :
                                Vinit[if_i_plus_node[phase]] = 1e-10
                                Vinit[if_i_minus_node[phase]] = abs(if_i_old[phase])
                                Vinit[mu_i_lb] = 1e-6   
                            else:
                                Vinit[if_i_plus_node[phase]] = 1e-10
                                Vinit[if_i_minus_node[phase]] = 1e-10
                                Vinit[mu_i_lb] = 1e-6  
                                
                else:
                    Vinit[node[i].node1_Vr] = old_V[old_node[i].node1_Vr]
                    Vinit[node[i].node2_Vr] = old_V[old_node[i].node2_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].node1_Vi] = old_V[old_node[i].node1_Vi]
                    Vinit[node[i].node2_Vi] = old_V[old_node[i].node2_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                        if_r_old = [old_V[old_node[i].node1_if_r], old_V[old_node[i].node2_if_r], old_V[old_node[i].nodeN_if_r]]
                        if_i_old = [old_V[old_node[i].node1_if_i], old_V[old_node[i].node2_if_i], old_V[old_node[i].nodeN_if_i]]
                        
                        if_r_plus_node =  [node[i].node1_if_r_plus, node[i].node2_if_r_plus, node[i].nodeN_if_r_plus]
                        if_i_plus_node =  [node[i].node1_if_i_plus, node[i].node2_if_i_plus, node[i].nodeN_if_i_plus]
                        if_r_minus_node = [node[i].node1_if_r_minus, node[i].node2_if_r_minus, node[i].nodeN_if_r_minus]
                        if_i_minus_node = [node[i].node1_if_i_minus, node[i].node2_if_i_minus, node[i].nodeN_if_i_minus]
                        
                        mu_r_lb = [node[i].node1_mu_r_lb, node[i].node2_mu_r_lb, node[i].nodeN_mu_r_lb]
                        mu_r_ub = [node[i].node1_mu_r_ub, node[i].node2_mu_r_ub, node[i].nodeN_mu_r_ub]
                        mu_i_lb = [node[i].node1_mu_i_lb, node[i].node2_mu_i_lb, node[i].nodeN_mu_i_lb]
                        mu_i_ub = [node[i].node1_mu_i_ub, node[i].node2_mu_i_ub, node[i].nodeN_mu_i_ub]
                        
                        
                        Lr = [node[i].node1_dual_eq_var_r, node[i].node2_dual_eq_var_r, node[i].nodeN_dual_eq_var_r]
                        Li = [node[i].node1_dual_eq_var_i, node[i].node2_dual_eq_var_i, node[i].nodeN_dual_eq_var_i] 
                        
                        old_Lr = [old_node[i].node1_dual_eq_var_r, old_node[i].node2_dual_eq_var_r, old_node[i].nodeN_dual_eq_var_r]
                        old_Li = [old_node[i].node1_dual_eq_var_i, old_node[i].node2_dual_eq_var_i, old_node[i].nodeN_dual_eq_var_i]
                        
                        for phase in range(0,3):
                            Vinit[Lr[phase]] = old_V[old_Lr[phase]]
                            Vinit[Li[phase]] = old_V[old_Li[phase]]
                            
                            if if_r_old[phase] > 1e-6:
                                Vinit[if_r_plus_node[phase]] = if_r_old[phase]
                                Vinit[if_r_minus_node[phase]] = 1e-10
                                Vinit[mu_r_ub] = 1e-6
                            elif -if_r_old[phase] > 1e-6:
                                Vinit[if_r_plus_node[phase]] = 1e-10
                                Vinit[if_r_minus_node[phase]] = abs(if_r_old[phase])  
                                Vinit[mu_r_lb] = 1e-6
                            else:
                                Vinit[if_r_plus_node[phase]] = 1e-10
                                Vinit[if_r_minus_node[phase]] = 1e-10
                                Vinit[mu_r_lb] = 1e-6      
                            if if_i_old[phase] > 1e-6:
                                Vinit[if_i_plus_node[phase]] = if_i_old[phase]
                                Vinit[if_i_minus_node[phase]] = 1e-10
                                Vinit[mu_i_ub] = 1e-6
                            elif -if_i_old[phase] > 1e-6 :
                                Vinit[if_i_plus_node[phase]] = 1e-10
                                Vinit[if_i_minus_node[phase]] = abs(if_i_old[phase])
                                Vinit[mu_i_lb] = 1e-6   
                            else:
                                Vinit[if_i_plus_node[phase]] = 1e-10
                                Vinit[if_i_minus_node[phase]] = 1e-10
                                Vinit[mu_i_lb] = 1e-6  
                            
    if obj_old == 0:
        name_node = old_filename + "_pf_node_" + float(lf_old) + ".pkl"
        name_V = old_filename + "_pf_V_" + float(lf_old) + ".npy"

        if obj_new == 0:
            Vinit = np.load(name_V)

        if obj_new == 1:
            with open(name_node, 'rb') as f:
                old_node = pickle.load(f)
            
            old_V = np.load(name_V)
            for i in range(1,len(node)):
                if node[i].isTriplex == 0:
                    Vinit[node[i].nodeA_Vr] = old_V[old_node[i].nodeA_Vr]
                    Vinit[node[i].nodeB_Vr] = old_V[old_node[i].nodeB_Vr]
                    Vinit[node[i].nodeC_Vr] = old_V[old_node[i].nodeC_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].nodeA_Vi] = old_V[old_node[i].nodeA_Vi]
                    Vinit[node[i].nodeB_Vi] = old_V[old_node[i].nodeB_Vi]
                    Vinit[node[i].nodeC_Vi] = old_V[old_node[i].nodeC_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                        if_r_plus_node = [node[i].nodeA_if_r_plus, node[i].nodeB_if_r_plus, node[i].nodeC_if_r_plus]
                        if_i_plus_node = [node[i].nodeA_if_i_plus, node[i].nodeB_if_i_plus, node[i].nodeC_if_i_plus]
                        if_r_minus_node = [node[i].nodeA_if_r_minus, node[i].nodeB_if_r_minus, node[i].nodeC_if_r_minus]
                        if_i_minus_node = [node[i].nodeA_if_i_minus, node[i].nodeB_if_i_minus, node[i].nodeC_if_i_minus]
                        
                        mu_r_lb = [node[i].nodeA_mu_r_lb, node[i].nodeB_mu_r_lb, node[i].nodeC_mu_r_lb]
                        mu_r_ub = [node[i].nodeA_mu_r_ub, node[i].nodeB_mu_r_ub, node[i].nodeC_mu_r_ub]
                        mu_i_lb = [node[i].nodeA_mu_i_lb, node[i].nodeB_mu_i_lb, node[i].nodeC_mu_i_lb]
                        mu_i_ub = [node[i].nodeA_mu_i_ub, node[i].nodeB_mu_i_ub, node[i].nodeC_mu_i_ub]
                        
                        Lr = [node[i].nodeA_dual_eq_var_r, node[i].nodeB_dual_eq_var_r, node[i].nodeC_dual_eq_var_r]
                        Li = [node[i].nodeA_dual_eq_var_i, node[i].nodeB_dual_eq_var_i, node[i].nodeC_dual_eq_var_i] 
                        
                        Vinit[if_r_plus_node] = 1e-9
                        Vinit[if_i_plus_node] = 1e-9
                        Vinit[if_r_minus_node] = 1e-9
                        Vinit[if_i_minus_node] = 1e-9
                        Vinit[Lr] = 1
                        Vinit[Li] = 1
                        Vinit[mu_r_lb] = 1e-6
                        Vinit[mu_i_lb] = 1e-6
                        Vinit[mu_r_ub] = 1e-6
                        Vinit[mu_i_ub] = 1e-6
                        
                                
                else:
                    Vinit[node[i].node1_Vr] = old_V[old_node[i].node1_Vr]
                    Vinit[node[i].node2_Vr] = old_V[old_node[i].node2_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].node1_Vi] = old_V[old_node[i].node1_Vi]
                    Vinit[node[i].node2_Vi] = old_V[old_node[i].node2_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                    
                        if_r_plus_node =  [node[i].node1_if_r_plus, node[i].node2_if_r_plus, node[i].nodeN_if_r_plus]
                        if_i_plus_node =  [node[i].node1_if_i_plus, node[i].node2_if_i_plus, node[i].nodeN_if_i_plus]
                        if_r_minus_node = [node[i].node1_if_r_minus, node[i].node2_if_r_minus, node[i].nodeN_if_r_minus]
                        if_i_minus_node = [node[i].node1_if_i_minus, node[i].node2_if_i_minus, node[i].nodeN_if_i_minus]
                        
                        mu_r_lb = [node[i].node1_mu_r_lb, node[i].node2_mu_r_lb, node[i].nodeN_mu_r_lb]
                        mu_r_ub = [node[i].node1_mu_r_ub, node[i].node2_mu_r_ub, node[i].nodeN_mu_r_ub]
                        mu_i_lb = [node[i].node1_mu_i_lb, node[i].node2_mu_i_lb, node[i].nodeN_mu_i_lb]
                        mu_i_ub = [node[i].node1_mu_i_ub, node[i].node2_mu_i_ub, node[i].nodeN_mu_i_ub]
                        
                        Lr = [node[i].node1_dual_eq_var_r, node[i].node2_dual_eq_var_r, node[i].nodeN_dual_eq_var_r]
                        Li = [node[i].node1_dual_eq_var_i, node[i].node2_dual_eq_var_i, node[i].nodeN_dual_eq_var_i] 
                        
                        Vinit[if_r_plus_node] = 1e-9
                        Vinit[if_i_plus_node] = 1e-9
                        Vinit[if_r_minus_node] = 1e-9
                        Vinit[if_i_minus_node] = 1e-9
                        Vinit[Lr] = 1
                        Vinit[Li] = 1
                        Vinit[mu_r_lb] = 1e-6
                        Vinit[mu_i_lb] = 1e-6
                        Vinit[mu_r_ub] = 1e-6
                        Vinit[mu_i_ub] = 1e-6

        if obj_new == 2:
            with open(name_node, 'rb') as f:
                old_node = pickle.load(f)
            
            old_V = np.load(name_V)
            for i in range(1,len(node)):
                if node[i].isTriplex == 0:
                    Vinit[node[i].nodeA_Vr] = old_V[old_node[i].nodeA_Vr]
                    Vinit[node[i].nodeB_Vr] = old_V[old_node[i].nodeB_Vr]
                    Vinit[node[i].nodeC_Vr] = old_V[old_node[i].nodeC_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].nodeA_Vi] = old_V[old_node[i].nodeA_Vi]
                    Vinit[node[i].nodeB_Vi] = old_V[old_node[i].nodeB_Vi]
                    Vinit[node[i].nodeC_Vi] = old_V[old_node[i].nodeC_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                        if_r_node = [node[i].nodeA_if_r, node[i].nodeB_if_r, node[i].nodeC_if_r]
                        if_i_node = [node[i].nodeA_if_i, node[i].nodeB_if_i, node[i].nodeC_if_i]
                        
                        Lr = [node[i].nodeA_dual_eq_var_r, node[i].nodeB_dual_eq_var_r, node[i].nodeC_dual_eq_var_r]
                        Li = [node[i].nodeA_dual_eq_var_i, node[i].nodeB_dual_eq_var_i, node[i].nodeC_dual_eq_var_i]  
                        
                        Vinit[if_r_node] = 1e-9
                        Vinit[if_i_node] = 1e-9
                        Vinit[Lr] = 1
                        Vinit[Li] = 1
                        
                                
                else:
                    Vinit[node[i].node1_Vr] = old_V[old_node[i].node1_Vr]
                    Vinit[node[i].node2_Vr] = old_V[old_node[i].node2_Vr]
                    Vinit[node[i].nodeN_Vr] = old_V[old_node[i].nodeN_Vr]
                    Vinit[node[i].node1_Vi] = old_V[old_node[i].node1_Vi]
                    Vinit[node[i].node2_Vi] = old_V[old_node[i].node2_Vi]
                    Vinit[node[i].nodeN_Vi] = old_V[old_node[i].nodeN_Vi]
                    
                    if node[i].bustype != 3:                
                    
                        if_r_node =  [node[i].node1_if_r, node[i].node2_if_r, node[i].nodeN_if_r]
                        if_i_node =  [node[i].node1_if_i, node[i].node2_if_i, node[i].nodeN_if_i]
                        
                        Lr = [node[i].node1_dual_eq_var_r, node[i].node2_dual_eq_var_r, node[i].nodeN_dual_eq_var_r]
                        Li = [node[i].node1_dual_eq_var_i, node[i].node2_dual_eq_var_i, node[i].nodeN_dual_eq_var_i] 
                        
                        Vinit[if_r_node] = 1e-9
                        Vinit[if_i_node] = 1e-9
                        Vinit[Lr] = 1
                        Vinit[Li] = 1
    return Vinit