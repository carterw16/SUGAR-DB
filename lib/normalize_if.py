#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 14:40:47 2021

@author: emfoster
"""
import numpy as np
from classes.Nodes import Nodes

def normalize_if(node, V, obj):
    normalized_val = np.zeros(len(V))
    max_normalized_val = np.zeros((len(node)))
    
    Sbase = 10000000
    
    
    if obj == 2:        
        for i in range(len(node)):
            if node[i].isTriplex == 0 and node[i].bustype != 3:
                if_r_index = [node[i].nodeA_if_r, node[i].nodeB_if_r, node[i].nodeC_if_r, node[i].nodeN_if_r]
                if_i_index = [node[i].nodeA_if_i, node[i].nodeB_if_i, node[i].nodeC_if_i, node[i].nodeN_if_i]   
                
                base_V = node[i].Vbase_LL*1000
                I_base = Sbase/(base_V * np.sqrt(3))
                
                normalized_val[if_r_index] = np.reshape(V[if_r_index]/I_base, -1)
                normalized_val[if_i_index] = np.reshape(V[if_i_index]/I_base, -1)           
                
                if max(abs(V[if_r_index])) > max(abs(V[if_i_index])):
                    max_normalized_val[i] = max(abs(V[if_r_index]))/I_base
                else:
                    max_normalized_val[i] = max(abs(V[if_i_index]))/I_base 
            elif node[i].isTriplex == 1 and node[i].bustype != 3:
                if_r_index = [node[i].node1_if_r, node[i].node2_if_r, node[i].nodeN_if_r]
                if_i_index = [node[i].node1_if_i, node[i].node2_if_i, node[i].nodeN_if_i]   
                
                base_V = node[i].Vbase_LL*1000
                I_base = Sbase/(base_V * np.sqrt(3))
                
                normalized_val[if_r_index] = np.reshape(V[if_r_index]/I_base, -1)
                normalized_val[if_i_index] = np.reshape(V[if_i_index]/I_base, -1)           
                
                if max(abs(V[if_r_index])) > max(abs(V[if_i_index])):
                    max_normalized_val[i] = max(abs(V[if_r_index]))/I_base
                else:
                    max_normalized_val[i] = max(abs(V[if_i_index]))/I_base                  
    
        normalized_if = abs(normalized_val[Nodes.if_index])
        
    if obj == 1:
        for i in range(len(node)):
            if node[i].isTriplex == 0 and node[i].bustype != 3:
                if_r_plus_index = [node[i].nodeA_if_r_plus,  node[i].nodeB_if_r_plus,  node[i].nodeC_if_r_plus,  node[i].nodeN_if_r_plus]        
                if_r_minus_index =[node[i].nodeA_if_r_minus, node[i].nodeB_if_r_minus, node[i].nodeC_if_r_minus, node[i].nodeN_if_r_minus]
                if_i_plus_index = [node[i].nodeA_if_i_plus,  node[i].nodeB_if_i_plus,  node[i].nodeC_if_i_plus,  node[i].nodeN_if_i_plus]   
                if_i_minus_index =[node[i].nodeA_if_i_minus, node[i].nodeB_if_i_minus, node[i].nodeC_if_i_minus, node[i].nodeN_if_i_minus] 
                
                base_V = node[i].Vbase_LL*1000
                I_base = Sbase/(base_V * np.sqrt(3))
                
                normalized_val[if_r_plus_index] = np.reshape(V[if_r_plus_index]/I_base, -1)
                normalized_val[if_i_plus_index] = np.reshape(V[if_i_plus_index]/I_base, -1)           
                normalized_val[if_r_minus_index] = np.reshape(V[if_r_minus_index]/I_base, -1)
                normalized_val[if_i_minus_index] = np.reshape(V[if_i_minus_index]/I_base, -1)   
                
                infeas_curr = np.append(V[if_r_plus_index] - V[if_r_minus_index], V[if_i_plus_index] - V[if_i_minus_index])
                max_normalized_val[i] = max(abs(infeas_curr))/I_base

            elif node[i].isTriplex == 1 and node[i].bustype != 3:
                if_r_plus_index = [node[i].node1_if_r_plus,  node[i].node2_if_r_plus,  node[i].nodeN_if_r_plus]        
                if_r_minus_index =[node[i].node1_if_r_minus, node[i].node2_if_r_minus, node[i].nodeN_if_r_minus]
                if_i_plus_index = [node[i].node1_if_i_plus,  node[i].node2_if_i_plus,  node[i].nodeN_if_i_plus]   
                if_i_minus_index =[node[i].node1_if_i_minus, node[i].node2_if_i_minus, node[i].nodeN_if_i_minus] 
               
                base_V = node[i].Vbase_LL*1000
                I_base = Sbase/(base_V * np.sqrt(3))
                
                normalized_val[if_r_plus_index] = np.reshape(V[if_r_plus_index]/I_base, -1)
                normalized_val[if_i_plus_index] = np.reshape(V[if_i_plus_index]/I_base, -1)           
                normalized_val[if_r_minus_index] = np.reshape(V[if_r_minus_index]/I_base, -1)
                normalized_val[if_i_minus_index] = np.reshape(V[if_i_minus_index]/I_base, -1)   
                
                infeas_curr = np.append(V[if_r_plus_index] - V[if_r_minus_index], V[if_i_plus_index] - V[if_i_minus_index])
                max_normalized_val[i] = max(abs(infeas_curr))/I_base                
    
        normalized_if = abs(normalized_val[Nodes.if_index])
        max_normalized_val[abs(max_normalized_val) < 1e-6] = 0
    return normalized_if, max_normalized_val
            # We can get base_power rating from xfmr
            # We can get primary/secondary voltage rating from xfmr
            # Base Current = Power Base/V Base
            # P.U. current = I/ Base Current
    