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

class InfeasibilitySources:
    
    def __init__(self, node):
        '''

        Parameters
        ----------
        stamp_dual : boolean variable that flags whether or not to incorporate dual variables
        obj : integer value of either 1 for the L1-norm or 2 for the L2-norm that denotes which sets of equations to stamp

        Returns
        -------
        None.

        '''
        
        self.connected_node_name = node.name
        self.assign_dual_equality_nodes(node)
    
    def assign_dual_equality_nodes(self, node):
        if node.isTriplex:
            self.node1_dual_eq_var_r = node.node1_dual_eq_var_r
            self.node2_dual_eq_var_r = node.node2_dual_eq_var_r
            self.nodeN_dual_eq_var_r = node.nodeN_dual_eq_var_r

            self.node1_dual_eq_var_i = node.node1_dual_eq_var_i
            self.node2_dual_eq_var_i = node.node2_dual_eq_var_i
            self.nodeN_dual_eq_var_i = node.nodeN_dual_eq_var_i

            self.dual_eq_var_r_nodes = [self.node1_dual_eq_var_r, self.node2_dual_eq_var_r, self.nodeN_dual_eq_var_r]
            self.dual_eq_var_i_nodes = [self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i]
        else:
            self.nodeA_dual_eq_var_r = node.nodeA_dual_eq_var_r
            self.nodeB_dual_eq_var_r = node.nodeB_dual_eq_var_r
            self.nodeC_dual_eq_var_r = node.nodeC_dual_eq_var_r
            self.nodeN_dual_eq_var_r = node.nodeN_dual_eq_var_r

            self.nodeA_dual_eq_var_i = node.nodeA_dual_eq_var_i
            self.nodeB_dual_eq_var_i = node.nodeB_dual_eq_var_i
            self.nodeC_dual_eq_var_i = node.nodeC_dual_eq_var_i
            self.nodeN_dual_eq_var_i = node.nodeN_dual_eq_var_i

            self.dual_eq_var_r_nodes = [self.nodeA_dual_eq_var_r, self.nodeB_dual_eq_var_r, self.nodeC_dual_eq_var_r, self.nodeN_dual_eq_var_r]
            self.dual_eq_var_i_nodes = [self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_i, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_i]

        

    
    
    

    
    
    
