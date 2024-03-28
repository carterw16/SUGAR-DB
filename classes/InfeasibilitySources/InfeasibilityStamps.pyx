#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 09:52:55 2021

@author: emfoster

  Implements stamps for the infeasibility class.

  Author(s): Elizabeth Foster
  Created Date: 10-11-2021
  Updated Date: 10-11-2021
  Email: emfoster@andrew.cmu.edu 
  Status: Development

"""
#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3
cpdef stampY(i, j, val, Y_val, Y_row, Y_col, idx):
    Y_val[idx] = val
    Y_row[idx] = i
    Y_col[idx] = j
    idx += 1

    return idx


cpdef stampJ(i, val, J_val, J_row, idx):
    J_val[idx] = val
    J_row[idx] = i
    idx += 1

    return idx 

cpdef stamp_equality_constraint(self, node, Y_val, Y_row, Y_col, idx_Y):
    if_r = self.if_r_nodes
    if_i = self.if_i_nodes
    
    if self.isTriplex:
        node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
        node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
    else:
        node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
        node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
    
    num_of_phases = len(if_r)
    
    for i in range(0, num_of_phases):
        # Infeasibility Currents are subtracted per the equality constraints in V rows
        idx_Y = stampY(node_Vr[i], if_r[i], -1, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(node_Vi[i], if_i[i], -1, Y_val, Y_row, Y_col, idx_Y)
    
    return Y_val, Y_row, Y_col, idx_Y

cpdef stamp_stationarity_constraints(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
    if_r_nodes = self.if_r_nodes
    if_i_nodes = self.if_i_nodes

    Lr_nodes = node.dual_eq_var_r_nodes
    Li_nodes = node.dual_eq_var_i_nodes
    
    if node.isTriplex:
        node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
        node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
    else:
        node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
        node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
    
    num_of_phases = len(if_r_nodes)
    
    alpha = self.obj_scaling

    if self.obj_type == 'power':
        for i in range(0, num_of_phases):
            Vr = V[node_Vr[i]]
            Vi = V[node_Vi[i]]
            if_r = V[if_r_nodes[i]]
            if_i = V[if_i_nodes[i]]
            # Equations for primal variable If_R
            # dL/difr = 2*ifr*(Vr**2 + Vi**2) - lambda_r = 0
            idx_Y = stampY(if_r_nodes[i], if_r_nodes[i], 2*alpha*(Vr**2 + Vi**2), Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_r_nodes[i], node_Vr[i], 4*alpha*if_r*Vr, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_r_nodes[i], node_Vi[i], 4*alpha*if_r*Vi, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_r_nodes[i], Lr_nodes[i], -1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_r_nodes[i], 4*alpha*if_r*(Vr+Vi), J_val, J_row, idx_J)

            # Equations for primal variable If_R
            # dL/difr = 2*ifr*(Vr**2 + Vi**2) - lambda_r = 0
            idx_Y = stampY(if_i_nodes[i], if_i_nodes[i], 2*alpha*(Vr**2 + Vi**2), Y_val, Y_row, Y_col, idx_Y)        
            idx_Y = stampY(if_i_nodes[i], node_Vr[i], 4*alpha*if_i*Vr, Y_val, Y_row, Y_col, idx_Y) 
            idx_Y = stampY(if_i_nodes[i], node_Vi[i], 4*alpha*if_i*Vi, Y_val, Y_row, Y_col, idx_Y) 
            idx_Y = stampY(if_i_nodes[i], Li_nodes[i], -1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_i_nodes[i], 4*alpha*if_i*(Vr+Vi), J_val, J_row, idx_J)

            # Equations for primal variable Vr
            # dL/dVr = .... + 2Vr*(ifr**2 + ifi**2) = 0
            idx_Y = stampY(Lr_nodes[i], if_r_nodes[i], 4*alpha*Vr*if_r, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(Lr_nodes[i], if_i_nodes[i], 4*alpha*Vr*if_i, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(Lr_nodes[i], node_Vr[i], 2*alpha*(if_r**2 + if_i**2), Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(Lr_nodes[i], 4*alpha*Vr*(if_r + if_i), J_val, J_row, idx_J)

            # Equations for primal variable Vi
            # dL/dVi = .... + 2Vi*(ifr**2 + ifi**2) = 0
            idx_Y = stampY(Li_nodes[i], if_r_nodes[i], 4*alpha*Vi*if_r, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(Li_nodes[i], if_i_nodes[i], 4*alpha*Vi*if_i, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(Li_nodes[i], Li_nodes[i], 2*alpha*(if_r**2 + if_i**2), Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(Li_nodes[i], 4*alpha*Vi*(if_r + if_i), J_val, J_row, idx_J)
    
    return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J  

# Nonlinear Function
cpdef stamp_ineq_duals(self, V, cs_eps, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
    mu_r_lb = self.dual_ineq_r_minus_nodes
    mu_i_lb = self.dual_ineq_i_minus_nodes
    mu_r_ub = self.dual_ineq_r_plus_nodes
    mu_i_ub = self.dual_ineq_i_plus_nodes

    if_r_plus  = self.if_r_plus_nodes
    if_r_minus = self.if_r_minus_nodes
    if_i_plus  = self.if_i_plus_nodes
    if_i_minus = self.if_i_minus_nodes 

    num_of_phases = len(if_r_plus)
    
    for i in range(0, num_of_phases):
        idx_Y = stampY(mu_r_lb[i], if_r_minus[i], V[mu_r_lb[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(mu_r_lb[i], mu_r_lb[i], V[if_r_minus[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(mu_r_lb[i], V[mu_r_lb[i]]* V[if_r_minus[i]] + cs_eps, J_val, J_row, idx_J)
        
        idx_Y = stampY(mu_r_ub[i], if_r_plus[i], V[mu_r_ub[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(mu_r_ub[i], mu_r_ub[i], V[if_r_plus[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(mu_r_ub[i], V[mu_r_ub[i]]* V[if_r_plus[i]] + cs_eps, J_val, J_row, idx_J)
        
        idx_Y = stampY(mu_i_lb[i], if_i_minus[i], V[mu_i_lb[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(mu_i_lb[i], mu_i_lb[i], V[if_i_minus[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(mu_i_lb[i], V[mu_i_lb[i]]* V[if_i_minus[i]] + cs_eps, J_val, J_row, idx_J)
        
        idx_Y = stampY(mu_i_ub[i], if_i_plus[i], V[mu_i_ub[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(mu_i_ub[i], mu_i_ub[i], V[if_i_plus[i]], Y_val, Y_row, Y_col, idx_Y)
        idx_J = stampJ(mu_i_ub[i], V[mu_i_ub[i]]* V[if_i_plus[i]] + cs_eps, J_val, J_row, idx_J)

    return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J  

cpdef stamp_infeas_current_L2(self, node, Y_val, Y_row, Y_col, idx_Y):
    if_r = self.if_r_nodes
    if_i = self.if_i_nodes
    
    if node.isTriplex:
        node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
        node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
    else:
        node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
        node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
    
    num_of_phases = len(if_r)
    
    for i in range(0, num_of_phases):
        # Infeasibility Currents are subtracted per the equality constraints in V rows
        idx_Y = stampY(node_Vr[i], if_r[i], -1, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(node_Vi[i], if_i[i], -1, Y_val, Y_row, Y_col, idx_Y)
    
    return Y_val, Y_row, Y_col, idx_Y

cpdef stamp_infeas_current_L1(self, node, Y_val, Y_row, Y_col, idx_Y):
    if_r_plus  = self.if_r_plus_nodes
    if_r_minus = self.if_r_minus_nodes
    if_i_plus  = self.if_i_plus_nodes
    if_i_minus = self.if_i_minus_nodes
    
    if node.isTriplex:
        node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
        node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
    else:
        node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
        node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
    
    num_of_phases = len(if_r_plus)
    
    for i in range(0, num_of_phases):
        idx_Y = stampY(node_Vr[i], if_r_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(node_Vr[i], if_r_plus[i], -1, Y_val, Y_row, Y_col, idx_Y)  
        idx_Y = stampY(node_Vi[i], if_i_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
        idx_Y = stampY(node_Vi[i], if_i_plus[i], -1, Y_val, Y_row, Y_col, idx_Y) 
    
    return Y_val, Y_row, Y_col, idx_Y  
     
   
cpdef stamp_infeas_dual_relationship(self, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
    # This function stamps the relationship between infeasibility currents and dual variables from the dL/dif stationarity constraint    
    Lr = self.dual_eq_var_r_nodes
    Li = self.dual_eq_var_i_nodes
    
    if self.obj == 'L1':
        mu_r_lb = self.dual_ineq_r_minus_nodes
        mu_i_lb = self.dual_ineq_i_minus_nodes
        mu_r_ub = self.dual_ineq_r_plus_nodes
        mu_i_ub = self.dual_ineq_i_plus_nodes

        if_r_plus  = self.if_r_plus_nodes
        if_r_minus = self.if_r_minus_nodes
        if_i_plus  = self.if_i_plus_nodes
        if_i_minus = self.if_i_minus_nodes

        num_of_phases = len(mu_r_lb)

        for i in range(0, num_of_phases):
            idx_Y = stampY(if_r_plus[i], Lr[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_r_plus[i], mu_r_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_r_plus[i], 1, J_val, J_row, idx_J)
 
            idx_Y = stampY(if_r_minus[i], Lr[i], -1, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_r_minus[i], mu_r_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_r_minus[i], 1, J_val, J_row, idx_J)
            
            idx_Y = stampY(if_i_plus[i], Li[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_i_plus[i], mu_i_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_i_plus[i], 1, J_val, J_row, idx_J)

            idx_Y = stampY(if_i_minus[i], Li[i], -1, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_i_minus[i], mu_i_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
            idx_J = stampJ(if_i_minus[i], 1, J_val, J_row, idx_J)  
    
    elif self.obj == 'L2' and self.obj_type == 'current':
        if_r = self.if_r_nodes
        if_i = self.if_i_nodes

        num_of_phases = len(if_r)

        for i in range(0, num_of_phases):
            # dL/dif = if - lambda = 0
            idx_Y = stampY(if_r[i], if_r[i], 2*self.obj_scaling, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_i[i], if_i[i], 2*self.obj_scaling, Y_val, Y_row, Y_col, idx_Y)        

            idx_Y = stampY(if_r[i], Lr[i], -1, Y_val, Y_row, Y_col, idx_Y)
            idx_Y = stampY(if_i[i], Li[i], -1, Y_val, Y_row, Y_col, idx_Y)
    
    return Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J

cpdef stamp_linear_infeasibility(self, node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
    # this function stamps all of the linear elements; it is called from the InfeasibilityAnalysis class
    if self.obj_type == 'current':
        Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J = self.stamp_infeas_dual_relationship(Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)                
  
        if self.obj == 'L2':
            Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L2(node, Y_val, Y_row, Y_col, idx_Y)
        
        if self.obj == 'L1':
            Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L1(node, Y_val, Y_row, Y_col, idx_Y)
            
    return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

cpdef stamp_nonlinear_infeasibility(self, node, V, cs_eps, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
    # this function stamps all of the nonlinear elements (from complementary slackness) and is called from the InfeasibilityAnalysis class
    if self.obj == 'L1' and self.source_type == 'current':
        Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_ineq_duals(V, cs_eps, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)

    if self.obj_type == 'power':
        Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_equality_constraint(node, Y_val, Y_row, Y_col, idx_Y)
        Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_stationarity_constraints(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
    
    return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

    
    
    
