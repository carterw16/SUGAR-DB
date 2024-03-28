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

	return Y_val, Y_row, Y_col, idx


cpdef stampJ(i, val, J_val, J_row, idx):
	J_val[idx] = val
	J_row[idx] = i
	idx += 1

	return J_val, J_row, idx
	
# Nonlinear Function
cpdef stamp_ineq_duals(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	mu_r_lb = [node.nodeA_mu_r_lb, node.nodeB_mu_r_lb, node.nodeC_mu_r_lb, node.nodeN_mu_r_lb]
	mu_i_lb = [node.nodeA_mu_i_lb, node.nodeB_mu_i_lb, node.nodeC_mu_i_lb, node.nodeN_mu_i_lb]  
	mu_r_ub = [node.nodeA_mu_r_ub, node.nodeB_mu_r_ub, node.nodeC_mu_r_ub, node.nodeN_mu_r_ub]
	mu_i_ub = [node.nodeA_mu_i_ub, node.nodeB_mu_i_ub, node.nodeC_mu_i_ub, node.nodeN_mu_i_ub]
	
	if_r_plus =  [node.nodeA_if_r_plus,  node.nodeB_if_r_plus,  node.nodeC_if_r_plus,  node.nodeN_if_r_plus]        
	if_r_minus = [node.nodeA_if_r_minus, node.nodeB_if_r_minus, node.nodeC_if_r_minus, node.nodeN_if_r_minus]
	if_i_plus =  [node.nodeA_if_i_plus,  node.nodeB_if_i_plus,  node.nodeC_if_i_plus,  node.nodeN_if_i_plus]   
	if_i_minus = [node.nodeA_if_i_minus, node.nodeB_if_i_minus, node.nodeC_if_i_minus, node.nodeN_if_i_minus]  
	
	for i in range(0,4):
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_lb[i], if_r_minus[i], V[mu_r_lb[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_lb[i], mu_r_lb[i], V[if_r_minus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_r_lb[i], V[mu_r_lb[i]]* V[if_r_minus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_ub[i], if_r_plus[i], V[mu_r_ub[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_ub[i], mu_r_ub[i], V[if_r_plus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_r_ub[i], V[mu_r_ub[i]]* V[if_r_plus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_lb[i], if_i_minus[i], V[mu_i_lb[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_lb[i], mu_i_lb[i], V[if_i_minus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_i_lb[i], V[mu_i_lb[i]]* V[if_i_minus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_ub[i], if_i_plus[i], V[mu_i_ub[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_ub[i], mu_i_ub[i], V[if_i_plus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_i_ub[i], V[mu_i_ub[i]]* V[if_i_plus[i]] + self.cs_eps, J_val, J_row, idx_J)


	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J  

# Nonlinear Function
cpdef stamp_ineq_duals_tplx(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	mu_r_lb = [node.node1_mu_r_lb, node.node2_mu_r_lb, node.nodeN_mu_r_lb]
	mu_i_lb = [node.node1_mu_i_lb, node.node2_mu_i_lb, node.nodeN_mu_i_lb]  
	mu_r_ub = [node.node1_mu_r_ub, node.node2_mu_r_ub, node.nodeN_mu_r_ub]
	mu_i_ub = [node.node1_mu_i_ub, node.node2_mu_i_ub, node.nodeN_mu_i_ub]
	
	if_r_plus =  [node.node1_if_r_plus,  node.node2_if_r_plus,  node.nodeN_if_r_plus]        
	if_r_minus = [node.node1_if_r_minus, node.node2_if_r_minus, node.nodeN_if_r_minus]
	if_i_plus =  [node.node1_if_i_plus,  node.node2_if_i_plus,  node.nodeN_if_i_plus]   
	if_i_minus = [node.node1_if_i_minus, node.node2_if_i_minus, node.nodeN_if_i_minus]  
	
	for i in range(0,3):
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_lb[i], if_r_minus[i], V[mu_r_lb[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_lb[i], mu_r_lb[i], V[if_r_minus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_r_lb[i], V[mu_r_lb[i]]* V[if_r_minus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_ub[i], if_r_plus[i], V[mu_r_ub[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_r_ub[i], mu_r_ub[i], V[if_r_plus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_r_ub[i], V[mu_r_ub[i]]* V[if_r_plus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_lb[i], if_i_minus[i], V[mu_i_lb[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_lb[i], mu_i_lb[i], V[if_i_minus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_i_lb[i], V[mu_i_lb[i]]* V[if_i_minus[i]] + self.cs_eps, J_val, J_row, idx_J)
		
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_ub[i], if_i_plus[i], V[mu_i_ub[i]], Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(mu_i_ub[i], mu_i_ub[i], V[if_i_plus[i]], Y_val, Y_row, Y_col, idx_Y)
		J_val, J_row, idx_J = stampJ(mu_i_ub[i], V[mu_i_ub[i]]* V[if_i_plus[i]] + self.cs_eps, J_val, J_row, idx_J)


	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J  

cpdef stamp_infeas_current_L2(self, node, Y_val, Y_row, Y_col, idx_Y):
	# Infeasibility Currents are subtracted per the equality constraints in V rows
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_Vr, node.nodeA_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_Vi, node.nodeA_if_i, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_Vr, node.nodeB_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_Vi, node.nodeB_if_i, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_Vr, node.nodeC_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_Vi, node.nodeC_if_i, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_Vr, node.nodeN_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_Vi, node.nodeN_if_i, -1, Y_val, Y_row, Y_col, idx_Y) 
	
	return Y_val, Y_row, Y_col, idx_Y

cpdef stamp_infeas_current_L2_tplx(self, node, Y_val, Y_row, Y_col, idx_Y):
	# Infeasibility Currents are subtracted per the equality constraints in V rows for triplex nodes
	Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_Vr, node.node1_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_Vi, node.node1_if_i, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_Vr, node.node2_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_Vi, node.node2_if_i, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_Vr, node.nodeN_if_r, -1, Y_val, Y_row, Y_col, idx_Y)
	Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_Vi, node.nodeN_if_i, -1, Y_val, Y_row, Y_col, idx_Y) 
	
	return Y_val, Y_row, Y_col, idx_Y
   
cpdef stamp_infeas_current_L1(self, node, Y_val, Y_row, Y_col, idx_Y):
	if_r_plus =  [node.nodeA_if_r_plus,  node.nodeB_if_r_plus,  node.nodeC_if_r_plus,  node.nodeN_if_r_plus]        
	if_r_minus = [node.nodeA_if_r_minus, node.nodeB_if_r_minus, node.nodeC_if_r_minus, node.nodeN_if_r_minus]
	if_i_plus =  [node.nodeA_if_i_plus,  node.nodeB_if_i_plus,  node.nodeC_if_i_plus,  node.nodeN_if_i_plus]   
	if_i_minus = [node.nodeA_if_i_minus, node.nodeB_if_i_minus, node.nodeC_if_i_minus, node.nodeN_if_i_minus] 
	
	node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
	node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
	
	for i in range(0,4):
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vr[i], if_r_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vr[i], if_r_plus[i], -1, Y_val, Y_row, Y_col, idx_Y)  
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vi[i], if_i_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vi[i], if_i_plus[i], -1, Y_val, Y_row, Y_col, idx_Y) 
	
	return Y_val, Y_row, Y_col, idx_Y      

cpdef stamp_infeas_current_L1_tplx(self, node, Y_val, Y_row, Y_col, idx_Y):
	if_r_plus =  [node.node1_if_r_plus,  node.node2_if_r_plus,  node.nodeN_if_r_plus]        
	if_r_minus = [node.node1_if_r_minus, node.node2_if_r_minus, node.nodeN_if_r_minus]
	if_i_plus =  [node.node1_if_i_plus,  node.node2_if_i_plus,  node.nodeN_if_i_plus]   
	if_i_minus = [node.node1_if_i_minus, node.node2_if_i_minus, node.nodeN_if_i_minus] 
	
	node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
	node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
	
	for i in range(0,3):
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vr[i], if_r_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vr[i], if_r_plus[i], -1, Y_val, Y_row, Y_col, idx_Y)  
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vi[i], if_i_minus[i], 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node_Vi[i], if_i_plus[i], -1, Y_val, Y_row, Y_col, idx_Y) 
	
	return Y_val, Y_row, Y_col, idx_Y  
	 
   
cpdef stamp_infeas_dual_relationship(self, node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	# This function stamps the relationship between infeasibility currents and dual variables from the dL/dif stationarity constraint
	if self.obj == 1:
		# L1-norm
		mu_r_lb = [node.nodeA_mu_r_lb, node.nodeB_mu_r_lb, node.nodeC_mu_r_lb, node.nodeN_mu_r_lb]
		mu_i_lb = [node.nodeA_mu_i_lb, node.nodeB_mu_i_lb, node.nodeC_mu_i_lb, node.nodeN_mu_i_lb]  
		mu_r_ub = [node.nodeA_mu_r_ub, node.nodeB_mu_r_ub, node.nodeC_mu_r_ub, node.nodeN_mu_r_ub]
		mu_i_ub = [node.nodeA_mu_i_ub, node.nodeB_mu_i_ub, node.nodeC_mu_i_ub, node.nodeN_mu_i_ub]
		
		if_r_plus =  [node.nodeA_if_r_plus,  node.nodeB_if_r_plus,  node.nodeC_if_r_plus,  node.nodeN_if_r_plus]        
		if_r_minus = [node.nodeA_if_r_minus, node.nodeB_if_r_minus, node.nodeC_if_r_minus, node.nodeN_if_r_minus]
		if_i_plus =  [node.nodeA_if_i_plus,  node.nodeB_if_i_plus,  node.nodeC_if_i_plus,  node.nodeN_if_i_plus]   
		if_i_minus = [node.nodeA_if_i_minus, node.nodeB_if_i_minus, node.nodeC_if_i_minus, node.nodeN_if_i_minus]  
		
		Lr = [node.nodeA_Lr, node.nodeB_Lr, node.nodeC_Lr, node.nodeN_Lr]
		Li = [node.nodeA_Li, node.nodeB_Li, node.nodeC_Li, node.nodeN_Li]
		
		for i in range(0,4):
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_plus[i], Lr[i], 1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_plus[i], mu_r_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_r_plus[i], 1, J_val, J_row, idx_J)
 
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_minus[i], Lr[i], -1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_minus[i], mu_r_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_r_minus[i], 1, J_val, J_row, idx_J)
			
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_plus[i], Li[i], 1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_plus[i], mu_i_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_i_plus[i], 1, J_val, J_row, idx_J)

			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_minus[i], Li[i], -1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_minus[i], mu_i_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_i_minus[i], 1, J_val, J_row, idx_J)  
	
	elif self.obj == 2:
		# L2-norm
	   
		# Relationship between infeasibility currents and lambda
		# dL/dif = if - lambda = 0
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_if_r, node.nodeA_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_if_i, node.nodeA_if_i, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_if_r, node.nodeB_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_if_i, node.nodeB_if_i, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_if_r, node.nodeC_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_if_i, node.nodeC_if_i, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_r, node.nodeN_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_i, node.nodeN_if_i, 1, Y_val, Y_row, Y_col, idx_Y) 
		
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_if_r, node.nodeA_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeA_if_i, node.nodeA_Li, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_if_r, node.nodeB_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeB_if_i, node.nodeB_Li, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_if_r, node.nodeC_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeC_if_i, node.nodeC_Li, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_r, node.nodeN_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_i, node.nodeN_Li, -1, Y_val, Y_row, Y_col, idx_Y)  
	
	return Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J

cpdef stamp_infeas_dual_relationship_tplx(self, node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	 # This function stamps the relationship between infeasibility currents and dual variables from the dL/dif stationarity constraint for triplex nodes
	if self.obj == 1:
		# L1-norm
		mu_r_lb = [node.node1_mu_r_lb, node.node2_mu_r_lb, node.nodeN_mu_r_lb]
		mu_i_lb = [node.node1_mu_i_lb, node.node2_mu_i_lb, node.nodeN_mu_i_lb]  
		mu_r_ub = [node.node1_mu_r_ub, node.node2_mu_r_ub, node.nodeN_mu_r_ub]
		mu_i_ub = [node.node1_mu_i_ub, node.node2_mu_i_ub, node.nodeN_mu_i_ub]
		
		if_r_plus =  [node.node1_if_r_plus,  node.node2_if_r_plus,  node.nodeN_if_r_plus]        
		if_r_minus = [node.node1_if_r_minus, node.node2_if_r_minus, node.nodeN_if_r_minus]
		if_i_plus =  [node.node1_if_i_plus,  node.node2_if_i_plus,  node.nodeN_if_i_plus]   
		if_i_minus = [node.node1_if_i_minus, node.node2_if_i_minus, node.nodeN_if_i_minus]  
		
		Lr = [node.node1_Lr, node.node2_Lr, node.nodeN_Lr]
		Li = [node.node1_Li, node.node2_Li, node.nodeN_Li]
		
		for i in range(0,3):
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_plus[i], Lr[i], 1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_plus[i], mu_r_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_r_plus[i], 1, J_val, J_row, idx_J)
 
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_minus[i], Lr[i], -1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_r_minus[i], mu_r_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_r_minus[i], 1, J_val, J_row, idx_J)
			
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_plus[i], Li[i], 1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_plus[i], mu_i_ub[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_i_plus[i], 1, J_val, J_row, idx_J)

			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_minus[i], Li[i], -1, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y = stampY(if_i_minus[i], mu_i_lb[i], 1, Y_val, Y_row, Y_col, idx_Y)
			J_val, J_row, idx_J = stampJ(if_i_minus[i], 1, J_val, J_row, idx_J) 
	
	elif self.obj == 2:
		# L2-norm
		
		# Relationship between infeasibility currents and lambda
		# dL/dif = if - lambda = 0
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_if_r, node.node1_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_if_i, node.node1_if_i, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_if_r, node.node2_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_if_i, node.node2_if_i, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_r, node.nodeN_if_r, 1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_i, node.nodeN_if_i, 1, Y_val, Y_row, Y_col, idx_Y) 
		
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_if_r, node.node1_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node1_if_i, node.node1_Li, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_if_r, node.node2_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.node2_if_i, node.node2_Li, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_r, node.nodeN_Lr, -1, Y_val, Y_row, Y_col, idx_Y)
		Y_val, Y_row, Y_col, idx_Y = stampY(node.nodeN_if_i, node.nodeN_Li, -1, Y_val, Y_row, Y_col, idx_Y)  
	
	return Y_val, Y_row, Y_col, idx_Y,  J_val, J_row, idx_Y, idx_J

cpdef stamp_linear_infeasibility(self, node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	# this function stamps all of the linear elements; it is called from stamp_linear_test
	if node.isTriplex:
		if self.obj == 2:
			Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L2_tplx(node, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J = self.stamp_infeas_dual_relationship_tplx(node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)            
		
		if self.obj == 1:
			Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L1_tplx(node, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J = self.stamp_infeas_dual_relationship_tplx(node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)                
  
	else:
		if self.obj == 2:
			Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L2(node, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J = self.stamp_infeas_dual_relationship(node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
		
		if self.obj == 1:
			Y_val, Y_row, Y_col, idx_Y = self.stamp_infeas_current_L1(node, Y_val, Y_row, Y_col, idx_Y)
			Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_Y, idx_J = self.stamp_infeas_dual_relationship(node, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)                
		
	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

cpdef stamp_nonlinear_infeasibility(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J):
	# this function stamps all of the nonlinear elements (from complementary slackness) and is called from stamp_nonlinear_test
	if self.obj == 1:
		if node.isTriplex:
			Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_ineq_duals_tplx(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
		else:
			Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_ineq_duals(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
		
	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

	
	
	
