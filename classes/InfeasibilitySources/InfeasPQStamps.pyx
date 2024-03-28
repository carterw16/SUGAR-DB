#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 05:59:55 2023

@author: emfoster

  Implements stamps for the PQ source infeasibility class.

  Author(s): Elizabeth Foster
  Created Date: 10-31-2023
  Updated Date: 10-31-2023
  Email: emfoster@andrew.cmu.edu 
  Status: Development

"""
#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3
import ipdb

cpdef stampY(self, i, j, val, Y_val, Y_row, Y_col, idx):
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return idx


cpdef stampJ(self, i, val, J_val, J_row, idx):
	J_val[idx] = val
	J_row[idx] = i
	idx += 1

	return idx 

cpdef stamp_nonlinear_infeasibility_PQ(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J, cs_tol):
	# Stamping the addition of the infeasibility source in the network equations in Vr, Vi nodes
	Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_equality_constraints(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
	# Stamping equations from (d-Lagrange)/(d-primal variable)
	Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_stationarity_constraints(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
	
	if self.obj == 'L1':
		Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_complementary_slackness(V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J, cs_tol)
	
	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

cpdef stamp_complementary_slackness(self, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J, cs_tol):
	if self.source_type == 'PQ':
		P_plus_nodes = self.P_plus_nodes
		P_minus_nodes = self.P_minus_nodes
		
		dual_ineq_r_plus_nodes = self.dual_ineq_r_plus_nodes
		dual_ineq_r_minus_nodes = self.dual_ineq_r_minus_nodes		
			
	#Q_plus_nodes = self.Q_plus_nodes
	#Q_minus_nodes = self.Q_minus_nodes

	dual_ineq_i_plus_nodes = self.dual_ineq_i_plus_nodes
	dual_ineq_i_minus_nodes = self.dual_ineq_i_minus_nodes

	num_of_phases = len(P_plus_nodes)
		
	for index in range(num_of_phases):
		mu_q_ub = V[dual_ineq_i_plus_nodes[index]]
		mu_q_lb = V[dual_ineq_i_minus_nodes[index]]
		
		Q_plus = 0
		Q_minus = 0
		
		P_plus = V[P_plus_nodes[index]] if self.source_type == 'PQ' else 0
		P_minus = V[P_minus_nodes[index]] if self.source_type == 'PQ' else 0

		if self.source_type == 'PQ':
			mu_p_ub = V[dual_ineq_r_plus_nodes[index]]
			mu_p_lb = V[dual_ineq_r_minus_nodes[index]]
			
			# -P_plus*mu_ub + cs_eps = 0
			hist_mu_p_ub = -P_plus*mu_p_ub + cs_tol
			idx_Y = self.stampY(dual_ineq_r_plus_nodes[index], dual_ineq_r_plus_nodes[index], -P_plus, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(dual_ineq_r_plus_nodes[index], P_plus_nodes[index], -mu_p_ub, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			idx_J = self.stampJ(dual_ineq_r_plus_nodes[index], -hist_mu_p_ub - 2*P_plus*mu_p_ub, Jnlin_val, Jnlin_row, idx_J)
				
			# -P_minus*mu_lb + cs_eps
			hist_mu_p_lb = -P_minus*mu_p_lb + cs_tol
			idx_Y = self.stampY(dual_ineq_r_minus_nodes[index], dual_ineq_r_minus_nodes[index], -P_minus, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(dual_ineq_r_minus_nodes[index], P_minus_nodes[index], -mu_p_lb, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			idx_J = self.stampJ(dual_ineq_r_minus_nodes[index], -hist_mu_p_lb - 2*P_minus*mu_p_lb, Jnlin_val, Jnlin_row, idx_J)

			
			# (debug) upper bounds

			#P_max = 1800000
			P_max = 3000000
			dual_ineq_r_plus_nodes_upper = self.dual_ineq_r_plus_nodes_upper
			dual_ineq_r_minus_nodes_upper = self.dual_ineq_r_minus_nodes_upper	

			#print("stamping lower bounds", dual_ineq_r_plus_nodes, dual_ineq_r_minus_nodes)

			#print("stamping upper bounds", dual_ineq_r_plus_nodes_upper, dual_ineq_r_minus_nodes_upper)

			mu_p_ub = V[dual_ineq_r_plus_nodes_upper[index]]
			mu_p_lb = V[dual_ineq_r_minus_nodes_upper[index]]
			
			#print("mu [upper] Pch {} mu [upper] Pd {}".format(mu_p_ub, mu_p_lb))
			#print("mu p ub {} mu p lb {} P_plus {} P minus {}".format(mu_p_ub, mu_p_lb, P_plus, P_minus ))

			# (P_plus-P_max)*mu_ub + cs_eps = 0

			hist_mu_p_ub = (P_plus-P_max)*mu_p_ub + cs_tol
			idx_Y = self.stampY(dual_ineq_r_plus_nodes_upper[index], dual_ineq_r_plus_nodes_upper[index], P_plus-P_max, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(dual_ineq_r_plus_nodes_upper[index], P_plus_nodes[index], mu_p_ub, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			idx_J = self.stampJ(dual_ineq_r_plus_nodes_upper[index], -(hist_mu_p_ub - (P_plus - P_max)*mu_p_ub - mu_p_ub*P_plus), Jnlin_val, Jnlin_row, idx_J)
				
			# (P_minus-P_max)*mu_lb + cs_eps
			hist_mu_p_lb = (P_minus-P_max)*mu_p_lb + cs_tol
			idx_Y = self.stampY(dual_ineq_r_minus_nodes_upper[index], dual_ineq_r_minus_nodes_upper[index], P_minus-P_max, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(dual_ineq_r_minus_nodes_upper[index], P_minus_nodes[index], mu_p_lb, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			idx_J = self.stampJ(dual_ineq_r_minus_nodes_upper[index], -(hist_mu_p_lb - (P_minus - P_max)*mu_p_lb - mu_p_lb*P_minus), Jnlin_val, Jnlin_row, idx_J)


			
	return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J

cpdef stamp_equality_constraints(self, node, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J):	
	if self.isTriplex:
		node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
		node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
	else:
		node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
		node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
	
	if self.obj == 'L2':
		if self.source_type == 'PQ':
			P_nodes = self.P_nodes
		Q_nodes = self.Q_nodes

		num_of_phases = len(Q_nodes)
		
		for index in range(num_of_phases):
			Vr = V[node_Vr[index]]
			Vi = V[node_Vi[index]]
			P = V[P_nodes[index]] if self.source_type == 'PQ' else 0
			Q = V[Q_nodes[index]]
			partials = self.calculate_partial_derivatives(Vr, Vi, P, Q, calc_hessian = False)
			
			# Equality Constraint - Real
			# Real Network Equation: ... + (P*Vr + Q*Vi)/(Vr**2 + Vi**2) = 0
			# Taylor Expansion: ... + 
			idx_Y = self.stampY(node_Vr[index], node_Vr[index], partials['dIr_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vr[index], node_Vi[index], partials['dIr_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vr[index], Q_nodes[index], partials['dIr_dQ'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Vr[index], P_nodes[index], partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Vr[index], (partials['dIr_dVr']*Vr + partials['dIr_dVi']*Vi + partials['dIr_dP']*P + partials['dIr_dQ']*Q - partials['Ir']), Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Vr[index], (partials['dIr_dVr']*Vr + partials['dIr_dVi']*Vi + 
									partials['dIr_dQ']*Q - partials['Ir']), Jnlin_val, Jnlin_row, idx_J)

			# Equality Constraint - Imag
			# Imag Network Equation: ... + (P*Vi - Q*Vr)/(Vr**2 + Vi**2) = 0
			idx_Y = self.stampY(node_Vi[index], node_Vr[index], partials['dIi_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vi[index], node_Vi[index], partials['dIi_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vi[index], Q_nodes[index], partials['dIi_dQ'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Vi[index], P_nodes[index], partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Vi[index], (partials['dIi_dVr']*Vr + partials['dIi_dVi']*Vi + 
								partials['dIi_dP']*P +  partials['dIi_dQ']*Q - partials['Ii']), Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Vi[index], (partials['dIi_dVr']*Vr + partials['dIi_dVi']*Vi + 
									partials['dIi_dQ']*Q - partials['Ii']), Jnlin_val, Jnlin_row, idx_J)
	
	if self.obj == 'L1':
		if self.source_type == 'PQ':
			P_plus_nodes = self.P_plus_nodes
			P_minus_nodes = self.P_minus_nodes

		num_of_phases = len(P_plus_nodes)
		
		for index in range(num_of_phases):
			Vr = V[node_Vr[index]]
			Vi = V[node_Vi[index]]
			P_plus = V[P_plus_nodes[index]] if self.source_type == 'PQ' else 0
			P_minus = V[P_minus_nodes[index]] if self.source_type == 'PQ' else 0

			Q_plus = 0 # placeholders for debugging
			Q_minus = 0 # placeholders for debugging

			partials = self.calculate_partial_derivatives(Vr, Vi, (P_plus - P_minus), (Q_plus - Q_minus), calc_hessian = False)
			
			# Equality Constraint - Real
			# Real Network Equation: ... + (P*Vr + Q*Vi)/(Vr**2 + Vi**2) = 0
			# Taylor Expansion: ... + 
			idx_Y = self.stampY(node_Vr[index], node_Vr[index], partials['dIr_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vr[index], node_Vi[index], partials['dIr_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Vr[index], P_plus_nodes[index], partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(node_Vr[index], P_minus_nodes[index], -1*partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Vr[index], (partials['dIr_dVr']*Vr + partials['dIr_dVi']*Vi +
													partials['dIr_dP']*(P_plus - P_minus) + partials['dIr_dQ']*(Q_plus - Q_minus) - partials['Ir']), 
													Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Vr[index], (partials['dIr_dVr']*Vr + partials['dIr_dVi']*Vi + 
									partials['dIr_dQ']*(Q_plus - Q_minus) - partials['Ir']), Jnlin_val, Jnlin_row, idx_J)

			# Equality Constraint - Imag
			# Imag Network Equation: ... + (P*Vi - Q*Vr)/(Vr**2 + Vi**2) = 0
			idx_Y = self.stampY(node_Vi[index], node_Vr[index], partials['dIi_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stampY(node_Vi[index], node_Vi[index], partials['dIi_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Vi[index], P_plus_nodes[index], partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(node_Vi[index], P_minus_nodes[index], -1*partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Vi[index], (partials['dIi_dVr']*Vr + partials['dIi_dVi']*Vi + 
								partials['dIi_dP']*(P_plus - P_minus) +  partials['dIi_dQ']*(Q_plus - Q_minus) - partials['Ii']), Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Vi[index], (partials['dIi_dVr']*Vr + partials['dIi_dVi']*Vi + 
									partials['dIi_dQ']*(Q_plus - Q_minus) - partials['Ii']), Jnlin_val, Jnlin_row, idx_J)	
	
	return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J

cpdef stamp_stationarity_constraints(self, node, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J):
	node_Lr = node.dual_eq_var_r_nodes
	node_Li = node.dual_eq_var_i_nodes
	
	if self.isTriplex:
		node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
		node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
	else:
		node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
		node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]
	
	if self.obj == 'L2':
		if self.source_type == 'PQ':
			P_nodes = self.P_nodes
		Q_nodes = self.Q_nodes

		num_of_phases = len(Q_nodes)
	
		for index in range(num_of_phases):
			Vr = V[node_Vr[index]]
			Vi = V[node_Vi[index]]
			P = V[P_nodes[index]] if self.source_type == 'PQ' else 0
			Q = V[Q_nodes[index]]
			Lr = V[node_Lr[index]]
			Li = V[node_Li[index]]

			partials = self.calculate_partial_derivatives(Vr, Vi, P, Q, calc_hessian = True)
			###################### STATIONARITY CONSTRAINTS ############################################
			# Equations for primal variable Vr
			# dL/dVr [stamped to Lr nodes]
			# dL/dVr = ... + Lr*dIpqr/dVr + Li*dIpqi/dVr
			idx_Y= self.stampY(node_Lr[index], node_Vr[index], partials['d2Ir_dVr2']*Lr + 
							partials['d2Ii_dVr2']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Vi[index], partials['d2Ir_dVrdVi']*Lr + 
							partials['d2Ii_dVrdVi']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], Q_nodes[index], partials['d2Ir_dVrdQ']*Lr + 
							partials['d2Ii_dVrdQ']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Lr[index], partials['dIr_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Li[index], partials['dIi_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Lr[index], P_nodes[index],  partials['d2Ir_dVrdP']*Lr + 
								partials['d2Ii_dVrdP']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

				idx_J= self.stampJ(node_Lr[index], (partials['d2Ir_dVr2']*Lr + partials['d2Ii_dVr2']*Li)*Vr +
									(partials['d2Ir_dVrdVi']*Lr + partials['d2Ii_dVrdVi']*Li)*Vi + 
									(partials['d2Ir_dVrdP']*Lr + partials['d2Ii_dVrdP']*Li)*P + 
									(partials['d2Ir_dVrdQ']*Lr + partials['d2Ii_dVrdQ']*Li)*Q, Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Lr[index], (partials['d2Ir_dVr2']*Lr + partials['d2Ii_dVr2']*Li)*Vr +
									(partials['d2Ir_dVrdVi']*Lr + partials['d2Ii_dVrdVi']*Li)*Vi + 
									(partials['d2Ir_dVrdQ']*Lr + partials['d2Ii_dVrdQ']*Li)*Q, Jnlin_val, Jnlin_row, idx_J)

			# Equations for primal variable Vi
			# dL/dVi [stamped to Li nodes]
			# dL/dVi = ... + Lr*dIpqr/dVi + Li*dIpqi/dVi
			idx_Y= self.stampY(node_Li[index], node_Vr[index], (partials['d2Ir_dVidVr']*Lr + 
							partials['d2Ii_dVidVr']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Vi[index], (partials['d2Ir_dVi2']*Lr + 
							partials['d2Ii_dVi2']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], Q_nodes[index], (partials['d2Ir_dVidQ']*Lr 
							+ partials['d2Ii_dVidQ']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Lr[index], partials['dIr_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Li[index], partials['dIi_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Li[index], P_nodes[index], (partials['d2Ir_dVidP']*Lr 
								+ partials['d2Ii_dVidP']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Li[index],(partials['d2Ir_dVidVr']*Lr + partials['d2Ii_dVidVr']*Li)*Vr
								+ (partials['d2Ir_dVi2']*Lr + partials['d2Ii_dVi2']*Li)*Vi +
								(partials['d2Ir_dVidP']*Lr + partials['d2Ii_dVidP']*Li)*P + 
								(partials['d2Ir_dVidQ']*Lr + partials['d2Ii_dVidQ']*Li)*Q, Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Li[index],(partials['d2Ir_dVidVr']*Lr + partials['d2Ii_dVidVr']*Li)*Vr 
								+ (partials['d2Ir_dVi2']*Lr + partials['d2Ii_dVi2']*Li)*Vi +
								(partials['d2Ir_dVidQ']*Lr + partials['d2Ii_dVidQ']*Li)*Q, Jnlin_val, Jnlin_row, idx_J)

			# Equations for primal variable P
			# dL/dP
			# dL/dP = ... + Lr*dIpqr/dP + Li*dIpqi/dP + df_obj/dP
			# f_obj = cost_val*(P^2 + Q^2)
			cost_val = self.obj_scaling#1/(self.Vnom_ln**2)#1#min(.001, 1-h_factor)
			if self.source_type == 'PQ':
				histP = Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] + 2*cost_val*P
				idx_Y= self.stampY(P_nodes[index], node_Vr[index], (partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_nodes[index], node_Vi[index], (partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_nodes[index], P_nodes[index], cost_val*2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_nodes[index], node_Lr[index], partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_nodes[index], node_Li[index], partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(P_nodes[index], -histP + (partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
								(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi + 
								Lr*partials['dIr_dP'] + Li*partials['dIi_dP']+ P*cost_val*2, Jnlin_val, Jnlin_row, idx_J)
			
			# Equations for primal variable Q
			# dL/dQ
			# dL/dQ = ... + Lr*dIpqr/dQ + Li*dIpqi/dQ + df_obj/dQ
			# f_obj = cost_val*(P^2 + Q^2)
			histQ = partials['dIr_dQ']*Lr + partials['dIi_dQ']*Li + 2*Q*cost_val
			idx_Y= self.stampY(Q_nodes[index], node_Vr[index], (partials['d2Ir_dQdVr']*Lr + partials['d2Ii_dQdVr']*Li), 
							Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(Q_nodes[index], node_Vi[index], (partials['d2Ir_dQdVi']*Lr + partials['d2Ii_dQdVi']*Li), 
							Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(Q_nodes[index], Q_nodes[index], cost_val*2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(Q_nodes[index], node_Lr[index], 
							partials['dIr_dQ'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(Q_nodes[index], node_Li[index], 
							partials['dIi_dQ'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

			idx_J= self.stampJ(Q_nodes[index], -histQ - (partials['d2Ir_dQdVr']*Lr + partials['d2Ii_dQdVr']*Li)*Vr + 
							(partials['d2Ir_dQdVi']*Lr + partials['d2Ii_dQdVi']*Li)*Vi + 
							partials['dIr_dQ']*Lr + partials['dIi_dQ']*Li + 2*Q*cost_val, Jnlin_val, Jnlin_row, idx_J)

	if self.obj == 'L1':
		if self.source_type == 'PQ':
			P_plus_nodes = self.P_plus_nodes
			P_minus_nodes = self.P_minus_nodes
			dual_ineq_r_plus_nodes = self.dual_ineq_r_plus_nodes
			dual_ineq_r_minus_nodes = self.dual_ineq_r_minus_nodes

		dual_ineq_i_plus_nodes = self.dual_ineq_i_plus_nodes
		dual_ineq_i_minus_nodes = self.dual_ineq_i_minus_nodes

		num_of_phases = len(P_plus_nodes)
	
		for index in range(num_of_phases):
			Vr = V[node_Vr[index]]
			Vi = V[node_Vi[index]]
			P_plus = V[P_plus_nodes[index]] if self.source_type == 'PQ' else 0
			P_minus = V[P_minus_nodes[index]] if self.source_type == 'PQ' else 0
			Q_plus = 0 #
			Q_minus = 0 #
			Lr = V[node_Lr[index]]
			Li = V[node_Li[index]]

			partials = self.calculate_partial_derivatives(Vr, Vi, (P_plus - P_minus), (Q_plus - Q_minus), calc_hessian = True)
			###################### STATIONARITY CONSTRAINTS ############################################
			# Equations for primal variable Vr
			# dL/dVr [stamped to Lr nodes]
			# dL/dVr = ... + Lr*dIpqr/dVr + Li*dIpqi/dVr
			hist_dLdVr = Lr*partials['dIr_dVr'] + Li*partials['dIi_dVr']
			idx_Y= self.stampY(node_Lr[index], node_Vr[index], partials['d2Ir_dVr2']*Lr + 
							partials['d2Ii_dVr2']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Vi[index], partials['d2Ir_dVrdVi']*Lr + 
							partials['d2Ii_dVrdVi']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Lr[index], partials['dIr_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], node_Li[index], partials['dIi_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Lr[index], P_plus_nodes[index],  partials['d2Ir_dVrdP']*Lr + 
								partials['d2Ii_dVrdP']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(node_Lr[index], P_minus_nodes[index],  -1*partials['d2Ir_dVrdP']*Lr + 
								-1*partials['d2Ii_dVrdP']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Lr[index], (-hist_dLdVr +  partials['dIr_dVr']*Lr + partials['dIi_dVr']*Li +
									partials['d2Ir_dVr2']*Lr + partials['d2Ii_dVr2']*Li)*Vr +
									(partials['d2Ir_dVrdVi']*Lr + partials['d2Ii_dVrdVi']*Li)*Vi + 
									(partials['d2Ir_dVrdP']*Lr + partials['d2Ii_dVrdP']*Li)*(P_plus - P_minus), Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Lr[index], (-hist_dLdVr +  partials['dIr_dVr']*Lr + partials['dIi_dVr']*Li +
									partials['d2Ir_dVr2']*Lr + partials['d2Ii_dVr2']*Li)*Vr +
									(partials['d2Ir_dVrdVi']*Lr + partials['d2Ii_dVrdVi']*Li)*Vi +  
									(partials['d2Ir_dVrdQ']*Lr + partials['d2Ii_dVrdQ']*Li)*(Q_plus - Q_minus), Jnlin_val, Jnlin_row, idx_J)
			# Equations for primal variable Vi
			# dL/dVi [stamped to Li nodes]
			# dL/dVi = ... + Lr*dIpqr/dVi + Li*dIpqi/dVi
			hist_dLdVi = partials['dIr_dVi']*Lr + partials['dIi_dVi']*Li
			idx_Y= self.stampY(node_Li[index], node_Vr[index], (partials['d2Ir_dVidVr']*Lr + 
							partials['d2Ii_dVidVr']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Vi[index], (partials['d2Ir_dVi2']*Lr + 
							partials['d2Ii_dVi2']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Lr[index], partials['dIr_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], node_Li[index], partials['dIi_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			if self.source_type == 'PQ':
				idx_Y= self.stampY(node_Li[index], P_plus_nodes[index], (partials['d2Ir_dVidP']*Lr 
								+ partials['d2Ii_dVidP']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(node_Li[index], P_minus_nodes[index], -1*(partials['d2Ir_dVidP']*Lr 
								+ partials['d2Ii_dVidP']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_J= self.stampJ(node_Li[index], -hist_dLdVi + partials['dIr_dVi']*Lr + partials['dIi_dVi']*Li +
								(partials['d2Ir_dVidVr']*Lr + partials['d2Ii_dVidVr']*Li)*Vr
								+ (partials['d2Ir_dVi2']*Lr + partials['d2Ii_dVi2']*Li)*Vi +
								(partials['d2Ir_dVidP']*Lr + partials['d2Ii_dVidP']*Li)*(P_plus - P_minus), Jnlin_val, Jnlin_row, idx_J)
			else:
				idx_J= self.stampJ(node_Li[index],-hist_dLdVi + partials['dIr_dVi']*Lr + partials['dIi_dVi']*Li +
								(partials['d2Ir_dVidVr']*Lr + partials['d2Ii_dVidVr']*Li)*Vr 
								+ (partials['d2Ir_dVi2']*Lr + partials['d2Ii_dVi2']*Li)*Vi +
								(partials['d2Ir_dVidQ']*Lr + partials['d2Ii_dVidQ']*Li)*(Q_plus - Q_minus), Jnlin_val, Jnlin_row, idx_J)

			if self.source_type == "PQ":

				pch_l2 = False
				pd_l2 = False
				c1 = 500
				c2 = -0.5
				# --- Pch
				# Equations for primal variable P+
				# dL/dP = Lr*dIpqr/dP + Li*dIpqi/dP + df_obj/dP + mu
				# new objective: f_obj = 0.5*C1*P_plus^2 - 0.5*C2*P_minus^2

				hist_dLdPplus = Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] - V[dual_ineq_r_plus_nodes[index]] + c1
				idx_Y= self.stampY(P_plus_nodes[index], node_Vr[index], (partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_plus_nodes[index], node_Vi[index], (partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_plus_nodes[index], node_Lr[index], 
								partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_plus_nodes[index], node_Li[index], 
								partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_plus_nodes[index], dual_ineq_r_plus_nodes[index], 
								-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

				if (pch_l2 == True):
					idx_Y = self.stampY(P_plus_nodes[index], P_plus_nodes[index], c1,
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
					hist_dLdPplus = Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] - V[dual_ineq_r_plus_nodes[index]] + c1 * P_plus
					idx_J= self.stampJ(P_plus_nodes[index], -hist_dLdPplus + Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] +
						(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
						(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi + c1 * P_plus
						- V[dual_ineq_r_plus_nodes[index]], Jnlin_val, Jnlin_row, idx_J)
				else:
					idx_J= self.stampJ(P_plus_nodes[index], -hist_dLdPplus + Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] +
									(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
									(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi
									- V[dual_ineq_r_plus_nodes[index]], Jnlin_val, Jnlin_row, idx_J)

				# --- Pd
				# Equations for primal variable P-
				# dL/dP = -Lr*dIpqr/dP - Li*dIpqi/dP - df_obj/dP + mu
				# new objective: f_obj = 0.5*C1*P_plus^2 - 0.5*C2*P_minus^2

				hist_dLdPminus = -Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] - V[dual_ineq_r_minus_nodes[index]] + c2
				idx_Y= self.stampY(P_minus_nodes[index], node_Vr[index], -1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_minus_nodes[index], node_Vi[index], -1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li), 
								Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_minus_nodes[index], node_Lr[index], 
								-1*partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_minus_nodes[index], node_Li[index], 
								-1*partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y= self.stampY(P_minus_nodes[index], dual_ineq_r_minus_nodes[index], 
								-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				
				if (pd_l2 == True):
					idx_Y = self.stampY(P_minus_nodes[index], P_minus_nodes[index], c2,
									Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
					hist_dLdPminus = -Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] - V[dual_ineq_r_minus_nodes[index]] + c2 * P_minus
					idx_J= self.stampJ(P_minus_nodes[index], -hist_dLdPminus - Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] +
						-1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
						-1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi + c2 * P_minus
						- V[dual_ineq_r_minus_nodes[index]], Jnlin_val, Jnlin_row, idx_J)
				else:
					idx_J= self.stampJ(P_minus_nodes[index], -hist_dLdPminus - Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] +
									-1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
									-1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi
									- V[dual_ineq_r_minus_nodes[index]], Jnlin_val, Jnlin_row, idx_J)
				
				# (debug) upper bounds to bound Pch and Pd values < some constant
				#print("dual ineq r plus nodes", dual_ineq_r_plus_nodes)
				#print("dual ineq r minus nodes", dual_ineq_r_minus_nodes)


				dual_ineq_r_plus_nodes_upper = self.dual_ineq_r_plus_nodes_upper
				dual_ineq_r_minus_nodes_upper = self.dual_ineq_r_minus_nodes_upper		
				#print("dual ineq r plus nodes upper", dual_ineq_r_plus_nodes_upper)
				#print("dual ineq r minus nodes upper", dual_ineq_r_minus_nodes_upper)		
				
				# upper bound matrix stamps
				# V[dual] and J stamp addition not needed as mu is a linear term

				# Pch
				idx_Y= self.stampY(P_plus_nodes[index], dual_ineq_r_plus_nodes_upper[index], 
				1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

				# Pd
				idx_Y= self.stampY(P_minus_nodes[index], dual_ineq_r_minus_nodes_upper[index], 
				1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)


				#print("phase {} Pch {} Pd {} dLdPch {} dLdPd {}".format(index, P_plus, P_minus, hist_dLdPplus, hist_dLdPminus))
			
	return(Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J)


cpdef calculate_partial_derivatives(self, Vr, Vi, P, Q, calc_hessian):
	# Inputs:
	# calc_hessian is a T/F flag that is false for power flow stamps or equality constraint stamps
	# and true for dL/dprimal stationarity constraints
	partials = {}
	v_mag = Vr**2 + Vi**2
	v_mag2 = v_mag**2
	v_mag3 = v_mag**3

	partials['Ir'] = (P*Vr + Q*Vi)/v_mag if v_mag != 0 else 0
	partials['Ii'] = (P*Vi - Q*Vr)/v_mag if v_mag != 0 else 0

	Ir_diff_Ii_same = (2*Q*(Vr**3) - 6*P*(Vr**2)*Vi
						- 6*Q*Vr*(Vi**2) + 2*P*(Vi**3))/v_mag3 if v_mag !=0 else 0
	Ir_same_Ii_diff = (2*P*(Vr**3) + 6*Q*(Vr**2)*Vi 
						- 6*P*Vr*(Vi**2) - 2*Q*(Vi**3))/v_mag3 if v_mag !=0 else 0

	numerator_2VrVi = 2*Vr*Vi/v_mag2 if v_mag !=0 else 0
	numerator_Vr2_Vi2 = (Vr**2 - Vi**2)/v_mag2 if v_mag !=0 else 0

	partials['dIr_dP'] = Vr/v_mag if v_mag != 0 else 0
	if calc_hessian:
		partials['d2Ir_dPdVr'] = -1*numerator_Vr2_Vi2
		partials['d2Ir_dPdVi'] = -1*numerator_2VrVi

	partials['dIi_dP'] = Vi/v_mag if v_mag != 0 else 0
	if calc_hessian:
		partials['d2Ii_dPdVr'] = -1*numerator_2VrVi
		partials['d2Ii_dPdVi'] = numerator_Vr2_Vi2

	partials['dIr_dQ'] = partials['dIi_dP']
	if calc_hessian:
		partials['d2Ir_dQdVr'] = -1*numerator_2VrVi
		partials['d2Ir_dQdVi'] = numerator_Vr2_Vi2

	partials['dIi_dQ'] = -partials['dIr_dP']
	if calc_hessian:
		partials['d2Ii_dQdVr'] = numerator_Vr2_Vi2
		partials['d2Ii_dQdVi'] = numerator_2VrVi

	partials['dIr_dVr'] = (-P*(Vr**2 - Vi**2) - 2*Q*Vr*Vi)/v_mag2 if v_mag !=0 else 0
	if calc_hessian:
		partials['d2Ir_dVr2'] = Ir_same_Ii_diff
		partials['d2Ir_dVrdP'] = -1*numerator_Vr2_Vi2
		partials['d2Ir_dVrdQ'] = -1*numerator_2VrVi
		partials['d2Ir_dVrdVi'] = -1*Ir_diff_Ii_same

	partials['dIi_dVr'] = (Q*(Vr**2 - Vi**2) - 2*P*Vr*Vi)/v_mag2 if v_mag !=0 else 0
	if calc_hessian:
		partials['d2Ii_dVr2'] = -1*Ir_diff_Ii_same
		partials['d2Ii_dVrdP'] = -numerator_2VrVi 
		partials['d2Ii_dVrdQ'] = numerator_Vr2_Vi2
		partials['d2Ii_dVrdVi'] = -1*Ir_same_Ii_diff

	partials['dIr_dVi'] = (Q*(Vr**2 - Vi**2) - 2*P*Vr*Vi)/v_mag2 if v_mag !=0 else 0
	if calc_hessian:
		partials['d2Ir_dVidVr'] = -1*Ir_diff_Ii_same
		partials['d2Ir_dVidP'] = -1*numerator_2VrVi
		partials['d2Ir_dVidQ'] = numerator_Vr2_Vi2
		partials['d2Ir_dVi2'] = -1*Ir_same_Ii_diff

	partials['dIi_dVi'] = (P*(Vr**2 - Vi**2) + 2*Q*Vr*Vi)/v_mag2 if v_mag !=0 else 0
	if calc_hessian:
		partials['d2Ii_dVidVr'] = -1*Ir_same_Ii_diff
		partials['d2Ii_dVidP'] = numerator_Vr2_Vi2
		partials['d2Ii_dVidQ'] = numerator_2VrVi
		partials['d2Ii_dVi2'] = Ir_diff_Ii_same
	
	return partials