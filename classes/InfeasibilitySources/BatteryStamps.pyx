#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 05:59:55 2023

@author: mhopkins

  Stamps for Battery object

  Author(s): Meshach Hopkins
  Created Date: 2-15-2024
  Updated Date: 2-15-2024
  Email: mhopkins@andrew.cmu.edu
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

cpdef stamp_nonlinear_infeasibility(self, node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J, cs_tol):
	# Stamping the addition of the infeasibility source in the network equations in Vr, Vi nodes
	Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_equality_constraints(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
	# Stamping equations from (d-Lagrange)/(d-primal variable)
	Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_stationarity_constraints(node, V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J)
	# complementary slackness
	Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_complementary_slackness(V, Y_val, Y_row, Y_col, idx_Y, J_val, J_row, idx_J, cs_tol)

	return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

cpdef stamp_complementary_slackness(self, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J, cs_tol):
	P_plus_nodes = self.P_plus_nodes
	P_minus_nodes = self.P_minus_nodes

	dual_ineq_r_plus_nodes = self.dual_ineq_r_plus_nodes
	dual_ineq_r_minus_nodes = self.dual_ineq_r_minus_nodes

	num_of_phases = len(P_plus_nodes)
	
	phases = ["A","B","C"]
	iter_phases = range(num_of_phases)
	if (self.single_phase):
		iter_phases = [phases.index(self.single_phase)]


	for index,v_index in enumerate(iter_phases):
		
		P_plus = V[P_plus_nodes[index]] 
		P_minus = V[P_minus_nodes[index]]


		mu_pch = V[dual_ineq_r_plus_nodes[index]]
		mu_pd = V[dual_ineq_r_minus_nodes[index]]

		# Stamping Power Distribution Lower Bound Stamps

		# -P_plus*mu_ub + cs_eps = 0
		hist_mu_pch = -P_plus*mu_pch + cs_tol
		idx_Y = self.stampY(dual_ineq_r_plus_nodes[index], dual_ineq_r_plus_nodes[index], -P_plus, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(dual_ineq_r_plus_nodes[index], P_plus_nodes[index], -mu_pch, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stampJ(dual_ineq_r_plus_nodes[index], -hist_mu_pch - 2*P_plus*mu_pch, Jnlin_val, Jnlin_row, idx_J)

		# -P_minus*mu_lb + cs_eps
		hist_mu_pd = -P_minus*mu_pd + cs_tol
		idx_Y = self.stampY(dual_ineq_r_minus_nodes[index], dual_ineq_r_minus_nodes[index], -P_minus, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(dual_ineq_r_minus_nodes[index], P_minus_nodes[index], -mu_pd, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stampJ(dual_ineq_r_minus_nodes[index], -hist_mu_pd - 2*P_minus*mu_pd, Jnlin_val, Jnlin_row, idx_J)


		# Stamping Power Distribution Upper Bound Stamps
		P_max = self.P_max
		dual_ineq_r_plus_nodes_upper = self.dual_ineq_r_plus_nodes_upper
		dual_ineq_r_minus_nodes_upper = self.dual_ineq_r_minus_nodes_upper

		mu_pch_upper = V[dual_ineq_r_plus_nodes_upper[index]]
		mu_pd_upper = V[dual_ineq_r_minus_nodes_upper[index]]

		# (P_plus-P_max)*mu_ub + cs_eps = 0
		hist_mu_p_ub = (P_plus-P_max)*mu_pch_upper + cs_tol
		idx_Y = self.stampY(dual_ineq_r_plus_nodes_upper[index], dual_ineq_r_plus_nodes_upper[index], P_plus-P_max, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(dual_ineq_r_plus_nodes_upper[index], P_plus_nodes[index], mu_pch_upper, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stampJ(dual_ineq_r_plus_nodes_upper[index], -(hist_mu_p_ub - (P_plus - P_max)*mu_pch_upper - mu_pch_upper*P_plus), Jnlin_val, Jnlin_row, idx_J)

		# (P_minus-P_max)*mu_lb + cs_eps
		hist_mu_p_lb = (P_minus-P_max)*mu_pd_upper + cs_tol
		idx_Y = self.stampY(dual_ineq_r_minus_nodes_upper[index], dual_ineq_r_minus_nodes_upper[index], P_minus-P_max, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(dual_ineq_r_minus_nodes_upper[index], P_minus_nodes[index], mu_pd_upper, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stampJ(dual_ineq_r_minus_nodes_upper[index], -(hist_mu_p_lb - (P_minus - P_max)*mu_pd_upper - mu_pd_upper*P_minus), Jnlin_val, Jnlin_row, idx_J)

		if (self.verbose):
			print("phase {}".format(v_index))
			print("[upper bounds] mu Pch {} mu Pd {}".format(mu_pch_upper, mu_pd_upper))
			print("[lower bounds] mu Pch {} mu Pd {}".format(mu_pch, mu_pd))


		# TODO: Add lower and upper bounds for reactive power output

		if self.source_type == 'PQ':
			Q_plus_nodes = self.Q_plus_nodes
			Q_minus_nodes = self.Q_minus_nodes

			dual_ineq_i_nodes = self.dual_ineq_i_nodes


	# Stamping SOC Lower and Upper Bounds
	Bt_nodes = self.Bt_nodes[0]
	dual_Bt_nodes = self.mu_index_Bt[0]
	dual_Bt_nodes_upper = self.mu_index_Bt_upper[0]
	Bt = V[Bt_nodes]
	mu_Bt = V[dual_Bt_nodes]
	mu_Bt_upper = V[dual_Bt_nodes_upper]
	Bt_max = 1
	Bt_min = 0

	# -Bt * mu_Bt + cs_eps = 0 (lower)
	hist_mu_Bt = -Bt*mu_Bt + cs_tol
	idx_Y = self.stampY(dual_Bt_nodes, dual_Bt_nodes, -Bt, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = self.stampY(dual_Bt_nodes, Bt_nodes, -mu_Bt, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
	idx_J = self.stampJ(dual_Bt_nodes, -hist_mu_Bt- 2*Bt*mu_Bt, Jnlin_val, Jnlin_row, idx_J)

	if (self.verbose):
		print("[upper bound] Bt {}".format(mu_Bt_upper))
		print("[lower bound] Bt {}".format(mu_Bt))
	
	# (Bt-0) * mu_Bt + cs_eps = 0 (upper)
	hist_mu_Bt_upper = (Bt-1)*mu_Bt_upper + cs_tol
	idx_Y = self.stampY(dual_Bt_nodes_upper, dual_Bt_nodes_upper, (Bt-1), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = self.stampY(dual_Bt_nodes_upper, Bt_nodes, mu_Bt_upper, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
	idx_J = self.stampJ(dual_Bt_nodes_upper, -(hist_mu_Bt_upper- (Bt-1)*mu_Bt_upper - mu_Bt_upper*Bt), Jnlin_val, Jnlin_row, idx_J)

			
	return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J

cpdef stamp_equality_constraints(self, node, V, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y, Jnlin_val, Jnlin_row, idx_J):
	if self.isTriplex:
		node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
		node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
	else:
		node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
		node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]


	P_plus_nodes = self.P_plus_nodes
	P_minus_nodes = self.P_minus_nodes
	Bt_nodes = self.Bt_nodes

	num_of_phases = len(P_plus_nodes)

	phases = ["A","B","C"]
	iter_phases = range(num_of_phases)
	if (self.single_phase):
		iter_phases = [phases.index(self.single_phase)]


	for index,v_index in enumerate(iter_phases):
		Vr = V[node_Vr[v_index]]
		Vi = V[node_Vi[v_index]]
		P_plus = V[P_plus_nodes[index]]
		P_minus = V[P_minus_nodes[index]]

		Q_plus = 0 # placeholders for debugging
		Q_minus = 0 # placeholders for debugging

		partials = self.calculate_partial_derivatives(Vr, Vi, (P_plus - P_minus), (Q_plus - Q_minus), calc_hessian = False)

		# Equality Constraint - Real
		# Real Network Equation: ... + (P*Vr + Q*Vi)/(Vr**2 + Vi**2) = 0
		# Taylor Expansion: ... + 
		idx_Y = self.stampY(node_Vr[v_index], node_Vr[v_index], partials['dIr_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Vr[v_index], node_Vi[v_index], partials['dIr_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)	
		idx_Y= self.stampY(node_Vr[v_index], P_plus_nodes[index], partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(node_Vr[v_index], P_minus_nodes[index], -1*partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J= self.stampJ(node_Vr[v_index], (partials['dIr_dVr']*Vr + partials['dIr_dVi']*Vi +
											partials['dIr_dP']*(P_plus - P_minus) + partials['dIr_dQ']*(Q_plus - Q_minus) - partials['Ir']), 
											Jnlin_val, Jnlin_row, idx_J)

		# Equality Constraint - Imag
		# Imag Network Equation: ... + (P*Vi - Q*Vr)/(Vr**2 + Vi**2) = 0
		idx_Y = self.stampY(node_Vi[v_index], node_Vr[v_index], partials['dIi_dVr'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Vi[v_index], node_Vi[v_index], partials['dIi_dVi'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(node_Vi[v_index], P_plus_nodes[index], partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(node_Vi[v_index], P_minus_nodes[index], -1*partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J= self.stampJ(node_Vi[v_index], (partials['dIi_dVr']*Vr + partials['dIi_dVi']*Vi + 
						partials['dIi_dP']*(P_plus - P_minus) +  partials['dIi_dQ']*(Q_plus - Q_minus) - partials['Ii']), Jnlin_val, Jnlin_row, idx_J)


	# SOC Equality Constraint
	node_Bt = Bt_nodes[0]
	Bt = V[node_Bt]

	# SOC constraint
	if (self.single_phase == ""):

		# Equality Stamp for SOC is Based on the Summation of All Powers

		# TODO: this will need to be changed when triplex nodes are added

		P_plusA = V[P_plus_nodes[0]]
		P_minusA = V[P_plus_nodes[0]]
		P_plusB = V[P_plus_nodes[1]]
		P_minusB = V[P_plus_nodes[1]]
		P_plusC = V[P_plus_nodes[2]]
		P_minusC = V[P_minus_nodes[2]]
		node_PchA = P_plus_nodes[0]
		node_PdA = P_minus_nodes[0]
		node_PchB = P_plus_nodes[1]
		node_PdB = P_minus_nodes[1]
		node_PchC = P_plus_nodes[2]
		node_PdC = P_minus_nodes[2]

		total_Pch = (P_plusA + P_plusB + P_plusC)
		total_Pd = (P_minusA + P_minusB + P_minusC)

		Bt_hist = self.Bt_prev - Bt + (total_Pch) * self.Mch - (total_Pd) * self.Md
		idx_Y = self.stampY(node_Bt, node_Bt, partials["dBtprevdBt"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PchA, partials["dBtprevdPch"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PdA, partials["dBtprevdPd"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PchB, partials["dBtprevdPch"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PdB, partials["dBtprevdPd"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PchC, partials["dBtprevdPch"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_PdC, partials["dBtprevdPd"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		idx_J = self.stampJ(
			node_Bt,
			(-Bt_hist + partials["dBtprevdBt"] * Bt + partials["dBtprevdPch"] * (total_Pch) + partials["dBtprevdPd"] * (total_Pd))
			, Jnlin_val, Jnlin_row, idx_J)


	else:

		# Equality Stamp for SOC is Based on One Phase

		Bt_hist = self.Bt_prev - Bt + P_plus * self.Mch - P_minus * self.Md
		node_Pch = P_plus_nodes[0]
		node_Pd = P_minus_nodes[0]
		idx_Y = self.stampY(node_Bt, node_Bt, partials["dBtprevdBt"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_Pch, partials["dBtprevdPch"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(node_Bt, node_Pd, partials["dBtprevdPd"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		idx_J = self.stampJ(
			node_Bt,
			(-Bt_hist + partials["dBtprevdBt"] * Bt + partials["dBtprevdPch"] * P_plus + partials["dBtprevdPd"] * P_minus)
			, Jnlin_val, Jnlin_row, idx_J)




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


	P_plus_nodes = self.P_plus_nodes
	P_minus_nodes = self.P_minus_nodes
	dual_ineq_r_plus_nodes = self.dual_ineq_r_plus_nodes
	dual_ineq_r_minus_nodes = self.dual_ineq_r_minus_nodes

	if (self.source_type == "PQ"):
		dual_ineq_i_nodes = self.dual_ineq_i_nodes

	num_of_phases = len(P_plus_nodes)

	#print("Battery Power Distribution and SOC Variables:")

	phases = ["A","B","C"]
	iter_phases = range(num_of_phases)
	if (self.single_phase):
		iter_phases = [phases.index(self.single_phase)]


	for i,index in enumerate(iter_phases):
		# i is the index within the battery phase list
		# index is the index within the phase voltages

		Vr = V[node_Vr[index]]
		Vi = V[node_Vi[index]]
		P_plus = V[P_plus_nodes[i]] 
		P_minus = V[P_minus_nodes[i]]
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
		if self.source_type == 'P':
			idx_Y= self.stampY(node_Lr[index], P_plus_nodes[i],  partials['d2Ir_dVrdP']*Lr + 
							partials['d2Ii_dVrdP']*Li, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Lr[index], P_minus_nodes[i],  -1*partials['d2Ir_dVrdP']*Lr + 
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
		if self.source_type == 'P':
			idx_Y= self.stampY(node_Li[index], P_plus_nodes[i], (partials['d2Ir_dVidP']*Lr 
							+ partials['d2Ii_dVidP']*Li), Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y= self.stampY(node_Li[index], P_minus_nodes[i], -1*(partials['d2Ir_dVidP']*Lr 
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


		pch_l2 = True
		pd_l2 = False
		c1 = self.C_ch
		c2 = self.C_d
		# --- Pch
		# Equations for primal variable P+
		# dL/dP = Lr*dIpqr/dP + Li*dIpqi/dP + df_obj/dP + mu
		# new objective: f_obj = 0.5*C1*P_plus^2 - 0.5*C2*P_minus^2

		hist_dLdPplus = Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] - V[dual_ineq_r_plus_nodes[i]] + c1
		idx_Y= self.stampY(P_plus_nodes[i], node_Vr[index], (partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li), 
						Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_plus_nodes[i], node_Vi[index], (partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li), 
						Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_plus_nodes[i], node_Lr[index], 
						partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_plus_nodes[i], node_Li[index], 
						partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_plus_nodes[i], dual_ineq_r_plus_nodes[i], 
						-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		if (pch_l2 == True):
			idx_Y = self.stampY(P_plus_nodes[i], P_plus_nodes[i], c1,
						Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			hist_dLdPplus = Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] - V[dual_ineq_r_plus_nodes[i]] + c1 * P_plus
			idx_J= self.stampJ(P_plus_nodes[i], -hist_dLdPplus + Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] +
				(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
				(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi + c1 * P_plus
				- V[dual_ineq_r_plus_nodes[i]], Jnlin_val, Jnlin_row, idx_J)
		else:
			idx_J= self.stampJ(P_plus_nodes[i], -hist_dLdPplus + Lr*partials['dIr_dP'] + Li*partials['dIi_dP'] +
							(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
							(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi
							- V[dual_ineq_r_plus_nodes[i]], Jnlin_val, Jnlin_row, idx_J)

		# --- Pd
		# Equations for primal variable P-
		# dL/dP = -Lr*dIpqr/dP - Li*dIpqi/dP - df_obj/dP + mu
		# new objective: f_obj = 0.5*C1*P_plus^2 - 0.5*C2*P_minus^2

		hist_dLdPminus = -Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] - V[dual_ineq_r_minus_nodes[i]] + c2
		idx_Y= self.stampY(P_minus_nodes[i], node_Vr[index], -1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li), 
						Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_minus_nodes[i], node_Vi[index], -1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li), 
						Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_minus_nodes[i], node_Lr[index], 
						-1*partials['dIr_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_minus_nodes[i], node_Li[index], 
						-1*partials['dIi_dP'], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y= self.stampY(P_minus_nodes[i], dual_ineq_r_minus_nodes[i], 
						-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		if (pd_l2 == True):
			idx_Y = self.stampY(P_minus_nodes[i], P_minus_nodes[i], c2,
							Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			hist_dLdPminus = -Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] - V[dual_ineq_r_minus_nodes[i]] + c2 * P_minus
			idx_J= self.stampJ(P_minus_nodes[i], -hist_dLdPminus - Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] +
				-1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
				-1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi + c2 * P_minus
				- V[dual_ineq_r_minus_nodes[i]], Jnlin_val, Jnlin_row, idx_J)
		else:
			idx_J= self.stampJ(P_minus_nodes[i], -hist_dLdPminus - Lr*partials['dIr_dP'] - Li*partials['dIi_dP'] +
							-1*(partials['d2Ir_dPdVr']*Lr + partials['d2Ii_dPdVr']*Li)*Vr + 
							-1*(partials['d2Ir_dPdVi']*Lr + partials['d2Ii_dPdVi']*Li)*Vi
							- V[dual_ineq_r_minus_nodes[i]], Jnlin_val, Jnlin_row, idx_J)
		
		dual_ineq_r_plus_nodes_upper = self.dual_ineq_r_plus_nodes_upper
		dual_ineq_r_minus_nodes_upper = self.dual_ineq_r_minus_nodes_upper

		# Pch
		idx_Y= self.stampY(P_plus_nodes[i], dual_ineq_r_plus_nodes_upper[i], 
		1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Pd
		idx_Y = self.stampY(P_minus_nodes[i], dual_ineq_r_minus_nodes_upper[i], 
		1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Battery Equality constraint Pch/Pd sensitivities

		node_Lb = self.dual_Bt_nodes[0]
		idx_Y = self.stampY(P_plus_nodes[i], node_Lb, self.Mch, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stampY(P_minus_nodes[i], node_Lb, -self.Md, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		if (self.verbose):
			print("phase {} Pch {} Pd {}".format(index, P_plus, P_minus))


	# Battery SOC
	node_Lb = self.dual_Bt_nodes[0]
	node_Bt = self.Bt_nodes[0]
	Lb = V[node_Lb]
	hist_Lb = -Lb
	Bt = V[self.Bt_nodes[0]]

	idx_Y = self.stampY(node_Lb, node_Lb, partials["dLbdlb"], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	idx_J = self.stampJ(node_Lb, hist_Lb - Lb * partials["dLbdlb"], Jnlin_val, Jnlin_row, idx_J)

	# Stationarity Dual Inequality Constraints
	Bt_nodes = self.Bt_nodes
	dual_Bt_nodes = self.mu_index_Bt
	dual_Bt_nodes_upper = self.mu_index_Bt_upper

	# SOC lower and upper bounds
	# lower inequality constraint
	idx_Y= self.stampY(node_Lb, dual_Bt_nodes[0], 
	-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	# upper inequality constraint
	idx_Y= self.stampY(node_Lb, dual_Bt_nodes_upper[0], 
	1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	if (self.verbose):
		print("Bt: {}".format(Bt))


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


	# --- BATTERY-EXCLUSIVE PARTIALS ---
	partials["dBtprevdBt"] = -1
	partials["dBtprevdPch"] = self.Mch
	partials["dBtprevdPd"] = -self.Md
	partials["dLbdlb"] = -1

	return partials
