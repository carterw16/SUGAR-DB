""" Implements tx stepping homotopy.

  Author(s): Elizabeth Foster
  Created Date: 04-05-2021
  Updated Date: 06-14-2021
  Email: emfoster@andrew.cmu.edu
  Status: Development

"""

import logging
import numpy as np
from classes.Regulators import VariableRegulator


class Homotopy:

	def __init__(self, tol, g_homotopy, b_homotopy, h_factor_init,
				 h_step_size_init):
		'''

		Parameters
		----------
		tol : scalar value; represents the tolerance for h_factor convergence
		g_homotopy : scalar value; represents the eta value multiplied by all conductances in tx_stepping
		b_homotopy : scalar value; represents the eta value multiplied by all susceptances in tx_stepping
		h_factor_init : intial h_factor, default is 1
		h_step_size_init : initial step size, default is 0.1

		Returns
		-------
		None.

		'''
		self.tolerance = tol

		self.g_value = g_homotopy
		self.b_value = b_homotopy

		self.h_factor = h_factor_init
		self.h_factor_last_converged = []

		self.step_size = h_step_size_init
		self.step_size_last_converged = h_step_size_init
		self.num_of_steps = 0

		self.last_V_converged = []
		self.Verr_prev = 1e10
		self.Vinit = []
		self.perturb_conditions = False

		self.tx_stepping_successful = False
		self.tx_stepping_failed = False

		self.flag_zero_fails = False
		self.rate_flag_fail = False
		self.rates = [0.11, 0.251, 0.51, 0.751, 0.91, 0.99, 0.999]

		# Stuck point settings
		self.consecutive_successes = 10
		self.stuck_point = False
		self.stuck_point_high_factor = 0.9999 #0.9999 # 0.99999
		self.stuck_point_medium_factor = 0.999
		self.stuck_point_low_factor = 0.9
		self.stuck_point_high_steps = 2
		self.stuck_point_low_steps = 7


	def increase_rate(self):
		''' This function is called to change the amount that the homotopy factor 
		decreases after failing '''

		# The homotopy factor that last converged
		current_rate = self.h_factor/self.h_factor_last_converged
		current_rate = round(current_rate, 8)
		
		for i in range(len(self.rates)):
			if abs(self.rates[i] - current_rate) < 1e-4:
				current_rate = self.rates[i]
				
			if i == 0:
				if current_rate < self.rates[i]:
					if current_rate * 1.5 > self.rates[i] and current_rate != 0:
						new_rate = current_rate*1.5
					else:
						new_rate = self.rates[i]
					break
				elif current_rate == self.rates[i]:
					new_rate = self.rates[i+1]
					break
					
			if i == len(self.rates) - 1:
				if current_rate * 1.5 > 1:
					self.rate_flag_fail = True
					new_rate = 1
				else:
					new_rate = current_rate * 1.5
				break
				
			if current_rate == self.rates[i]:
				new_rate = self.rates[i+1]
				break
			
			if i != len(self.rates) - 1:
				if current_rate < self.rates[i+1] and current_rate > self.rates[i]:
					if current_rate * 1.5 <= self.rates[i+1]:
						new_rate = current_rate*1.5
						break
					elif (self.rates[i+1] - current_rate) < 1e-4:
						continue
					else:
						new_rate = self.rates[i+1]
						break
		
		self.new_increase_rate = new_rate
					
	def decrease_rate(self):
		''' This function is called to change the amount that the homotopy factor 
		decreases after failing '''

		# The homotopy factor that last converged
		current_rate = self.h_factor*0.5
		new_rate = round(current_rate, 8)
		
		self.new_decrease_rate = new_rate

	def increase_h_factor(self):
		'''
		Purpose
		--------
		This function increases the homotopy factor if the N-R hasn't converged for that
		value of the homotopy factor within a set number of iterations

		Returns
		-------
		None.

		'''
		self.increase_rate()
		h_difference = self.h_factor_last_converged - self.h_factor_last_converged*self.new_increase_rate

		if h_difference < 1e-7:
			self.tx_stepping_failed = True

		else:
			self.h_factor = self.h_factor_last_converged * self.new_increase_rate
			self.h_factor = round(self.h_factor, 8)

			if self.h_factor > 1:
				self.h_factor = 1
			if self.h_factor < 0:
				self.h_factor = 0

			logging.info("Increasing homotopy factor to %0.8f", round(self.h_factor,8))            

	def decrease_h_factor(self):
		'''
		Purpose
		--------
		This function decreases the homotopy factor if the N-R has converged for that
		value of the homotopy factor

		Returns
		-------
		None.

		'''
		# The N-R loop for a given homotopy factor has converged. Now we
		# want to step the homotopy factor down.	
		if self.flag_zero_fails == False and self.stuck_point == False:
			self.h_factor_last_converged = self.h_factor
			self.step_size_last_converged = self.step_size
			self.h_factor -= self.step_size
			self.h_factor = round(self.h_factor, 8)
			if abs(self.h_factor) <= 1e-7 or self.h_factor < 0:
				self.h_factor = 0
		
		elif self.stuck_point == True:
			if self.consecutive_successes < self.stuck_point_high_steps:
				self.step_size_last_converged = self.h_factor_last_converged - self.h_factor
				self.h_factor_last_converged = self.h_factor
				self.h_factor = self.h_factor*self.stuck_point_high_factor # 0.99999 
				self.step_size = round(self.h_factor_last_converged - self.h_factor, 12)
			
			elif self.consecutive_successes >= self.stuck_point_high_steps and self.consecutive_successes < self.stuck_point_low_steps:
				self.step_size_last_converged = self.h_factor_last_converged - self.h_factor
				self.h_factor_last_converged = self.h_factor
				self.h_factor = self.h_factor*self.stuck_point_medium_factor 
				self.step_size = round(self.h_factor_last_converged - self.h_factor, 12)
			
			else:
				self.step_size_last_converged = self.h_factor_last_converged - self.h_factor
				self.h_factor_last_converged = self.h_factor
				self.h_factor = self.h_factor*self.stuck_point_low_factor 
				self.step_size = round(self.h_factor_last_converged - self.h_factor, 12)
				if self.consecutive_successes > self.stuck_point_low_steps:
					self.stuck_point == False
			
		else:
			self.decrease_rate()
			self.step_size_last_converged = self.h_factor_last_converged - self.h_factor
			self.h_factor_last_converged = self.h_factor
			self.h_factor = self.h_factor - self.step_size_last_converged #self.new_decrease_rate
			self.step_size = round(self.h_factor_last_converged - self.h_factor, 8)
			self.flag_zero_fails = False

			if self.h_factor > 1:
				self.h_factor = 1
			if self.h_factor < 0:
				self.h_factor = 0


		logging.info('Decreasing homotopy factor to %0.12f' %
					 round(self.h_factor, 12))

	def run_tx_stepping(self, Verr_max, Vsol):
		'''
		Purpose
		--------
		This function increments the homotopy factor as necessary during the
		Newton-Raphson process.

		Parameters
		----------
		Verr_max: scalar
		Vsol: solution vector with values calculated from Y and J for this iteration

		Returns
		-------
		None.

		'''

		# Cases in tx-stepping:
		# 1. No convergence in a set number of iterations for a given value of the homotopy factor
		# 2. Newton-Raphson is oscillating between different error values [no convergence]
		# 3. Error has converged and the homotopy factor can be decreased
		# 4. Error has converved and the homotopy factor was 0, so tx stepping was successful!
		# 5. No convergence and homotopy factor is so small that convergence is unlikely

		if Verr_max <= self.tolerance:
			# Case 4
			if self.h_factor == 0:
				self.tx_stepping_successful = True
			# Case 3
			else:
				self.decrease_h_factor()
				self.num_of_steps = 0
				self.Verr_prev = 1e10
				self.last_V_converged = Vsol
				self.consecutive_successes += 1
		
		elif self.num_of_steps >= 1:
			if self.num_of_steps >= 8 and self.Verr_prev < Verr_max and self.h_factor < 1:
				if self.h_factor == 0:
					self.flag_zero_fails = True
				self.num_of_steps = 0
				self.increase_h_factor()
				Vsol = self.last_V_converged
				self.consecutive_successes = 0

			# Case 1 or 5
			elif self.num_of_steps >= 20 and Verr_max > self.tolerance and self.h_factor < 1:
				self.num_of_steps = 0
				self.increase_h_factor()
				Vsol = self.last_V_converged
				self.consecutive_successes = 0
				if self.h_factor == 0:
					self.flag_zero_fails = True

			elif self.num_of_steps >= 8 and self.h_factor == 1:
				# For when there is no convergence for h = 1, try increasing the
				# B and G values
				self.g_value = self.g_value * 1.5
				self.b_value = self.b_value * 1.5
				logging.info('Increasing B and G values')
				self.num_of_steps = 0
			
			elif abs(self.Verr_prev - Verr_max) < 1e-6:
				self.num_of_steps = 0
				self.increase_h_factor()
				Vsol = self.last_V_converged
				self.consecutive_successes = 0
				if self.h_factor == 0:
					self.flag_zero_fails = True

			else:
				self.num_of_steps += 1
				self.Verr_prev = Verr_max
		else:
			self.num_of_steps += 1
			self.Verr_prev = Verr_max

		if self.consecutive_successes > 9:
			self.stuck_point = False

		if self.tx_stepping_failed == True:
			self.stuck_point = True
			self.perturb_conditions = True
			self.consecutive_successes = 0
			#self.step_size = 0.1
			#self.step_size_converged = 0.1
			self.tx_stepping_failed = False
			self.h_factor = self.h_factor_last_converged * 0.999999
			Vsol = self.last_V_converged
			self.Verr_prev = 1e10
			self.num_of_steps = 0
			self.step_size = self.h_factor * 0.5
			logging.info("Homotopy failed. Perturbing initial conditions")
			logging.info('Restarting homotopy factor at %0.12f' %
					 round(self.h_factor, 12))
			self.consecutive_successes = 0
			#self.h_factor = 0.5
			# elif self.perturbGB_flag == True:
			#     self.g_value  = self.g_value * 1.5
			#     self.b_value = self.b_value * 1.5
			#     if self.g_value > 1000 or self.b_value > 1000:
			#         logging.info("Homotopy has failed at homotopy factor %0.8f" %self.h_factor)  
			#     else:          
			#         self.step_size = 1
			#         self.step_size_last_converged = 1
			#         self.tx_stepping_failed = False
			#         self.h_factor = 1
			#         Vsol = self.last_V_converged

		
		return Vsol

	def stamp_tx_homotopy(self, nodes, node_key, overheadline, undergroundline,
						  triplexline, xfmr, stamped_ground, Y_val, Y_row,
						  Y_col, obj, J_val, J_row, V):
		'''
		Purpose
		--------
		Stamp all of the values associated with tx_stepping

		Parameters
		----------
		nodes : from Nodes class; size is number of nodes in original system
		overheadline : from Lines class; size is number of lines in original system
		undergroundline : from Lines class; size is number of lines in original system
		triplexline : from Lines class; size is number of lines in original system
		xfmr : from Transformer class; size is number of transformers in original system
		stamped_ground : indices of ground locations
		Y_val : arbitrarily large empty column vector.
		Y_row : arbitrarily large empty column vector.
		Y_col : arbitrarily large empty column vector.

		Returns
		-------
		Y_row : column vector the same size as input Y_row but it now has stamped row indices
		Y_col : column vector the same size as input Y_col but it now has stamped column indices
		Y_val : column vector the same size as input Y_rval but it now has stamped values
		Y_idx : index location for the next stamp

		'''

		Y_idx = 0
		for ele in range(len(overheadline)):
			Y_val, Y_row, Y_col, Y_idx = overheadline[ele].stamp_tx_stepping(node_key, nodes, Y_val, Y_row, Y_col,
																			 Y_idx, self.h_factor, self.g_value,
																			 self.b_value)

		for ele in range(len(undergroundline)):
			Y_val, Y_row, Y_col, Y_idx = undergroundline[ele].stamp_tx_stepping(node_key, nodes, Y_val, Y_row, Y_col,
																				Y_idx, self.h_factor, self.g_value,
																				self.b_value)

		for ele in range(len(triplexline)):
			Y_val, Y_row, Y_col, Y_idx = triplexline[ele].stamp_tx_stepping(node_key, nodes, Y_val, Y_row, Y_col, Y_idx,
																			self.h_factor, self.g_value, self.b_value)

		for ele in range(len(xfmr)):
			Y_val, Y_row, Y_col, Y_idx, stamped_ground = xfmr[ele].stamp_linear(
				Y_val, Y_row, Y_col, Y_idx, stamped_ground, True, self.h_factor,
				self.g_value, self.b_value)
		
		if obj == 'L1':
			alpha = -1e-6/2 
			idx_J = 0
			for ele in range(len(nodes)):
				if nodes[ele].isTriplex and nodes[ele].bustype != 3:
					Vr = [nodes[ele].node1_Vr,  nodes[ele].node2_Vr,  nodes[ele].nodeN_Vr] 
					Vi = [nodes[ele].node1_Vi,  nodes[ele].node2_Vi,  nodes[ele].nodeN_Vi] 

					if_r_plus =  [nodes[ele].node1_if_r_plus,  nodes[ele].node2_if_r_plus,  nodes[ele].nodeN_if_r_plus]        
					if_r_minus = [nodes[ele].node1_if_r_minus, nodes[ele].node2_if_r_minus, nodes[ele].nodeN_if_r_minus]
					if_i_plus =  [nodes[ele].node1_if_i_plus,  nodes[ele].node2_if_i_plus,  nodes[ele].nodeN_if_i_plus]   
					if_i_minus = [nodes[ele].node1_if_i_minus, nodes[ele].node2_if_i_minus, nodes[ele].nodeN_if_i_minus]  

					Lr = [nodes[ele].node1_dual_eq_var_r,  nodes[ele].node2_dual_eq_var_r,  nodes[ele].nodeN_dual_eq_var_r] 
					Li = [nodes[ele].node1_dual_eq_var_i,  nodes[ele].node2_dual_eq_var_i,  nodes[ele].nodeN_dual_eq_var_i] 
					 
					mu_r_lb = nodes.dual_ineq_r_lb_nodes
					mu_i_lb = nodes.dual_ineq_i_lb_nodes
					mu_r_ub = nodes.dual_ineq_r_ub_nodes
					mu_i_ub = nodes.dual_ineq_i_ub_nodes     
					 
					for i in range(0,3):
# # # 						J_val, J_row, idx_J = stampJ(if_i_minus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_r_minus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_i_plus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_r_plus[i], -self.h_factor, J_val, J_row, idx_J) 
						
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], if_r_plus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], if_i_plus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], if_r_minus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], if_i_minus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], Lr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], Li[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], Lr[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], Li[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], if_r_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], if_r_minus[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], if_i_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], if_i_minus[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						J_val, J_row, idx_J = stampJ(mu_i_lb[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_r_lb[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_i_ub[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_r_ub[i], self.h_factor, J_val, J_row, idx_J)
						
						# Y_val, Y_row, Y_col, Y_idx = stampY(sr_plus[i], if_r_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si_plus[i], if_i_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(sr_minus[i], if_r_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si_minus[i], if_i_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], sr_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], sr_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], si_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], si_minus[i],  -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
						# J_val, J_row, idx_J = stampJ(sr_plus[i], self.h_factor, J_val, J_row, idx_J)
						# J_val, J_row, idx_J = stampJ(si_plus[i], self.h_factor, J_val, J_row, idx_J)
						# J_val, J_row, idx_J = stampJ(sr_minus[i], self.h_factor, J_val, J_row, idx_J)
						# J_val, J_row, idx_J = stampJ(si_minus[i], self.h_factor, J_val, J_row, idx_J)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_ub[i], mu_r_ub[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_lb[i], mu_r_lb[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_i_ub[i], mu_i_ub[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_lb[i], mu_i_lb[i], -alpha, Y_val, Y_row, Y_col, Y_idx)

						# Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], sr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], si[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						# Y_val, Y_row, Y_col, Y_idx = stampY(sr[i], sr[i], 2, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si[i], si[i], 2, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(sr[i], Lr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si[i], Li[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)

				
				elif nodes[ele].bustype != 3:   
					Vr = [nodes[ele].nodeA_Vr,  nodes[ele].nodeB_Vr, nodes[ele].nodeC_Vr,  nodes[ele].nodeN_Vr] 
					Vi = [nodes[ele].nodeA_Vi,  nodes[ele].nodeB_Vi, nodes[ele].nodeC_Vi, nodes[ele].nodeN_Vi] 

					sr = [nodes[ele].nodeA_sr,  nodes[ele].nodeB_sr, nodes[ele].nodeC_sr, nodes[ele].nodeN_sr] 
					si = [nodes[ele].nodeA_si,  nodes[ele].nodeB_si, nodes[ele].nodeC_si, nodes[ele].nodeN_si] 
					 
					if_r_plus =  [nodes[ele].nodeA_if_r_plus,  nodes[ele].nodeB_if_r_plus,  nodes[ele].nodeC_if_r_plus,  nodes[ele].nodeN_if_r_plus]        
					if_r_minus = [nodes[ele].nodeA_if_r_minus, nodes[ele].nodeB_if_r_minus, nodes[ele].nodeC_if_r_minus, nodes[ele].nodeN_if_r_minus]
					if_i_plus =  [nodes[ele].nodeA_if_i_plus,  nodes[ele].nodeB_if_i_plus,  nodes[ele].nodeC_if_i_plus,  nodes[ele].nodeN_if_i_plus]   
					if_i_minus = [nodes[ele].nodeA_if_i_minus, nodes[ele].nodeB_if_i_minus, nodes[ele].nodeC_if_i_minus, nodes[ele].nodeN_if_i_minus]  

					Lr = [nodes[ele].nodeA_dual_eq_var_r,  nodes[ele].nodeB_dual_eq_var_r, nodes[ele].nodeC_dual_eq_var_r, nodes[ele].nodeN_dual_eq_var_r] 
					Li = [nodes[ele].nodeA_dual_eq_var_i,  nodes[ele].nodeB_dual_eq_var_i, nodes[ele].nodeC_dual_eq_var_i, nodes[ele].nodeN_dual_eq_var_i]
					 
					mu_r_lb = nodes.dual_ineq_r_lb_nodes
					mu_i_lb = nodes.dual_ineq_i_lb_nodes
					mu_r_ub = nodes.dual_ineq_r_ub_nodes
					mu_i_ub = nodes.dual_ineq_i_ub_nodes     
					 
					for i in range(0,4):
# # # 						J_val, J_row, idx_J = stampJ(if_i_minus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_r_minus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_i_plus[i], -self.h_factor, J_val, J_row, idx_J) 
# # # 						J_val, J_row, idx_J = stampJ(if_r_plus[i], -self.h_factor, J_val, J_row, idx_J)

# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], if_r_plus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], if_i_plus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], if_r_minus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# # 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], if_i_minus[i], -2*self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], Lr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], Li[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], Lr[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], Li[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
						Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], if_r_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], if_r_minus[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], if_i_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], if_i_minus[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_ub[i], mu_r_ub[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_lb[i], mu_r_lb[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_i_ub[i], mu_i_ub[i], -alpha, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(mu_r_lb[i], mu_i_lb[i], -alpha, Y_val, Y_row, Y_col, Y_idx)

						J_val, J_row, idx_J = stampJ(mu_i_lb[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_r_lb[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_i_ub[i], self.h_factor, J_val, J_row, idx_J) 
						J_val, J_row, idx_J = stampJ(mu_r_ub[i], self.h_factor, J_val, J_row, idx_J)
						
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(sr_plus[i], if_r_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(si_plus[i], if_i_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(sr_minus[i], if_r_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(si_minus[i], if_i_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_plus[i], sr_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_r_minus[i], sr_minus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_plus[i], si_plus[i], -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
# 						# Y_val, Y_row, Y_col, Y_idx = stampY(if_i_minus[i], si_minus[i],  -self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						
# 						# J_val, J_row, idx_J = stampJ(sr_plus[i], self.h_factor, J_val, J_row, idx_J)
# 						# J_val, J_row, idx_J = stampJ(si_plus[i], self.h_factor, J_val, J_row, idx_J)
# 						# J_val, J_row, idx_J = stampJ(sr_minus[i], self.h_factor, J_val, J_row, idx_J)
# 						# J_val, J_row, idx_J = stampJ(si_minus[i], self.h_factor, J_val, J_row, idx_J)
						# Y_val, Y_row, Y_col, Y_idx = stampY(Vr[i], sr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(Vi[i], si[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)

						# Y_val, Y_row, Y_col, Y_idx = stampY(sr[i], sr[i], 2, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si[i], si[i], 2, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(sr[i], Lr[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
						# Y_val, Y_row, Y_col, Y_idx = stampY(si[i], Li[i], self.h_factor, Y_val, Y_row, Y_col, Y_idx)
		
		return (Y_row, Y_col, Y_val, Y_idx, stamped_ground)

def stampY(i, j, val, Y_val, Y_row, Y_col, idx):
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return Y_val, Y_row, Y_col, idx

def stampJ(i, val, J_val, J_row, idx):
	J_val[idx] = val
	J_row[idx] = i
	idx += 1

	return J_val, J_row, idx
