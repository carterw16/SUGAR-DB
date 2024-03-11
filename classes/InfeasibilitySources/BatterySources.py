"""
Three-phase battery system classes

Author(s): Meshach Hopkins
Created Date: 02-13-2024
Updated Date: 02-13-2024
Email: mhopkins@andrew.cmu.edu
Status: Development

Parses the power grid elements, initializes their respective classes, and then returns a list for each element.

"""

from types import MethodType
from classes.InfeasibilitySources.base import InfeasibilitySources
from .BatteryStamps import (stampY, stampJ, stamp_nonlinear_infeasibility,
									stamp_complementary_slackness, stamp_equality_constraints,
									stamp_stationarity_constraints, calculate_partial_derivatives)

import ipdb


class BatterySources(InfeasibilitySources):
	def __init__(self, node, node_index, obj_scalar, source_type, stamp_neutral_infeas_source_flag,
			  ID="bat", P_max = 3000000, P_min = 0, Mch = 0.0000005, Md = 0.0000005, Bt_prev = 0,
			  C_ch = 1, C_d = -0.5, single_phase = "", verbose = True):
		
		InfeasibilitySources.__init__(self, node)
		self.source_type = source_type
		# the source type for batteries are "P" and "PQ"

		# initialize battery parameters
		self.ID = ID
		self.P_max = P_max
		self.P_min = P_min
		self.Mch = Mch
		self.Md = Md
		self.C_ch = C_ch
		self.C_d = C_d
		self.Bt_prev = Bt_prev
		self.single_phase = single_phase
		self.verbose = verbose

		# stamping functions
		self.stamp_neutral = stamp_neutral_infeas_source_flag
		self.stampY = MethodType(stampY, self)
		self.stampJ = MethodType(stampJ, self)
		self.stamp_nonlinear_infeasibility = MethodType(stamp_nonlinear_infeasibility, self)
		self.stamp_complementary_slackness = MethodType(stamp_complementary_slackness, self)
		self.stamp_equality_constraints = MethodType(stamp_equality_constraints, self)
		self.stamp_stationarity_constraints = MethodType(stamp_stationarity_constraints, self)
		self.calculate_partial_derivatives = MethodType(calculate_partial_derivatives, self)

		# assign and update nodes
		node_index = self.assign_nodes(node, node_index)
		self.updated_node_index = node_index

		self.obj_scaling = obj_scalar #1/(node.Vnom**2)


	def initialize(self, Vinit):
		if (self.single_phase == ""):
			Vinit[self.nodeA_p_plus] = 0.0001
			Vinit[self.nodeB_p_plus] = 0.0001
			Vinit[self.nodeC_p_plus] = 0.0001

			Vinit[self.nodeA_p_minus] = 500
			Vinit[self.nodeB_p_minus] = 500
			Vinit[self.nodeC_p_minus] = 500

			Vinit[self.nodeA_dual_ineq_var_p_plus] = 0.0001
			Vinit[self.nodeB_dual_ineq_var_p_plus] = 0.0001
			Vinit[self.nodeC_dual_ineq_var_p_plus] = 0.0001

			Vinit[self.nodeA_dual_ineq_var_p_plus_upper] = 0.0001
			Vinit[self.nodeB_dual_ineq_var_p_plus_upper] = 0.0001
			Vinit[self.nodeC_dual_ineq_var_p_plus_upper] = 0.0001

			Vinit[self.nodeA_dual_ineq_var_p_minus] = 0.0001
			Vinit[self.nodeB_dual_ineq_var_p_minus] = 0.0001
			Vinit[self.nodeC_dual_ineq_var_p_minus] = 0.0001

			Vinit[self.nodeA_dual_ineq_var_p_minus_upper] = 0.0001
			Vinit[self.nodeB_dual_ineq_var_p_minus_upper] = 0.0001
			Vinit[self.nodeC_dual_ineq_var_p_minus_upper] = 0.0001
		else:
			Vinit[self.node_p_plus] = 0.0001

			Vinit[self.node_p_minus] = 500

			Vinit[self.node_dual_ineq_var_p_plus] = 0.0001

			Vinit[self.node_dual_ineq_var_p_plus_upper] = 0.0001

			Vinit[self.node_dual_ineq_var_p_minus] = 0.0001

			Vinit[self.node_dual_ineq_var_p_minus_upper] = 0.0001

		# SOC Initialization

		Vinit[self.Bt] = self.Bt_prev
		Vinit[self.lambda_Bt] = 0.0001
		Vinit[self.node_dual_ineq_Bt] = 0.0001
		Vinit[self.node_dual_ineq_Bt_upper] = 0.0001
		
		return Vinit


	def assign_nodes(self, node, node_index):
		# So that batteries are isolated from Infeasibility currents and several can be placed at an individual node
		self.isTriplex = node.isTriplex

		if (self.isTriplex == False):
			if (self.single_phase == ""):
				# Three-Phase Battery

				self.nodeA_p_plus = node_index.__next__()
				self.nodeB_p_plus = node_index.__next__()
				self.nodeC_p_plus = node_index.__next__()

				self.nodeA_p_minus = node_index.__next__()
				self.nodeB_p_minus = node_index.__next__()
				self.nodeC_p_minus = node_index.__next__()

				self.nodeA_dual_ineq_var_p_plus = node_index.__next__()
				self.nodeB_dual_ineq_var_p_plus = node_index.__next__()
				self.nodeC_dual_ineq_var_p_plus = node_index.__next__()

				self.nodeA_dual_ineq_var_p_plus_upper = node_index.__next__()
				self.nodeB_dual_ineq_var_p_plus_upper = node_index.__next__()
				self.nodeC_dual_ineq_var_p_plus_upper = node_index.__next__()

				self.nodeA_dual_ineq_var_p_minus = node_index.__next__()
				self.nodeB_dual_ineq_var_p_minus = node_index.__next__()
				self.nodeC_dual_ineq_var_p_minus = node_index.__next__()

				self.nodeA_dual_ineq_var_p_minus_upper = node_index.__next__()
				self.nodeB_dual_ineq_var_p_minus_upper = node_index.__next__()
				self.nodeC_dual_ineq_var_p_minus_upper = node_index.__next__()


				# Power distribution Nodes +/-
				self.P_plus_nodes = [self.nodeA_p_plus, self.nodeB_p_plus,
					self.nodeC_p_plus]
				self.P_minus_nodes = [self.nodeA_p_minus, self.nodeB_p_minus,
						self.nodeC_p_minus]
				
				# Lower Inequality Nodes set
				self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_var_p_plus, self.nodeB_dual_ineq_var_p_plus,
						self.nodeC_dual_ineq_var_p_plus]
				self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_var_p_minus, self.nodeB_dual_ineq_var_p_minus,
						self.nodeC_dual_ineq_var_p_minus]
				
				# Upper Inequality Nodes set
				self.dual_ineq_r_plus_nodes_upper = [self.nodeA_dual_ineq_var_p_plus_upper, 
									self.nodeB_dual_ineq_var_p_plus_upper,self.nodeC_dual_ineq_var_p_plus_upper]
				self.dual_ineq_r_minus_nodes_upper = [self.nodeA_dual_ineq_var_p_minus_upper, 
									self.nodeB_dual_ineq_var_p_minus_upper,self.nodeC_dual_ineq_var_p_minus_upper]
			else:
				# Single-Phase Battery

				self.node_p_plus = node_index.__next__()

				self.node_p_minus = node_index.__next__()

				self.node_dual_ineq_var_p_plus = node_index.__next__()

				self.node_dual_ineq_var_p_plus_upper = node_index.__next__()

				self.node_dual_ineq_var_p_minus = node_index.__next__()

				self.node_dual_ineq_var_p_minus_upper = node_index.__next__()


				# Power distribution Nodes +/-
				self.P_plus_nodes = [self.node_p_plus]
				self.P_minus_nodes = [self.node_p_minus]
				
				# Lower Inequality Nodes set
				self.dual_ineq_r_plus_nodes = [self.node_dual_ineq_var_p_plus]
				self.dual_ineq_r_minus_nodes = [self.node_dual_ineq_var_p_minus]
				
				# Upper Inequality Nodes set
				self.dual_ineq_r_plus_nodes_upper = [self.node_dual_ineq_var_p_plus_upper]
				self.dual_ineq_r_minus_nodes_upper = [self.node_dual_ineq_var_p_minus_upper]


			# SOC equality terms
			# SOC evolution equality constraints
			self.Bt = node_index.__next__()
			self.lambda_Bt = node_index.__next__()

			# Upper and Lower Inequality bounds for SOC
			# lowerer inequality bounds
			self.node_dual_ineq_Bt = node_index.__next__()

			# upper inequality bounds
			self.node_dual_ineq_Bt_upper = node_index.__next__()

			self.Bt_nodes = [self.Bt]
			self.dual_Bt_nodes = [self.lambda_Bt]

			# SOC dynamics full set of inequalities
			self.mu_index_Bt = [self.node_dual_ineq_Bt]
			self.mu_index_Bt_upper = [self.node_dual_ineq_Bt_upper]


			# Power Distribution full set
			self.infeas_real_var_index = self.P_plus_nodes + self.P_minus_nodes
			self.mu_index = self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes
			self.mu_index_upper = self.dual_ineq_r_plus_nodes_upper + self.dual_ineq_r_minus_nodes_upper


			"""
			self.P_plus_nodes = node.P_plus_nodes
			self.P_minus_nodes = node.P_minus_nodes
			# lower-bound nodes
			self.dual_ineq_r_plus_nodes = node.dual_ineq_r_plus_nodes
			self.dual_ineq_r_minus_nodes = node.dual_ineq_r_minus_nodes
			# upper-bound nodes
			self.dual_ineq_r_plus_nodes_upper = [node.nodeA_dual_ineq_var_p_plus_upper, node.nodeB_dual_ineq_var_p_plus_upper, node.nodeC_dual_ineq_var_p_plus_upper]
			self.dual_ineq_r_minus_nodes_upper = [node.nodeA_dual_ineq_var_p_minus_upper, node.nodeB_dual_ineq_var_p_minus_upper, node.nodeC_dual_ineq_var_p_minus_upper]
			"""


			if (self.source_type == 'PQ'):
				# TODO: define Triplex node assignments
				"""
				#self.Q_nodes = node.Q_nodes
				self.Q_plus_nodes = node.Q_plus_nodes
				self.Q_minus_nodes = node.Q_minus_nodes
				# lower-bound nodes
				self.dual_ineq_i_nodes = node.dual_ineq_i_nodes
				# upper-bound nodes
				# TODO: have to add these nodes and store them internally in battery object
				self.dual_ineq_i_nodes_upper = [node.nodeA_dual_ineq_var_p_plus_upper, node.nodeB_dual_ineq_var_p_plus_upper, node.nodeC_dual_ineq_var_p_plus_upper]
				"""
				raise Exception("Battery placed at Triplex node [Functionality Not Implemented]")
		else:
			# TODO: define Triplex node assignments
			raise Exception("Battery placed at Triplex node [Functionality Not Implemented]")

		return node_index

	def calc_residuals(self, node, V, res_eqn, cs_tol):
		#
		#
		# TODO: calc_residuals function needs to be adjusted
		#
		#
		if node.isTriplex:
			node_Vr = [node.node1_Vr, node.node2_Vr, node.nodeN_Vr]
			node_Vi = [node.node1_Vi, node.node2_Vi, node.nodeN_Vi]
		else:
			node_Vr = [node.nodeA_Vr, node.nodeB_Vr, node.nodeC_Vr, node.nodeN_Vr]
			node_Vi = [node.nodeA_Vi, node.nodeB_Vi, node.nodeC_Vi, node.nodeN_Vi]

		node_Lr = node.dual_eq_var_r_nodes
		node_Li = node.dual_eq_var_i_nodes


		if self.obj == 'L1':
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
				Lr = V[node_Lr[index]]
				Li = V[node_Li[index]]

				P_plus = V[P_plus_nodes[index]] if self.source_type == 'PQ' else 0
				P_minus = V[P_minus_nodes[index]] if self.source_type == 'PQ' else 0
				Q_plus = 0 #V[Q_plus_nodes[index]]
				Q_minus = 0 #V[Q_minus_nodes[index]]

				mu_q_ub = V[dual_ineq_i_plus_nodes[index]]
				mu_q_lb = V[dual_ineq_i_minus_nodes[index]]

				partials = self.calculate_partial_derivatives(Vr, Vi, (P_plus - P_minus), (Q_plus - Q_minus), calc_hessian = False)

				res_eqn[node_Vr[index]] += partials['Ir']
				res_eqn[node_Vi[index]] += partials['Ii']

				res_eqn[node_Lr[index]] += partials['dIr_dVr']*Lr + partials['dIi_dVr']*Li
				res_eqn[node_Li[index]] += partials['dIr_dVi']*Lr + partials['dIi_dVi']*Li

				#res_eqn[Q_plus_nodes[index]] += partials['dIr_dQ']*Lr + partials['dIi_dQ']*Li + 1 - mu_q_ub
				#res_eqn[Q_minus_nodes[index]] += -partials['dIr_dQ']*Lr - partials['dIi_dQ']*Li + 1 - mu_q_lb

				res_eqn[dual_ineq_i_plus_nodes[index]] += -mu_q_ub*Q_plus + cs_tol
				res_eqn[dual_ineq_i_minus_nodes[index]] += -mu_q_lb*Q_minus + cs_tol

				mu_p_ub = V[dual_ineq_r_plus_nodes[index]]
				mu_p_lb = V[dual_ineq_r_minus_nodes[index]]

				res_eqn[P_plus_nodes[index]] += partials['dIr_dP']*Lr + partials['dIi_dP']*Li + 1 - mu_p_ub
				res_eqn[P_minus_nodes[index]] += -partials['dIr_dP']*Lr - partials['dIi_dP']*Li + 1 - mu_p_lb

				res_eqn[dual_ineq_r_plus_nodes[index]] += -mu_p_ub*P_plus + cs_tol
				res_eqn[dual_ineq_r_minus_nodes[index]] += -mu_p_lb*P_minus + cs_tol

		return res_eqn