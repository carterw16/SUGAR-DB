""" Connection points in three-phase power flow

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-10-2017
  Updated Date: 10-14-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""

from __future__ import division

import logging
import operator
from itertools import count
from classes.GlobalVars import _PHASE
import numpy as np
import math
import ipdb

def calcVrange(nodes, Vsol):
	Va_min = []
	Vb_min = []
	Vc_min = []
	# calculate the min and max angles and voltages on all nodes
	for nodeEle in nodes:
		nodeEle.calcMagAng(Vsol, True)

	for nodeEle in nodes:
		Vmin_idx, Va = min(enumerate([nodeEle.V_mag[0]]))
		Vmin_idx, Vb = min(enumerate([nodeEle.V_mag[1]]))
		try:
			Vmin_idx, Vc = min(enumerate([nodeEle.V_mag[2]]))
		except IndexError:
			Vc = 0
		Va_min.append(Va)
		Vb_min.append(Vb)
		Vc_min.append(Vc)

	Va_min = np.asarray(Va_min)

	Vmax_node, Vmax = max(enumerate([nodeEle.V_mag for nodeEle in nodes]),
						  key=operator.itemgetter(1))
	Vmin = min(Va_min[np.nonzero(Va_min)])
	Vmax_ang_node, Vmax_ang = max(enumerate(
		[nodeEle.V_ang for nodeEle in nodes]),
								  key=operator.itemgetter(1))
	Vmin_ang_node, Vmin_ang = min(enumerate(
		[nodeEle.V_ang for nodeEle in nodes]),
								  key=operator.itemgetter(1))
	return Vmax_node, max(Vmax), Vmin_ang_node, Vmin, Vmax_ang_node, max(
		Vmax_ang), Vmin_ang_node, min(Vmin_ang)


class Nodes:
	ExistTriplex = False
	voltage_index = []
	Vr_index = []
	Vi_index = []
	Vr_mag = []
	Vi_mag = []
	V_mag = []
	Q_index = []
	I_index = []
	Tap_index = []
	Q_mag = []
	# Transformer Ground Nodes
	xfmr_primal_ground = []
	xfmr_dual_ground = []
	# Infeasibility Variable Indices
	if_index = []
	infeas_real_var_index = []
	infeas_imag_var_index = []
	if_r_index = []
	if_i_index = []
	if_r_plus_index = []
	if_r_minus_index = []
	if_i_plus_index = []
	if_i_minus_index = []
	# Equality Dual Variable Indices
	L_index = []
	Lr_index = []
	Li_index = []
	# Inequality Dual Variable Indices
	mu_index = []
	mu_r_index = []
	mu_i_index = []
	mu_index_upper = []
	mu_r_index_upper = []
	# Voltange Unbalance Variable Indices
	voltage_unbalance_primal_nodes = []
	voltage_unbalance_dual_nodes = []
	voltage_unbalance_mu_nodes = []
	voltage_unbalance_tracking_nodes = []
	# Voltage Bound Variable Indices
	vmag2_index = []
	Lvmag2_index = []
	umax_vmag2_index = []
	umin_vmag2_index = []
	vi_bounding_inds = []
	vr_bounding_inds = []


	@staticmethod
	def almostEqual(d1, d2, epsilon=10**-7):
		return abs(d2 - d1) < epsilon

	@staticmethod
	def _totuple(a):
		try:
			return tuple(i for i in a)
		except TypeError:
			return a

	@staticmethod
	def _convertline2phase(Vab, Vbc, Vca):
		_Veq = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
		_VLL = np.array([Vab, Vbc, Vca])
		_Vabc = np.linalg.solve(_Veq, _VLL)
		return Nodes._totuple(_Vabc)

	@staticmethod
	def sind(x):
		return np.sin(float(x) * np.pi / 180)

	@staticmethod
	def cosd(x):
		return np.cos(float(x) * np.pi / 180)

	def update_Vnom(self, Vsol):

		try:
			Va = abs(Vsol[self.nodeA_Vr] + 1j * Vsol[self.nodeA_Vi])
			err_nom = abs(Va - self.Vnom / np.sqrt(3))
			Vb = abs(Vsol[self.nodeB_Vr] + 1j * Vsol[self.nodeB_Vi])
			err_nom = abs(Vb - self.Vnom / np.sqrt(3))
		except AttributeError:
			Vc = abs(Vsol[self.nodeC_Vr] + 1j * Vsol[self.nodeC_Vi])
			err_nom = abs(Vc - self.Vnom / np.sqrt(3))

		try:
			self.Vnom = self.Vnom / np.sqrt(3)
		except err_nom > 0.10 * self.Vnom:
			self.Vnom = self.Vnom

	def __init__(self,
				 ID,
				 name,
				 phases,
				 nominal_voltage,
				 voltage_setpoint,
				 node_connect_type,
				 isTriplex=False,
				 bustype=0,
				 maxvoltage_error=1e-2,
				 busflags=0,
				 referencebus=None,
				 meanrepairtime=None,
				 degree=None,
				 stamp_dual= False,
				 voltage_type='LL'):
		self.ID = ID
		self.name = name
		self.phases = int(phases)
		self.voltage_type = voltage_type
		# Nominal voltage is coming as LN
		self.Vbase_LN = nominal_voltage * 1e-3 #if voltage_type != 'LL' else (
			# nominal_voltage * 1e-3) / np.sqrt(3)
		if voltage_type == 'LL':
			self.Vbase_LL = nominal_voltage * 1e-3
		else:
			if isTriplex:
				self.Vbase_LL = nominal_voltage * 1e-3
			else:
				self.Vbase_LL = nominal_voltage * np.sqrt(3) * 1e-3
		self.Vnom = float(nominal_voltage)
		self.Vset = float(voltage_setpoint)
		# Unknown - 0 | PQ - 1 | PV -2 | Swing - 3
		self.bustype = int(bustype)
		self.stamp_dual = stamp_dual
		# Seperate triplex nodes from nodes (Assign ABCN to nodes)
		self.isTriplex = isTriplex
		Nodes.ExistTriplex = self.isTriplex
		self.print_slack = False
		self.V_unb = 0  # Voltage Unbalance
		# Check the node connect type
		# 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
		self.node_connect_type = node_connect_type
		# Check the degree of the nodes
		# self.degree = int(degree)
		if not self.isTriplex:
			self.Va = (self.Vnom * self.Vset) + 0 * 1j
			self.Vb = (self.Vnom * self.Vset) * (self.cosd(-120) +
												 self.sind(-120) * 1j)
			self.Vc = (self.Vnom * self.Vset) * (self.cosd(120) +
												 self.sind(120) * 1j)
			self.Vn = complex(0, 0)

			if self.node_connect_type == 0:
				self.Vab = self.Va - self.Vb
				self.Vbc = self.Vb - self.Vc  # [volts]
				self.Vca = self.Vc - self.Va
				self.Vn = complex(0, 0)

			Nodes.Vr_mag.extend(4 * [np.absolute(self.Va)])
			Nodes.Vi_mag.extend(4 * [np.absolute(self.Va)])

			Nodes.V_mag.extend(8 * [np.absolute(self.Va)])

		else:
			# Define V1, V2 and Vn [In triplex node Va,
			# Vb, Vc corresponds to V1, V2 and VN]
			if self.phases & 0x01 == 1:
				self.V1 = self.Vnom + 0 * 1j
			elif self.phases & 0x02 == 2:
				self.V1 = self.Vnom * self.cosd(120) + self.Vnom * self.sind(
					120) * 1j
			elif self.phases & 0x04 == 4:
				self.V1 = self.Vnom * self.cosd(240) + self.Vnom * self.sind(
					240) * 1j
			# Initialize hot wire 2
			self.V2 = -self.V1  # Is the negative of V1
			# Initialize neutral wire
			self.Vn = complex(0, 0)
			Nodes.Vr_mag.extend(3 * [np.absolute(self.V1)])
			Nodes.Vi_mag.extend(3 * [np.absolute(self.V1)])
			Nodes.V_mag.extend(6 * [np.absolute(self.V1)])

		# Additional Parameters
		self.Ia = 0
		self.Ib = 0
		self.Ic = 0
		self.maxvoltage_error = maxvoltage_error
		self.busflags = int(
			busflags) if busflags else None  # Represents a PV bus in future
		self.referencebus = referencebus
		self.meanrepairtime = meanrepairtime
		self.obj_scalar = 1e-3


	def __str__(self):
		return 'Node name is :' + str(self.name)

	def assign_nodes(self, node_index_, obj, source_type, stamp_slack_bus_flag, stamp_infeas_neutral_source, casetype):
		# Assign indices initially for A, B, C, N
		self.nodeA_Vr = node_index_.__next__()
		self.nodeA_Vi = node_index_.__next__()
		self.nodeB_Vr = node_index_.__next__()
		self.nodeB_Vi = node_index_.__next__()
		self.nodeC_Vr = node_index_.__next__()
		self.nodeC_Vi = node_index_.__next__()
		self.nodeN_Vr = node_index_.__next__()
		self.nodeN_Vi = node_index_.__next__()

		self.Vr_nodes = [self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr]
		self.Vi_nodes = [self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi]

		# create dual voltage nodes if performing an optimization
		if self.stamp_dual:
			node_index_ = self.assign_dual_eq_var_nodes(node_index_, obj, stamp_slack_bus_flag)
			if self.bustype == 3 and stamp_slack_bus_flag == False:
				self.node_set = [
						self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
						self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
						self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
						self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
					] 
			else:
				if source_type == 'current':
					node_index_ = self.assign_infeas_current_L2_nodes(node_index_, obj, stamp_slack_bus_flag, stamp_infeas_neutral_source)
					
					if obj == 'L2':
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_if_r, self.nodeB_if_r, self.nodeC_if_r, self.nodeN_if_r,
								self.nodeA_if_i, self.nodeB_if_i, self.nodeC_if_i, self.nodeN_if_i
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_if_r, self.nodeB_if_r, self.nodeC_if_r,
								self.nodeA_if_i, self.nodeB_if_i, self.nodeC_if_i
							]		
					
					# if using the L1-norm, assign the inequality dual variables
					elif obj == 'L1':
						node_index_ = self.assign_ineq_dual_nodes(node_index_, source_type, stamp_infeas_neutral_source)
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_if_r_plus, self.nodeB_if_r_plus, self.nodeC_if_r_plus, self.nodeN_if_r_plus,
								self.nodeA_if_r_minus, self.nodeB_if_r_minus, self.nodeC_if_r_minus, self.nodeN_if_r_minus,
								self.nodeA_if_i_plus, self.nodeB_if_i_plus, self.nodeC_if_i_plus, self.nodeN_if_i_plus,
								self.nodeA_if_i_minus, self.nodeB_if_i_minus, self.nodeC_if_i_minus, self.nodeN_if_i_minus
							] 
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_if_r_plus, self.nodeB_if_r_plus, self.nodeC_if_r_plus, 
								self.nodeA_if_r_minus, self.nodeB_if_r_minus, self.nodeC_if_r_minus,
								self.nodeA_if_i_plus, self.nodeB_if_i_plus, self.nodeC_if_i_plus,
								self.nodeA_if_i_minus, self.nodeB_if_i_minus, self.nodeC_if_i_minus
							] 

				elif source_type == 'GB' or source_type == 'B':
					node_index_ = self.assign_infeas_GB_nodes(node_index_, source_type, obj, stamp_infeas_neutral_source)
					if source_type == 'GB':
						if stamp_infeas_neutral_source:
								self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_B, self.nodeB_B, self.nodeC_B, self.nodeN_B,
								self.nodeA_G, self.nodeB_G, self.nodeC_G, self.nodeN_G
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_B, self.nodeB_B, self.nodeC_B,
								self.nodeA_G, self.nodeB_G, self.nodeC_G
							]
					else:
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_B, self.nodeB_B, self.nodeC_B, self.nodeN_B
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_B, self.nodeB_B, self.nodeC_B
							]
				
				elif source_type == 'PQ' or source_type == 'Q':
					node_index_ = self.assign_infeas_PQ_nodes(node_index_, source_type, obj, stamp_infeas_neutral_source)
					if source_type == 'PQ' and obj == 'L2':
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_Q, self.nodeB_Q, self.nodeC_Q, self.nodeN_Q,
								self.nodeA_P, self.nodeB_P, self.nodeC_P, self.nodeN_P
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_Q, self.nodeB_Q, self.nodeC_Q, 
								self.nodeA_P, self.nodeB_P, self.nodeC_P
							]
					elif source_type == 'Q' and obj == 'L2':
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_Q, self.nodeB_Q, self.nodeC_Q, self.nodeN_Q
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_Q, self.nodeB_Q, self.nodeC_Q
							]
					elif source_type == 'PQ' and obj == 'L1':
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_p_plus, self.nodeB_p_plus, self.nodeC_p_plus, self.nodeN_p_plus,
								self.nodeA_p_minus, self.nodeB_p_minus, self.nodeC_p_minus, self.nodeN_p_minus,
								self.nodeA_q_plus, self.nodeB_q_plus, self.nodeC_q_plus, self.nodeN_q_plus,
								self.nodeA_q_minus, self.nodeB_q_minus, self.nodeC_q_minus, self.nodeN_q_minus
							] 
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_p_plus, self.nodeB_p_plus, self.nodeC_p_plus, 
								self.nodeA_p_minus, self.nodeB_p_minus, self.nodeC_p_minus, 
								self.nodeA_q_plus, self.nodeB_q_plus, self.nodeC_q_plus, 
								self.nodeA_q_minus, self.nodeB_q_minus, self.nodeC_q_minus
							] 

					elif source_type == 'Q' and obj == 'L1':
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_q_plus, self.nodeB_q_plus, self.nodeC_q_plus, self.nodeN_q_plus,
								self.nodeA_q_minus, self.nodeB_q_minus, self.nodeC_q_minus, self.noden_q_minus
							]
						else:
							self.node_set = [
								self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
								self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi,
								self.nodeA_dual_eq_var_r, self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_r, self.nodeB_dual_eq_var_i,
								self.nodeC_dual_eq_var_r, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_r, self.nodeN_dual_eq_var_i,
								self.nodeA_q_plus, self.nodeB_q_plus, self.nodeC_q_plus, 
								self.nodeA_q_minus, self.nodeB_q_minus, self.nodeC_q_minus
							]
		else:
			self.node_set = [
				self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
				self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi
			]

		# Check for PV nodes
		if self.bustype == 2:
			self.nodeA_Veq = node_index_.__next__()
			self.nodeB_Veq = node_index_.__next__()
			self.nodeC_Veq = node_index_.__next__()
		else:
			self.nodeA_Veq = None
			self.nodeB_Veq = None
			self.nodeC_Veq = None
		# Add extra equations for slack bus
		if self.bustype == 3:
			self.slack_nodeA_Vr = node_index_.__next__()
			self.slack_nodeA_Vi = node_index_.__next__()
			self.slack_nodeB_Vr = node_index_.__next__()
			self.slack_nodeB_Vi = node_index_.__next__()
			self.slack_nodeC_Vr = node_index_.__next__()
			self.slack_nodeC_Vi = node_index_.__next__()

			# Create dual slack nodes if performing optimization
			if self.stamp_dual:
				self.slack_nodeA_Lr = node_index_.__next__()
				self.slack_nodeA_Li = node_index_.__next__()
				self.slack_nodeB_Lr = node_index_.__next__()
				self.slack_nodeB_Li = node_index_.__next__()
				self.slack_nodeC_Lr = node_index_.__next__()
				self.slack_nodeC_Li = node_index_.__next__()
			
			# if casetype == 1:
			# 	self.nodeA_Vr_prime = node_index_.__next__()
			# 	self.nodeA_Vi_prime = node_index_.__next__()
			# 	self.nodeB_Vr_prime = node_index_.__next__()
			# 	self.nodeB_Vi_prime = node_index_.__next__()
			# 	self.nodeC_Vr_prime = node_index_.__next__()
			# 	self.nodeC_Vi_prime = node_index_.__next__()


		# Collect the indices of the voltage vector
		Nodes.voltage_index.extend([
			self.nodeA_Vr, self.nodeA_Vi, self.nodeB_Vr, self.nodeB_Vi,
			self.nodeC_Vr, self.nodeC_Vi, self.nodeN_Vr, self.nodeN_Vi
		])

		# Make sure the order is the same for the indices and normalizing magnitude
		Nodes.Vr_index.extend(
			[self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr])
		Nodes.Vi_index.extend(
			[self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi])
		
		if self.stamp_dual:
			Nodes.Lr_index.extend(
				[self.nodeA_dual_eq_var_r, self.nodeB_dual_eq_var_r, self.nodeC_dual_eq_var_r, self.nodeN_dual_eq_var_r])
			self.dual_eq_var_r_nodes = [self.nodeA_dual_eq_var_r, self.nodeB_dual_eq_var_r, self.nodeC_dual_eq_var_r, self.nodeN_dual_eq_var_r]
			Nodes.Li_index.extend(
				[self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_i, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_i])
			self.dual_eq_var_i_nodes = [self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_i, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_i]		

		return node_index_

	def assign_dual_eq_var_nodes(self, node_index_, obj, stamp_slack_bus_flag):
		# Function that assigns dual variables for each equality constraint 
		# and additionally assigns infeasibility current sources at each node but
		# the slack node
		self.nodeA_dual_eq_var_r = node_index_.__next__()
		self.nodeA_dual_eq_var_i = node_index_.__next__()
		self.nodeB_dual_eq_var_r = node_index_.__next__()
		self.nodeB_dual_eq_var_i = node_index_.__next__()
		self.nodeC_dual_eq_var_r = node_index_.__next__()
		self.nodeC_dual_eq_var_i = node_index_.__next__()
		self.nodeN_dual_eq_var_r = node_index_.__next__()
		self.nodeN_dual_eq_var_i = node_index_.__next__()

		
		Nodes.L_index.extend([self.nodeA_dual_eq_var_r, self.nodeB_dual_eq_var_r, self.nodeC_dual_eq_var_r, self.nodeN_dual_eq_var_r,
							  self.nodeA_dual_eq_var_i, self.nodeB_dual_eq_var_i, self.nodeC_dual_eq_var_i, self.nodeN_dual_eq_var_i])
		
		return node_index_

	def assign_infeas_current_L2_nodes(self, node_index_, obj, stamp_slack_bus_flag, stamp_infeas_neutral_source):
		if obj != 'L1' and (self.bustype!= 3 or stamp_slack_bus_flag):
			self.nodeA_if_r = node_index_.__next__()
			self.nodeA_if_i = node_index_.__next__()
			self.nodeB_if_r = node_index_.__next__()
			self.nodeB_if_i = node_index_.__next__()
			self.nodeC_if_r = node_index_.__next__()
			self.nodeC_if_i = node_index_.__next__()

			Nodes.if_r_index.extend([self.nodeA_if_r, self.nodeB_if_r,
									 self.nodeC_if_r])
			
			Nodes.if_i_index.extend([self.nodeA_if_i, self.nodeB_if_i,
									 self.nodeC_if_i])   
			
			Nodes.if_index.extend([self.nodeA_if_r, self.nodeA_if_i,
								   self.nodeB_if_r, self.nodeB_if_i, 
								   self.nodeC_if_r, self.nodeC_if_i])
			
			if stamp_infeas_neutral_source:
				self.nodeN_if_r = node_index_.__next__()
				self.nodeN_if_i = node_index_.__next__() 

				Nodes.if_r_index.extend([self.nodeN_if_r])
			
				Nodes.if_i_index.extend([self.nodeN_if_i])   
			
				Nodes.if_index.extend([self.nodeN_if_r, self.nodeN_if_i]) 
			
				self.if_r_nodes = [self.nodeA_if_r, self.nodeB_if_r,
										self.nodeC_if_r, self.nodeN_if_r]
				self.if_i_nodes = [self.nodeA_if_i, self.nodeB_if_i,
										self.nodeC_if_i, self.nodeN_if_i]
			else:
				self.if_r_nodes = [self.nodeA_if_r, self.nodeB_if_r,
										self.nodeC_if_r]
				self.if_i_nodes = [self.nodeA_if_i, self.nodeB_if_i,
										self.nodeC_if_i]
		
		return node_index_
	
	def assign_infeas_GB_nodes(self, node_index_, source_type, obj, stamp_infeas_neutral_source):
		if obj == 'L2':
			self.nodeA_B = node_index_.__next__()
			self.nodeB_B = node_index_.__next__()
			self.nodeC_B = node_index_.__next__()
			Nodes.infeas_imag_var_index.extend([self.nodeA_B, self.nodeB_B,
									self.nodeC_B])
			if source_type == 'GB':
				self.nodeA_G = node_index_.__next__()
				self.nodeB_G = node_index_.__next__()
				self.nodeC_G = node_index_.__next__()
				Nodes.infeas_real_var_index.extend([self.nodeA_G, self.nodeB_G,
									self.nodeC_G]) 
			
			if stamp_infeas_neutral_source:
				self.nodeN_B = node_index_.__next__()
				Nodes.infeas_imag_var_index.extend([self.nodeN_B])
				self.B_nodes = [self.nodeA_B, self.nodeB_B,
									self.nodeC_B, self.nodeN_B]
				if source_type == 'GB':
					self.nodeN_G = node_index_.__next__() 
					Nodes.infeas_real_var_index.extend([self.nodeN_G])
					self.G_nodes = [self.nodeA_G, self.nodeB_G,
									self.nodeC_G, self.nodeN_G] 
			else:
				self.B_nodes = [self.nodeA_B, self.nodeB_B,
									self.nodeC_B]
				if source_type == 'GB':
					self.G_nodes = [self.nodeA_G, self.nodeB_G,
									self.nodeC_G] 				

		return node_index_
	
	def assign_infeas_PQ_nodes(self, node_index_, source_type, obj, stamp_infeas_neutral_source):
		if self.isTriplex:
			if obj == 'L2' and (source_type == 'PQ' or source_type == 'Q'):
				self.node1_Q = node_index_.__next__()
				self.node2_Q = node_index_.__next__()
				if stamp_infeas_neutral_source:
					self.nodeN_Q = node_index_.__next__()
					self.Q_nodes = [self.node1_Q, self.node2_Q, self.nodeN_Q]
				else:
					self.Q_nodes = [self.node1_Q, self.node2_Q]
				Nodes.infeas_imag_var_index += self.Q_nodes
				
				if source_type == 'PQ':
					self.node1_P = node_index_.__next__()
					self.node2_P = node_index_.__next__()
					if stamp_infeas_neutral_source:
						self.nodeN_P = node_index_.__next__() 
						self.P_nodes = [self.node1_P, self.node2_P, self.nodeN_P]
					else:
						self.P_nodes = [self.node1_P, self.node2_P]
					Nodes.infeas_real_var_index += self.P_nodes
				
			if obj == 'L1' and (source_type == 'PQ' or source_type == 'Q'):
				self.node1_q_plus = node_index_.__next__()
				self.node1_q_minus = node_index_.__next__()
				self.node2_q_plus = node_index_.__next__()
				self.node2_q_minus = node_index_.__next__()

				self.node1_dual_ineq_var_q_plus = node_index_.__next__()
				self.node1_dual_ineq_var_q_minus = node_index_.__next__()
				self.node2_dual_ineq_var_q_plus = node_index_.__next__()
				self.node2_dual_ineq_var_q_minus = node_index_.__next__()

				if stamp_infeas_neutral_source:
					self.nodeN_q_plus = node_index_.__next__()
					self.nodeN_q_minus = node_index_.__next__()

					self.nodeN_dual_ineq_var_q_plus = node_index_.__next__()
					self.nodeN_dual_ineq_var_q_minus = node_index_.__next__()

					self.Q_plus_nodes = [self.node1_q_plus, self.node2_q_plus, self.nodeN_q_plus]
					self.Q_minus_nodes = [self.node1_q_minus, self.node2_q_minus, self.nodeN_q_minus]
					self.dual_ineq_i_plus_nodes = [self.node1_dual_ineq_var_q_plus, self.node2_dual_ineq_var_q_plus,
							self.nodeN_dual_ineq_var_q_plus]
					self.dual_ineq_i_minus_nodes = [self.node1_dual_ineq_var_q_minus, self.node2_dual_ineq_var_q_minus,
							self.nodeN_dual_ineq_var_q_minus]
				else:
					self.Q_plus_nodes = [self.node1_q_plus, self.node2_q_plus]
					self.Q_minus_nodes = [self.node1_q_minus, self.node2_q_minus]
					self.dual_ineq_i_plus_nodes = [self.node1_dual_ineq_var_q_plus, self.node2_dual_ineq_var_q_plus]
					self.dual_ineq_i_minus_nodes = [self.node1_dual_ineq_var_q_minus, self.node2_dual_ineq_var_q_minus]	
							
				if source_type == 'PQ':
					self.node1_p_plus = node_index_.__next__()
					self.node1_p_minus = node_index_.__next__()
					self.node2_p_plus = node_index_.__next__()
					self.node2_p_minus = node_index_.__next__()

					self.node1_dual_ineq_var_p_plus = node_index_.__next__()
					self.node1_dual_ineq_var_p_minus = node_index_.__next__()
					self.node2_dual_ineq_var_p_plus = node_index_.__next__()
					self.node2_dual_ineq_var_p_minus = node_index_.__next__()

					if stamp_infeas_neutral_source:
						self.nodeN_p_plus = node_index_.__next__()
						self.nodeN_p_minus = node_index_.__next__()

						self.nodeN_dual_ineq_var_p_plus = node_index_.__next__()
						self.nodeN_dual_ineq_var_p_minus = node_index_.__next__()
						
						self.P_plus_nodes = [self.node1_p_plus, self.node2_p_plus, self.nodeN_p_plus]
						self.P_minus_nodes = [self.node1_p_minus, self.node2_p_minus, self.nodeN_p_minus]
						self.dual_ineq_r_plus_nodes = [self.node1_dual_ineq_var_p_plus, self.node2_dual_ineq_var_p_plus,
								self.nodeN_dual_ineq_var_p_plus]
						self.dual_ineq_r_minus_nodes = [self.node1_dual_ineq_var_q_minus, self.node2_dual_ineq_var_p_minus,
								self.nodeN_dual_ineq_var_p_minus]
					else:
						self.P_plus_nodes = [self.node1_p_plus, self.node2_p_plus]
						self.P_minus_nodes = [self.node1_p_minus, self.node2_p_minus]
						self.dual_ineq_r_plus_nodes = [self.node1_dual_ineq_var_p_plus, self.node2_dual_ineq_var_p_plus]
						self.dual_ineq_r_minus_nodes = [self.node1_dual_ineq_var_p_minus, self.node2_dual_ineq_var_p_minus]
					
				Nodes.infeas_imag_var_index += self.Q_plus_nodes + self.Q_minus_nodes
				Nodes.mu_index += self.dual_ineq_i_plus_nodes + self.dual_ineq_i_minus_nodes
				Nodes.mu_i_index += self.dual_ineq_i_plus_nodes + self.dual_ineq_i_minus_nodes

				if source_type == 'PQ':
					Nodes.infeas_real_var_index += self.P_plus_nodes + self.P_minus_nodes
					Nodes.mu_index_upper += self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes
					Nodes.mu_r_index_upper += self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes
			
		else:
			if obj == 'L2':
				self.nodeA_Q = node_index_.__next__()
				self.nodeB_Q = node_index_.__next__()
				self.nodeC_Q = node_index_.__next__()
				Nodes.infeas_imag_var_index.extend([self.nodeA_Q, self.nodeB_Q,
										self.nodeC_Q])
				
				if source_type == 'PQ':
					self.nodeA_P = node_index_.__next__()
					self.nodeB_P = node_index_.__next__()
					self.nodeC_P = node_index_.__next__()
					Nodes.infeas_real_var_index.extend([self.nodeA_P, self.nodeB_P,
										self.nodeC_P]) 

				if stamp_infeas_neutral_source:
					self.nodeN_Q = node_index_.__next__()
					Nodes.infeas_imag_var_index.extend([self.nodeN_Q])
					self.Q_nodes = [self.nodeA_Q, self.nodeB_Q,
										self.nodeC_Q, self.nodeN_Q]
					if source_type == 'PQ':
						self.nodeN_P = node_index_.__next__() 
						Nodes.infeas_real_var_index.extend([self.nodeN_P]) 
						self.P_nodes = [self.nodeA_P, self.nodeB_P,
											self.nodeC_P, self.nodeN_P]
				else:
					self.Q_nodes = [self.nodeA_Q, self.nodeB_Q,
										self.nodeC_Q]
					if source_type == 'PQ': 
						self.P_nodes = [self.nodeA_P, self.nodeB_P,
											self.nodeC_P]

			if obj == 'L1' and (source_type == 'PQ' or source_type == 'Q'):
				if source_type == 'PQ':
					self.nodeA_p_plus = node_index_.__next__()
					self.nodeA_p_minus = node_index_.__next__()
					self.nodeB_p_plus = node_index_.__next__()
					self.nodeB_p_minus = node_index_.__next__()
					self.nodeC_p_plus = node_index_.__next__()
					self.nodeC_p_minus = node_index_.__next__()

					self.nodeA_dual_ineq_var_p_plus = node_index_.__next__()
					self.nodeA_dual_ineq_var_p_minus = node_index_.__next__()
					self.nodeB_dual_ineq_var_p_plus = node_index_.__next__()
					self.nodeB_dual_ineq_var_p_minus = node_index_.__next__()
					self.nodeC_dual_ineq_var_p_plus = node_index_.__next__()
					self.nodeC_dual_ineq_var_p_minus = node_index_.__next__()

					# (debug) upper bounds
					self.nodeA_dual_ineq_var_p_plus_upper = node_index_.__next__()
					self.nodeA_dual_ineq_var_p_minus_upper = node_index_.__next__()
					self.nodeB_dual_ineq_var_p_plus_upper = node_index_.__next__()
					self.nodeB_dual_ineq_var_p_minus_upper = node_index_.__next__()
					self.nodeC_dual_ineq_var_p_plus_upper = node_index_.__next__()
					self.nodeC_dual_ineq_var_p_minus_upper = node_index_.__next__()

				self.nodeA_q_plus = node_index_.__next__()
				self.nodeA_q_minus = node_index_.__next__()
				self.nodeB_q_plus = node_index_.__next__()
				self.nodeB_q_minus = node_index_.__next__()
				self.nodeC_q_plus = node_index_.__next__()
				self.nodeC_q_minus = node_index_.__next__()

				self.nodeA_dual_ineq_var_q_plus = node_index_.__next__()
				self.nodeA_dual_ineq_var_q_minus = node_index_.__next__()
				self.nodeB_dual_ineq_var_q_plus = node_index_.__next__()
				self.nodeB_dual_ineq_var_q_minus = node_index_.__next__()
				self.nodeC_dual_ineq_var_q_plus = node_index_.__next__()
				self.nodeC_dual_ineq_var_q_minus = node_index_.__next__()

				if stamp_infeas_neutral_source:
					if source_type == 'PQ':
						self.nodeN_p_plus = node_index_.__next__()
						self.nodeN_p_minus = node_index_.__next__()

						self.nodeN_dual_ineq_var_p_plus = node_index_.__next__()
						self.nodeN_dual_ineq_var_p_minus = node_index_.__next__()

						# (debug) upper bounds
						self.nodeN_dual_ineq_var_p_plus_upper = node_index_.__next__()
						self.nodeN_dual_ineq_var_p_minus_upper = node_index_.__next__()
						
						self.P_plus_nodes = [self.nodeA_p_plus, self.nodeB_p_plus,
							self.nodeC_p_plus, self.nodeN_p_plus]
						self.P_minus_nodes = [self.nodeA_p_minus, self.nodeB_p_minus,
								self.nodeC_p_minus, self.nodeN_p_minus]
						self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_var_p_plus, self.nodeB_dual_ineq_var_p_plus,
								self.nodeC_dual_ineq_var_p_plus,self.nodeN_dual_ineq_var_p_plus]
						self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_var_q_minus, self.nodeB_dual_ineq_var_p_minus,
								self.nodeC_dual_ineq_var_p_minus,self.nodeN_dual_ineq_var_p_minus]
						
						# (debug) upper nodes
						self.dual_ineq_r_plus_nodes_upper = [self.nodeA_dual_ineq_var_p_plus_upper, 
										   self.nodeB_dual_ineq_var_p_plus_upper,self.nodeC_dual_ineq_var_p_plus_upper,
											self.nodeN_dual_ineq_var_p_plus_upper]
						self.dual_ineq_r_minus_nodes_upper = [self.nodeA_dual_ineq_var_p_minus_upper, 
										   self.nodeB_dual_ineq_var_p_minus_upper,self.nodeC_dual_ineq_var_p_minus_upper,
											self.nodeN_dual_ineq_var_p_minus_upper]

					self.nodeN_q_plus = node_index_.__next__()
					self.nodeN_q_minus = node_index_.__next__()

					self.nodeN_dual_ineq_var_q_plus = node_index_.__next__()
					self.nodeN_dual_ineq_var_q_minus = node_index_.__next__()

					self.Q_plus_nodes = [self.nodeA_q_plus, self.nodeB_q_plus,
							self.nodeC_q_plus, self.nodeN_q_plus]
					self.Q_minus_nodes = [self.nodeA_q_minus, self.nodeB_q_minus,
							self.nodeC_q_minus, self.nodeN_q_minus]
					self.dual_ineq_i_plus_nodes = [self.nodeA_dual_ineq_var_q_plus, self.nodeB_dual_ineq_var_q_plus,
							self.nodeC_dual_ineq_var_q_plus,self.nodeN_dual_ineq_var_q_plus]
					self.dual_ineq_i_minus_nodes = [self.nodeA_dual_ineq_var_q_minus, self.nodeB_dual_ineq_var_q_minus,
							self.nodeC_dual_ineq_var_q_minus,self.nodeN_dual_ineq_var_q_minus]
				else:
					if source_type == 'PQ':
						self.P_plus_nodes = [self.nodeA_p_plus, self.nodeB_p_plus,
							self.nodeC_p_plus]
						self.P_minus_nodes = [self.nodeA_p_minus, self.nodeB_p_minus,
								self.nodeC_p_minus]
						self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_var_p_plus, self.nodeB_dual_ineq_var_p_plus,
								self.nodeC_dual_ineq_var_p_plus]
						self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_var_p_minus, self.nodeB_dual_ineq_var_p_minus,
								self.nodeC_dual_ineq_var_p_minus]
						
						# (debug) upper nodes
						self.dual_ineq_r_plus_nodes_upper = [self.nodeA_dual_ineq_var_p_plus_upper, 
										   self.nodeB_dual_ineq_var_p_plus_upper,self.nodeC_dual_ineq_var_p_plus_upper]
						self.dual_ineq_r_minus_nodes_upper = [self.nodeA_dual_ineq_var_p_minus_upper, 
										   self.nodeB_dual_ineq_var_p_minus_upper,self.nodeC_dual_ineq_var_p_minus_upper]

					self.Q_plus_nodes = [self.nodeA_q_plus, self.nodeB_q_plus,
							self.nodeC_q_plus]
					self.Q_minus_nodes = [self.nodeA_q_minus, self.nodeB_q_minus,
							self.nodeC_q_minus]
					self.dual_ineq_i_plus_nodes = [self.nodeA_dual_ineq_var_q_plus, self.nodeB_dual_ineq_var_q_plus,
							self.nodeC_dual_ineq_var_q_plus]
					self.dual_ineq_i_minus_nodes = [self.nodeA_dual_ineq_var_q_minus, self.nodeB_dual_ineq_var_q_minus,
							self.nodeC_dual_ineq_var_q_minus]
				
				Nodes.infeas_imag_var_index += self.Q_plus_nodes + self.Q_minus_nodes
				Nodes.mu_index += self.dual_ineq_i_plus_nodes + self.dual_ineq_i_minus_nodes
				Nodes.mu_i_index += self.dual_ineq_i_plus_nodes + self.dual_ineq_i_minus_nodes

				if source_type == 'PQ':
					Nodes.infeas_real_var_index += self.P_plus_nodes + self.P_minus_nodes
					Nodes.mu_index += self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes
					Nodes.mu_r_index += self.dual_ineq_r_plus_nodes + self.dual_ineq_r_minus_nodes

					# debug nodes
					Nodes.mu_index_upper += self.dual_ineq_r_plus_nodes_upper + self.dual_ineq_r_minus_nodes_upper
					Nodes.mu_r_index_upper += self.dual_ineq_r_plus_nodes_upper + self.dual_ineq_r_minus_nodes_upper

		return node_index_
	
	def assign_ineq_dual_nodes(self, node_index_, source_type, stamp_infeas_neutral_source):
		# Function that assigns dual inquality variables for each inequality constraint
		# at each node, in addition to the slack variables for the L1 norm formulation  
		if source_type == 'current':
			self.nodeA_if_r_plus = node_index_.__next__()
			self.nodeA_if_r_minus = node_index_.__next__()
			self.nodeA_if_i_plus = node_index_.__next__()
			self.nodeA_if_i_minus = node_index_.__next__()
			
			self.nodeB_if_r_plus = node_index_.__next__()
			self.nodeB_if_r_minus = node_index_.__next__()               
			self.nodeB_if_i_plus = node_index_.__next__()
			self.nodeB_if_i_minus = node_index_.__next__() 
			
			self.nodeC_if_r_plus = node_index_.__next__()
			self.nodeC_if_r_minus = node_index_.__next__()               
			self.nodeC_if_i_plus = node_index_.__next__()
			self.nodeC_if_i_minus = node_index_.__next__()  

			if stamp_infeas_neutral_source:
				self.nodeN_if_r_plus = node_index_.__next__()
				self.nodeN_if_r_minus = node_index_.__next__()               
				self.nodeN_if_i_plus = node_index_.__next__()
				self.nodeN_if_i_minus = node_index_.__next__()

				self.if_r_plus_nodes = [self.nodeA_if_r_plus, self.nodeB_if_r_plus,
										self.nodeC_if_r_plus, self.nodeN_if_r_plus]
				self.if_r_minus_nodes = [self.nodeA_if_r_minus, self.nodeB_if_r_minus,
											self.nodeC_if_r_minus, self.nodeN_if_r_minus]
				self.if_i_plus_nodes = [self.nodeA_if_i_plus, self.nodeB_if_i_plus,
											self.nodeC_if_i_plus, self.nodeN_if_i_plus]
				self.if_i_minus_nodes = [self.nodeA_if_i_minus, self.nodeB_if_i_minus,
											self.nodeC_if_i_minus, self.nodeN_if_i_minus] 
			else:
				self.if_r_plus_nodes = [self.nodeA_if_r_plus, self.nodeB_if_r_plus,
										self.nodeC_if_r_plus]
				self.if_r_minus_nodes = [self.nodeA_if_r_minus, self.nodeB_if_r_minus,
											self.nodeC_if_r_minus]
				self.if_i_plus_nodes = [self.nodeA_if_i_plus, self.nodeB_if_i_plus,
											self.nodeC_if_i_plus]
				self.if_i_minus_nodes = [self.nodeA_if_i_minus, self.nodeB_if_i_minus,
											self.nodeC_if_i_minus]
			
			Nodes.if_r_plus_index += self.if_r_plus_nodes
			Nodes.if_r_minus_index += self.if_r_minus_nodes
				
			Nodes.if_i_plus_index += self.if_i_plus_nodes  
			Nodes.if_i_minus_index += self.if_i_minus_nodes 
			
			self.nodeA_dual_ineq_r_minus = node_index_.__next__()
			self.nodeA_dual_ineq_i_minus = node_index_.__next__()
			self.nodeB_dual_ineq_r_minus = node_index_.__next__()
			self.nodeB_dual_ineq_i_minus = node_index_.__next__()
			self.nodeC_dual_ineq_r_minus = node_index_.__next__()
			self.nodeC_dual_ineq_i_minus = node_index_.__next__()
			
			self.nodeA_dual_ineq_r_plus = node_index_.__next__()
			self.nodeA_dual_ineq_i_plus = node_index_.__next__()
			self.nodeB_dual_ineq_r_plus = node_index_.__next__()
			self.nodeB_dual_ineq_i_plus = node_index_.__next__()
			self.nodeC_dual_ineq_r_plus = node_index_.__next__()
			self.nodeC_dual_ineq_i_plus = node_index_.__next__()

			if stamp_infeas_neutral_source:
				self.nodeN_dual_ineq_r_minus = node_index_.__next__()
				self.nodeN_dual_ineq_i_minus = node_index_.__next__()
			
				self.nodeN_dual_ineq_r_plus = node_index_.__next__()
				self.nodeN_dual_ineq_i_plus = node_index_.__next__()
				
				self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_r_minus, self.nodeB_dual_ineq_r_minus, 
									self.nodeC_dual_ineq_r_minus, self.nodeN_dual_ineq_r_minus]
			
				self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_r_plus, self.nodeB_dual_ineq_r_plus,
										self.nodeC_dual_ineq_r_plus, self.nodeN_dual_ineq_r_plus]
				
				self.dual_ineq_i_minus_nodes = [self.nodeA_dual_ineq_i_minus, self.nodeB_dual_ineq_i_minus, 
									self.nodeC_dual_ineq_i_minus, self.nodeN_dual_ineq_i_minus]
			
				self.dual_ineq_i_plus_nodes = [self.nodeA_dual_ineq_i_plus, self.nodeB_dual_ineq_i_plus,
										self.nodeC_dual_ineq_i_plus, self.nodeN_dual_ineq_i_plus]
			else:
				self.dual_ineq_r_minus_nodes = [self.nodeA_dual_ineq_r_minus, self.nodeB_dual_ineq_r_minus, 
									self.nodeC_dual_ineq_r_minus]
			
				self.dual_ineq_r_plus_nodes = [self.nodeA_dual_ineq_r_plus, self.nodeB_dual_ineq_r_plus,
										self.nodeC_dual_ineq_r_plus]
				
				self.dual_ineq_i_minus_nodes = [self.nodeA_dual_ineq_i_minus, self.nodeB_dual_ineq_i_minus, 
									self.nodeC_dual_ineq_i_minus]
			
				self.dual_ineq_i_plus_nodes = [self.nodeA_dual_ineq_i_plus, self.nodeB_dual_ineq_i_plus,
										self.nodeC_dual_ineq_i_plus]	
		
			Nodes.mu_index += self.dual_ineq_r_minus_nodes + self.dual_ineq_r_plus_nodes + self.dual_ineq_i_minus_nodes + self.dual_ineq_i_plus_nodes
			
			Nodes.if_index += self.if_r_plus_nodes + self.if_r_minus_nodes + self.if_i_plus_nodes + self.if_i_minus_nodes
			
			Nodes.mu_r_index += self.dual_ineq_r_minus_nodes + self.dual_ineq_r_plus_nodes
								
			Nodes.mu_i_index += self.dual_ineq_i_minus_nodes + self.dual_ineq_i_plus_nodes
		
		return node_index_
	
	def assign_triplex_nodes(self, node_index_, obj, source_type,  stamp_tplx_infeas_sources, stamp_infeas_neutral_source):
		self.node1_Vr = node_index_.__next__()
		self.node1_Vi = node_index_.__next__()
		self.node2_Vr = node_index_.__next__()
		self.node2_Vi = node_index_.__next__()
		self.nodeN_Vr = node_index_.__next__()
		self.nodeN_Vi = node_index_.__next__()

		Nodes.voltage_index.extend([
				self.node1_Vr, self.node1_Vi, self.node2_Vr, self.node2_Vi,
				self.nodeN_Vr, self.nodeN_Vi
			])

		if self.stamp_dual:
			self.node1_dual_eq_var_r = node_index_.__next__()
			self.node1_dual_eq_var_i = node_index_.__next__()
			self.node2_dual_eq_var_r = node_index_.__next__()
			self.node2_dual_eq_var_i = node_index_.__next__()
			self.nodeN_dual_eq_var_r = node_index_.__next__()
			self.nodeN_dual_eq_var_i = node_index_.__next__()
			
			Nodes.L_index.extend([self.node1_dual_eq_var_r, self.node2_dual_eq_var_r, self.nodeN_dual_eq_var_r,
								 self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i])
			
			Nodes.Lr_index.extend(
					[self.node1_dual_eq_var_r, self.node2_dual_eq_var_r, self.nodeN_dual_eq_var_r])
			Nodes.Li_index.extend(
					[self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i])
			
			self.dual_eq_var_r_nodes = [self.node1_dual_eq_var_r, self.node2_dual_eq_var_r, self.nodeN_dual_eq_var_r]
			self.dual_eq_var_i_nodes = [self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i]
			
			if stamp_tplx_infeas_sources:
				if obj == 'L2' and source_type == 'current':
					self.node1_if_r = node_index_.__next__()
					self.node1_if_i = node_index_.__next__()
					self.node2_if_r = node_index_.__next__()
					self.node2_if_i = node_index_.__next__()

					if stamp_infeas_neutral_source:
						self.nodeN_if_r = node_index_.__next__()
						self.nodeN_if_i = node_index_.__next__()  
						
						self.if_r_nodes = [self.node1_if_r, self.node2_if_r,
											self.nodeN_if_r]
						self.if_i_nodes = [self.node1_if_i, self.node2_if_i,
												self.nodeN_if_i]
											
						self.node_set = [
							self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
							self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
							self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
							self.node1_if_r, self.node2_if_r, self.nodeN_if_r, 
							self.node1_if_i, self.node2_if_i, self.nodeN_if_i    
						] 
					else:
						self.if_r_nodes = [self.node1_if_r, self.node2_if_r]
						self.if_i_nodes = [self.node1_if_i, self.node2_if_i]
						
						self.node_set = [
							self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
							self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
							self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
							self.node1_if_r, self.node2_if_r, self.node1_if_i, self.node2_if_i    
						] 
						
					Nodes.if_r_index += self.if_r_nodes
					
					Nodes.if_r_index += self.if_i_nodes
				
				if obj == 'L1' and source_type == 'current':
					node_index_ = self.assign_ineq_dual_nodes_tplx(node_index_, source_type, stamp_infeas_neutral_source)

					self.node_set = [
						self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
						self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
						self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
						self.node1_if_r_plus, self.node2_if_r_plus, self.nodeN_if_r_plus, 
						self.node1_if_r_minus, self.node2_if_r_minus, self.nodeN_if_r_minus, 
						self.node1_if_i_plus, self.node2_if_i_plus, self.nodeN_if_i_plus,
						self.node1_if_i_minus, self.node2_if_i_minus, self.nodeN_if_i_minus
				] 
  				
				if obj == 'L2' and (source_type == 'GB' or source_type == 'B'):
					self.node1_B = node_index_.__next__()
					self.node2_B = node_index_.__next__()
					if stamp_infeas_neutral_source:
						self.nodeN_B = node_index_.__next__()
						self.B_nodes = [self.node1_B, self.node2_B, self.nodeN_B]
					else:
						self.B_nodes = [self.node1_B, self.node2_B]
					Nodes.infeas_imag_var_index += self.B_nodes

					if source_type == 'GB':
						self.node1_G = node_index_.__next__()
						self.node2_G = node_index_.__next__()
						if stamp_infeas_neutral_source:
							self.nodeN_G = node_index_.__next__() 
							self.G_nodes = [self.node1_G, self.node2_G, self.nodeN_G] 
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_B, self.node2_B, self.nodeN_B,
								self.node1_G, self.node2_G, self.nodeN_G   
							]  
						else:
							self.G_nodes = [self.node1_G, self.node2_G] 
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_B, self.node2_B,
								self.node1_G, self.node2_G 
							] 
						Nodes.infeas_real_var_index += self.G_nodes  
					
					else: 
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_B, self.node2_B, self.nodeN_B   
							]
						else:
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_B, self.node2_B   
							]

				if obj == 'L2' and (source_type == 'PQ' or source_type == 'Q'):
					node_index_ = self.assign_infeas_PQ_nodes(node_index_, source_type, obj, stamp_infeas_neutral_source)
					
					if source_type == 'PQ':
						if stamp_infeas_neutral_source:
							self.node_set = [
									self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
									self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
									self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
									self.node1_Q, self.node2_Q, self.nodeN_Q,
									self.node1_P, self.node2_P, self.nodeN_P
								]
						else:
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_Q, self.node2_Q, self.node1_P, self.node2_P
							]
					else:
						if stamp_infeas_neutral_source:
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_Q, self.node2_Q, self.nodeN_Q
							] 
						else:
							self.node_set = [
								self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
								self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
								self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
								self.node1_Q, self.node2_Q
							] 

				if obj == 'L1' and (source_type == 'PQ' or source_type == 'Q'):
					node_index_ = self.assign_infeas_PQ_nodes(node_index_, source_type, obj, stamp_infeas_neutral_source)
					if source_type == 'PQ':
						if stamp_infeas_neutral_source:
							self.node_set = [
									self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
									self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
									self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
									self.node1_p_plus, self.node2_p_plus, self.nodeN_p_plus,
									self.node1_p_minus, self.node2_p_minus, self.nodeN_p_minus,
									self.node1_q_plus, self.node2_q_plus, self.nodeN_q_plus,
									self.node1_q_minus, self.node2_q_minus, self.nodeN_q_minus
								]
						else:
							self.node_set = [
									self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
									self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
									self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
									self.node1_p_plus, self.node2_p_plus, 
									self.node1_p_minus, self.node2_p_minus,
									self.node1_q_plus, self.node2_q_plus, 
									self.node1_q_minus, self.node2_q_minus
								]
					else:
						if stamp_infeas_neutral_source:
							self.node_set = [
									self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
									self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
									self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
									self.node1_q_plus, self.node2_q_plus, self.nodeN_q_plus,
									self.node1_q_minus, self.node2_q_minus, self.nodeN_q_minus
								]
						else:
							self.node_set = [
									self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
									self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
									self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i,
									self.node1_q_plus, self.node2_q_plus, 
									self.node1_q_minus, self.node2_q_minus
								] 

			else:
				self.node_set = [
						self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
						self.node2_Vi, self.nodeN_Vi, self.node1_dual_eq_var_r, self.node2_dual_eq_var_r,
						self.nodeN_dual_eq_var_r, self.node1_dual_eq_var_i, self.node2_dual_eq_var_i, self.nodeN_dual_eq_var_i
				]                               
		else:
			self.node_set = [
				self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
				self.node2_Vi, self.nodeN_Vi
			]
			
		Nodes.Vr_index.extend([self.node1_Vr, self.node2_Vr, self.nodeN_Vr])
		Nodes.Vi_index.extend([self.node1_Vi, self.node2_Vi, self.nodeN_Vi])
			
		return node_index_
			
	def assign_ineq_dual_nodes_tplx(self, node_index_, source_type, stamp_infeas_neutral_source):
		# Function that assigns dual inquality variables for each inequality constraint
		# at each triplex node, in addition to the slack variables for the L1 norm formulation	
		if source_type == 'current':
			self.node1_if_r_plus = node_index_.__next__()
			self.node1_if_r_minus = node_index_.__next__()
			self.node1_if_i_plus = node_index_.__next__()
			self.node1_if_i_minus = node_index_.__next__()
			
			self.node2_if_r_plus = node_index_.__next__()
			self.node2_if_r_minus = node_index_.__next__()               
			self.node2_if_i_plus = node_index_.__next__()
			self.node2_if_i_minus = node_index_.__next__() 

			if stamp_infeas_neutral_source:
				self.nodeN_if_r_plus = node_index_.__next__()
				self.nodeN_if_r_minus = node_index_.__next__()               
				self.nodeN_if_i_plus = node_index_.__next__()
				self.nodeN_if_i_minus = node_index_.__next__() 
				self.if_r_plus_nodes = [self.node1_if_r_plus, self.node2_if_r_plus]
				self.if_r_minus_nodes = [self.node1_if_r_minus, self.node2_if_r_minus]
				self.if_i_plus_nodes = [self.node1_if_i_plus, self.node2_if_i_plus]
				self.if_i_minus_nodes = [self.node1_if_i_minus, self.node2_if_i_minus]
			else:
				self.if_r_plus_nodes = [self.node1_if_r_plus, self.node2_if_r_plus]
				self.if_r_minus_nodes = [self.node1_if_r_minus, self.node2_if_r_minus]
				self.if_i_plus_nodes = [self.node1_if_i_plus, self.node2_if_i_plus]
				self.if_i_minus_nodes = [self.node1_if_i_minus, self.node2_if_i_minus]		
			
			Nodes.if_r_plus_index += self.if_r_plus_nodes
			Nodes.if_r_minus_index += self.if_r_minus_nodes
				
			Nodes.if_i_plus_index += self.if_i_plus_nodes
			Nodes.if_i_minus_index += self.if_i_minus_nodes 
		
			self.node1_dual_ineq_r_minus = node_index_.__next__()
			self.node1_dual_ineq_i_minus = node_index_.__next__()
			self.node2_dual_ineq_r_minus = node_index_.__next__()
			self.node2_dual_ineq_i_minus = node_index_.__next__()
			
			self.node1_dual_ineq_r_plus = node_index_.__next__()
			self.node1_dual_ineq_i_plus = node_index_.__next__()
			self.node2_dual_ineq_r_plus = node_index_.__next__()
			self.node2_dual_ineq_i_plus = node_index_.__next__()

			if stamp_infeas_neutral_source:
				self.nodeN_dual_ineq_r_minus = node_index_.__next__()
				self.nodeN_dual_ineq_i_minus = node_index_.__next__()
				self.nodeN_dual_ineq_r_plus = node_index_.__next__()
				self.nodeN_dual_ineq_i_plus = node_index_.__next__()
				self.dual_ineq_r_minus_nodes = [self.node1_dual_ineq_r_minus, self.node2_dual_ineq_r_minus, 
				   				self.nodeN_dual_ineq_r_minus]
				self.dual_ineq_i_minus_nodes = [self.node1_dual_ineq_i_minus, self.node2_dual_ineq_i_minus, 
										self.nodeN_dual_ineq_i_minus]
				self.dual_ineq_r_plus_nodes = [self.node1_dual_ineq_r_plus, self.node2_dual_ineq_r_plus,
										self.nodeN_dual_ineq_r_plus]
				self.dual_ineq_i_plus_nodes = [self.node1_dual_ineq_i_plus, self.node2_dual_ineq_i_plus,
										self.nodeN_dual_ineq_i_plus]
			else:
				self.dual_ineq_r_minus_nodes = [self.node1_dual_ineq_r_minus, self.node2_dual_ineq_r_minus]
				self.dual_ineq_i_minus_nodes = [self.node1_dual_ineq_i_minus, self.node2_dual_ineq_i_minus]
				self.dual_ineq_r_plus_nodes = [self.node1_dual_ineq_r_plus, self.node2_dual_ineq_r_plus]
				self.dual_ineq_i_plus_nodes = [self.node1_dual_ineq_i_plus, self.node2_dual_ineq_i_plus]
		
			Nodes.mu_index += self.dual_ineq_r_minus_nodes + self.dual_ineq_r_plus_nodes + self.dual_ineq_i_minus_nodes + self.dual_ineq_i_plus_nodes
			
			Nodes.if_index += self.if_r_plus_nodes + self.if_r_minus_nodes + self.if_i_plus_nodes + self.if_i_minus_nodes
			
			Nodes.mu_r_index += self.dual_ineq_r_minus_nodes + self.dual_ineq_r_plus_nodes
								
			Nodes.mu_i_index += self.dual_ineq_i_minus_nodes + self.dual_ineq_i_plus_nodes
		
		return node_index_

	def assign_voltage_limit_nodes(self, node_index_):
		if self.bustype == 3:
			# try not limiting voltage on slack bus for now?
			return
		# if len(Nodes.bounded_bus_list) > 0 and self.name not in Nodes.bounded_bus_list:
		# 	return
		self.vmag2_inds = []
		self.Lvmag2_inds = []
		self.Lvmag2_hard_bounds = []
		self.umax_vmag2_inds = []
		self.umin_vmag2_inds = []
		self.vr_bounding_inds = []
		self.vi_bounding_inds = []
		self.Lr_bounding_inds = []
		self.Li_bounding_inds = []

		if self.isTriplex:
			self.node1_vmag2_index = node_index_.__next__()
			self.node1_Lvmag2_index = node_index_.__next__()
			self.node1_umax_vmag2_index = node_index_.__next__()
			self.node1_umin_vmag2_index = node_index_.__next__()
			self.node1_dual_eq_var_hard_bound = node_index_.__next__()

			self.node2_vmag2_index = node_index_.__next__()
			self.node2_Lvmag2_index = node_index_.__next__()
			self.node2_umax_vmag2_index = node_index_.__next__()
			self.node2_umin_vmag2_index = node_index_.__next__()
			self.node2_dual_eq_var_hard_bound = node_index_.__next__()

			self.vr_bounding_inds = [self.node1_Vr, self.node2_Vr]
			self.vi_bounding_inds = [self.node1_Vi, self.node2_Vi]
			self.Lr_bounding_inds = [self.node1_dual_eq_var_r, self.node2_dual_eq_var_r]
			self.Li_bounding_inds = [self.node1_dual_eq_var_i, self.node2_dual_eq_var_i]

			self.vmag2_inds.extend([self.node1_vmag2_index, self.node2_vmag2_index])			
			self.Lvmag2_inds.extend([self.node1_Lvmag2_index, self.node2_Lvmag2_index])
			self.umax_vmag2_inds.extend([self.node1_umax_vmag2_index, self.node2_umax_vmag2_index])
			self.umin_vmag2_inds.extend([self.node1_umin_vmag2_index, self.node2_umin_vmag2_index])
			self.Lvmag2_hard_bounds.extend([self.node1_dual_eq_var_hard_bound, self.node2_dual_eq_var_hard_bound])

			# self.node_set.extend([
			# 	self.node1_vmag2_index, self.node2_vmag2_index,
			# 	self.node1_Lvmag2_index, self.node2_Lvmag2_index,
			# 	self.node1_umax_vmag2_index, self.node2_umax_vmag2_index,
			# 	self.node1_umin_vmag2_index, self.node2_umin_vmag2_index
			# ])
		else:
			phase_cxn_str = list(_PHASE.keys())[list(_PHASE.values()).index(self.phases)]
			if 'A' in phase_cxn_str:
				self.nodeA_vmag2_index = node_index_.__next__()
				self.nodeA_Lvmag2_index = node_index_.__next__()
				self.nodeA_umax_vmag2_index = node_index_.__next__()
				self.nodeA_umin_vmag2_index = node_index_.__next__()
				self.nodeA_dual_eq_var_hard_bound = node_index_.__next__()
				self.vmag2_inds.append(self.nodeA_vmag2_index)
				self.Lvmag2_inds.append(self.nodeA_Lvmag2_index)
				self.umax_vmag2_inds.append(self.nodeA_umax_vmag2_index)
				self.umin_vmag2_inds.append(self.nodeA_umin_vmag2_index)
				self.vr_bounding_inds.append(self.nodeA_Vr)
				self.vi_bounding_inds.append(self.nodeA_Vi)
				self.Lr_bounding_inds.append(self.nodeA_dual_eq_var_r)
				self.Li_bounding_inds.append(self.nodeA_dual_eq_var_i)
				self.Lvmag2_hard_bounds.append(self.nodeA_dual_eq_var_hard_bound)
			if 'B' in phase_cxn_str:
				self.nodeB_vmag2_index = node_index_.__next__()
				self.nodeB_Lvmag2_index = node_index_.__next__()
				self.nodeB_umax_vmag2_index = node_index_.__next__()
				self.nodeB_umin_vmag2_index = node_index_.__next__()
				self.nodeB_dual_eq_var_hard_bound = node_index_.__next__()				
				self.vmag2_inds.append(self.nodeB_vmag2_index)
				self.Lvmag2_inds.append(self.nodeB_Lvmag2_index)
				self.umax_vmag2_inds.append(self.nodeB_umax_vmag2_index)
				self.umin_vmag2_inds.append(self.nodeB_umin_vmag2_index)
				self.vr_bounding_inds.append(self.nodeB_Vr)
				self.vi_bounding_inds.append(self.nodeB_Vi)
				self.Lr_bounding_inds.append(self.nodeB_dual_eq_var_r)
				self.Li_bounding_inds.append(self.nodeB_dual_eq_var_i)
				self.Lvmag2_hard_bounds.append(self.nodeB_dual_eq_var_hard_bound)
			if 'C' in phase_cxn_str:
				self.nodeC_vmag2_index = node_index_.__next__()
				self.nodeC_Lvmag2_index = node_index_.__next__()
				self.nodeC_umax_vmag2_index = node_index_.__next__()
				self.nodeC_umin_vmag2_index = node_index_.__next__()
				self.nodeC_dual_eq_var_hard_bound = node_index_.__next__()
				self.vmag2_inds.append(self.nodeC_vmag2_index)
				self.Lvmag2_inds.append(self.nodeC_Lvmag2_index)
				self.umax_vmag2_inds.append(self.nodeC_umax_vmag2_index)
				self.umin_vmag2_inds.append(self.nodeC_umin_vmag2_index)
				self.vr_bounding_inds.append(self.nodeC_Vr)
				self.vi_bounding_inds.append(self.nodeC_Vi)
				self.Lr_bounding_inds.append(self.nodeC_dual_eq_var_r)
				self.Li_bounding_inds.append(self.nodeC_dual_eq_var_i)
				self.Lvmag2_hard_bounds.append(self.nodeC_dual_eq_var_hard_bound)
			
		Nodes.vmag2_index.extend(self.vmag2_inds)
		Nodes.Lvmag2_index.extend(self.Lvmag2_inds)
		Nodes.umax_vmag2_index.extend(self.umax_vmag2_inds)
		Nodes.umin_vmag2_index.extend(self.umin_vmag2_inds)
		Nodes.vr_bounding_inds.extend(self.vr_bounding_inds)
		Nodes.vi_bounding_inds.extend(self.vi_bounding_inds)

		
		self.node_set.extend((self.vmag2_inds + self.Lvmag2_inds + self.umax_vmag2_inds + self.umin_vmag2_inds))

	def assign_voltage_unbalance_nodes(self, node_index_):
		# From Global Variables class, phases = 7 is ABC and self.phases is ABCN
		if self.phases == 7 or self.phases == 15:
			# New primal variables: intermediate nodes translating ABC to positive and negative sequence components
			self.extra_unbalance_node_pos_real = node_index_.__next__()
			self.extra_unbalance_node_pos_imag = node_index_.__next__()
			self.extra_unbalance_node_neg_real = node_index_.__next__()
			self.extra_unbalance_node_neg_imag = node_index_.__next__()

			# New dual equality variables (lambda)
			self.dual_unb_eq_var_pr = node_index_.__next__()
			self.dual_unb_eq_var_pi = node_index_.__next__()
			self.dual_unb_eq_var_nr = node_index_.__next__()
			self.dual_unb_eq_var_ni = node_index_.__next__()

			# New primal variable: intermediate variable tracking voltage unbalance 
			self.extra_voltage_unbalance_node_sqd = node_index_.__next__()
			self.dual_unb_eq_var_total = node_index_.__next__()

			# New dual inequality variables (mu)
			# Upper bound
			self.dual_unb_ineq_var_ub = node_index_.__next__()
			
			Nodes.voltage_unbalance_primal_nodes.extend([self.extra_unbalance_node_pos_real, self.extra_unbalance_node_pos_imag, 
												self.extra_unbalance_node_neg_real, self.extra_unbalance_node_neg_imag, self.extra_voltage_unbalance_node_sqd])
			Nodes.voltage_unbalance_dual_nodes.extend([self.dual_unb_eq_var_pr, self.dual_unb_eq_var_pi, 
											  	self.dual_unb_eq_var_nr, self.dual_unb_eq_var_ni, 
												self.dual_unb_eq_var_total, self.dual_unb_ineq_var_ub])
			Nodes.voltage_unbalance_tracking_nodes.extend([self.extra_voltage_unbalance_node_sqd])
			Nodes.voltage_unbalance_mu_nodes.extend([self.dual_unb_ineq_var_ub])
			self.node_set.extend([self.extra_unbalance_node_pos_real, self.extra_unbalance_node_pos_imag, 
												self.extra_unbalance_node_neg_real, self.extra_unbalance_node_neg_imag, 
			 									self.dual_unb_eq_var_pr, self.dual_unb_eq_var_pi, 
											  	self.dual_unb_eq_var_nr, self.dual_unb_eq_var_ni, self.dual_unb_ineq_var_ub])#, self.dual_unb_ineq_var_ub, self.dual_unb_ineq_var_lb])
		
		return node_index_

	def calcMagAng(self, V, flag_Vrange=False):
		if not self.isTriplex:
			# Line to ground voltage
			_Vag = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi])
			_Vbg = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi])
			_Vcg = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi])
			# Line to neutral voltage
			_Van = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi]) - (
				V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
			_Vbn = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi]) - (
				V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
			_Vcn = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi]) - (
				V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
			# Line to Line voltage
			_Vab = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi]) - (
				V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi])
			_Vbc = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi]) - (
				V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi])
			_Vca = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi]) - (
				V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi])

			# Mag and Ang Line to Gnd
			Vag_mag = float(np.absolute(_Vag))
			Vbg_mag = float(np.absolute(_Vbg))
			Vcg_mag = float(np.absolute(_Vcg))
			Vag_ang = float(np.angle(_Vag, deg=True)) if Vag_mag != 0 else 0
			Vbg_ang = float(np.angle(_Vbg, deg=True)) if Vbg_mag != 0 else 0
			Vcg_ang = float(np.angle(_Vcg, deg=True)) if Vcg_mag != 0 else 0
			Vag_ang_rad = float(np.angle(_Vag,
										 deg=False)) if Vag_mag != 0 else 0
			Vbg_ang_rad = float(np.angle(_Vbg,
										 deg=False)) if Vbg_mag != 0 else 0
			Vcg_ang_rad = float(np.angle(_Vcg,
										 deg=False)) if Vcg_mag != 0 else 0
			# Mag and Ang Line to Line
			Vab_mag = float(np.absolute(_Vab))
			Vbc_mag = float(np.absolute(_Vbc))
			Vca_mag = float(np.absolute(_Vca))
			Vab_ang = float(np.angle(_Vab, deg=True))
			Vbc_ang = float(np.angle(_Vbc, deg=True))
			Vca_ang = float(np.angle(_Vca, deg=True))
			#
			if self.phases & 0x08 == int(0x08):
				Vab_mag_pu = Vab_mag / (self.Vbase_LL * 1e3)
				Vbc_mag_pu = Vbc_mag / (self.Vbase_LL * 1e3)
				Vca_mag_pu = Vca_mag / (self.Vbase_LL * 1e3)
			else:
				Vab_mag_pu = Vab_mag / (self.Vbase_LL * 1e3)
				Vbc_mag_pu = Vbc_mag / (self.Vbase_LL * 1e3)
				Vca_mag_pu = Vca_mag / (self.Vbase_LL * 1e3)
			# Mag and Ang Line to Neutral
			Van_mag = float(np.absolute(_Van))
			Vbn_mag = float(np.absolute(_Vbn))
			Vcn_mag = float(np.absolute(_Vcn))
			Van_ang = float(np.angle(_Van, deg=True))
			Vbn_ang = float(np.angle(_Vbn, deg=True))
			Vcn_ang = float(np.angle(_Vcn, deg=True))
			Van_mag_pu = Van_mag / (self.Vbase_LN * 1e3)
			Vbn_mag_pu = Vbn_mag / (self.Vbase_LN * 1e3)
			Vcn_mag_pu = Vcn_mag / (self.Vbase_LN * 1e3)
			Vag_mag_pu = Vag_mag / self.Vnom
			Vbg_mag_pu = Vbg_mag / self.Vnom
			Vcg_mag_pu = Vcg_mag / self.Vnom

			# Store results
			# Line to neutral
			self.Vln = (np.around(Van_mag,
								  2), np.around(Vbn_mag,
												2), np.around(Vcn_mag, 2))
			self.Vln_pu = (np.around(Van_mag_pu,
									 5), np.around(Vbn_mag_pu,
												   5), np.around(Vcn_mag_pu, 5))
			self.Vll = (np.around(Vab_mag * 1e-3,
								  5), np.around(Vbc_mag * 1e-3, 5),
						np.around(Vca_mag * 1e-3, 5))
			self.Vll_pu = (np.around(Vab_mag_pu,
									 5), np.around(Vbc_mag_pu,
												   5), np.around(Vca_mag_pu, 5))
			self.Vlg = (np.around(Vag_mag * 1e-3,
								  5), np.around(Vbg_mag * 1e-3, 5),
						np.around(Vcg_mag * 1e-3, 5))
			self.Vlg_pu = (np.around(Vag_mag_pu,
									 5), np.around(Vbg_mag_pu,
												   5), np.around(Vcg_mag_pu, 5))
			self.Vbase_LL = np.around(self.Vbase_LL, 5)
			self.Vnom = np.around(self.Vnom, 5)

			if not flag_Vrange:
				self.V_mag = self.Vlg
				self.Vll_mag = self.Vll
			else:
				self.V_mag = self.Vlg_pu
				self.Vll_mag = self.Vll_pu

			self.V_ang = (np.around(Vag_ang,
									1), np.around(Vbg_ang,
												  1), np.around(Vcg_ang, 1))
			self.Vll_ang = (np.around(Vab_ang,
									  5), np.around(Vbc_ang,
													5), np.around(Vca_ang, 5))
			Vmag = np.asarray(self.V_mag)
			Vmax = np.max(Vmag) if Vmag.any() != 0 else 0
			Vavg = Vmag[Vmag != 0].mean() if Vmag.any() != 0 else 0

			if self.phases == 7 or self.phases == 15:
				try:
					self.V_unb = math.sqrt(float(V[self.extra_voltage_unbalance_node_sqd]))# np.abs(Vmax - Vavg) / Vavg if Vmag.any() != 0 else 0.0
				except:
					Var = float(V[self.nodeA_Vr])
					Vai = float(V[self.nodeA_Vi])
					Vbr = float(V[self.nodeB_Vr])
					Vbi = float(V[self.nodeB_Vi])
					Vcr = float(V[self.nodeC_Vr]) 
					Vci = float(V[self.nodeC_Vi])

					Vpr = 1/3*(Var - 0.5*Vbr - math.sqrt(3)/2*Vbi - 0.5*Vcr + math.sqrt(3)/2*Vci)
					Vpi = 1/3*(Vai + math.sqrt(3)/2*Vbr - 0.5*Vbi - math.sqrt(3)/2*Vcr -0.5*Vci)
					Vnr = 1/3*(Var - 0.5*Vbr + math.sqrt(3)/2*Vbi - 0.5*Vcr - math.sqrt(3)/2*Vci)
					Vni = 1/3*(Vai - math.sqrt(3)/2*Vbr - 0.5*Vbi + math.sqrt(3)/2*Vcr -0.5*Vci)
					self.V_unb = math.sqrt((Vnr**2 + Vni**2)/(Vpr**2 + Vpi**2))
				# self.V_unb = np.abs(Vmax - Vavg) / Vavg 
				self.V_unb = np.around(self.V_unb * 100, 6)
			else:
				self.V_unb = 'N/A'
			# Calculate swing current
			if self.bustype == 3:
				self.Type = "Slack"
				currentA = -V[self.slack_nodeA_Vr] - 1j * V[self.slack_nodeA_Vi]
				currentB = -V[self.slack_nodeB_Vr] - 1j * V[self.slack_nodeB_Vi]
				currentC = -V[self.slack_nodeC_Vr] - 1j * V[self.slack_nodeC_Vi]

				self.Ia = currentA
				self.Ib = currentB
				self.Ic = currentC

				self.Ia_mag = float(np.absolute(currentA))
				self.Ib_mag = float(np.absolute(currentB))
				self.Ic_mag = float(np.absolute(currentC))
				self.Ia_ang = float(np.angle(currentA, deg=True))
				self.Ib_ang = float(np.angle(currentB, deg=True))
				self.Ic_ang = float(np.angle(currentC, deg=True))
				# Calculate Slack Bus Power
				self.Sa = currentA * self.Va
				self.Sb = currentB * self.Vb
				self.Sc = currentC * self.Vc
				self.S = np.absolute(self.Sa.real) + np.absolute(
					self.Sb.real) + np.absolute(self.Sc.real) + 1j * (
						np.absolute(self.Sa.imag) + np.absolute(self.Sb.imag) +
						np.absolute(self.Sc.imag))
				if self.bustype == 3 and self.print_slack:
					logging.debug('Slack bus power consumption is %f + %fi' %
								  (self.S.real, self.S.imag))
					logging.debug(
						'Slack Bus Voltages are : (%f, %f, %f, %f, %f, %f)' %
						(Van_mag, Vbn_mag, Vcn_mag, Van_ang, Vbn_ang, Vcn_ang))
					logging.debug(
						'Slack Bus Currents are : (%f, %f, %f, %f, %f, %f)' %
						(self.Ia_mag, self.Ib_mag, self.Ic_mag, self.Ia_ang,
						 self.Ib_ang, self.Ic_ang))
		else:
			# Line to ground voltage
			_V1g = (V[self.node1_Vr] + 1j * V[self.node1_Vi])
			_V2g = -(V[self.node2_Vr] + 1j * V[self.node2_Vi])
			_V12 = (V[self.node1_Vr] + 1j * V[self.node1_Vi]) - (
				V[self.node2_Vr] + 1j * V[self.node2_Vi])
			# Mag and Ang Line to Gnd
			V1g_mag = float(np.absolute(_V1g))
			V2g_mag = float(np.absolute(_V2g))
			V12_mag = float(np.absolute(_V12))
			V1g_pu = V1g_mag / self.Vnom
			V2g_pu = V2g_mag / self.Vnom
			V12_pu = V12_mag / self.Vbase_LL
			V1g_ang = float(np.angle(_V1g, deg=True)) if V1g_mag else 0
			V1g_ang_rad = float(np.angle(_V1g, deg=False)) if V1g_mag else 0
			V2g_ang = float(np.angle(_V2g, deg=True)) if V2g_mag else 0
			V2g_ang_rad = float(np.angle(_V2g, deg=False)) if V2g_mag else 0
			V12_ang = float(np.angle(_V12, deg=True)) if V12_mag else float(0)
			self.Vln = (np.around(V1g_mag,
								  5), np.around(V2g_mag, 5))
			self.Vll = np.around(V12_mag * 1e-3, 5)
			self.Vln_pu = (V1g_pu, V2g_pu)
			self.Vll_pu = np.around(V12_pu * 1e-3, 5)
			self.V_mag = (V1g_mag, V2g_mag)
			self.V_ang = (np.around(V1g_ang, 5), np.around(V2g_ang, 5))
			self.Vll_ang = np.around(V12_ang, 5)
			self.V_unb = 0


class TriplexNodes:

	def __init__(self, name, ptr, bustype, current1, current12, current2,
				 currentN, impedance1, impedance12, impedance2,
				 maximum_voltage_error, mean_repair_time, nominal_voltage,
				 parent, phases, power1, power12, power2, reference_bus, shunt1,
				 shunt12, shunt2, voltage1, voltage12, voltage1N, voltage2,
				 voltage2N, voltageN):
		self.name = name
		self.ptr = int(ptr)
		self.bustype = int(bustype)
		self.parent = int(parent) if parent else None
		self.phases = int(phases)
		self.nominal_voltage = int(nominal_voltage)
		self.voltage1 = (voltage1) if voltage1 else None
		self.voltage2 = (voltage2) if voltage2 else None
		self.voltage12 = (voltage12) if voltage12 else None
		self.voltage1N = (voltage1N) if voltage1N else None
		self.voltage2N = (voltage2N) if voltage2N else None
		self.voltageN = (voltageN) if voltage12 else None
		self.reference_bus = reference_bus
		self.maximum_voltage_error = maximum_voltage_error
		self.mean_repair_time = mean_repair_time
		# Triplex Node Load Parameters
		self.current1 = (current1) if current1 else None
		self.current2 = (current2) if current2 else None
		self.current12 = (current12) if current12 else None
		self.currentN = (currentN) if currentN else None
		self.impedance1 = (impedance1) if impedance1 else None
		self.impedance2 = (impedance2) if impedance2 else None
		self.impedance12 = (impedance12) if impedance12 else None
		self.power1 = (power1) if power1 else None
		self.power2 = (power2) if power2 else None
		self.power12 = (power12) if power12 else None
		self.shunt1 = (shunt1) if shunt1 else None
		self.shunt2 = (shunt2) if shunt2 else None
		self.shunt12 = (shunt12) if shunt12 else None
