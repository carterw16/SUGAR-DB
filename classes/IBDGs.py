"""Provides various Inverter-Based Distributed Generator models for distribution system analysis.

    Author(s): Naeem Turner-Bandele
    Created Date: 09-25-2019
    Updated Date: 10-14-2021
    Email: nturnerb@cmu.edu
    Status: Development

"""

# import built-in modules
from __future__ import division

from types import SimpleNamespace

# import third-party modules
import numpy as np
import logging
import warnings

# import libraries
from lib.parse_setpoints import parse_setpoints

# import custom modules
from classes.Elements import NonlinearElement
from classes.GlobalVars import MVAbase, connected_nodes, _PHASE
from classes.IBDGVoltageControl import QVSplineControl, QVPatchingControl, QVLoopControl
from classes.Nodes import Nodes


class IBDG(NonlinearElement):

	def __init__(self,
	             ID,
	             name,
	             phases,
	             nominal_voltage,
	             voltage_setpoint,
	             connection_type,
	             rated_power,
	             P,
	             Q,
	             active_rating,
	             reactive_rating,
	             power_factor,
	             irradiance,
	             Pmpp,
	             inverter_settings,
	             is_triplex,
	             sensor="GVCC",
	             fixed_Q=True):

		super(IBDG, self).__init__()

		self.ID = ID
		self.name = name  # Name of the DG
		self.phases = int(phases) if phases else None
		self.ibdg_phases = [0x01, 0x02, 0x04]
		self.Vnom = nominal_voltage
		self.Vset = voltage_setpoint if voltage_setpoint != 0 else 1.0
		self.connection_type = connection_type
		self.is_triplex = is_triplex
		self.base_power = MVAbase * 1e6
		self.rated_power = rated_power if rated_power else np.abs(complex(active_rating, reactive_rating))
		self.P3 = Pmpp * irradiance if Pmpp else P
		self.Q3 = Q
		self.Q3_prev = 0.0
		self.pf = power_factor
		self.Sensor = sensor
		self.fixed_Q = fixed_Q

		self.Q3_final = 0
		self.Qmax = self.Q3 * 2 if self.Q3 != 0 else self.P3 / 2
		self.Qmin = -self.Q3 * 2 if self.Q3 != 0 else -self.P3 / 2
		self.Qscale = self.Qmin

		# Parse Inverter Settings
		settings = parse_setpoints(inverter_settings) if inverter_settings else None
		if inverter_settings:
			self.inverter_settings = next(setting for setting in settings if setting["Name"] == self.name)
		else:
			self.inverter_settings = None

		self.Qlim_type = self.inverter_settings["QControl"] if self.inverter_settings else "Patching"
		self.P3_max = active_rating
		self.ibdg_counter = 0

		# starting voltage setpoints for the IBDGs
		self.Va_pu = float(self.inverter_settings["Va_pu"]) if self.inverter_settings else 1.0
		self.Vb_pu = float(self.inverter_settings["Vb_pu"]) if self.inverter_settings else 1.0
		self.Vc_pu = float(self.inverter_settings["Vc_pu"]) if self.inverter_settings else 1.0

		# # Initialize Inverter Control Parameters # #
		self.V1 = 0
		self.V2 = 0
		self.V3 = 0
		self.V4 = 0

		self.perc_Q = 0

		self.Q2A = 0
		self.Q3A = 0
		self.Q4A = 0

		self.p1 = 0
		self.p2 = 0
		self.p3 = 0
		self.p4 = 0

		self.Vph = self.Vnom / np.sqrt(3)
		self.flag_Qlim = False
		self.Calc_Q = True
		self.flag_Qmax = False
		self.flag_Qmin = False
		self.h_factor_init = -1
		self.Qphase = 'C'

		# # Calculate and initialize IBDG params # #
		self.alpha_r = 0
		self.alpha_i = 0
		self.beta_r = 0
		self.beta_i = 0
		self.set_ibdg_params()

		# # Initialize Transformation Matrices # #
		self.a = (-1 / 2) + 1j * np.sqrt(3) / 2
		self.a2 = (-1 / 2) - 1j * np.sqrt(3) / 2

		sqrt3_2 = np.sqrt(3) / 2

		# Pr = Real Positive Sequence Transformation Matrix
		self.Pr = (1 / 3) * np.array([[1, 0, -1 / 2, -sqrt3_2, -1 / 2, sqrt3_2],
		                              [-1 / 2, sqrt3_2, 1, 0, -1 / 2, -sqrt3_2],
		                              [-1 / 2, -sqrt3_2, -1 / 2, sqrt3_2, 1, 0]
		                              ])

		# Pi = Imaginary Positive Sequence Transformation Matrix
		self.Pi = (1 / 3) * np.array([[0, 1, sqrt3_2, -1 / 2, -sqrt3_2, -1 / 2],
		                              [-sqrt3_2, -1 / 2, 0, 1, sqrt3_2, -1 / 2],
		                              [sqrt3_2, -1 / 2, -sqrt3_2, -1 / 2, 0, 1]
		                              ])

		# Nr = Real Negative Sequence Transformation Matrix
		self.Nr = (1 / 3) * np.array([[1, 0, -1 / 2, sqrt3_2, -1 / 2, -sqrt3_2],
		                              [-1 / 2, -sqrt3_2, 1, 0, -1 / 2, sqrt3_2],
		                              [-1 / 2, sqrt3_2, -1 / 2, -sqrt3_2, 1, 0]
		                              ])

		# Ni = Imaginary Negative Sequence Transformation Matrix
		self.Ni = (1 / 3) * np.array([[0, 1, -sqrt3_2, -1 / 2, sqrt3_2, -1 / 2],
		                              [sqrt3_2, -1 / 2, 0, 1, -sqrt3_2, -1 / 2],
		                              [-sqrt3_2, -1 / 2, sqrt3_2, -1 / 2, 0, 1]
		                              ])

		# # Initialize IBDG Currents # #
		self.Ia_mag = 0.0
		self.Ia_ang = 0.0
		self.Ib_mag = 0.0
		self.Ib_ang = 0.0
		self.Ic_mag = 0.0
		self.Ic_ang = 0.0

		# # Initialize IBDG Powers # #
		self.Pa = 0.0
		self.Pb = 0.0
		self.Pc = 0.0
		self.Qa = 0.0
		self.Qb = 0.0
		self.Qc = 0.0

		# # Assign and Initialize Nodes # #
		# Create Reactive Power Nodes and current measurement nodes if enabled
		self.node_Q3 = None



	@staticmethod
	def calc_transforms(Var,
	                    Vai,
	                    Vbr,
	                    Vbi,
	                    Vcr,
	                    Vci):
		Va = Var + 1j * Vai
		Vb = Vbr + 1j * Vbi
		Vc = Vcr + 1j * Vci

		a = (-1 / 2) + 1j * np.sqrt(3) / 2
		a2 = (-1 / 2) - 1j * np.sqrt(3) / 2

		Vpos = (1 / 3) * (Va + Vb * a + Vc * a2)
		Vneg = (1 / 3) * (Va + Vb * a2 + Vc * a)

		return Vpos, Vneg

	def assign_nodes(self,
	                 node_index_):
		if not self.fixed_Q:
			self.node_Q3 = node_index_.__next__()
			# Assign the new Q nodes to the the Nodes Q_index
			Nodes.Q_index.extend([self.node_Q3])
			Nodes.Q_mag.extend(([self.Q3]))

		# Add the DG id to connected nodes
		connected_nodes.add(self.ID)

		return node_index_

	def symmetrical_components(self,
	                           Va,
	                           Vb,
	                           Vc):
		a = complex((-1 / 2), np.sqrt(3) / 2)
		a2 = complex((-1 / 2), -np.sqrt(3) / 2)

		T = np.array([
			[1, 1, 1],
			[1, a, a2],
			[1, a2, a]
		])
		Vabc = np.array([Va, Vb, Vc])
		Va_seq = (1 / 3) * T.dot(Vabc)
		Vb_seq = np.array([Va_seq[0], a2 * Va_seq[1], a * Va_seq[2]])
		Vc_seq = np.array([Va_seq[0], a * Va_seq[1], a2 * Va_seq[2]])

		self.Vzero = np.array([Va_seq[0], Vb_seq[1], Vc_seq[2]])
		self.Vpos = np.array([Va_seq[1], Vb_seq[1], Vc_seq[1]])
		self.Vneg = np.array([Va_seq[2], Vb_seq[2], Vc_seq[2]])

	def set_ibdg_params(self):
		if self.Sensor == "GVCC":
			self.alpha_r = 1
		elif self.Sensor == "CVCC":
			self.alpha_r = 1
			self.alpha_i = 0
			self.beta_r = 0
		elif self.Sensor == "GVGC":
			self.alpha_r = 1
			self.alpha_i = 0
			self.beta_r = 0
			self.beta_i = 0
		else:
			self.alpha_r = 1
			self.alpha_i = 0
			self.beta_r = 0
			self.beta_i = 0

	def get_nodes_states(self,
	                     node,
						 node_key,
	                     V):
		"""
		Find the indices of the solution vector and the values at those indices.
		:param node: The vector of all system nodes.
		:param V: The solution vector.
		:return: Three dictionaries which hold the phase and sequence nodes.
		"""
		Var = 0.0
		Vai = 0.0
		Vbr = 0.0
		Vbi = 0.0
		Vcr = 0.0
		Vci = 0.0

		Va_mag2 = 0.0
		Vb_mag2 = 0.0
		Vc_mag2 = 0.0

		Va_mag = 0.0
		Vb_mag = 0.0
		Vc_mag = 0.0

		nodeA_Vr = -1
		nodeA_Vi = -1
		nodeB_Vr = -1
		nodeB_Vi = -1
		nodeC_Vr = -1
		nodeC_Vi = -1

		if self.phases & 0x1 == 1:  # Check for phase A
			if not self.is_triplex:
				nodeA_Vr = node[node_key[self.ID]].nodeA_Vr
				nodeA_Vi = node[node_key[self.ID]].nodeA_Vi
			else:
				nodeA_Vr = node[node_key[self.ID]].node2_Vr
				nodeA_Vi = node[node_key[self.ID]].node2_Vi
			Var = V[nodeA_Vr]
			Vai = V[nodeA_Vi]
		if self.phases & 0x2 == 2:  # Check for phase B
			if not self.is_triplex:
				nodeB_Vr = node[node_key[self.ID]].nodeB_Vr
				nodeB_Vi = node[node_key[self.ID]].nodeB_Vi
			else:
				nodeB_Vr = node[node_key[self.ID]].node2_Vr
				nodeB_Vi = node[node_key[self.ID]].node2_Vi
			Vbr = V[nodeB_Vr]
			Vbi = V[nodeB_Vi]

		if self.phases & 0x4 == 4:  # Check for phase C
			if not self.is_triplex:
				nodeC_Vr = node[node_key[self.ID]].nodeC_Vr
				nodeC_Vi = node[node_key[self.ID]].nodeC_Vi
			else:
				nodeC_Vr = node[node_key[self.ID]].node2_Vr
				nodeC_Vi = node[node_key[self.ID]].node2_Vi
			Vcr = V[nodeC_Vr]
			Vci = V[nodeC_Vi]
		if not self.fixed_Q:
			node_Q3 = self.node_Q3
			node[node_key[self.ID]].node_set.extend([self.node_Q3])

		# Find the node voltages and reactive power

		if not self.fixed_Q:
			(Q3) = (V[self.node_Q3])
		else:
			Q3 = self.Q3
			node_Q3 = -1

		Va_ang = np.angle(Var + 1j * Vai)
		Vb_ang = np.angle(Vbr + 1j * Vbi)
		Vc_ang = np.angle(Vcr + 1j * Vci)

		Vabc = {
			'VR': [Var, Vbr, Vcr],
			'VI': [Vai, Vbi, Vci],
			'nodeVR': [nodeA_Vr, nodeB_Vr, nodeC_Vr],
			'nodeVI': [nodeA_Vi, nodeB_Vi, nodeC_Vi],
			'Q': [Q3],
			'nodeQ3': [node_Q3],
			'ang': [Va_ang, Vb_ang, Vc_ang],
			'Vmag2': [Va_mag2, Vb_mag2, Vc_mag2],
			'Vmag': [Va_mag, Vb_mag, Vc_mag]
		}

		Vabc = SimpleNamespace(**Vabc)

		return Vabc

	def determine_setpoint_profile(self):
		"""Set the inverter control setpoints V1, V2, V3, and V4 based on base voltage setpoint.

		Returns: None

		"""
		# Profile 1
		if 0.995 < self.Vset < 1.005:
			self.V1 = 0.98
			self.V2 = 0.995
			self.V3 = 1.005
			self.V4 = 1.02
		# Profile 2
		elif 1.005 < self.Vset < 1.015:
			self.V1 = 0.99
			self.V2 = 1.005
			self.V3 = 1.015
			self.V4 = 1.03
		elif 1.01 < self.Vset < 1.02:
			self.V1 = 0.995
			self.V2 = 1.01
			self.V3 = 1.02
			self.V4 = 1.035
		# Profile 3
		elif 1.015 < self.Vset < 1.025:
			self.V1 = 1.0
			self.V2 = 1.015
			self.V3 = 1.025
			self.V4 = 1.04
		# Profile 4
		elif 1.02 < self.Vset < 1.03:
			self.V1 = 1.005
			self.V2 = 1.02
			self.V3 = 1.03
			self.V4 = 1.045
		elif 1.025 < self.Vset < 1.035:
			self.V1 = 1.01
			self.V2 = 1.025
			self.V3 = 1.035
			self.V4 = 1.05
		elif 1.03 < self.Vset < 1.04:
			self.V1 = 1.015
			self.V2 = 1.03
			self.V3 = 1.04
			self.V4 = 1.055
		# Profile 5
		elif 1.035 < self.Vset < 1.045:
			self.V1 = 1.02
			self.V2 = 1.035
			self.V3 = 1.045
			self.V4 = 1.06
		# Profile 6
		elif 1.05 < self.Vset < 1.06:
			self.V1 = 1.035
			self.V2 = 1.05
			self.V3 = 1.06
			self.V4 = 1.075
		# Profile 7
		elif 1.065 < self.Vset < 1.075:
			self.V1 = 1.05
			self.V2 = 1.065
			self.V3 = 1.075
			self.V4 = 1.09
		# Profile 8
		elif 1.08 < self.Vset < 1.09:
			self.V1 = 1.065
			self.V2 = 1.08
			self.V3 = 1.09
			self.V4 = 1.15
		# Profile 9
		elif 0.98 < self.Vset < 0.99:
			self.V1 = 0.965
			self.V2 = 0.98
			self.V3 = 0.99
			self.V4 = 1.005
		elif 0.985 < self.Vset < 0.995:
			self.V1 = 0.97
			self.V2 = 0.985
			self.V3 = 0.995
			self.V4 = 1.01
		# Profile 10
		elif 0.965 < self.Vset < 0.975:
			self.V1 = 0.95
			self.V2 = 0.965
			self.V3 = 0.975
			self.V4 = 0.99
		elif 1.045 < self.Vset < 1.055:
			self.V1 = 1.03
			self.V2 = 1.045
			self.V3 = 1.055
			self.V4 = 1.07
		else:
			self.V2 = self.Vset - 0.005
			self.V3 = self.Vset + 0.005
			self.V1 = self.V2 - 0.015
			self.V4 = self.V3 + 0.015

	def initialize_setpoints(self):
		"""Create voltage setpoints for the inverter reactive power controls.

		Args:
			inverter_settings (Dict):  Stores the voltage setpoints and any other relevant parameters for reactive
			power control.

		Returns: None but initializes class attributes V1, V2, V3, and V4.

		"""
		inverter_settings = self.inverter_settings
		if inverter_settings is not None and "V1" in inverter_settings:
			self.V1 = float(inverter_settings["V1"])
			self.V2 = float(inverter_settings["V2"])
			self.V3 = float(inverter_settings["V3"])
			self.V4 = float(inverter_settings["V4"])
		else:
			self.determine_setpoint_profile()

	def initialize_ibdg(self,
	                    node,
						node_key,
	                    V,
	                    lf):

		node_id = node_key[self.ID]

		Var = 0.0
		Vai = 0.0
		Vbr = 0.0
		Vbi = 0.0
		Vcr = 0.0
		Vci = 0.0

		if self.phases & 0x1 == 1:  # Check for phase A
			if not self.is_triplex:
				Var = V[node[node_id].nodeA_Vr] * self.Va_pu
				Vai = V[node[node_id].nodeA_Vi] * self.Va_pu
			else:
				Var = V[node[node_id].node2_Vr] * self.Va_pu
				Vai = V[node[node_id].node2_Vi] * self.Va_pu
		if self.phases & 0x2 == 2:  # Check for phase B
			if not self.is_triplex:
				Vbr = V[node[node_id].nodeB_Vr] * self.Vb_pu
				Vbi = V[node[node_id].nodeB_Vi] * self.Vb_pu
			else:
				Vbr = V[node[node_id].node2_Vr] * self.Vb_pu
				Vbi = V[node[node_id].node2_Vi] * self.Vb_pu
		if self.phases & 0x4 == 4:  # Check for phase C
			if not self.is_triplex:
				Vcr = V[node[node_id].nodeC_Vr] * self.Vc_pu
				Vci = V[node[node_id].nodeC_Vi] * self.Vc_pu
			else:
				Vcr = V[node[node_id].node2_Vr] * self.Vc_pu
				Vci = V[node[node_id].node2_Vi] * self.Vc_pu

		Vabc_mag = np.array([
			Var, Vai, Vbr, Vbi, Vcr,
			Vci
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		self.V_pos = Vposr + 1j * Vposi
		self.V_neg = Vnegr + 1j * Vnegi

		# # Calculate Initial Positive and Negative Sequence Voltages # #
		Va = complex(Var, Vai)
		Vb = complex(Vbr, Vbi)
		Vc = complex(Vcr, Vci)
		self.symmetrical_components(Va, Vb, Vc)

		# # calculate k1 and k2 for FPNSC type DG
		if isinstance(self, FPNSC):
			self.determine_k1_k2()

		# # Overcurrent Control # #
		if not self.fixed_Q:
			Q3_ref = self.Q3
			Vpos, Vneg = self.calc_transforms(Var, Vai, Vbr, Vbi, Vcr, Vci)
			Imax = self.calc_max_current(Q3_ref, lf)
			self.calc_Q_limits(Imax)
			self.Qscale = self.Qmin
			Vnom = self.Vnom

			Va = float(np.abs(Var + 1j * Vai))
			Vb = float(np.abs(Vbr + 1j * Vbi))
			Vc = float(np.abs(Vcr + 1j * Vci))

			if np.abs(self.Vnom - 120.0) > 10:
				if np.abs(Va - self.Vph) < 30:
					Vnom_ph = self.Vph
				else:
					Vnom_ph = self.Vph * np.sqrt(3)
			else:
				Vnom_ph = self.Vnom

			self.Vnom_ph = Vnom_ph

			if self.Qphase == 'A':
				if self.phases & 0x1 == 1:  # Check for phase A
					self.Vset = Va / Vnom_ph
				elif self.phases & 0x2 == 2:  # Check for phase B
					self.Vset = Vb / Vnom_ph
				elif self.phases & 0x4 == 4:  # Check for phase C
					self.Vset = Vc / Vnom_ph

			elif self.Qphase == 'B':
				if self.phases & 0x1 == 1:  # Check for phase A
					self.Vset = Va / Vnom_ph
				elif self.phases & 0x2 == 2:  # Check for phase B
					self.Vset = Vb / Vnom_ph
				elif self.phases & 0x4 == 4:  # Check for phase C
					self.Vset = Vc / Vnom_ph

			else:
				if self.phases & 0x1 == 1:  # Check for phase A
					self.Vset = Va / Vnom_ph
				elif self.phases & 0x2 == 2:  # Check for phase B
					self.Vset = Vb / Vnom_ph
				elif self.phases & 0x4 == 4:  # Check for phase C
					self.Vset = Vc / Vnom_ph
			# self.Vset = Vc / Vnom_ph

			# # Correct Voltage Setpoint if Incorrect # #
			if self.Vset > 1.1:
				if self.phases & 0x1 == 1:  # Check for phase A
					self.Vset = Va / (Vnom_ph * np.sqrt(3))
				elif self.phases & 0x2 == 2:  # Check for phase B
					self.Vset = Vb / (Vnom_ph * np.sqrt(3))
				elif self.phases & 0x4 == 4:  # Check for phase C
					self.Vset = Vc / (Vnom_ph * np.sqrt(3))

			self.initialize_setpoints()

			# #  Initialize Q-V Control # #
			if self.Qlim_type == "Patching":
				self.QVControl = QVPatchingControl(self.V1, self.V2, self.V3, self.V4,
				                                   self.Vnom_ph, self.Qscale)
			elif self.Qlim_type == "Spline":
				self.QVControl = QVSplineControl(self.V1, self.V2, self.V3, self.V4,
				                                 self.Vnom, self.Qscale)
			else:
				self.QVControl = QVLoopControl(self.Qmax, self.Qmin, self.rated_power, self.Vnom_ph, self.Va_pu,
				                               self.Vb_pu, self.Vc_pu, self.V1, self.V2, self.V3, self.V4)

		if not self.fixed_Q:
			V[self.node_Q3] = self.QVControl.get_Q(self.Vset * self.Vnom)

		Q = V[self.node_Q3]
		if Q > self.Qmax:
			self.h_factor_init = Q / self.Qmax
		elif Q < self.Qmin:
			self.h_factor_init = Q / self.Qmin
		else:
			self.h_factor_init = 1

		# # Update Nodal Voltages for IBDGs # #
		if self.phases & 0x1 == 1:  # Check for phase A
			if not self.is_triplex:
				V[node[node_id].nodeA_Vr] = Var
				V[node[node_id].nodeA_Vi] = Vai
			else:
				V[node[node_id].node2_Vr] = Var
				V[node[node_id].node2_Vi] = Vai
		if self.phases & 0x2 == 2:  # Check for phase B
			if not self.is_triplex:
				V[node[node_id].nodeB_Vr] = Vbr
				V[node[node_id].nodeB_Vi] = Vbi
			else:
				V[node[node_id].node2_Vr] = Vbr
				V[node[node_id].node2_Vi] = Vbi
		if self.phases & 0x4 == 4:  # Check for phase C
			if not self.is_triplex:
				V[node[node_id].nodeC_Vr] = Vcr
				V[node[node_id].nodeC_Vi] = Vci
			else:
				V[node[node_id].node2_Vr] = Vcr
				V[node[node_id].node2_Vi] = Vci

		return V

	def calc_Q3(self,
	            Vabc):
		Var = Vabc.VR[0]
		Vbr = Vabc.VR[1]
		Vcr = Vabc.VR[2]

		Vai = Vabc.VI[0]
		Vbi = Vabc.VI[1]
		Vci = Vabc.VI[2]

		Q3 = Vabc.Q[0]
		# self.Q3_prev = Q3

		Va = np.abs(Var + 1j * Vai)
		Vb = np.abs(Vbr + 1j * Vbi)
		Vc = np.abs(Vcr + 1j * Vci)


		Va_set = self.Vset
		Vb_set = self.Vset
		Vc_set = self.Vset
		Vnom = self.Vnom

		Var_pu = Var / Vnom
		Vai_pu = Vai / Vnom
		Vbr_pu = Vbr / Vnom
		Vbi_pu = Vbi / Vnom
		Vcr_pu = Vcr / Vnom
		Vci_pu = Vci / Vnom

		np.abs(Var_pu + 1j * Vai_pu)
		np.abs(Vbr_pu + 1j * Vbi_pu)
		np.abs(Vcr_pu + 1j * Vci_pu)

		if self.Qlim_type == 'Discontinuous':
			FQ3 = Va_set ** 2 - Var ** 2 - Vai ** 2 \
			      + Vb_set ** 2 - Vbr ** 2 - Vbi ** 2 + Vc_set ** 2 - Vcr ** 2 - Vci ** 2
			dFQ3_dVar = -2 * Var
			dFQ3_dVai = -2 * Vai

			dFQ3_dVbr = -2 * Vbr
			dFQ3_dVbi = -2 * Vbi

			dFQ3_dVcr = -2 * Vcr
			dFQ3_dVci = -2 * Vci

			dFQ3_dQ3 = 0

			dFQ3 = {
				'dVr': [dFQ3_dVar, dFQ3_dVbr, dFQ3_dVcr],
				'dVi': [dFQ3_dVai, dFQ3_dVbi, dFQ3_dVci],
				'dQ': dFQ3_dQ3
			}

		else:
			# # Piecewise Continuous Q-V Control # #

			if self.Qphase == 'A':
				if self.phases & 0x1 == 1:  # Check for phase A
					V_k = Va
					Vr_k = Var
					Vi_k = Vai
				elif self.phases & 0x2 == 2:  # Check for phase B
					V_k = Vb
					Vr_k = Vbr
					Vi_k = Vbi
				elif self.phases & 0x4 == 4:  # Check for phase C
					V_k = Vc
					Vr_k = Vcr
					Vi_k = Vci
			elif self.Qphase == 'B':
				if self.phases & 0x1 == 1:  # Check for phase A
					V_k = Va
					Vr_k = Var
					Vi_k = Vai
				elif self.phases & 0x2 == 2:  # Check for phase B
					V_k = Vb
					Vr_k = Vbr
					Vi_k = Vbi
				elif self.phases & 0x4 == 4:  # Check for phase C
					V_k = Vc
					Vr_k = Vcr
					Vi_k = Vci
			else:
				if self.phases & 0x1 == 1:  # Check for phase A
					V_k = Va
					Vr_k = Var
					Vi_k = Vai
				elif self.phases & 0x2 == 2:  # Check for phase B
					V_k = Vb
					Vr_k = Vbr
					Vi_k = Vbi
				elif self.phases & 0x4 == 4:  # Check for phase C
					V_k = Vc
					Vr_k = Vcr
					Vi_k = Vci

			phase_val_list = list(_PHASE.values())
			phase_key_list = list(_PHASE.keys())
			phase_idx = phase_val_list.index(self.phases)

			phase = phase_key_list[phase_idx] if self.phases != 'ABC' or self.phases != 'ABCN' else self.Qphase

			if self.Qlim_type != 'Q-V_Switch':
				self.flag_Qmax = False
				self.flag_Qlim = False
				self.flag_Qmin = False
				dFQ3, FQ3, flag_Qmax, flag_Qmin, flag_Qlim = self.QVControl.control_QV(
					V_k, Vr_k, Vi_k, Q3, phase)
			else:
				dFQ3, FQ3, flag_Qmax, flag_Qmin, flag_Qlim = self.QVControl.control_QV(
					V_k, Vr_k, Vi_k, Q3, phase, self.flag_Qlim, self.flag_Qmin, self.flag_Qmax)

			self.flag_Qmax = flag_Qmax
			self.flag_Qmin = flag_Qmin
			self.flag_Qlim = flag_Qlim

		dFQ3 = SimpleNamespace(**dFQ3)

		return dFQ3, FQ3

	def stamp_Q3_partials(self,
	                      Vabc,
	                      dFQ3,
	                      FQ3_hist,
	                      Ynlin_val,
	                      Ynlin_row,
	                      Ynlin_col,
	                      Jnlin_val,
	                      Jnlin_row,
	                      idx_Y,
	                      idx_J,
	                      homotopy_enabled=False,
	                      h_factor=None):

		nodeA_Vr = Vabc.nodeVR[0]
		nodeB_Vr = Vabc.nodeVR[1]
		nodeC_Vr = Vabc.nodeVR[2]

		nodeA_Vi = Vabc.nodeVI[0]
		nodeB_Vi = Vabc.nodeVI[1]
		nodeC_Vi = Vabc.nodeVI[2]

		node_Q3 = Vabc.nodeQ3[0]

		Var = Vabc.VR[0]
		Vbr = Vabc.VR[1]
		Vcr = Vabc.VR[2]

		Vai = Vabc.VI[0]
		Vbi = Vabc.VI[1]
		Vci = Vabc.VI[2]

		Q3 = Vabc.Q[0]

		if homotopy_enabled and self.Qlim_type != 'Q-V_Switch':
			h_bound = (1 + (self.h_factor_init + 1) * h_factor)
		else:
			h_bound = 1

		if self.flag_Qmax and self.flag_Qlim:
			Qmax = self.Qmax * h_bound
			idx_Y = self.stamp_Y(node_Q3, node_Q3, 1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_J = self.stamp_J(node_Q3, Qmax, Jnlin_val, Jnlin_row, idx_J)
		elif self.flag_Qmin and self.flag_Qlim:
			Qmin = self.Qmin * h_bound
			idx_Y = self.stamp_Y(node_Q3, node_Q3, 1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_J = self.stamp_J(node_Q3, Qmin, Jnlin_val, Jnlin_row, idx_J)
		else:

			if self.phases & 0x01:
				# # Stamp Voltage A Magnitude # #
				Vr = Var
				Vi = Vai
				Vmag2 = Vabc.Vmag2[0]

				Vmag2_hist = Vr ** 2 + Vi ** 2 - Vmag2

				dVmag2_dVr = 2 * Vr
				dVmag2_dVi = 2 * Vi
				dVmag2_dVmag2 = -1

				_Vmag2 = Vmag2_hist - dVmag2_dVr * Vr - dVmag2_dVi * Vi - dVmag2_dVmag2 * Vmag2

				# # Stamp Q3 Va partials # #
				idx_Y = self.stamp_Y(node_Q3, nodeA_Vr, dFQ3.dVr[0], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y = self.stamp_Y(node_Q3, nodeA_Vi, dFQ3.dVi[0], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

			if self.phases & 0x02:
				# # Stamp Voltage B Magnitude # #
				Vr = Vbr
				Vi = Vbi
				Vmag2 = Vabc.Vmag2[1]

				Vmag2_hist = Vr ** 2 + Vi ** 2 - Vmag2

				dVmag2_dVr = 2 * Vr
				dVmag2_dVi = 2 * Vi
				dVmag2_dVmag2 = -1

				_Vmag2 = Vmag2_hist - dVmag2_dVr * Vr - dVmag2_dVi * Vi - dVmag2_dVmag2 * Vmag2

				# # Stamp Q3 Vb partials # #
				idx_Y = self.stamp_Y(node_Q3, nodeB_Vr, dFQ3.dVr[1], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y = self.stamp_Y(node_Q3, nodeB_Vi, dFQ3.dVi[1], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

			if self.phases & 0x04:
				# # Stamp Voltage C Magnitude # #
				Vr = Vcr
				Vi = Vci
				Vmag2 = Vabc.Vmag2[2]

				Vmag2_hist = Vr ** 2 + Vi ** 2 - Vmag2

				dVmag2_dVr = 2 * Vr
				dVmag2_dVi = 2 * Vi
				dVmag2_dVmag2 = -1

				_Vmag2 = Vmag2_hist - dVmag2_dVr * Vr - dVmag2_dVi * Vi - dVmag2_dVmag2 * Vmag2

				# # Stamp Q3 Vc partials # #
				idx_Y = self.stamp_Y(node_Q3, nodeC_Vr, dFQ3.dVr[2], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
				idx_Y = self.stamp_Y(node_Q3, nodeC_Vi, dFQ3.dVi[2], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		idx_Y = self.stamp_Y(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		FQ3 = FQ3_hist - dFQ3.dVr[0] * Var - dFQ3.dVi[0] * Vai - dFQ3.dVr[1] * Vbr - \
		      dFQ3.dVi[1] * Vbi - dFQ3.dVr[2] * Vcr - dFQ3.dVi[2] * Vci - dFQ3.dQ * Q3

		idx_J = self.stamp_J(node_Q3, -FQ3, Jnlin_val, Jnlin_row, idx_J)

		return idx_Y, idx_J

	def stamp_ibdg_partials(self,
	                        Vseq,
	                        phase,
	                        dIdgrPos,
	                        dIdgiPos,
	                        dIdgrNeg,
	                        dIdgiNeg,
	                        IdgrPos_hist,
	                        IdgiPos_hist,
	                        IdgrNeg_hist,
	                        IdgiNeg_hist,
	                        Ynlin_val,
	                        Ynlin_row,
	                        Ynlin_col,
	                        Jnlin_val,
	                        Jnlin_row,
	                        idx_Y,
	                        idx_J):

		Vposr = Vseq["VR"][0][phase]
		Vnegr = Vseq["VR"][1][phase]
		Vposi = Vseq["VI"][0][phase]
		Vnegi = Vseq["VI"][1][phase]
		Q3 = Vseq["Q3"][0]

		nodePos_Vr = Vseq["nodeVR"][phase]
		nodeNeg_Vr = Vseq["nodeVR"][phase]
		nodePos_Vi = Vseq["nodeVI"][phase]
		nodeNeg_Vi = Vseq["nodeVI"][phase]
		node_Q3 = Vseq["nodeQ3"][0]

		# Idgr Pos Stamps
		idx_Y = self.stamp_Y(nodePos_Vr, nodePos_Vr, dIdgrPos.dVposr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vr, nodePos_Vi, dIdgrPos.dVposi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vr, nodeNeg_Vr, dIdgrPos.dVnegr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vr, nodeNeg_Vi, dIdgrPos.dVnegi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Idgi Pos Stamps
		idx_Y = self.stamp_Y(nodePos_Vi, nodePos_Vr, dIdgiPos.dVposr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vi, nodePos_Vi, dIdgiPos.dVposi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vi, nodeNeg_Vr, dIdgiPos.dVnegr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodePos_Vi, nodeNeg_Vi, dIdgiPos.dVnegi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Idgr Neg Stamps
		idx_Y = self.stamp_Y(nodeNeg_Vr, nodePos_Vr, dIdgrNeg.dVposr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vr, nodePos_Vi, dIdgrNeg.dVposi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vr, nodeNeg_Vr, dIdgrNeg.dVnegr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vr, nodeNeg_Vi, dIdgrNeg.dVnegi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Idgi Pos Stamps
		idx_Y = self.stamp_Y(nodeNeg_Vi, nodePos_Vr, dIdgiNeg.dVposr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vi, nodePos_Vi, dIdgiNeg.dVposi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vi, nodeNeg_Vr, dIdgiNeg.dVnegr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = self.stamp_Y(nodeNeg_Vi, nodeNeg_Vi, dIdgiNeg.dVnegi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		#  Q Stamps
		if not self.fixed_Q:
			idx_Y = self.stamp_Y(nodePos_Vr, node_Q3, dIdgrPos.dQ3, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(nodePos_Vi, node_Q3, dIdgiPos.dQ3, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(nodeNeg_Vr, node_Q3, dIdgrNeg.dQ3, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(nodeNeg_Vi, node_Q3, dIdgiNeg.dQ3, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		else:
			dIdgrPos.dQ3 = 0
			dIdgiPos.dQ3 = 0
			dIdgrNeg.dQ3 = 0
			dIdgiNeg.dQ3 = 0

		IdgrPos = IdgrPos_hist - dIdgrPos.dVposr * Vposr - dIdgrPos.dVposi * Vposi - dIdgrPos.dVnegr * Vnegr - \
		          dIdgrPos.dVnegi * Vnegi - dIdgrPos.dQ3 * Q3

		IdgiPos = IdgiPos_hist - dIdgiPos.dVposr * Vposr - dIdgiPos.dVposi * Vposi - dIdgiPos.dVnegr * Vnegr - \
		          dIdgiPos.dVnegi * Vnegi - dIdgiPos.dQ3 * Q3

		IdgrNeg = IdgrNeg_hist - dIdgrNeg.dVposr * Vposr - dIdgrNeg.dVposi * Vposi - dIdgrNeg.dVnegr * Vnegr \
		          - dIdgrNeg.dVnegi * Vnegi - dIdgrNeg.dQ3 * Q3

		IdgiNeg = IdgiNeg_hist - dIdgiNeg.dVposr * Vposr - dIdgiNeg.dVposi * Vposi - dIdgiNeg.dVnegr * Vnegr - \
		          dIdgiNeg.dVnegi * Vnegi - dIdgiNeg.dQ3 * Q3

		idx_J = self.stamp_J(nodePos_Vr, -IdgrPos, Jnlin_val, Jnlin_row, idx_J)
		idx_J = self.stamp_J(nodePos_Vi, -IdgiPos, Jnlin_val, Jnlin_row, idx_J)
		idx_J = self.stamp_J(nodeNeg_Vr, -IdgrNeg, Jnlin_val, Jnlin_row, idx_J)
		idx_J = self.stamp_J(nodeNeg_Vi, -IdgiNeg, Jnlin_val, Jnlin_row, idx_J)

		return idx_Y, idx_J


class PNSC(IBDG):

	def __init__(self,
	             ID,
	             name,
	             phases,
	             nominal_voltage,
	             voltage_setpoint,
	             connection_type,
	             rated_power,
	             P,
	             Q,
	             active_rating,
	             reactive_rating,
	             power_factor,
	             irradiance,
	             Pmpp,
	             inverter_settings,
	             is_triplex,
	             sensor="GVCC",
	             fixed_Q=True):

		super(PNSC, self).__init__(ID, name, phases,
		                           nominal_voltage, voltage_setpoint,
		                           connection_type, rated_power, P,
		                           Q, active_rating, reactive_rating, power_factor, irradiance,
		                           Pmpp, inverter_settings, is_triplex, sensor, fixed_Q)

	@staticmethod
	def calc_hist_currents(Vpos,
	                       Vneg,
	                       Q3,
	                       P3,
	                       alphar,
	                       alphai,
	                       betar,
	                       betai):

		Idgr_Pos = np.zeros((3, 1))
		Idgi_Pos = np.zeros((3, 1))
		Idgr_Neg = np.zeros((3, 1))
		Idgi_Neg = np.zeros((3, 1))

		for _phase in range(3):
			Vposr = Vpos.real[_phase]
			Vposi = Vpos.imag[_phase]
			Vnegr = Vneg.real[_phase]
			Vnegi = Vneg.imag[_phase]
			vmag = Vposr ** 2 + Vposi ** 2 - Vnegr ** 2 - Vnegi ** 2

			if vmag != 0:
				Idgr_Pos[_phase] = ((Q3 * (Vposi * alphar + alphai * Vposr) + P3 * (
						-alphai * Vposi + Vposr * alphar)) / vmag) \
				                   - Vposi * betai + Vposr * betar
				Idgi_Pos[_phase] = ((P3 * (Vposi * alphar + alphai * Vposr) + Q3 * (
						alphai * Vposi - Vposr * alphar)) / vmag) \
				                   + Vposi * betar + Vposr * betai if vmag != 0 else 0
				Idgr_Neg[_phase] = ((-Q3 * (Vnegi * alphar + alphai * Vnegr) + P3 * (
						alphai * Vnegi - Vnegr * alphar)) / vmag) \
				                   - Vnegi * betai + Vnegr * betar
				Idgi_Neg[_phase] = ((-P3 * (Vnegi * alphar + alphai * Vnegr) + Q3 * (
						-alphai * Vnegi + Vnegr * alphar)) / vmag) \
				                   + Vnegi * betar + Vnegr * betai

		return Idgr_Pos, Idgi_Pos, Idgr_Neg, Idgi_Neg

	@staticmethod
	def calc_partials(Vseq,
	                  phase,
	                  P3,
	                  Q3,
	                  alphar,
	                  alphai,
	                  betar,
	                  betai):

		Vposr = Vseq["VR"][0][phase]
		Vnegr = Vseq["VR"][1][phase]
		Vposi = Vseq["VI"][0][phase]
		Vnegi = Vseq["VI"][1][phase]

		vmag = Vposr ** 2 + Vposi ** 2 - Vnegr ** 2 - Vnegi ** 2
		vmag2 = vmag ** 2

		# Real Partial Derivatives of Idg Pos
		dIdgrPos_dVposr = Q3 * (2 * (-alphar * Vposi * Vposr - alphai * Vposr ** 2) / vmag2 + (alphai / vmag)) \
		                  + P3 * (2 * (alphai * Vposi * Vposr - alphar * Vposr ** 2) / vmag2 + (alphar / vmag)) - betar
		dIdgrPos_dVposi = Q3 * (2 * (-alphai * Vposi * Vposr - alphar * Vposi ** 2) / vmag2 + (alphar / vmag)) \
		                  + P3 * (2 * (-alphar * Vposi * Vposr + alphai * Vposi ** 2) / vmag2 - (alphai / vmag)) + \
		                  betai
		dIdgrPos_dVnegr = (
				                  2 * Q3 *
				                  (alphar * Vnegr * Vposi + alphai * Vnegr * Vposr) + 2 * P3 *
				                  (-alphai * Vnegr * Vposi + alphar * Vnegr * Vposr)) / vmag2
		dIdgrPos_dVnegi = (
				                  2 * Q3 *
				                  (alphar * Vnegi * Vposi + alphai * Vnegi * Vposr) + 2 * P3 *
				                  (-alphai * Vnegi * Vposi + alphar * Vnegi * Vposr)) / vmag2
		dIdgrPos_dQ3 = (Vposi * alphar + alphai * Vposr) / vmag

		# Imaginary Partial Derivatives of Idg Pos
		dIdgiPos_dVposr = P3 * (2 * (-alphar * Vposi * Vposr - alphai * Vposr ** 2) / vmag2 + (alphai / vmag)) \
		                  + Q3 * (2 * (-alphai * Vposi * Vposr + alphar * Vposr ** 2) / vmag2 - (alphar / vmag)) - \
		                  betai
		dIdgiPos_dVposi = P3 * (2 * (-alphai * Vposi * Vposr - alphar * Vposi ** 2) / vmag2 + (alphar / vmag)) \
		                  + Q3 * (2 * (alphar * Vposi * Vposr - alphai * Vposi ** 2) / vmag2 + (alphai / vmag)) - betar
		dIdgiPos_dVnegr = (
				                  2 * P3 *
				                  (alphar * Vnegr * Vposi + alphai * Vnegr * Vposr) + 2 * Q3 *
				                  (alphai * Vnegr * Vposi - alphar * Vnegr * Vposr)) / vmag2
		dIdgiPos_dVnegi = (
				                  2 * P3 *
				                  (alphar * Vnegi * Vposi + alphai * Vnegi * Vposr) + 2 * Q3 *
				                  (alphai * Vnegi * Vposi - alphar * Vnegi * Vposr)) / vmag2
		dIdgiPos_dQ3 = (alphai * Vposi - Vposr * alphar) / vmag

		# Real Partial Derivatives of Idg Neg
		dIdgrNeg_dVposr = (
				                  2 * Q3 *
				                  (alphar * Vnegi * Vposr + alphai * Vnegr * Vposr) + 2 * P3 *
				                  (-alphai * Vnegi * Vposr + alphar * Vnegr * Vposr)) / vmag2
		dIdgrNeg_dVposi = (
				                  2 * Q3 *
				                  (alphar * Vnegi * Vposi + alphai * Vnegr * Vposi) + 2 * P3 *
				                  (-alphai * Vnegi * Vposi + alphar * Vnegr * Vposi)) / vmag2
		dIdgrNeg_dVnegr = Q3 * (2 * (-alphar * Vnegi * Vnegr - alphai * Vnegr ** 2) / vmag2 - (alphai / vmag)) \
		                  + P3 * (2 * (alphai * Vnegi * Vnegr - alphar * Vnegr ** 2) / vmag2 - (alphar / vmag)) - betar
		dIdgrNeg_dVnegi = Q3 * (2 * (-alphai * Vnegi * Vnegr - alphar * Vnegi ** 2) / vmag2 - (alphar / vmag)) \
		                  + P3 * (2 * (-alphar * Vnegi * Vnegr + alphai * Vnegi ** 2) / vmag2 + (alphai / vmag)) + \
		                  betai
		dIdgrNeg_dQ3 = -(Vnegi * alphar + alphai * Vnegr) / vmag

		# Imaginary Partial Derivatives of Idg Neg
		dIdgiNeg_dVposr = (
				                  2 * P3 *
				                  (alphar * Vnegi * Vposr + alphai * Vnegr * Vposr) + 2 * Q3 *
				                  (alphai * Vnegi * Vposr - alphar * Vnegr * Vposr)) / vmag2
		dIdgiNeg_dVposi = (
				                  2 * P3 *
				                  (alphar * Vnegi * Vposi + alphai * Vnegr * Vposi) + 2 * Q3 *
				                  (alphai * Vnegi * Vposi - alphar * Vnegr * Vposi)) / vmag2
		dIdgiNeg_dVnegr = P3 * (2 * (-alphar * Vnegi * Vnegr - alphai * Vnegr ** 2) / vmag2 - (alphai / vmag)) \
		                  + Q3 * (2 * (-alphai * Vnegi * Vnegr + alphar * Vnegr ** 2) / vmag2 + (alphar / vmag)) - \
		                  betai
		dIdgiNeg_dVnegi = P3 * (2 * (-alphai * Vnegi * Vnegr - alphar * Vnegi ** 2) / vmag2 - (alphar / vmag)) \
		                  + Q3 * (2 * (alphar * Vnegi * Vnegr - alphai * Vnegi ** 2) / vmag2 - (alphai / vmag)) - betar
		dIdgiNeg_dQ3 = (-alphai * Vnegi + Vnegr * alphar) / vmag

		# Add partial derivatives and currents at iteration k to dictionaries
		dIdgrPos = {
			'dVposr': dIdgrPos_dVposr,
			'dVposi': dIdgrPos_dVposi,
			'dVnegr': dIdgrPos_dVnegr,
			'dVnegi': dIdgrPos_dVnegi,
			'dQ3': dIdgrPos_dQ3
		}

		dIdgiPos = {
			'dVposr': dIdgiPos_dVposr,
			'dVposi': dIdgiPos_dVposi,
			'dVnegr': dIdgiPos_dVnegr,
			'dVnegi': dIdgiPos_dVnegi,
			'dQ3': dIdgiPos_dQ3
		}

		dIdgrNeg = {
			'dVposr': dIdgrNeg_dVposr,
			'dVposi': dIdgrNeg_dVposi,
			'dVnegr': dIdgrNeg_dVnegr,
			'dVnegi': dIdgrNeg_dVnegi,
			'dQ3': dIdgrNeg_dQ3
		}

		dIdgiNeg = {
			'dVposr': dIdgiNeg_dVposr,
			'dVposi': dIdgiNeg_dVposi,
			'dVnegr': dIdgiNeg_dVnegr,
			'dVnegi': dIdgiNeg_dVnegi,
			'dQ3': dIdgiNeg_dQ3
		}

		dIdgrPos = SimpleNamespace(**dIdgrPos)
		dIdgiPos = SimpleNamespace(**dIdgiPos)
		dIdgrNeg = SimpleNamespace(**dIdgrNeg)
		dIdgiNeg = SimpleNamespace(**dIdgiNeg)

		return dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg

	def calc_ibdg_outputs(self,
	                      node,
						  node_key,
	                      Vsol,
	                      flag_print):

		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		Vabc = self.get_nodes_states(node, node_key, Vsol)

		Var = Vabc.VR[0]
		Vbr = Vabc.VR[1]
		Vcr = Vabc.VR[2]

		Vai = Vabc.VI[0]
		Vbi = Vabc.VI[1]
		Vci = Vabc.VI[2]

		Va = (Var + Vai * 1j)
		Vb = (Vbr + Vbi * 1j)
		Vc = (Vcr + Vci * 1j)

		Vabc_mag = np.array([Var, Vai, Vbr, Vbi, Vcr, Vci], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		IdgrPos = np.zeros((3, 1))
		IdgiPos = np.zeros((3, 1))
		IdgrNeg = np.zeros((3, 1))
		IdgiNeg = np.zeros((3, 1))

		Q3 = Vabc.Q[0]
		self.Q3_final = Q3
		P3 = self.P3

		for v in range(3):
			vmag = Vposr[v] ** 2 + Vposi[v] ** 2 - Vnegr[v] ** 2 - Vnegi[v] ** 2

			IdgrPos[v] = ((Q3 * (Vposi[v] * alphar + alphai * Vposr[v]) + P3 *
			               (Vposr[v] * alphar - alphai * Vposi[v])) /
			              vmag) - Vposi[v] * betai + Vposr[v] * betar

			IdgiPos[v] = ((Q3 * (Vposi[v] * alphai - Vposr[v] * alphar) + P3 *
			               (Vposi[v] * alphar + alphai * Vposr[v])) /
			              vmag) + Vposi[v] * betar + Vposr[v] * betai

			IdgrNeg[v] = ((-Q3 * (Vnegi[v] * alphar + alphai * Vnegr[v]) + P3 *
			               (alphai * Vnegi[v] - Vnegr[v] * alphar)) /
			              vmag) - Vnegi[v] * betai + Vnegr[v] * betar

			IdgiNeg[v] = ((-P3 * (Vnegi[v] * alphar + alphai * Vnegr[v]) + Q3 *
			               (-alphai * Vnegi[v] + Vnegr[v] * alphar)) /
			              vmag) + Vnegi[v] * betar + Vnegr[v] * betai

		Iar = IdgrPos[0] + IdgrNeg[0]
		Iai = IdgiPos[0] + IdgiNeg[0]
		Ibr = IdgrPos[1] + IdgrNeg[1]
		Ibi = IdgiPos[1] + IdgiNeg[1]
		Icr = IdgrPos[2] + IdgrNeg[2]
		Ici = IdgiPos[2] + IdgiNeg[2]

		Ia = (Iar + 1j * Iai) / 3
		Ib = (Ibr + 1j * Ibi) / 3
		Ic = (Icr + 1j * Ici) / 3

		self.Ia_mag = float(np.abs(Ia))
		self.Ia_ang = float(np.angle(Ia, True))

		self.Ib_mag = float(np.abs(Ib))
		self.Ib_ang = float(np.angle(Ib, True))

		self.Ic_mag = float(np.abs(Ic))
		self.Ic_ang = float(np.angle(Ic, True))

		# self.perc_Q = abs(self.Q3_final[0] / self.Qscale)

		S3 = Va * np.conj(Ia) + Vb * np.conj(Ib) + Vc * np.conj(Ic)
		if flag_print:
			logging.info("---------------------------------------------------------")
			logging.info(self.name)
			logging.info("Currents are: Ia = %f < %f | Ib = %f < %f | Ic = %f < %f" %
			             (self.Ia_mag, self.Ia_ang, self.Ib_mag, self.Ib_ang,
			              self.Ic_mag, self.Ic_ang))
			logging.info("Three-Phase Powers are: P = %f kW | Q = %f kVAR" %
			             (-S3.real * 1e-3, S3.imag * 1e-3))

	def calc_max_current(self,
	                     Vposr,
	                     Vnegr,
	                     Vposi,
	                     Vnegi,
	                     Q3,
	                     lf):
		Qcurr = Q3
		Pref = -self.P3_max * lf

		Vpos = Vposr + 1j * Vposi
		Vneg = Vnegr + 1j * Vnegi

		Vpos_mag = np.abs(Vpos)
		Vneg_mag = np.abs(Vneg)

		Ipeak = np.zeros((3, 1))
		delta = (np.angle(Vpos) ** 2 - np.angle(Vneg) ** 2) / 2
		shift = np.array([0, np.pi / 3, -np.pi / 3])
		gamma = np.array(
			([delta[0] + shift[0], delta[1] + shift[1], delta[2] + shift[2]]))

		for v in range(3):
			A1 = ((Vpos_mag[v] - Vneg_mag[v]) /
			      (Vpos_mag[v] ** 2 - Vneg_mag[v] ** 2))
			A2 = ((Vpos_mag[v] + Vneg_mag[v]) /
			      (Vpos_mag[v] ** 2 - Vneg_mag[v] ** 2))

			IpL = Pref * A1
			IpS = Pref * A2
			IqL = Qcurr * A1
			IqS = Qcurr * A2

			C1 = IpL * np.cos(gamma[v]) - IqL * np.sin(gamma[v])
			C2 = -IqS * np.cos(gamma[v]) - IpS * np.sin(gamma[v])

			Ipeak[v] = np.sqrt(C1 ** 2 + C2 ** 2)

		return Ipeak

	def calc_Q_limits(self,
	                  Vpos,
	                  Vneg,
	                  Imax):
		Qalt = np.zeros((3, 1))

		Vpos_mag = np.abs(Vpos) / np.sqrt(3)
		Vneg_mag = np.abs(Vneg) / np.sqrt(3)
		delta = (np.angle(Vpos) ** 2 - np.angle(Vneg) ** 2) / 2
		shift = np.array([0, np.pi / 3, -np.pi / 3])
		gamma = np.array(
			([delta[0] + shift[0], delta[1] + shift[1], delta[2] + shift[2]]))

		Imax = Imax.reshape((3, 1))

		alt1 = ((Imax ** 2) * (((Vneg_mag - Vpos_mag) ** 2) *
		                       ((Vneg_mag + Vpos_mag) ** 2)))
		alt2 = (Vneg_mag ** 2 + Vpos_mag ** 2 +
		        2 * Vpos_mag * Vneg_mag * np.cos(2 * gamma))
		for phase in range(3):
			Qalt[phase] = np.sqrt(alt1.item((phase, 0)) / alt2.item((phase, 0)))
		Q3 = np.min(Qalt)
		Q3_loc = np.argmin(Qalt)

		if Q3_loc == 0:
			self.Qphase = 'A'
		elif Q3_loc == 1:
			self.Qphase = 'B'
		else:
			self.Qphase = 'C'

		if Q3 < 0:
			self.Qmin = -Q3
			self.Qmax = Q3
		elif Q3 > 0:
			self.Qmax = -Q3
			self.Qmin = Q3

	def stamp_nonlinear(self, node_key, node, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    Jnlin_val, Jnlin_row, idx_Y, idx_J, homotopy_enabled=False, h_factor=None):

		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		Vabc = self.get_nodes_states(node, node_key, V)

		Vabc_mag = np.array([
			Vabc.VR[0], Vabc.VI[0], Vabc.VR[1], Vabc.VI[1], Vabc.VR[2],
			Vabc.VI[2]
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		Vseq = {
			'VR': [Vposr, Vnegr],
			'VI': [Vposi, Vnegi],
			'Q3': Vabc.Q,
			'nodeVR': Vabc.nodeVR,
			'nodeVI': Vabc.nodeVI,
			'nodeQ3': Vabc.nodeQ3
		}

		Vpos = Vposr + 1j * Vposi
		Vneg = Vnegr + 1j * Vnegi
		Q3_ref = Vabc.Q[0]

		P3_ref = self.P3

		# # real and imaginary stamps # #

		IdgrPos_hist, IdgiPos_hist, IdgrNeg_hist, IdgiNeg_hist = self.calc_hist_currents(
			Vpos, Vneg, Q3_ref, P3_ref, alphar, alphai, betar, betai)
		self.IdgrPos = IdgrPos_hist
		self.IdgiPos = IdgiPos_hist
		self.IdgrNeg = IdgrNeg_hist
		self.IdgiNeg = IdgiNeg_hist

		for _phase_index in range(3):
			if self.ibdg_phases[_phase_index] & self.phases == int(
					self.ibdg_phases[_phase_index]):
				IdgrPos_phase = self.IdgrPos[_phase_index]
				IdgiPos_phase = self.IdgiPos[_phase_index]
				IdgrNeg_phase = self.IdgrNeg[_phase_index]
				IdgiNeg_phase = self.IdgiNeg[_phase_index]

				node_Vr_from = Vabc.nodeVR[_phase_index]
				node_Vi_from = Vabc.nodeVI[_phase_index]

				dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg = self.calc_partials(
					Vseq, _phase_index, P3_ref, Q3_ref, alphar, alphai, betar,
					betai)

				idx_Y, idx_J = super(PNSC, self).stamp_ibdg_partials(
					Vseq, _phase_index, dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg,
					IdgrPos_phase, IdgiPos_phase, IdgrNeg_phase, IdgiNeg_phase,
					Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
					idx_Y, idx_J)

		if not self.fixed_Q:
			dFQ3, FQ3_hist = self.calc_Q3(Vabc)

			idx_Y, idx_J = super(PNSC,
			                     self).stamp_Q3_partials(Vabc, dFQ3, FQ3_hist,
			                                             Ynlin_val, Ynlin_row,
			                                             Ynlin_col, Jnlin_val,
			                                             Jnlin_row, idx_Y, idx_J, homotopy_enabled, h_factor)

		return idx_Y, idx_J


class BPSC(IBDG):

	def __init__(self,
	             ID,
	             name,
	             phases,
	             nominal_voltage,
	             voltage_setpoint,
	             connection_type,
	             rated_power,
	             P,
	             Q,
	             active_rating,
	             reactive_rating,
	             power_factor,
	             irradiance,
	             Pmpp,
	             inverter_settings,
	             is_triplex,
	             sensor="GVCC",
	             fixed_Q=True):

		super(BPSC, self).__init__(ID, name, phases,
		                           nominal_voltage, voltage_setpoint,
		                           connection_type, rated_power, P,
		                           Q, active_rating, reactive_rating, power_factor, irradiance,
		                           Pmpp, inverter_settings, is_triplex, sensor, fixed_Q)

	@staticmethod
	def calc_hist_currents(Vpos, Vneg, Q3, P3, alphar, alphai, betar, betai):
		Vposr = Vpos.real
		Vposi = Vpos.imag
		Vnegr = Vneg.real
		Vnegi = Vneg.imag

		vmag = Vposr ** 2 + Vposi ** 2

		Idgr_Pos = ((Q3 * (Vposi * alphar + alphai * Vposr) + P3 * (Vposr * alphar - alphai * Vposi)) / vmag) \
		           + Vposi * betai - betar * Vposr
		Idgi_Pos = ((P3 * (Vposi * alphar + alphai * Vposr) + Q3 * (-alphar * Vposr + alphai * Vposi)) / vmag) \
		           - Vposr * betai - betar * Vposi
		Idgr_Neg = Vnegr * betar - betai * Vnegi
		Idgi_Neg = betar * Vnegi + betai * Vnegr

		return Idgr_Pos, Idgi_Pos, Idgr_Neg, Idgi_Neg

	@staticmethod
	def calc_partials(Vseq, phase, P3, Q3, alphar, alphai, betar, betai):
		Vposr = Vseq["VR"][0][phase]
		Vposi = Vseq["VI"][0][phase]

		vmag = Vposr ** 2 + Vposi ** 2
		vmag2 = vmag ** 2

		# Real Partial Derivatives of Idg Pos
		dIdgrPos_dVposr = Q3 * (-2 * (alphar * Vposi * Vposr + alphai * Vposr ** 2) / vmag2 + (alphai / vmag)) \
		                  + P3 * (2 * ((-alphar * Vposr ** 2 + alphai * Vposi * Vposr) / vmag2) + (
				alphar / vmag)) - betar
		dIdgrPos_dVposi = P3 * (2 * (-alphar * Vposi * Vposr + alphai * Vposi ** 2) / vmag2 - (alphai / vmag)) \
		                  + Q3 * (2 * ((-alphar * Vposi ** 2 - alphai * Vposi * Vposr) / vmag2) + (
				alphar / vmag)) + betai
		dIdgrPos_dQ3 = (Vposi * alphar + alphai * Vposr) / vmag

		# Imaginary Partial Derivatives of Idg Pos
		dIdgiPos_dVposr = P3 * (2 * (-alphar * Vposi * Vposr - alphai * Vposr ** 2) / vmag2 + (alphai / vmag)) \
		                  + Q3 * (2 * ((alphar * Vposr ** 2 - alphai * Vposi * Vposr) / vmag2) - (
				alphar / vmag)) - betai
		dIdgiPos_dVposi = Q3 * (2 * (alphar * Vposi * Vposr - alphai * Vposi ** 2) / vmag2 + (alphai / vmag)) \
		                  + P3 * (2 * ((-alphar * Vposi ** 2 - alphai * Vposi * Vposr) / vmag2) + (
				alphar / vmag)) - betar
		dIdgiPos_dQ3 = (-alphar * Vposr + alphai * Vposi) / vmag

		# Real Partial Derivatives of Idg Neg
		dIdgrNeg_dVnegi = -betai
		dIdgrNeg_dVnegr = betar

		# Imaginary Partial Derivatives of Idg Neg
		dIdgiNeg_dVnegr = betai
		dIdgiNeg_dVnegi = betar

		# Add partial derivatives and currents at iteration k to dictionaries
		dIdgrPos = {
			'dVposr': dIdgrPos_dVposr,
			'dVposi': dIdgrPos_dVposi,
			'dVnegr': 0,
			'dVnegi': 0,
			'dQ3': dIdgrPos_dQ3
		}

		dIdgiPos = {
			'dVposr': dIdgiPos_dVposr,
			'dVposi': dIdgiPos_dVposi,
			'dVnegr': 0,
			'dVnegi': 0,
			'dQ3': dIdgiPos_dQ3
		}

		dIdgrNeg = {
			'dVposr': 0,
			'dVposi': 0,
			'dVnegr': dIdgrNeg_dVnegr,
			'dVnegi': dIdgrNeg_dVnegi,
			'dQ3': 0
		}

		dIdgiNeg = {
			'dVposr': 0,
			'dVposi': 0,
			'dVnegr': dIdgiNeg_dVnegr,
			'dVnegi': dIdgiNeg_dVnegi,
			'dQ3': 0
		}

		dIdgrPos = SimpleNamespace(**dIdgrPos)
		dIdgiPos = SimpleNamespace(**dIdgiPos)
		dIdgrNeg = SimpleNamespace(**dIdgrNeg)
		dIdgiNeg = SimpleNamespace(**dIdgiNeg)

		return dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg

	def stamp_nonlinear(self, node_key, node, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    Jnlin_val, Jnlin_row, idx_Y, idx_J, homotopy_enabled=False, h_factor=None):
		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		P3_ref = self.P3

		Vabc = self.get_nodes_states(node, node_key, V)

		Vabc_mag = np.array([
			Vabc.VR[0], Vabc.VI[0], Vabc.VR[1], Vabc.VI[1], Vabc.VR[2],
			Vabc.VI[2]
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		Vseq = {
			'VR': [Vposr, Vnegr],
			'VI': [Vposi, Vnegi],
			'Q3': Vabc.Q,
			'nodeVR': Vabc.nodeVR,
			'nodeVI': Vabc.nodeVI,
			'nodeQ3': Vabc.nodeQ3
		}
		''' real and imaginary stamps '''

		Vpos = Vposr + 1j * Vposi
		Vneg = Vnegr + 1j * Vnegi
		Q3_ref = Vabc.Q[0]
		# Q3_ref *= lf

		IdgrPos_hist, IdgiPos_hist, IdgrNeg_hist, IdgiNeg_hist = self.calc_hist_currents(
			Vpos, Vneg, Q3_ref, P3_ref, alphar, alphai, betar, betai)

		self.IdgrPos = IdgrPos_hist
		self.IdgiPos = IdgiPos_hist
		self.IdgrNeg = IdgrNeg_hist
		self.IdgiNeg = IdgiNeg_hist

		for _phase_index in range(3):
			dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg = self.calc_partials(
				Vseq, _phase_index, P3_ref, Q3_ref, alphar, alphai, betar,
				betai)

			IdgrPos_phase = self.IdgrPos[_phase_index]
			IdgiPos_phase = self.IdgiPos[_phase_index]
			IdgrNeg_phase = self.IdgrNeg[_phase_index]
			IdgiNeg_phase = self.IdgiNeg[_phase_index]

			idx_Y, idx_J = super(BPSC, self).stamp_ibdg_partials(
				Vseq, _phase_index, dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg,
				IdgrPos_phase, IdgiPos_phase, IdgrNeg_phase, IdgiNeg_phase,
				Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y,
				idx_J)
		if not self.fixed_Q:
			dFQ3, FQ3_hist = self.calc_Q3(Vabc)

			idx_Y, idx_J = self.stamp_Q3_partials(Vabc, dFQ3, FQ3_hist, Ynlin_val,
			                                      Ynlin_row, Ynlin_col, Jnlin_val,
			                                      Jnlin_row, idx_Y, idx_J, homotopy_enabled, h_factor)

		return idx_Y, idx_J

	def calc_ibdg_outputs(self, node, node_key, Vsol, flag_print):

		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		Vabc = self.get_nodes_states(node, node_key, Vsol)

		Va = Vabc.VR[0] + Vabc.VI[0] * 1j
		Vb = Vabc.VR[1] + Vabc.VI[1] * 1j
		Vc = Vabc.VR[2] + Vabc.VI[2] * 1j

		Vabc_mag = np.array([
			Vabc.VR[0], Vabc.VI[0], Vabc.VR[1], Vabc.VI[1], Vabc.VR[2],
			Vabc.VI[2]
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		IdgrPos = np.zeros((3, 1))
		IdgiPos = np.zeros((3, 1))
		IdgrNeg = np.zeros((3, 1))
		IdgiNeg = np.zeros((3, 1))

		Q3 = Vabc.Q[0]
		self.Q3_final = Q3
		P3 = self.P3

		for v in range(3):
			vmag = Vposr[v] ** 2 + Vposi[v] ** 2

			IdgrPos[v] = ((Q3 * (Vposi[v] * alphar + alphai * Vposr[v]) + P3 * (
					Vposr[v] * alphar - alphai * Vposi[v])) / vmag) \
			             + Vposi[v] * betai - betar * Vposr[v]
			IdgiPos[v] = ((P3 * (Vposi[v] * alphar + alphai * Vposr[v]) + Q3 * (
					-alphar * Vposr[v] + alphai * Vposi[v])) / vmag) \
			             - Vposr[v] * betai - betar * Vposi[v]

			IdgrNeg[v] = Vnegr[v] * betar - betai * Vnegi[v]
			IdgiNeg[v] = betar * Vnegi[v] + betai * Vnegr[v]

		Iar = IdgrPos[0] + IdgrNeg[0]
		Iai = IdgiPos[0] + IdgiNeg[0]
		Ibr = IdgrPos[1] + IdgrNeg[1]
		Ibi = IdgiPos[1] + IdgiNeg[1]
		Icr = IdgrPos[2] + IdgrNeg[2]
		Ici = IdgiPos[2] + IdgiNeg[2]

		Ia = (Iar + 1j * Iai) / 3
		Ib = (Ibr + 1j * Ibi) / 3
		Ic = (Icr + 1j * Ici) / 3

		self.Ia_mag = float(np.abs(Ia))
		self.Ia_ang = float(np.angle(Ia, True))

		self.Ib_mag = float(np.abs(Ib))
		self.Ib_ang = float(np.angle(Ib, True))

		self.Ic_mag = float(np.abs(Ic))
		self.Ic_ang = float(np.angle(Ic, True))

		S3 = Va * np.conj(Ia) + Vb * np.conj(Ib) + Vc * np.conj(Ic)
		if flag_print:
			logging.info("---------------------------------------------------------")
			logging.info(self.name)
			logging.info("Currents are: Ia = %f < %f | Ib = %f < %f | Ic = %f < %f" %
			             (self.Ia_mag, self.Ia_ang, self.Ib_mag, self.Ib_ang, self.Ic_mag,
			              self.Ic_ang))
			logging.info("Three-Phase Powers are: P = %f kW | Q = %f kVAR" %
			             (-S3.real * 1e-3, S3.imag * 1e-3))


class FPNSC(IBDG):

	def __init__(self,
	             ID,
	             name,
	             phases,
	             nominal_voltage,
	             voltage_setpoint,
	             connection_type,
	             rated_power,
	             P,
	             Q,
	             active_rating,
	             reactive_rating,
	             power_factor,
	             irradiance,
	             Pmpp,
	             inverter_settings,
	             is_triplex,
	             sensor="GVCC",
	             fixed_Q=True):

		super(FPNSC, self).__init__(ID, name, phases,
		                            nominal_voltage, voltage_setpoint,
		                            connection_type, rated_power, P,
		                            Q, active_rating, reactive_rating, power_factor, irradiance,
		                            Pmpp, inverter_settings, is_triplex, sensor, fixed_Q)
		self.k1 = np.ones((3, 1))
		self.k2 = np.ones((3, 1))
		self.Vpos = np.ones((3, 1))
		self.Vneg = np.ones((3, 1))
		self.Vzero = np.ones((3, 1))

	def determine_k1_k2(self):
		Vpos2 = np.linalg.norm(self.Vpos) ** 2
		Vneg2 = np.linalg.norm(self.Vneg) ** 2

		np.seterr(all='warn')
		with warnings.catch_warnings():
			warnings.filterwarnings("error", category=RuntimeWarning)
			try:
				if np.isclose(Vpos2, Vneg2):
					k1 = Vpos2 / 0.0
				else:
					k1 = Vpos2 / (Vpos2 - Vneg2)
				k2 = Vpos2 / (Vpos2 + Vneg2)
				self.k1 = np.array([[k1], [k1], [k1]])
				self.k2 = np.array([[k2], [k2], [k2]])
			except RuntimeWarning:
				self.k1 = np.ones((3, 1))
				self.k2 = np.ones((3, 1))

	@staticmethod
	def calc_hist_currents(Vpos, Vneg, Q3, P3, alphar, alphai, betar, betai, k1,
	                       k2):
		Vposr = Vpos.real.reshape((3, 1))
		Vposi = Vpos.imag.reshape((3, 1))
		Vnegr = Vneg.real.reshape((3, 1))
		Vnegi = Vneg.imag.reshape((3, 1))

		vmag_pos = Vposr ** 2 + Vposi ** 2
		vmag_neg = Vnegr ** 2 + Vnegi ** 2

		Idgr_Pos = ((k2 * Q3 * (Vposi * alphar + alphai * Vposr) + k1 * P3 * (
				Vposr * alphar - alphai * Vposi)) / vmag_pos) \
		           - Vposi * betai + betar * Vposr
		Idgi_Pos = ((k1 * P3 * (Vposi * alphar + alphai * Vposr) + k2 * Q3 * (
				-Vposr * alphar + alphai * Vposi)) / vmag_pos) \
		           + Vposr * betai + betar * Vposi

		Idgr_Neg = Q3 * ((Vnegi * alphar - k2 * Vnegi * alphar + alphai * Vnegr - alphai * k2 * Vnegr) / vmag_neg) \
		           + P3 * ((-alphai * Vnegi + alphai * k1 * Vnegi + alphar * Vnegr - alphar * k1 * Vnegr) / vmag_neg) \
		           - betai * Vnegi + betar * Vnegr
		Idgi_Neg = Q3 * ((alphai * Vnegi - alphai * k2 * Vnegi - alphar * Vnegr + alphar * k2 * Vnegr) / vmag_neg) \
		           + P3 * ((Vnegi * alphar - k1 * Vnegi * alphar + alphai * Vnegr - alphai * k1 * Vnegr) / vmag_neg) \
		           + Vnegi * betar + betai * Vnegr

		return Idgr_Pos, Idgi_Pos, Idgr_Neg, Idgi_Neg

	@staticmethod
	def calc_partials(Vseq, phase, P3, Q3, alphar, alphai, betar, betai, k1,
	                  k2):

		Vposr = Vseq["VR"][0][phase]
		Vnegr = Vseq["VR"][1][phase]
		Vposi = Vseq["VI"][0][phase]
		Vnegi = Vseq["VI"][1][phase]

		k1 = k1[phase]
		k2 = k2[phase]

		vmag_pos = Vposr ** 2 + Vposi ** 2
		vmag_neg = Vnegr ** 2 + Vnegi ** 2

		vmag2_pos = vmag_pos ** 2
		vmag2_neg = vmag_neg ** 2

		dIdgrPos_dVposr = betar + Q3 * k2 * (
				(2 * (-alphar * Vposi * Vposr - alphai * Vposr ** 2) / vmag2_pos) + (alphai / vmag_pos)) \
		                  + P3 * k1 * ((2 * (alphai * Vposi * Vposr - alphar * Vposr ** 2) / vmag2_pos) + (
				alphar / vmag_pos))
		dIdgrPos_dVposi = -betai + Q3 * k2 * (
				(2 * (-alphar * Vposi ** 2 - alphai * Vposi * Vposr) / vmag2_pos) + (alphar / vmag_pos)) \
		                  + P3 * k1 * ((2 * (alphai * Vposi ** 2 - alphar * Vposi * Vposr) / vmag2_pos) - (
				alphai / vmag_pos))
		dIdgrPos_dQ3 = k2 * (Vposi * alphar + alphai * Vposr) / vmag_pos

		dIdgiPos_dVposr = betai + P3 * k1 * (
				(2 * (-alphar * Vposi * Vposr - alphai * Vposr ** 2) / vmag2_pos) + (alphai / vmag_pos)) \
		                  + Q3 * k2 * ((2 * (- alphai * Vposi * Vposr + alphar * Vposr ** 2) / vmag2_pos) - (
				alphar / vmag_pos))
		dIdgiPos_dVposi = betar + P3 * k1 * (
				(2 * (-alphar * Vposi ** 2 - alphai * Vposi * Vposr) / vmag2_pos) + (alphar / vmag_pos)) \
		                  + Q3 * k2 * ((2 * (-alphai * Vposi ** 2 + alphar * Vposi * Vposr) / vmag2_pos) + (
				alphai / vmag_pos))
		dIdgiPos_dQ3 = k2 * (-Vposr * alphar + alphai * Vposi) / vmag_pos

		dIdgrNeg_dVnegr = betar + Q3 * (
				((-2 * alphar * Vnegi * Vnegr - 2 * alphai * Vnegr ** 2 + k2 * (2 * alphar * Vnegi * Vnegr +
				                                                                2 * alphai * Vnegr ** 2)) /
				 vmag2_neg) + (
						alphai / vmag_neg) - (alphai * k2 / vmag_neg)) + \
		                  P3 * (((2 * alphai * Vnegi * Vnegr - 2 * alphar * Vnegr ** 2 + 2 * alphar * k1 * Vnegr ** 2
		                          - 2 * alphai * k1 * Vnegi * Vnegr) / vmag2_neg) + (alphar / vmag_neg) - (
				                        alphar * k1 / vmag_neg))
		dIdgrNeg_dVnegi = -betai + P3 * (
				((2 * alphai * Vnegi ** 2 - 2 * alphar * Vnegi * Vnegr + 2 * k1 * alphar * Vnegi * Vnegr -
				  2 * k1 * alphai * Vnegi ** 2) / vmag2_neg) - (alphai / vmag_neg) + (alphai * k1 / vmag_neg)) + \
		                  Q3 * (((-2 * alphar * Vnegi ** 2 - 2 * alphai * Vnegi * Vnegr + 2 * alphar * k2 * Vnegi ** 2
		                          + 2 * alphai * k2 * Vnegi * Vnegr) / vmag2_neg) + (alphar / vmag_neg) - (
				                        alphar * k2 / vmag_neg))
		dIdgrNeg_dQ3 = (Vnegi * alphar - k2 * Vnegi * alphar + alphai * Vnegr -
		                alphai * k2 * Vnegr) / vmag_neg

		dIdgiNeg_dVnegr = betai + P3 * (
				((-2 * alphar * Vnegi * Vnegr - 2 * alphai * Vnegr ** 2 + 2 * k1 * alphar * Vnegi * Vnegr +
				  2 * k1 * alphai * Vnegr ** 2) / vmag2_neg) + (alphai / vmag_neg) - (alphai * k1 / vmag_neg)) + \
		                  Q3 * (((- 2 * alphai * Vnegi * Vnegr + 2 * alphar * Vnegr ** 2 - 2 * alphar * k2 * Vnegr ** 2
		                          + 2 * alphai * k2 * Vnegi * Vnegr) / vmag2_neg) - (alphar / vmag_neg) + (
				                        alphar * k2 / vmag_neg))
		dIdgiNeg_dVnegi = betar + P3 * (
				((- 2 * alphar * Vnegi ** 2 - 2 * alphai * Vnegi * Vnegr + 2 * k1 * alphai * Vnegi * Vnegr +
				  2 * alphar * k1 * Vnegi ** 2) / vmag2_neg) + (alphar / vmag_neg) - (alphar * k1 / vmag_neg)) + \
		                  Q3 * (((-2 * alphai * Vnegi ** 2 + 2 * alphai * k2 * Vnegi ** 2 + 2 * alphar * Vnegi * Vnegr
		                          - 2 * alphar * k2 * Vnegi * Vnegr) / vmag2_neg) + (alphai / vmag_neg) - (
				                        alphai * k2 / vmag_neg))
		dIdgiNeg_dQ3 = (alphai * Vnegi - alphai * k2 * Vnegi - alphar * Vnegr +
		                alphar * k2 * Vnegr) / vmag_neg

		# Add partial derivatives and currents at iteration k to dictionaries
		dIdgrPos = {
			'dVposr': dIdgrPos_dVposr,
			'dVposi': dIdgrPos_dVposi,
			'dVnegr': 0,
			'dVnegi': 0,
			'dQ3': dIdgrPos_dQ3
		}

		dIdgiPos = {
			'dVposr': dIdgiPos_dVposr,
			'dVposi': dIdgiPos_dVposi,
			'dVnegr': 0,
			'dVnegi': 0,
			'dQ3': dIdgiPos_dQ3
		}

		dIdgrNeg = {
			'dVposr': 0,
			'dVposi': 0,
			'dVnegr': dIdgrNeg_dVnegr,
			'dVnegi': dIdgrNeg_dVnegi,
			'dQ3': dIdgrNeg_dQ3
		}

		dIdgiNeg = {
			'dVposr': 0,
			'dVposi': 0,
			'dVnegr': dIdgiNeg_dVnegr,
			'dVnegi': dIdgiNeg_dVnegi,
			'dQ3': dIdgiNeg_dQ3
		}

		dIdgrPos = SimpleNamespace(**dIdgrPos)
		dIdgiPos = SimpleNamespace(**dIdgiPos)
		dIdgrNeg = SimpleNamespace(**dIdgrNeg)
		dIdgiNeg = SimpleNamespace(**dIdgiNeg)

		return dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg

	def calc_ibdg_outputs(self, node, node_key, Vsol, print_flag=False):

		k1 = self.k1
		k2 = self.k2

		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		Vabc = self.get_nodes_states(node, node_key, Vsol)

		Va = Vabc.VR[0] + Vabc.VI[0] * 1j
		Vb = Vabc.VR[1] + Vabc.VI[1] * 1j
		Vc = Vabc.VR[2] + Vabc.VI[2] * 1j

		Vabc_mag = np.array([
			Vabc.VR[0], Vabc.VI[0], Vabc.VR[1], Vabc.VI[1], Vabc.VR[2],
			Vabc.VI[2]
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		IdgrPos = np.zeros((3, 1))
		IdgiPos = np.zeros((3, 1))
		IdgrNeg = np.zeros((3, 1))
		IdgiNeg = np.zeros((3, 1))

		Q3 = Vabc.Q[0]
		self.Q3_final = float(Q3)
		P3 = self.P3
		self.perc_Q = np.around(np.abs(self.Q3_final / self.Qmin), 5)

		for v in range(3):
			if self.ibdg_phases[v] & self.phases == int(self.ibdg_phases[v]):
				vmag_pos = Vposr[v] ** 2 + Vposi[v] ** 2
				vmag_neg = Vnegr[v] ** 2 + Vnegi[v] ** 2

				IdgrPos[v] = ((k2[v] * Q3 * (Vposi[v] * alphar + alphai * Vposr[v]) + k1[v] * P3 * (
						Vposr[v] * alphar - alphai * Vposi[v])) / vmag_pos) \
				             - Vposi[v] * betai + betar * Vposr[v]
				IdgiPos[v] = ((k1[v] * P3 * (Vposi[v] * alphar + alphai * Vposr[v]) + k2[v] * Q3 * (
						-Vposr[v] * alphar + alphai * Vposi[v])) / vmag_pos) \
				             + Vposr[v] * betai + betar * Vposi[v]

				IdgrNeg[v] = Q3 * ((Vnegi[v] * alphar - k2[v] * Vnegi[v] * alphar + alphai * Vnegr[v] - alphai * k2[
					v] *
				                    Vnegr[v]) / vmag_neg) \
				             + P3 * ((-alphai * Vnegi[v] + alphai * k1[v] * Vnegi[v] + alphar * Vnegr[v] - alphar * k1[
					v] *
				                      Vnegr[v]) / vmag_neg) \
				             - betai * Vnegi[v] + betar * Vnegr[v]
				IdgiNeg[v] = Q3 * ((alphai * Vnegi[v] - alphai * k2[v] * Vnegi[v] - alphar * Vnegr[v] + alphar * k2[
					v] *
				                    Vnegr[v]) / vmag_neg) \
				             + P3 * ((Vnegi[v] * alphar - k1[v] * Vnegi[v] * alphar + alphai * Vnegr[v] - alphai * k1[
					v] *
				                      Vnegr[v]) / vmag_neg) \
				             + Vnegi[v] * betar + betai * Vnegr[v]

		Iar = IdgrPos[0] + IdgrNeg[0]
		Iai = IdgiPos[0] + IdgiNeg[0]
		Ibr = IdgrPos[1] + IdgrNeg[1]
		Ibi = IdgiPos[1] + IdgiNeg[1]
		Icr = IdgrPos[2] + IdgrNeg[2]
		Ici = IdgiPos[2] + IdgiNeg[2]

		Ia = (Iar + 1j * Iai) / 3
		Ib = (Ibr + 1j * Ibi) / 3
		Ic = (Icr + 1j * Ici) / 3

		self.Ia_mag = np.around(float(np.abs(Ia)), 5)
		self.Ia_ang = np.around(float(np.angle(Ia, True)), 5)

		self.Ib_mag = np.around(float(np.abs(Ib)), 5)
		self.Ib_ang = np.around(float(np.angle(Ib, True)), 5)

		self.Ic_mag = np.around(float(np.abs(Ic)), 5)
		self.Ic_ang = np.around(float(np.angle(Ic, True)), 5)

		self.Qmax = np.around(self.Qmax, 5)
		self.Qmin = np.around(self.Qmin, 5)
		self.Q3_final = np.around(self.Q3_final, 5)

		# self.perc_Q = abs(self.Q3_final[0] / self.Qscale)

		S3 = Va * np.conj(Ia) + Vb * np.conj(Ib) + Vc * np.conj(Ic)

		if print_flag:
			logging.info("---------------------------------------------------------")
			logging.info(self.name)
			logging.info("Currents are: Ia = %f < %f | Ib = %f < %f | Ic = %f < %f" %
			             (self.Ia_mag, self.Ia_ang, self.Ib_mag, self.Ib_ang,
			              self.Ic_mag, self.Ic_ang))
			logging.info("Three-Phase Powers are: P = %f kW | Q = %f kVAR" %
			             (S3.real * 1e-3, S3.imag * 1e-3))

	def calc_max_current(self, Q3, lf):
		Qcurr = Q3
		Pref = -self.P3 * lf
		k1 = self.k1
		k2 = self.k2

		Vpos_mag = np.linalg.norm(self.Vpos)
		Vneg_mag = np.linalg.norm(self.Vneg)
		Vpos_ang = np.linalg.norm(np.angle(self.Vpos))
		Vneg_ang = np.linalg.norm(np.angle(self.Vneg))

		gamma_a = float((Vpos_ang - Vneg_ang) / 2)
		gamma_b = gamma_a + (np.pi / 3)
		gamma_c = gamma_a - (np.pi / 3)
		gamma = np.array([[gamma_a], [gamma_b], [gamma_c]])

		Vneg_mag = Vneg_mag if Vneg_mag != 0 else 1e-14

		IpL = Pref * ((k1 / Vpos_mag) + ((1 - k1) / Vneg_mag))
		IpS = Pref * ((k1 / Vpos_mag) - ((1 - k1) / Vneg_mag))
		IqL = Qcurr * ((k2 / Vpos_mag) + ((1 - k2) / Vneg_mag))
		IqS = Qcurr * ((k2 / Vpos_mag) - ((1 - k2) / Vneg_mag))

		A1 = IpL * np.cos(gamma) - IqL * np.sin(gamma)
		B1 = -IqS * np.cos(gamma) - IpS * np.sin(gamma)

		Ipeak = np.sqrt(A1 ** 2 + B1 ** 2)

		return Ipeak

	def calc_Q_limits(self, Imax):
		P3_ref = 0
		k2 = self.k2
		k1 = self.k1

		Q = np.zeros((3, 1))
		Q_max = np.zeros((3, 1))
		Q_min = np.zeros((3, 1))

		Vpos_mag = np.linalg.norm(self.Vpos)
		Vneg_mag = np.linalg.norm(self.Vneg)
		Vpos_ang = np.linalg.norm(np.angle(self.Vpos))
		Vneg_ang = np.linalg.norm(np.angle(self.Vneg))

		Vneg_mag = Vneg_mag if Vneg_mag != 0 else 1e-14

		gamma_a = (Vpos_ang - Vneg_ang) / 2
		gamma_b = gamma_a + (np.pi / 3)
		gamma_c = gamma_a - (np.pi / 3)
		gamma = np.array([[gamma_a], [gamma_b], [gamma_c]])

		Imax = Imax.reshape((3, 1))

		if P3_ref != 0:
			# simultaneous active and reactive power delivery
			a = (k2 * Vneg_mag ** 2) + ((1 - k2) ** 2 * (Vpos_mag ** 2)) - (
					2 * k2 * (1 - k2) * np.cos(2 * gamma) * Vpos_mag * Vneg_mag)
			b = P3_ref * ((2 * k1 + 2 * k2 - 4 * k1 * k2) * Vpos_mag * Vneg_mag * np.sin(2 * gamma))
			c = (P3_ref ** 2) * ((k1 * Vneg_mag ** 2) + (((1 - k1) ** 2) * Vpos_mag ** 2)
			                     + 2 * k1 * (1 - k1) * Vneg_mag * Vpos_mag * np.cos(2 * gamma)) - \
			    (Imax ** 2) * (Vneg_mag ** 2) * (Vpos_mag ** 2)
			for phase in range(3):
				Q_max[phase] = np.roots([a[phase].item(), b[phase].item(), c[phase].item()])[0]
				Q_min[phase] = np.roots([a[phase].item(), b[phase].item(), c[phase].item()])[1]
			Q = Q_max
		else:
			# injection of maximum reactive power
			alt1 = ((Imax ** 2) * (Vpos_mag ** 2) * (Vneg_mag ** 2))
			alt2 = (k2 ** 2) * (Vneg_mag ** 2) + ((1 - k2) ** 2 * (Vpos_mag ** 2)) - (
					2 * k2 * (1 - k2) * np.cos(2 * gamma)) \
			       * Vpos_mag * Vneg_mag
			Q = np.sqrt(alt1 / alt2)

		Q3 = np.min(Q)
		Q3_loc = np.argmin(Q)

		if Q3_loc == 0:
			self.Qphase = 'A'
		elif Q3_loc == 1:
			self.Qphase = 'B'
		else:
			self.Qphase = 'C'

		if Q3 < 0:
			self.Qmin = -Q3
			self.Qmax = Q3
		elif Q3 > 0:
			self.Qmax = -Q3
			self.Qmin = Q3

	def stamp_nonlinear(self, node_key, node, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    Jnlin_val, Jnlin_row, idx_Y, idx_J, homotopy_enabled=False, h_factor=None):

		alphar = self.alpha_r
		alphai = self.alpha_i
		betar = self.beta_r
		betai = self.beta_i

		Vabc = self.get_nodes_states(node, node_key, V)

		Vabc_mag = np.array([
			Vabc.VR[0], Vabc.VI[0], Vabc.VR[1], Vabc.VI[1], Vabc.VR[2],
			Vabc.VI[2]
		], dtype=object)

		Vposr = self.Pr.dot(Vabc_mag)
		Vnegr = self.Nr.dot(Vabc_mag)
		Vposi = self.Pi.dot(Vabc_mag)
		Vnegi = self.Ni.dot(Vabc_mag)

		Vseq = {
			'VR': [Vposr, Vnegr],
			'VI': [Vposi, Vnegi],
			'Q3': Vabc.Q,
			'nodeVR': Vabc.nodeVR,
			'nodeVI': Vabc.nodeVI,
			'nodeQ3': Vabc.nodeQ3
		}

		Vpos = Vposr + 1j * Vposi
		Vneg = Vnegr + 1j * Vnegi

		Q3_ref = Vabc.Q[0]

		P3_ref = self.P3

		k1 = self.k1
		k2 = self.k2

		IdgrPos_hist, IdgiPos_hist, IdgrNeg_hist, IdgiNeg_hist = self.calc_hist_currents(
			self.Vpos, self.Vneg, Q3_ref, P3_ref, alphar, alphai, betar, betai, k1, k2)
		self.IdgrPos = IdgrPos_hist
		self.IdgiPos = IdgiPos_hist
		self.IdgrNeg = IdgrNeg_hist
		self.IdgiNeg = IdgiNeg_hist

		for _phase_index in range(3):
			if self.ibdg_phases[_phase_index] & self.phases == int(
					self.ibdg_phases[_phase_index]):
				dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg = self.calc_partials(
					Vseq, _phase_index, P3_ref, Q3_ref, alphar, alphai, betar,
					betai, k1, k2)

				IdgrPos_phase = self.IdgrPos[_phase_index]
				IdgiPos_phase = self.IdgiPos[_phase_index]
				IdgrNeg_phase = self.IdgrNeg[_phase_index]
				IdgiNeg_phase = self.IdgiNeg[_phase_index]

				idx_Y, idx_J = super(FPNSC, self).stamp_ibdg_partials(
					Vseq, _phase_index, dIdgrPos, dIdgiPos, dIdgrNeg, dIdgiNeg,
					IdgrPos_phase, IdgiPos_phase, IdgrNeg_phase, IdgiNeg_phase,
					Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
					idx_Y, idx_J)

		if not self.fixed_Q:
			dFQ3, FQ3_hist = self.calc_Q3(Vabc)

			idx_Y, idx_J = super(FPNSC,
			                     self).stamp_Q3_partials(Vabc, dFQ3, FQ3_hist,
			                                             Ynlin_val, Ynlin_row,
			                                             Ynlin_col, Jnlin_val,
			                                             Jnlin_row, idx_Y, idx_J, homotopy_enabled, h_factor)
		self.ibdg_counter += 1

		return idx_Y, idx_J
