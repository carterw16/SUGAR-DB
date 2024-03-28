"""Power delivery element that steps-up or steps-down voltage.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 04-10-2017
  Updated Date: 10-14-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Development

  A power deliver element with two or more windings that steps-up or steps-down voltage in power system. Types of
  transformers currently implemented include Wye-Wye, Gwye-Gwye, Delta-Delta, Delta-Gwye, Single Phase, Center-Tap,
  Delta-Wye, Wye-Delta, and Gwye-Delta.

"""
from __future__ import division

from itertools import count
from types import MethodType, SimpleNamespace
import numpy as np

from .Elements import LinearElement
from .GlobalVars import _XFMR
from .TransformerStamp import  (stamp_tx_stepping, stamp_primary, stamp_secondary, stamp_ground_neutral,
								stamp_linear) 
from .Nodes import Nodes



class Transformer(LinearElement):
	_ids = count(0)

	# Create variable names
	_WYE_WYE = -1
	_GWYE_GWYE = 1
	_DELTA_DELTA = 2
	_DELTA_GWYE = 3
	_SINGLE_PHASE = 4
	_CENTER_TAP = 5
	# _GWYE_GWYE = 6
	_DELTA_WYE = 7
	_WYE_DELTA = 8
	_GWYE_DELTA = 9
	# Phase name
	phaseA = 0x01
	phaseB = 0x02
	phaseC = 0x04

	@staticmethod
	def sind(x):
		return np.sin(float(x) * np.pi / 180)

	@staticmethod
	def cosd(x):
		return np.cos(float(x) * np.pi / 180)

	@staticmethod
	def calc_G_B(r, x):
		if r == 0:
			G = None
		else:
			G = r / (r**2 + x**2)
		if x == 0:
			B = None
		else:
			B = -x / (r**2 + x**2)
		return G, B
		
	def __init__(self,
				 name,
				 ID,
				 casetype,
				 phases,
				 from_node,
				 to_node,
				 nominal_voltage,
				 connect_type,
				 install_type,
				 power_rating,
				 primary_voltage,
				 secondary_voltage,
				 resistance,
				 reactance,
				 shunt_impedance,
				 connected_regulator,
				 is_center_tap=False,
				 is_grounded=False,
				 full_load_loss=None,
				 no_load_loss=None,
				 phase_shift=False,
				 stamp_dual=False):
 		
		"""

		Args:
			node_index_ (iterator): tracking the total number of nodes
			node (list): contains all existing node objects
			name (str): name of the transformer
			ID (str): the ID of the transformer
			phases (int): the phases that the transformer has
			from_node (str): node at the sending end of the transformer
			to_node (str): node at the receiving end of the transformer
			nominal_voltage (float): nominal voltage of the transformer
			connect_type (int): the connection type of the transformer (e.g, Y-Y, D-Y, etc.)
			install_type (str): The mounting type of the transformer: one of {POLETOP, PADMOUNT, VAULT}
			power_rating (float): Rated power of the transformer. Could be three-phase or single-phase.
			primary_voltage (float): nominal voltage at the sending end of the transformer
			secondary_voltage (float): nominal voltage at the receiving end of the transformer
			resistance (list): Resistances between the transformer windings
			reactance (list): Reactances between the transformer windings
			shunt_impedance (complex): Shunt impedance of the transformer
			connected_regulator (Regulator): If the transformer is linked to a regulator, this is a Regulator instance.
			is_center_tap (bool): Identifies if the transformer is center tap or not.
			is_grounded (bool): Identifies if the transformer is grounded or not.
			full_load_loss (float): Percent Losses at rated load
			no_load_loss (float): The no-load loss for a zero sequence short-circuit test on the entire transformer
			phase_shift (float): the phase shift angle in degrees if applicable, otherwise None.
		"""
		super(Transformer, self).__init__()
		self.casetype = casetype
		self.name = name
		self.id = self._ids.__next__()
		self.ID = ID
		self.phases = int(phases)
		self.from_node = from_node
		self.to_node = to_node
		self.nominal_voltage = float(
			nominal_voltage) if nominal_voltage else None
		self.connect_type = connect_type
		self.install_type = install_type
		self.primary_voltage = primary_voltage
		self.secondary_voltage = secondary_voltage
		self.phase_shift = float(phase_shift) if phase_shift else None
		self.is_center_tap = is_center_tap
		self.is_grounded = is_grounded
		self.full_load_loss = float(full_load_loss) if full_load_loss else None
		self.no_load_loss = float(no_load_loss) if no_load_loss else None
		self.shunt_impedance = shunt_impedance if shunt_impedance else None
		self.connected_regulator = connected_regulator
		self.stamp_dual = stamp_dual

		# # Link Cython Methods to Class # #
		self.stamp_linear = MethodType(stamp_linear, self)
		self.stamp_primary = MethodType(stamp_primary, self)
		self.stamp_secondary = MethodType(stamp_secondary, self)
		self.stamp_ground_neutral = MethodType(stamp_ground_neutral, self)
		self.stamp_tx_stepping = MethodType(stamp_tx_stepping, self)


		# Calculated Parameters
		if self.connect_type == _XFMR['WYE_WYE'] or self.connect_type == _XFMR[
				'GWYE_GWYE']:  # Wye-Wye or Gwye-Gwye
			self.Vpri = self.primary_voltage
			self.Vsec = self.secondary_voltage
			self.ang = self.phase_shift
		elif self.connect_type == _XFMR[
				'DELTA_DELTA']:  # If primary voltage is given line to line
			self.Vpri = self.primary_voltage * np.sqrt(3)
			self.Vsec = self.secondary_voltage * np.sqrt(3)
			self.ang = self.phase_shift
		elif self.connect_type == _XFMR[
				'DELTA_GWYE']:  # Delta - wye or Delta - Gwye
			self.Vpri = self.primary_voltage * np.sqrt(3)
			self.Vsec = self.secondary_voltage
			self.ang = self.phase_shift
		elif self.connect_type == _XFMR['WYE_DELTA']:  # Wye-Delta or Gwye-Delta
			self.Vpri = self.primary_voltage
			self.Vsec = self.secondary_voltage * np.sqrt(3)
			self.ang = self.phase_shift
		elif self.connect_type == _XFMR[
				'SINGLE_PHASE_CENTER_TAPPED']:  # center-tap transformer
			self.Vpri = self.primary_voltage
			self.Vsec = self.secondary_voltage
			self.ang = None
		else:
			print('Invalid Transformer Type')

		# Assign Nodes
		self.nodeVr_meas = {'A': 0, 'B': 0, 'C': 0}
		self.nodeVi_meas = {'A': 0, 'B': 0, 'C': 0}
		if self.connect_type not in [self._SINGLE_PHASE, self._CENTER_TAP]:
			self.xfmr_type = 'Three-Phase'
			self.isSinglePhase = False
		   
			if self.ang:
				self.all_phases = {0: "A", 1: "B", 2: "C", 3: "N"}
				self.nodes_secondary_link = {'A':[], 'B':[], 'C': []}
				self.nodes_primary_link = {'A': [], 'B': [], 'C': []}
		 
		else:
			self.xfmr_type = 'Center Tap'
			self.isSinglePhase = True
			

		# Calculate turns ratio r and x (For single phase do not use sqrt(3))
		self.tr = self.Vpri / self.Vsec
		if self.tr > 1:
			self.is_stepdown = True
		else:
			self.is_stepdown = False

		#if self.is_center_tap and len(power_rating) > 1:
		if self.is_center_tap and isinstance(power_rating, list):
			self._r = resistance[0] / 0.5
			self._x = reactance[0] / 0.8
			if self.phases & 0x1 == 1:  # Check for phase A
				self.power_rating = power_rating[0]
			elif self.phases & 0x2 == 2:
				self.power_rating = power_rating[1]
			else:
				self.power_rating = power_rating[2]
		elif self.is_center_tap and isinstance(power_rating, float):
			self._r = resistance[0] / 0.5
			self._x = reactance[0] / 0.8
			self.power_rating = power_rating
		else:
			self._r = resistance
			self._x = reactance
			self.power_rating = power_rating

		self.r = (self._r * self.Vsec**2) / self.power_rating
		self.x = (self._x * self.Vsec**2) / self.power_rating

		if self.connect_type == 5:  # Reference Kersting (2007)
			_ZbaseH = self.Vpri**2 / self.power_rating
			_ZbaseL = self.Vsec**2 / self.power_rating
			# Correct this if impedance, impedance1 and impedance2 is given instead
			self._r0 = resistance[0]
			self._r1 = resistance[1]
			self._r2 = resistance[2]
			self._x0 = reactance[0]
			self._x1 = reactance[1]
			self._x2 = reactance[2]
			# Convert the impedance to ohms
			self.r0 = self._r0 * _ZbaseH
			self.x0 = self._x0 * _ZbaseH
			self.r1 = self._r1 * _ZbaseL
			self.x1 = self._x1 * _ZbaseL
			self.r2 = self._r2 * _ZbaseL
			self.x2 = self._x2 * _ZbaseL

			(self.G0, self.B0) = self.calc_G_B(self.r0, self.x0)
			(self.G1, self.B1) = self.calc_G_B(self.r1, self.x1)
			(self.G2, self.B2) = self.calc_G_B(self.r2, self.x2)

		(self.G, self.B) = self.calc_G_B(self.r, self.x)

		if self.shunt_impedance:
			if self.connect_type != 5:
				self.xShunt = (self.shunt_impedance *
							   self.Vsec**2) / self.power_rating
			else:
				self.xShunt = (self.shunt_impedance *
							   self.Vpri**2) / self.power_rating
			self.hasShunt = True
		else:
			self.hasShunt = False
		if self.hasShunt:
			(self.Gshunt, self.Bshunt) = self.calc_G_B(self.xShunt.real,
													   self.xShunt.imag)
			self.hasShunt = True
		else:
			self.Gshunt, self.Bshunt = 0.0, 0.0
			self.hasShunt = False

		# Initialize Transformer Output Currents (L-N)
		if self.connect_type != _XFMR['SINGLE_PHASE_CENTER_TAPPED']:
			self.Ia_mag = 0
			self.Ia_ang = 0
			self.Ib_mag = 0
			self.Ib_ang = 0
			self.Ic_mag = 0
			self.Ic_ang = 0
		else:
			self.I1_mag = 0
			self.I2_mag = 0
			self.I1_ang = 0
			self.I2_ang = 0
			
	def assign_nodes_3phase(self, node_key, node_index_, node):
		# There are likely to be 12 additional nodes for the Transformer
		# 6 additional rows for the voltage source eqns. (No Phase Shifters
		# 12 additional rows for the voltage sources eqn (w Phase Shifters
		self.nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr
		self.nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi
		self.nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr
		self.nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi
		self.nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
		self.nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi
		self.nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
		self.nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi
		self.nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
		self.nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi
		self.nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
		self.nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi
		self.nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
		self.nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi
		self.nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
		self.nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi
		
		if self.stamp_dual:
			self.nodeA_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
			self.nodeA_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i
			self.nodeB_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
			self.nodeB_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
			self.nodeC_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
			self.nodeC_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i
			self.nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
			self.nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i
			self.nodeA_Lr_to = node[node_key[self.to_node]].nodeA_dual_eq_var_r
			self.nodeA_Li_to = node[node_key[self.to_node]].nodeA_dual_eq_var_i
			self.nodeB_Lr_to = node[node_key[self.to_node]].nodeB_dual_eq_var_r
			self.nodeB_Li_to = node[node_key[self.to_node]].nodeB_dual_eq_var_i
			self.nodeC_Lr_to = node[node_key[self.to_node]].nodeC_dual_eq_var_r
			self.nodeC_Li_to = node[node_key[self.to_node]].nodeC_dual_eq_var_i
			self.nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
			self.nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i
		#############################################################
		# Add ground voltage source if needed
		if self.connect_type in [
				self._DELTA_DELTA, self._WYE_WYE, self._DELTA_GWYE,
				self._GWYE_GWYE, self._DELTA_WYE, self._WYE_DELTA,
				self._GWYE_DELTA
		]:
			# Add for the three phase extra node in the secondary
			# circuit (for both real and imaginary)
			self.nodeA_Vr_pos_secondary = node_index_.__next__()
			self.nodeA_Vi_pos_secondary = node_index_.__next__()
			self.nodeB_Vr_pos_secondary = node_index_.__next__()
			self.nodeB_Vi_pos_secondary = node_index_.__next__()
			self.nodeC_Vr_pos_secondary = node_index_.__next__()
			self.nodeC_Vi_pos_secondary = node_index_.__next__()

			# Add for the three phases extra node for the voltage sources
			self.nodeA_Vr_primary = node_index_.__next__()
			self.nodeA_Vi_primary = node_index_.__next__()
			self.nodeB_Vr_primary = node_index_.__next__()
			self.nodeB_Vi_primary = node_index_.__next__()
			self.nodeC_Vr_primary = node_index_.__next__()
			self.nodeC_Vi_primary = node_index_.__next__()
			
			if self.stamp_dual:
				# Add for the three phase extra node in the secondary
				# circuit (for both real and imaginary)
				self.nodeA_Lr_pos_secondary = node_index_.__next__()
				self.nodeA_Li_pos_secondary = node_index_.__next__()
				self.nodeB_Lr_pos_secondary = node_index_.__next__()
				self.nodeB_Li_pos_secondary = node_index_.__next__()
				self.nodeC_Lr_pos_secondary = node_index_.__next__()
				self.nodeC_Li_pos_secondary = node_index_.__next__()
	
				# Add for the three phases extra node for the voltage sources
				self.nodeA_Lr_primary = node_index_.__next__()
				self.nodeA_Li_primary = node_index_.__next__()
				self.nodeB_Lr_primary = node_index_.__next__()
				self.nodeB_Li_primary = node_index_.__next__()
				self.nodeC_Lr_primary = node_index_.__next__()
				self.nodeC_Li_primary = node_index_.__next__()

			# Add additional voltage sources for DELTA-WYE and WYE-DELTA
			if self.ang:
				self.nodeA_Vr_secondary_link = node_index_.__next__()
				self.nodeA_Vi_secondary_link = node_index_.__next__()
				self.nodeB_Vr_secondary_link = node_index_.__next__()
				self.nodeB_Vi_secondary_link = node_index_.__next__()
				self.nodeC_Vr_secondary_link = node_index_.__next__()
				self.nodeC_Vi_secondary_link = node_index_.__next__()

				self.nodeA_Vr_primary_link = node_index_.__next__()
				self.nodeA_Vi_primary_link = node_index_.__next__()
				self.nodeB_Vr_primary_link = node_index_.__next__()
				self.nodeB_Vi_primary_link = node_index_.__next__()
				self.nodeC_Vr_primary_link = node_index_.__next__()
				self.nodeC_Vi_primary_link = node_index_.__next__()
				
				if self.stamp_dual:
					self.nodeA_Lr_secondary_link = node_index_.__next__()
					self.nodeA_Li_secondary_link = node_index_.__next__()
					self.nodeB_Lr_secondary_link = node_index_.__next__()
					self.nodeB_Li_secondary_link = node_index_.__next__()
					self.nodeC_Lr_secondary_link = node_index_.__next__()
					self.nodeC_Li_secondary_link = node_index_.__next__()
	
					self.nodeA_Lr_primary_link = node_index_.__next__()
					self.nodeA_Li_primary_link = node_index_.__next__()
					self.nodeB_Lr_primary_link = node_index_.__next__()
					self.nodeB_Li_primary_link = node_index_.__next__()
					self.nodeC_Lr_primary_link = node_index_.__next__()
					self.nodeC_Li_primary_link = node_index_.__next__()
					
					self.nodes_secondary_link["A"] = [
						self.nodeA_Vr_secondary_link, self.nodeA_Vi_secondary_link,
						self.nodeA_Lr_secondary_link, self.nodeA_Li_secondary_link
					]
					self.nodes_secondary_link["B"] = [
						self.nodeB_Vr_secondary_link, self.nodeB_Vi_secondary_link,
						self.nodeB_Lr_secondary_link, self.nodeB_Li_secondary_link
					]
					self.nodes_secondary_link["C"] = [
						self.nodeC_Vr_secondary_link, self.nodeC_Vi_secondary_link,
						self.nodeC_Lr_secondary_link, self.nodeC_Li_secondary_link
					]
	
					self.nodes_primary_link["A"] = [
						self.nodeA_Vr_primary_link, self.nodeA_Vi_primary_link,
						self.nodeA_Lr_primary_link, self.nodeA_Li_primary_link
					]
					self.nodes_primary_link["B"] = [
						self.nodeB_Vr_primary_link, self.nodeB_Vi_primary_link,
						self.nodeB_Lr_primary_link, self.nodeB_Li_primary_link
					]
					self.nodes_primary_link["C"] = [
						self.nodeC_Vr_primary_link, self.nodeC_Vi_primary_link,
						self.nodeC_Lr_primary_link, self.nodeC_Li_primary_link
					]
				
				else:
					self.nodes_secondary_link["A"] = [
						self.nodeA_Vr_secondary_link, self.nodeA_Vi_secondary_link
					]
					self.nodes_secondary_link["B"] = [
						self.nodeB_Vr_secondary_link, self.nodeB_Vi_secondary_link
					]
					self.nodes_secondary_link["C"] = [
						self.nodeC_Vr_secondary_link, self.nodeC_Vi_secondary_link
					]
	
					self.nodes_primary_link["A"] = [
						self.nodeA_Vr_primary_link, self.nodeA_Vi_primary_link
					]
					self.nodes_primary_link["B"] = [
						self.nodeB_Vr_primary_link, self.nodeB_Vi_primary_link
					]
					self.nodes_primary_link["C"] = [
						self.nodeC_Vr_primary_link, self.nodeC_Vi_primary_link
					]
				

			# Check for ground and add further nodes
			if self.connect_type in [self._GWYE_DELTA, self._GWYE_GWYE]:
				self.nodeGnd_Vr_primary = node_index_.__next__()
				self.nodeGnd_Vi_primary = node_index_.__next__()
				
				if self.stamp_dual:
					self.nodeGnd_Lr_primary = node_index_.__next__()
					self.nodeGnd_Li_primary = node_index_.__next__()
				Nodes.xfmr_primal_ground.append([self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary])
				Nodes.xfmr_dual_ground.append([self.nodeGnd_Lr_primary, self.nodeGnd_Li_primary])
			
			if self.connect_type in [self._DELTA_GWYE, self._GWYE_GWYE]:
				self.nodeGnd_Vr_secondary = node_index_.__next__()
				self.nodeGnd_Vi_secondary = node_index_.__next__()                
				
				if self.stamp_dual:
					self.nodeGnd_Lr_secondary = node_index_.__next__()
					self.nodeGnd_Li_secondary = node_index_.__next__()     
				Nodes.xfmr_primal_ground.append([self.nodeGnd_Vr_secondary, self.nodeGnd_Vi_secondary])
				Nodes.xfmr_dual_ground.append([self.nodeGnd_Lr_secondary, self.nodeGnd_Li_secondary])
			   
		else:
			print('Incorrect transformer type for assigning nodes')
		return node_index_

	# A separate function to assign node for center-tap transformer
	# from node should be AN, BN (regular node)
	# to node should be triplex node (1, 2, N)
	def assign_nodes_centertap(self, node_key, node_index_, node):
		# First find the Vr and Vi of the single phase that is connected
		if self.phases & 0x01 == 1:  # Phase A
			self.node1_Vr_from = node[node_key[self.from_node]].nodeA_Vr
			self.node1_Vi_from = node[node_key[self.from_node]].nodeA_Vi
		elif self.phases & 0x02 == 2:  # Phase B
			self.node1_Vr_from = node[node_key[self.from_node]].nodeB_Vr
			self.node1_Vi_from = node[node_key[self.from_node]].nodeB_Vi
		elif self.phases & 0x04 == 4:  # Phase C
			self.node1_Vr_from = node[node_key[self.from_node]].nodeC_Vr
			self.node1_Vi_from = node[node_key[self.from_node]].nodeC_Vi
		self.nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
		self.nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi
		# Secondary Hot Wires
		self.node1_Vr_to = node[node_key[self.to_node]].node1_Vr
		self.node1_Vi_to = node[node_key[self.to_node]].node1_Vi
		self.node2_Vr_to = node[node_key[self.to_node]].node2_Vr
		self.node2_Vi_to = node[node_key[self.to_node]].node2_Vi
		self.nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
		self.nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi
		
		
		if self.stamp_dual:
			if self.phases & 0x01 == 1:  # Phase A
				self.node1_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
				self.node1_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i
			elif self.phases & 0x02 == 2:  # Phase B
				self.node1_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
				self.node1_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
			elif self.phases & 0x04 == 4:  # Phase C
				self.node1_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
				self.node1_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i
			self.nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
			self.nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i
			# Secondary Hot Wires
			self.node1_Lr_to = node[node_key[self.to_node]].node1_dual_eq_var_r
			self.node1_Li_to = node[node_key[self.to_node]].node1_dual_eq_var_i
			self.node2_Lr_to = node[node_key[self.to_node]].node2_dual_eq_var_r
			self.node2_Li_to = node[node_key[self.to_node]].node2_dual_eq_var_i
			self.nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
			self.nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i           
		# Assign new nodes
		try:
			if self.connect_type == self._CENTER_TAP:
				if not self.phase_shift:
					# 2 additional nodes on the primary side
					self.node_Vr_CT_primary = node_index_.__next__()
					self.node_Vi_CT_primary = node_index_.__next__()

					# 4 additional extra nodes on sending side
					self.node1_Vr_sending = node_index_.__next__()
					self.node1_Vi_sending = node_index_.__next__()
					self.node2_Vr_sending = node_index_.__next__()
					self.node2_Vi_sending = node_index_.__next__()

					# 4 additional rows for voltage sources
					self.node1_Vr_secondary = node_index_.__next__()
					self.node1_Vi_secondary = node_index_.__next__()
					self.node2_Vr_secondary = node_index_.__next__()
					self.node2_Vi_secondary = node_index_.__next__()

					# 2 additional variables from grounding neutrals
					self.nodeGnd_Vr_secondary = node_index_.__next__()
					self.nodeGnd_Vi_secondary = node_index_.__next__()
					self.nodeGnd_Vr_primary = node_index_.__next__()
					self.nodeGnd_Vi_primary = node_index_.__next__()
					
					if self.stamp_dual:
						# 2 additional nodes on the primary side
						self.node_Lr_CT_primary = node_index_.__next__()
						self.node_Li_CT_primary = node_index_.__next__()
	
						# 4 additional extra nodes on sending side
						self.node1_Lr_sending = node_index_.__next__()
						self.node1_Li_sending = node_index_.__next__()
						self.node2_Lr_sending = node_index_.__next__()
						self.node2_Li_sending = node_index_.__next__()
	
						# 4 additional rows for voltage sources
						self.node1_Lr_secondary = node_index_.__next__()
						self.node1_Li_secondary = node_index_.__next__()
						self.node2_Lr_secondary = node_index_.__next__()
						self.node2_Li_secondary = node_index_.__next__()
	
						# # 2 additional variables from grounding neutrals
						self.nodeGnd_Lr_secondary = node_index_.__next__()
						self.nodeGnd_Li_secondary = node_index_.__next__()
						self.nodeGnd_Lr_primary = node_index_.__next__()
						self.nodeGnd_Li_primary = node_index_.__next__()  
				else:
					pass
			else:
				raise Exception
		except Exception:
			print('Undefined Transformer Type')
		return node_index_

	def calc_currents(self, V):

		if self.connect_type == _XFMR['SINGLE_PHASE_CENTER_TAPPED']:
			# Record the secondary phase 1 current
			V1r_secondary = V[self.node1_Vr_sending]
			V1i_secondary = V[self.node1_Vi_sending]
			V1r_to = V[self.node1_Vr_to]
			V1i_to = V[self.node1_Vi_to]
			I1r = (V1r_to - V1r_secondary) * self.G1 - self.B1 * (V1i_to - V1i_secondary) - self.Bshunt * V1i_to + \
				  self.Gshunt * V1r_to
			I1i = (V1i_to - V1i_secondary) * self.G1 + self.B1 * (V1r_to - V1r_secondary) + self.Bshunt * V1r_to + \
				  self.Gshunt * V1i_to
			self.I1_mag = np.around(np.abs(complex(I1r, I1i)), 2)
			self.I1_ang = np.around(np.angle(complex(I1r, I1i), deg=True), 2)

			# Record the secondary phase 2 current
			V2r_secondary = V[self.node2_Vr_sending]
			V2i_secondary = V[self.node2_Vi_sending]
			V2r_to = V[self.node2_Vr_to]
			V2i_to = V[self.node2_Vi_to]
			I2r = (V2r_to - V2r_secondary) * self.G2 - self.B2 * (V2i_to - V2i_secondary) - self.Bshunt * V2i_to + \
				  self.Gshunt * V2r_to
			I2i = (V2i_to - V2i_secondary) * self.G2 + self.B2 * (V2r_to - V2r_secondary) + self.Bshunt * V2r_to + \
				  self.Gshunt * V2i_to
			self.I2_mag = np.around(np.abs(complex(I2r, I2i)), 2)
			self.I2_ang = np.around(np.angle(complex(I2r, I2i), deg=True), 2)
		else:
			# Record the secondary phase a current
			Var_secondary = V[self.nodeA_Vr_pos_secondary]
			Vai_secondary = V[self.nodeA_Vi_pos_secondary]
			Var_to = V[self.nodeA_Vr_to]
			Vai_to = V[self.nodeA_Vi_to]
			Iar = (Var_to - Var_secondary) * self.G - self.B * (Vai_to - Vai_secondary) - self.Bshunt * Vai_to + \
				  self.Gshunt * Var_to
			Iai = (Vai_to - Vai_secondary) * self.G + self.B * (Var_to - Var_secondary) + self.Bshunt * Var_to + \
				  self.Gshunt * Vai_to
			self.Ia_mag = np.around(np.abs(complex(Iar, Iai)), 2)
			self.Ia_ang = np.around(np.angle(complex(Iar, Iai), deg=True), 2)

			# Record the secondary phase b current
			Vbr_secondary = V[self.nodeB_Vr_pos_secondary]
			Vbi_secondary = V[self.nodeB_Vi_pos_secondary]
			Vbr_to = V[self.nodeB_Vr_to]
			Vbi_to = V[self.nodeB_Vi_to]
			Ibr = (Vbr_to - Vbr_secondary) * self.G - self.B * (Vbi_to - Vbi_secondary) - self.Bshunt * Vbi_to + \
				  self.Gshunt * Vbr_to
			Ibi = (Vbi_to - Vbi_secondary) * self.G + self.B * (Vbr_to - Vbr_secondary) + self.Bshunt * Vbr_to + \
				  self.Gshunt * Vbi_to
			self.Ib_mag = np.around(np.abs(complex(Ibr, Ibi)), 2)
			self.Ib_ang = np.around(np.angle(complex(Ibr, Ibi), deg=True), 2)

			# Record the secondary phase c current
			Vcr_secondary = V[self.nodeC_Vr_pos_secondary]
			Vci_secondary = V[self.nodeC_Vi_pos_secondary]
			Vcr_to = V[self.nodeC_Vr_to]
			Vci_to = V[self.nodeC_Vi_to]
			Icr = (Vcr_to - Vcr_secondary) * self.G - self.B * (Vci_to - Vci_secondary) - self.Bshunt * Vci_to + \
				  self.Gshunt * Vcr_to
			Ici = (Vci_to - Vci_secondary) * self.G + self.B * (Vcr_to - Vcr_secondary) + self.Bshunt * Vcr_to + \
				  self.Gshunt * Vci_to
			self.Ic_mag = np.around(np.abs(complex(Icr, Ici)), 2)
			self.Ic_ang = np.around(np.angle(complex(Icr, Ici), deg=True), 2)
	
	def get_nodes(self):
		"""
		Find the indices of the solution vector and the values at those indices.
		:param node: The vector of all system nodes.
		:param V: The solution vector.
		:return: Two dictionaries which hold the from/to phase nodes for the transmission lines.
		"""
		if self.casetype == 1:
			if self.stamp_dual:
				nodes_from = {
					'VR': [],
					'VI': [],
					'LR': [],
					'LI': []
				}
	
				nodes_to = {
					'VR': [],
					'VI': [],
					'LR': [],
					'LI': []
				}

				nodes_intermediate = {
					'VR POS SEC': [],
					'VI POS SEC': [],
					'VR PRIM': [],
					'VI PRIM': [],
					'LR POS SEC': [],
					'LI POS SEC': [],
					'LR PRIM': [],
					'LI PRIM': []
			  	}
			else:
				nodes_from = {
					'VR': [],
					'VI': []
				}
	
				nodes_to = {
					'VR': [],
					'VI': []
				} 

				nodes_intermediate = {
					'VR POS SEC': [],
					'VI POS SEC': [],
					'VR PRIM': [],
					'VI PRIM': []
			  	}
			print(self.phases)
			if 1 == self.phases:
				nodeA_Vr_from = self.nodeA_Vr_from	 
				nodeA_Vi_from = self.nodeA_Vi_from

				nodes_from['VR'].append(nodeA_Vr_from)
				nodes_from['VI'].append(nodeA_Vi_from)
				
				nodeA_Vr_to = self.nodeA_Vr_to
				nodeA_Vi_to = self.nodeA_Vi_to

				nodes_to['VR'].append(nodeA_Vr_to)
				nodes_to['VI'].append(nodeA_Vi_to)
				
				# Collect the secondary positive terminal nodes
				nodes_intermediate['VR POS SEC'].append(self.nodeA_Vr_pos_secondary)
				nodes_intermediate['VI POS SEC'].append(self.nodeA_Vi_pos_secondary)

				nodes_intermediate['VR PRIM'].append(self.nodeA_Vr_primary)
				nodes_intermediate['VI PRIM'].append(self.nodeA_Vi_primary)
		
				if self.stamp_dual:
					nodeA_Lr_from =  self.nodeA_Lr_from
					nodeA_Li_from =  self.nodeA_Li_from

					nodes_from['LR'].append(nodeA_Lr_from)
					nodes_from['LI'].append(nodeA_Li_from)

					nodeA_Lr_to =  self.nodeA_Lr_to
					nodeA_Li_to =  self.nodeA_Li_to

					nodes_to['LR'].append(nodeA_Lr_to)
					nodes_to['LI'].append(nodeA_Li_to)
					
					# Collect the secondary positive terminal nodes
					nodes_intermediate['LR POS SEC'].append(self.nodeA_Lr_pos_secondary)
					nodes_intermediate['LI POS SEC'].append(self.nodeA_Li_pos_secondary)

					nodes_intermediate['LR PRIM'].append(self.nodeA_Lr_primary)
					nodes_intermediate['LI PRIM'].append(self.nodeA_Li_primary)
			
			if 2 == self.phases:
				nodeB_Vr_from = self.nodeB_Vr_from
				nodeB_Vi_from = self.nodeB_Vi_from

				nodes_from['VR'].append(nodeB_Vr_from)
				nodes_from['VI'].append(nodeB_Vi_from)

				nodeB_Vr_to = self.nodeB_Vr_to
				nodeB_Vi_to = self.nodeB_Vi_to

				nodes_to['VR'].append(nodeB_Vr_to)
				nodes_to['VI'].append(nodeB_Vi_to)
				
				# Collect the secondary positive terminal nodes
				nodes_intermediate['VR POS SEC'].append(self.nodeB_Vr_pos_secondary)
				nodes_intermediate['VI POS SEC'].append(self.nodeB_Vi_pos_secondary)

				nodes_intermediate['VR PRIM'].append(self.nodeB_Vr_primary)
				nodes_intermediate['VI PRIM'].append(self.nodeB_Vi_primary)
				
				if self.stamp_dual:
					nodeB_Lr_from = self.nodeB_Lr_from
					nodeB_Li_from = self.nodeB_Li_from
					
					nodes_from['LR'].append(nodeB_Lr_from)
					nodes_from['LI'].append(nodeB_Li_from)

					nodeB_Lr_to = self.nodeB_Lr_to
					nodeB_Li_to = self.nodeB_Li_to

					nodes_to['LR'].append(nodeB_Lr_to)
					nodes_to['LI'].append(nodeB_Li_to)

					# Collect the secondary positive terminal nodes
					nodes_intermediate['LR POS SEC'].append(self.nodeB_Lr_pos_secondary)
					nodes_intermediate['LI POS SEC'].append(self.nodeB_Li_pos_secondary)

					nodes_intermediate['LR PRIM'].append(self.nodeB_Lr_primary)
					nodes_intermediate['LI PRIM'].append(self.nodeB_Li_primary)
			
			if 4 == self.phases:
				nodeC_Vr_from = self.nodeC_Vr_from
				nodeC_Vi_from = self.nodeC_Vi_from

				nodes_from['VR'].append(nodeC_Vr_from)
				nodes_from['VI'].append(nodeC_Vi_from)

				nodeC_Vr_to = self.nodeC_Vr_to
				nodeC_Vi_to = self.nodeC_Vi_to

				nodes_to['VR'].append(nodeC_Vr_to)
				nodes_to['VI'].append(nodeC_Vi_to)

				# Collect the secondary positive terminal nodes
				nodes_intermediate['VR POS SEC'].append(self.nodeC_Vr_pos_secondary)
				nodes_intermediate['VI POS SEC'].append(self.nodeC_Vi_pos_secondary)

				nodes_intermediate['VR PRIM'].append(self.nodeC_Vr_primary)
				nodes_intermediate['VI PRIM'].append(self.nodeC_Vi_primary)
				if self.stamp_dual:
					nodeC_Lr_from = self.nodeC_Lr_from
					nodeC_Li_from = self.nodeC_Li_from

					nodes_from['LR'].append(nodeC_Lr_from)
					nodes_from['LI'].append(nodeC_Li_from)

					nodeC_Lr_to = self.nodeC_Lr_to
					nodeC_Li_to = self.nodeC_Li_to

					nodes_to['LR'].append(nodeC_Lr_to)
					nodes_to['LI'].append(nodeC_Li_to)

					# Collect the secondary positive terminal nodes
					nodes_intermediate['LR POS SEC'].append(self.nodeC_Lr_pos_secondary)
					nodes_intermediate['LI POS SEC'].append(self.nodeC_Li_pos_secondary)

					nodes_intermediate['LR PRIM'].append(self.nodeC_Lr_primary)
					nodes_intermediate['LI PRIM'].append(self.nodeC_Li_primary)
			
			nodeN_Vr_from = self.nodeN_Vr_from
			nodeN_Vi_from = self.nodeN_Vi_from 
			
			nodes_from['VR'].append(nodeN_Vr_from)
			nodes_from['VI'].append(nodeN_Vi_from)

			if self.stamp_dual:
					nodeN_Lr_from = self.nodeN_Lr_from
					nodeN_Li_from = self.nodeN_Li_from

					nodes_from['LR'].append(nodeN_Lr_from)
					nodes_from['LI'].append(nodeN_Li_from)

			wyeSecondarySet = {
			self._WYE_WYE, self._DELTA_GWYE, self._GWYE_GWYE,
			self._DELTA_WYE
			}
			
			if self.connect_type in wyeSecondarySet:
				nodeN_Vr_to = self.nodeN_Vr_to
				nodeN_Vi_to = self.nodeN_Vi_to

				nodes_to['VR'].append(nodeN_Vr_to)
				nodes_to['VI'].append(nodeN_Vi_to)
				if self.stamp_dual:					
					nodeN_Lr_to = self.nodeN_Lr_to
					nodeN_Li_to = self.nodeN_Lr_to

					nodes_to['LR'].append(nodeN_Lr_to)
					nodes_to['LI'].append(nodeN_Li_to)

		return nodes_from, nodes_to, nodes_intermediate


