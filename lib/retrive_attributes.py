"""Retrieve the properties of power grid objects.

Author(s): Naeem Turner-Bandele
Created Date: 06-30-2020
Updated Date: 10-18-2020
Email: nturnerb@cmu.edu
Status: Development

Retrieves the properties of power grid objects. If a necessary property does not exist or is set to none,
then it uses default values.

"""

import re
import copy
import numpy as np

from classes.GlobalVars import _PHASE, _XFMR, _nodes_key


def convert_to_str(input_list, seperator):
	"""Convert list to string.

	   A helper method used to convert DiTTo's phase lists to a string SUGAR can use.

	   Args:
		   input_list:
			 A list.
		   seperator:
			 How the list elements are seperated when converted to a string.
	   Returns:
		   A string.
	"""
	final_str = seperator.join(filter(None, input_list))
	return final_str


def retrive_attributes(obj, obj_type, casetype=0, tap_controls=None):
	"""Retrieves the attributes of an object within the casedata dictionary.

		Retrieves properties of a power grid object to use in the initialization process.

		Args:
		  obj:
			A power grid element.Type dictionary.
		  obj_type:
			Specifies the type of power grid element that the object corresponds to. Type string.
		  casetype:
			Used to differentiate between GridLab-D and OpenDSS cases.

		Returns:
		  The object's relevant properties as variables.

		"""
	# Create filters to remove unecessary text
	remove_special = re.compile(r'[^A-Za-z0-9]+')
	remove_loadStr = re.compile(r'load\_')
	remove_node = re.compile(r'node\:')
	remove_load = re.compile(r'load\:')
	obj_strips = ('node:', 'load:', 'load_')

	# create phase index and specify all possible phases
	all_phases = ['A', 'B', 'C', 'N']
	_phase_idx = {'A': 0, 'B': 1, 'C': 2, 'D': 3}

	# check if the object has a parent element
	parent = getattr(obj, 'connecting_element', None)
	if parent:
		_id = parent
		_nodes_key[obj.name] = parent
	else:
		_id = getattr(obj, 'name')

	# set the objects id and strip any irrelevant text from the id
	_id = _id.strip('"')
	_id = remove_loadStr.sub('', _id)

	# determine bus type
	bustype = 3 if hasattr(
		obj, 'is_sourcebus') and getattr(obj, 'is_sourcebus') is True else 1

	# set nominal voltage
	nom_voltage = getattr(obj, 'nominal_voltage', None) if hasattr(
		obj, 'nominal_voltage') else None

	# determine if the object is a power delivery element that has a from node, if so, set the from property
	from_node = getattr(obj, 'from_element', None)
	if hasattr(obj, 'from_element'):

		if from_node in _nodes_key:
			from_node = _nodes_key[from_node]

		from_node = from_node.strip('"')
		from_node = remove_node.sub('', from_node)
		from_node = remove_load.sub('', from_node)
		if from_node in obj_strips:
			from_node = [
				from_node.strip(item)
				for item in obj_strips
				if item in from_node
			].pop()
	else:
		from_node = None

	# determine if the object is a power delivery element that has a to node, if so, set the to property
	to_node = getattr(obj, 'to_element', None)
	if hasattr(obj, 'to_element'):
		if to_node in _nodes_key:
			to_node = _nodes_key[to_node]
		to_node = to_node.strip('"')
		to_node = remove_node.sub('', to_node) if hasattr(
			obj, 'to_element') else None
		to_node = remove_load.sub('', to_node) if hasattr(
			obj, 'to_element') else None
		if to_node in obj_strips:
			to_node = [
				to_node.strip(item) for item in obj_strips if item in to_node
			].pop()
	else:
		to_node = None

	# # Retrieve the objects remaining properties depending on the object type # #
	if obj_type == 'load':
		phase_list = [load.phase for load in obj.phase_loads]
		phases = _PHASE[convert_to_str(phase_list, '')]
		is_delta = True if getattr(obj, 'connection_type',
								   None) == "D" else False
		load_class = 0
		powerA = complex(0, 0)
		powerB = complex(0, 0)
		powerC = complex(0, 0)
		impedanceA = complex(0, 0)
		impedanceB = complex(0, 0)
		impedanceC = complex(0, 0)
		currentA = complex(0, 0)
		currentB = complex(0, 0)
		currentC = complex(0, 0)

		# # Determine load type by model number # #
		# 1 - Constant PQ | 2 - Constant Impedance | 3 - Constant P and Quadratic Q | 4 - Exponential
		# | 5 - Constant I | 6 - Constant P and Fixed Q^2 | 7 - Constant P and Quadratic Q (varies with V^2) |
		# 8 - ZIP
		if obj.phase_loads.__len__() > 2:
			if obj.phase_loads[0].model == 8 or obj.phase_loads[1].model == 8 or obj.phase_loads[2].model == 8:
				load_type = "ZIP"
			elif obj.phase_loads[0].model == 4 or obj.phase_loads[1].model == 4 or obj.phase_loads[2].model == 4:
				load_type = "Exponential"
			else:
				load_type = "Constant"
		elif obj.phase_loads.__len__() > 1 and obj.phase_loads.__len__() < 3:
			if obj.phase_loads[0].model == 8 or obj.phase_loads[1].model == 8:
				load_type = "ZIP"
			elif obj.phase_loads[0].model == 4 or obj.phase_loads[1].model == 4:
				load_type = "Exponential"
			else:
				load_type = "Constant"
		else:
			if obj.phase_loads[0].model == 8:
				load_type = "ZIP"
			elif obj.phase_loads[0].model == 4:
				load_type = "Exponential"
			else:
				load_type = "Constant"

		percent_Z = [] if load_type == "ZIP" else None
		percent_I = [] if load_type == "ZIP" else None
		percent_P = [] if load_type == "ZIP" else None

		CVRwatts = []
		CVRvars = []

		CVRwatts_phaseA = 0.0
		CVRwatts_phaseB = 0.0
		CVRwatts_phaseC = 0.0

		CVRvars_phaseA = 0.0
		CVRvars_phaseB = 0.0
		CVRvars_phaseC = 0.0

		powerAB = getattr(obj, 'p_AB', complex(0, 0))
		powerBC = getattr(obj, 'p_BC', complex(0, 0))
		powerCA = getattr(obj, 'p_CA', complex(0, 0))

		# Check the node connect type
		# 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
		if is_delta:
			connection_type = 0
		else:
			connection_type = 1

		for load in obj.phase_loads:
			powerA = complex(load.p, load.q) if load.phase == 'A' else powerA
			powerB = complex(load.p, load.q) if load.phase == 'B' else powerB
			powerC = complex(load.p, load.q) if load.phase == 'C' else powerC
			impedanceA = load.z if load.phase == 'A' else impedanceA
			impedanceB = load.z if load.phase == 'B' else impedanceB
			impedanceC = load.z if load.phase == 'C' else impedanceC
			currentA = load.i if load.phase == 'A' else currentA
			currentB = load.i if load.phase == 'B' else currentB
			currentC = load.i if load.phase == 'C' else currentC

			if load.model == 4:
				CVRwatts_phaseA = load.cvrwatts if load.phase == 'A' else CVRwatts_phaseA
				CVRwatts_phaseB = load.cvrwatts if load.phase == 'B' else CVRwatts_phaseB
				CVRwatts_phaseC = load.cvrwatts if load.phase == 'C' else CVRwatts_phaseC

				CVRvars_phaseA = load.cvrvars if load.phase == 'A' else CVRvars_phaseA
				CVRvars_phaseB = load.cvrvars if load.phase == 'B' else CVRvars_phaseB
				CVRvars_phaseC = load.cvrvars if load.phase == 'C' else CVRvars_phaseC

				CVRwatts = [CVRwatts_phaseA, CVRwatts_phaseB, CVRwatts_phaseC]
				CVRvars = [CVRvars_phaseA, CVRvars_phaseB, CVRvars_phaseC]

			if load.model == 8:
				percent_Z.append(
						[load.ppercentimpedance, load.qpercentimpedance])
				percent_I.append([load.ppercentcurrent, load.qpercentcurrent])
				percent_P.append([load.ppercentpower, load.qpercentpower])

		return _id, phases, parent, nom_voltage, is_delta, load_class, powerA, powerB, powerC, \
			   impedanceA, impedanceB, impedanceC, currentA, currentB, currentC, load_type, percent_Z, \
			   percent_I, percent_P, powerAB, powerBC, powerCA, connection_type, CVRwatts, CVRvars

	elif obj_type == 'triplex_load':
		phase_list = [phase.default_value for phase in obj.phases]
		phases = _PHASE[convert_to_str(phase_list, '')]
		current1 = getattr(obj, "current1", None)
		current2 = getattr(obj, "current2", None)
		current12 = getattr(obj, "current12", None)
		currentN = getattr(obj, "currentN", None)
		impedance1 = getattr(obj, "impedance1", None)
		impedance2 = getattr(obj, "impedance2", None)
		impedance12 = getattr(obj, "impedance12", None)
		power1 = getattr(obj, "power1", None)
		power2 = getattr(obj, "power2", None)
		power12 = getattr(obj, "power12", None)
		shunt1 = getattr(obj, "shunt1", None)
		shunt12 = getattr(obj, "shunt12", None)
		shunt2 = getattr(obj, "shunt2", None)
		connection_type = 2

		return _id, parent, phases, current1, current12, current2, currentN, impedance1, impedance2, \
			   impedance12, nom_voltage, power1, power12, power2, shunt1, shunt12, shunt2, connection_type

	elif obj_type == 'SWING':
		phase_list = [phase.default_value for phase in obj.phases]
		phases = _PHASE[convert_to_str(phase_list, '')]
		is_triplex = getattr(obj, 'is_triplex', None) if hasattr(
			obj, 'is_triplex') else False

		phase_angle = getattr(obj, 'phase_angle', 0)

		voltage_setpoint = np.around(obj.per_unit, 3)
		voltage_A = getattr(obj, 'voltage_A', None)
		voltage_B = getattr(obj, 'voltage_B', None)
		voltage_C = getattr(obj, 'voltage_C', None)

		# Check the node connect type
		# 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
		if not is_triplex:
			if phases == _PHASE["ABCD"]:
				connection_type = 0
			else:
				connection_type = 1
		else:
			connection_type = 2
		
		slack_parameters = {}
		if casetype != 0:
			slack_parameters['file type'] = 'DSS'

			# Open DSS uses a thevinin equivalent, not an infinite bus, for the slack bus
			pos_sequence = obj.positive_sequence_impedance
			zero_sequence = obj.zero_sequence_impedance

			z_mutual = (zero_sequence - pos_sequence)/3
			z_series = (zero_sequence - 2*z_mutual)

			impedance_matrix = np.array([[z_series, z_mutual, z_mutual],[z_mutual, z_series, z_mutual],[z_mutual, z_mutual, z_series]])
			slack_parameters['impedance_matrix'] = impedance_matrix
		
		else:
			slack_parameters['file type'] = 'GLD'

		return _id, phases, nom_voltage, phase_angle, bustype, is_triplex, connection_type, voltage_setpoint, voltage_A, \
			   voltage_B, voltage_C, slack_parameters

	elif obj_type == 'transformer':
		winding1 = obj.windings[0]
		winding2 = obj.windings[1]
		if casetype == 0:
			phase_list = [phase.default_value for phase in obj.phases]
		else:
			phase_list = [phase.phase for phase in winding1.phase_windings]
		is_center_tap = obj.is_center_tap
		phases = _PHASE[convert_to_str(phase_list, '')]
		primary_voltage = winding1.nominal_voltage
		secondary_voltage = winding2.nominal_voltage
		shunt_impedance = getattr(obj, 'shunt_impedance', None)
		winding3 = obj.windings[2] if is_center_tap else None

		if winding1.connection_type == "Y" and winding2.connection_type == "Y":
# 			if phases in [_PHASE['AN'], _PHASE['BN'], _PHASE['CN']]:
# 				connection_type = _XFMR["SINGLE_PHASE"]
# 			else:
			connection_type = _XFMR["WYE_WYE"]
		elif winding1.connection_type == "D" and winding2.connection_type == "D":
			connection_type = _XFMR["DELTA_DELTA"]
		elif winding1.connection_type == "D" and winding2.connection_type == "Y":
			connection_type = _XFMR["DELTA_GWYE"]
		elif winding1.connection_type == "Y" and winding2.connection_type == "D":
			connection_type = _XFMR["WYE_DELTA"]
		elif winding1.connection_type == "SINGLE_PHASE":
			connection_type = _XFMR["SINGLE_PHASE"]
		elif is_center_tap:
			connection_type = _XFMR["SINGLE_PHASE_CENTER_TAPPED"]
		else:
			connection_type = _XFMR["UNKNOWN"]

		install_type = obj.install_type
		is_grounded = winding2.is_grounded
		if is_center_tap:
			reactance = obj.reactances
			resistance = [
				winding1.resistance, winding2.resistance, winding3.resistance
			]
			power_rating = winding1.rated_power
		else:
			reactance = obj.reactances.pop(
			) if casetype == 0 else obj.reactances.pop() / 100
			resistances = winding1.resistance + winding2.resistance
			resistance = resistances if casetype == 0 else resistances / 100
			if resistance is None:
				resistance = 0
			power_rating = winding1.rated_power

		loadloss = obj.loadloss
		noload_loss = obj.noload_loss
		phase_shift = obj.phase_shift if casetype == 0 else 0

		return _id, phases, from_node, to_node, nom_voltage, primary_voltage, secondary_voltage, connection_type, \
			   install_type, is_center_tap, is_grounded, reactance, resistance, shunt_impedance, loadloss, \
			   noload_loss, phase_shift, power_rating

	elif obj_type == 'line':
		line_parameters = {}
		if casetype == 0:
			phase_list = [phase.default_value for phase in obj.phases]
			conductor_phases = convert_to_str(
				[wire.phase for wire in obj.wires], '')
			line_parameters['file type'] = 'GLD'
		else:
			phase_list = [wire.phase for wire in obj.wires]
			conductor_phases = None
			line_parameters['file type'] = 'DSS'

		phases = _PHASE[convert_to_str(phase_list, '')]

		unused_phases = set(all_phases).difference(phase_list)
		unused_phases = [ele for ele in unused_phases
						] if casetype != 0 else None
		spacing = copy.deepcopy(obj.spacing) if hasattr(obj, 'spacing') else None
		line_type = getattr(obj, 'line_type') if getattr(
			obj, 'line_type') is not None else 'overhead'
		line_parameters["line type"] = line_type
		# EF Note - removed change to length
		length = obj.length / (5280) if casetype == 0 else obj.length  # convert to miles if glm

		# If DSS, the length is in meters. we want length in miles
		impedance_matrix = np.array(obj.impedance_matrix)
		capacitance_matrix = np.array(obj.capacitance_matrix)
		diameter = None
		outer_diameter = None
		neutral_diameter = None
		neutral_strands = None
		conductor_diameter = None
		conductor_resistance = None
		conductor_gmr = None
		neutral_resistance = None
		neutral_gmr = None
		insulation_thickness = None

		# if this is a gridlabd case, then retrieve the relevant line properties because ditto incorrectly calculates
		# the impedance and capacitance matrices for gridlabd cases
		if casetype == 0:
			line_parameters['given Z'] = False
			if spacing:
				if hasattr(spacing, '_distance_AN'):
					spacing._distance_AN = float(spacing._distance_AN)/12 
				if hasattr(spacing, '_distance_BN'):
					spacing._distance_BN = float(spacing._distance_BN)/12 
				if hasattr(spacing, '_distance_CN'):
					spacing._distance_CN = float(spacing._distance_CN)/12 
				if hasattr(spacing, '_distance_AB'):
					spacing._distance_AB = float(spacing._distance_AB)/12 
				if hasattr(spacing, '_distance_AC'):
					spacing._distance_AC = float(spacing._distance_AC)/12 
				if hasattr(spacing, '_distance_BC'):
					spacing._distance_BC = float(spacing._distance_BC)/12 
			
			line_parameters['spacing'] = spacing
			
			if line_type == 'overhead':
				diameter_dict = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'N': 0.0}
				for wire in obj.wires:
					if wire.diameter is not None:
						# EF note - removed unit converstion
						diameter_dict[
							'A'] = wire.diameter if wire.phase == 'A' else diameter_dict[
								'A']
						diameter_dict[
							'B'] = wire.diameter if wire.phase == 'B' else diameter_dict[
								'B']
						diameter_dict[
							'C'] = wire.diameter if wire.phase == 'C' else diameter_dict[
								'C']
						diameter_dict[
							'N'] = wire.diameter if wire.phase == 'N' else diameter_dict[
								'N']
						# diameter_dict[
						#     'A'] = wire.diameter * 39.3701 if wire.phase == 'A' else diameter_dict[
						#         'A']
						# diameter_dict[
						#     'B'] = wire.diameter * 39.3701 if wire.phase == 'B' else diameter_dict[
						#         'B']
						# diameter_dict[
						#     'C'] = wire.diameter * 39.3701 if wire.phase == 'C' else diameter_dict[
						#         'C']
						# diameter_dict[
						#     'N'] = wire.diameter * 39.3701 if wire.phase == 'N' else diameter_dict[
						#         'N']
					else:
						diameter_dict[
							'A'] = 0.5 if wire.phase == 'A' else diameter_dict[
								'A']
						diameter_dict[
							'B'] = 0.5 if wire.phase == 'B' else diameter_dict[
								'B']
						diameter_dict[
							'C'] = 0.5 if wire.phase == 'C' else diameter_dict[
								'C']
						diameter_dict[
							'N'] = 0.5 if wire.phase == 'N' else diameter_dict[
								'N']
				diameter = list(diameter_dict.values())
				line_parameters['diameter'] = diameter
				line_parameters['line type'] = 'overhead'
			
			if line_type == 'triplex':
				diameter_dict = {'A': 0.0, 'B': 0.0, 'N': 0.0}
				insulation_thickness_dict = {'A': 0.0, 'B': 0.0, 'N': 0.0}
				for wire in obj.wires:
					if wire.diameter is not None:
						# EF note - removed unit converstion
						diameter_dict[
							'A'] = wire.diameter if wire.phase == 'A' else diameter_dict[
								'A']
						diameter_dict[
							'B'] = wire.diameter if wire.phase == 'B' else diameter_dict[
								'B']
						diameter_dict[
							'N'] = wire.diameter if wire.phase == 'N' else diameter_dict[
								'N']
						# diameter_dict[
						#     'A'] = wire.diameter * 39.3701 if wire.phase == 'A' else diameter_dict[
						#         'A']
						# diameter_dict[
						#     'B'] = wire.diameter * 39.3701 if wire.phase == 'B' else diameter_dict[
						#         'B']
						# diameter_dict[
						#     'C'] = wire.diameter * 39.3701 if wire.phase == 'C' else diameter_dict[
						#         'C']
						# diameter_dict[
						#     'N'] = wire.diameter * 39.3701 if wire.phase == 'N' else diameter_dict[
						#         'N']
					else:
						diameter_dict[
							'A'] = 0.5 if wire.phase == 'A' else diameter_dict[
								'A']
						diameter_dict[
							'B'] = 0.5 if wire.phase == 'B' else diameter_dict[
								'B']
						diameter_dict[
							'N'] = 0.5 if wire.phase == 'N' else diameter_dict[
								'N']
					
					if wire.insulation_thickness is not None:
						# EF note - removed unit converstion
						insulation_thickness_dict[
							'A'] = wire.insulation_thickness if wire.phase == 'A' else insulation_thickness_dict[
								'A']
						insulation_thickness_dict[
							'B'] = wire.insulation_thickness if wire.phase == 'B' else insulation_thickness_dict[
								'B']
						insulation_thickness_dict[
							'N'] = wire.insulation_thickness if wire.phase == 'N' else insulation_thickness_dict[
								'N']
					else:
						insulation_thickness_dict[
							'A'] = 0.05 if wire.phase == 'A' else insulation_thickness_dict[
								'A']
						insulation_thickness_dict[
							'B'] = 0.05 if wire.phase == 'B' else insulation_thickness_dict[
								'B']
						insulation_thickness_dict[
							'N'] = 0.05 if wire.phase == 'N' else insulation_thickness_dict[
								'N']

				diameter = list(diameter_dict.values())
				insulation_thickness = list(insulation_thickness_dict.values())
				line_parameters['diameter'] = diameter
				line_parameters['insulation_thickness'] = insulation_thickness
				line_parameters['line type'] = 'triplex'

			if line_type == 'underground':
				impedance_matrix = None
				# set parameter values to defaults initially. defaults taken from old parser
				outer_diameter_dict = {
					'A': 1.56,
					'B': 1.56,
					'C': 1.56,
					'N': 0.0
				}
				neutral_diameter_dict = {
					'A': 0.0808,
					'B': 0.0808,
					'C': 0.0808,
					'N': 0.00
				}
				neutral_strands_dict = {
					'A': 16.0,
					'B': 16.0,
					'C': 16.0,
					'N': 0.0
				}
				conductor_diameter_dict = {
					'A': 0.814,
					'B': 0.814,
					'C': 0.814,
					'N': 0.0
				}
				conductor_resistance_dict = {
					'A': 0.0,
					'B': 0.0,
					'C': 0.0,
					'N': 0.0
				}
				conductor_gmr_dict = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'N': 0.0}
				neutral_resistance_dict = {
					'A': 9.3747,
					'B': 9.3747,
					'C': 9.3747,
					'N': 0.0
				}
				neutral_gmr_dict = {
					'A': 0.00262,
					'B': 0.00262,
					'C': 0.00262,
					'N': 0.0
				}

				for wire in obj.wires:
					if getattr(wire, 'outer_diameter', None):
						outer_diameter_val = getattr(wire, 'outer_diameter')
					else:
						outer_diameter_val = 1.56
					if getattr(wire, 'concentric_neutral_diameter', None):
						neutral_diameter_val = getattr(
							wire, 'concentric_neutral_diameter')
					else:
						neutral_diameter_val = 0.0808
					if getattr(wire, 'concentric_neutral_nstrand', None):
						neutral_strands_val = getattr(
							wire, 'concentric_neutral_nstrand')
					else:
						neutral_strands_val = int(16)
					if getattr(wire, 'conductor_diameter', None):
						conductor_diameter_val = getattr(
							wire, 'conductor_diameter')
					else:
						conductor_diameter_val = 0.814
					if getattr(wire, 'concentric_neutral_resistance', None):
						# EF note - removed unit conversion
						neutral_resistance_val = getattr(
							wire, 'concentric_neutral_resistance') 
						# neutral_resistance_val = getattr(
						#     wire, 'concentric_neutral_resistance') * 1609.34
					else:
						neutral_resistance_val = 9.3747
					if getattr(wire, 'concentric_neutral_gmr', None):
						# EF note - removed unit conversion
						neutral_gmr_val = getattr(
							wire, 'concentric_neutral_gmr')
						# neutral_gmr_val = getattr(
						#     wire, 'concentric_neutral_gmr') * 3.28084
					else:
						neutral_gmr_val = 0.00262
					outer_diameter_dict[
						'A'] = outer_diameter_val if wire.phase == 'A' else outer_diameter_dict[
							'A']
					outer_diameter_dict[
						'B'] = outer_diameter_val if wire.phase == 'B' else outer_diameter_dict[
							'B']
					outer_diameter_dict[
						'C'] = outer_diameter_val if wire.phase == 'C' else outer_diameter_dict[
							'C']
					outer_diameter_dict[
						'N'] = 0.0 if wire.phase == 'N' else outer_diameter_dict[
							'N']

					neutral_diameter_dict['A'] = neutral_diameter_val if wire.phase == 'A' else \
					 neutral_diameter_dict['A']
					neutral_diameter_dict['B'] = neutral_diameter_val if wire.phase == 'B' else \
					 neutral_diameter_dict['B']
					neutral_diameter_dict['C'] = neutral_diameter_val if wire.phase == 'C' else \
					 neutral_diameter_dict['C']
					neutral_diameter_dict['N'] = 0.0 if wire.phase == 'N' else \
					 neutral_diameter_dict['N']

					neutral_strands_dict['A'] = neutral_strands_val if wire.phase == 'A' else \
					 neutral_strands_dict['A']
					neutral_strands_dict['B'] = neutral_strands_val if wire.phase == 'B' else \
					 neutral_strands_dict['B']
					neutral_strands_dict['C'] = neutral_strands_val if wire.phase == 'C' else \
					 neutral_strands_dict['C']
					neutral_strands_dict['N'] = 0.0 if wire.phase == 'N' else \
					 neutral_strands_dict['N']

					conductor_diameter_dict['A'] = conductor_diameter_val if wire.phase == 'A' else \
					conductor_diameter_dict['A']
					conductor_diameter_dict['B'] = conductor_diameter_val if wire.phase == 'B' else \
					conductor_diameter_dict['B']
					conductor_diameter_dict['C'] = conductor_diameter_val if wire.phase == 'C' else \
					conductor_diameter_dict['C']
					conductor_diameter_dict['N'] = conductor_diameter_val if wire.phase == 'N' else \
					 conductor_diameter_dict['N']
					
					neutral_resistance_dict['A'] = neutral_resistance_val if wire.phase == 'A' else \
					 neutral_resistance_dict['A']
					neutral_resistance_dict['B'] = neutral_resistance_val if wire.phase == 'B' else \
					 neutral_resistance_dict['B']
					neutral_resistance_dict['C'] = neutral_resistance_val if wire.phase == 'C' else \
					 neutral_resistance_dict['C']
					neutral_resistance_dict['N'] = 0.0 if wire.phase == 'N' else \
					 neutral_resistance_dict['N']

					neutral_gmr_dict['A'] = neutral_gmr_val if wire.phase == 'A' else \
					 neutral_gmr_dict['A']
					neutral_gmr_dict['B'] = neutral_gmr_val if wire.phase == 'B' else \
					 neutral_gmr_dict['B']
					neutral_gmr_dict['C'] = neutral_gmr_val if wire.phase == 'C' else \
					 neutral_gmr_dict['C']
					neutral_gmr_dict['N'] = 0.0 if wire.phase == 'N' else \
					 neutral_gmr_dict['N']
				outer_diameter = list(outer_diameter_dict.values())
				neutral_diameter = list(neutral_diameter_dict.values())
				neutral_strands = list(neutral_strands_dict.values())
				conductor_diameter = list(conductor_diameter_dict.values())
				neutral_resistance = list(neutral_resistance_dict.values())
				neutral_gmr = list(neutral_gmr_dict.values())

				line_parameters['line type'] = 'underground'
				line_parameters['outer_diameter'] = outer_diameter
				line_parameters['neutral_diameter'] = neutral_diameter
				line_parameters['neutral_strands'] = neutral_strands
				line_parameters['conductor_diameter'] = conductor_diameter
				line_parameters['neutral_resistance'] = neutral_resistance
				line_parameters['neutral_gmr'] = neutral_gmr
			 
			conductor_resistance_dict = {
					'A': 0.0,
					'B': 0.0,
					'C': 0.0,
					'N': 0.0
				}  
			conductor_gmr_dict = {'A': 0.0, 'B': 0.0, 'C': 0.0, 'N': 0.0}
				
			for wire in obj.wires:     
				# EF note - removed unit conversion
				conductor_resistance_dict['A'] = wire.resistance if wire.phase == 'A' else \
				 conductor_resistance_dict['A']
				conductor_resistance_dict['B'] = wire.resistance if wire.phase == 'B' else \
				 conductor_resistance_dict['B']
				conductor_resistance_dict['C'] = wire.resistance if wire.phase == 'C' else \
				 conductor_resistance_dict['C']
				conductor_resistance_dict['N'] = wire.resistance if wire.phase == 'N' else \
				 conductor_resistance_dict['N']

					# conductor_resistance_dict['A'] = wire.resistance * 1609.34 if wire.phase == 'A' else \
					#  conductor_resistance_dict['A']
					# conductor_resistance_dict['B'] = wire.resistance * 1609.34 if wire.phase == 'B' else \
					#  conductor_resistance_dict['B']
					# conductor_resistance_dict['C'] = wire.resistance * 1609.34 if wire.phase == 'C' else \
					#  conductor_resistance_dict['C']
					# conductor_resistance_dict['N'] = wire.resistance * 1609.34 if wire.phase == 'N' else \
					#  conductor_resistance_dict['N']

					# EF note - removed unit conversion
				conductor_gmr_dict[
					'A'] = wire.gmr if wire.phase == 'A' else conductor_gmr_dict[
						'A']
				conductor_gmr_dict[
					'B'] = wire.gmr if wire.phase == 'B' else conductor_gmr_dict[
						'B']
				conductor_gmr_dict[
					'C'] = wire.gmr if wire.phase == 'C' else conductor_gmr_dict[
						'C']
				conductor_gmr_dict[
					'N'] = wire.gmr if wire.phase == 'N' else conductor_gmr_dict[
						'N']
					# conductor_gmr_dict[
					#     'A'] = wire.gmr * 3.28084 if wire.phase == 'A' else conductor_gmr_dict[
					#         'A']
					# conductor_gmr_dict[
					#     'B'] = wire.gmr * 3.28084 if wire.phase == 'B' else conductor_gmr_dict[
					#         'B']
					# conductor_gmr_dict[
					#     'C'] = wire.gmr * 3.28084 if wire.phase == 'C' else conductor_gmr_dict[
					#         'C']
					# conductor_gmr_dict[
					#     'N'] = wire.gmr * 3.28084 if wire.phase == 'N' else conductor_gmr_dict[
					#         'N']

				conductor_resistance = conductor_resistance_dict
				conductor_gmr = conductor_gmr_dict

				line_parameters['conductor_resistance'] = conductor_resistance
				line_parameters['conductor_gmr'] = conductor_gmr

		elif casetype == 1:
			line_parameters['length'] = obj.length / 1609.34 # length in m * 1 mi/1609.34 m
			if obj.capacitance_matrix is not None:
				line_parameters['given C'] = True
				line_parameters['capacitance matrix'] = np.array(obj.capacitance_matrix) * 1609.34 # Ohms/m * 1609 m/mi
			
			if obj.impedance_matrix is not None:
				line_parameters['given Z'] = True
				line_parameters['impedance matrix'] = np.array(obj.impedance_matrix) * 1609.34 # F/m * 1609 m/mi

		return _id, line_type, phases, nom_voltage, from_node, to_node, length, spacing, diameter, impedance_matrix, \
			   capacitance_matrix, outer_diameter, neutral_diameter, neutral_strands, conductor_diameter, \
			   conductor_resistance, conductor_gmr, neutral_resistance, neutral_gmr, unused_phases, conductor_phases, insulation_thickness, line_parameters

	elif obj_type == 'capacitor':
		if casetype == 0:
			phases = _PHASE[remove_special.sub('', obj.phases)]
		else:
			phase_list = [phase.phase for phase in obj.phase_capacitors]
			phases = _PHASE[convert_to_str(phase_list, '')]
		connect_type = 0 if obj.connection_type == 'WYE' or obj.connection_type == "Y" else 1
		capacitorA = 0
		capacitorB = 0
		capacitorC = 0
		delay = obj.delay
		pt_phase = obj.pt_phase
		upper_limit = obj.high
		lower_limit = obj.low

		for capacitor in obj.phase_capacitors:
			capacitorA = capacitor.var if capacitor.phase == 'A' else capacitorA
			capacitorB = capacitor.var if capacitor.phase == 'B' else capacitorB
			capacitorC = capacitor.var if capacitor.phase == 'C' else capacitorC

		return _id, parent, phases, nom_voltage, connect_type, capacitorA, capacitorB, capacitorC, pt_phase, delay, \
			   upper_limit, lower_limit

	elif obj_type == 'fuse':
		phase_A_state = 1
		phase_B_state = 1
		phase_C_state = 1

		try:
			phase_list = [phase.default_value for phase in obj.phases]
		except AttributeError:
			phase_list = [wire.phase for wire in obj.wires]

		phases = _PHASE[convert_to_str(phase_list, '')]
		for fuse in obj.wires:
			phase_A_state = 0 if fuse.is_open and fuse.phase == 'A' else phase_A_state
			phase_B_state = 0 if fuse.is_open and fuse.phase == 'B' else phase_B_state
			phase_C_state = 0 if fuse.is_open and fuse.phase == 'C' else phase_C_state

		return _id, from_node, to_node, phases, phase_A_state, phase_B_state, phase_C_state

	elif obj_type == 'switch':
		switch_parameters = {}
		phase_A_state = 1 if getattr(obj, 'status', None) == "CLOSED" else 0
		phase_B_state = 1 if getattr(obj, 'status', None) == "CLOSED" else 0
		phase_C_state = 1 if getattr(obj, 'status', None) == "CLOSED" else 0
		if casetype == 0:
			phase_list = [phase.default_value for phase in obj.phases]
		else:
			phase_list = [wire.phase for wire in obj.wires]
		phases = _PHASE[convert_to_str(phase_list, '')]
		for wire in obj.wires:
			phase_A_state = 1 if not wire.is_open and wire.phase == 'A' else phase_A_state
			phase_B_state = 1 if not wire.is_open and wire.phase == 'B' else phase_B_state
			phase_C_state = 1 if not wire.is_open and wire.phase == 'C' else phase_C_state
		
		switch_parameters['impedance matrix'] = obj.impedance_matrix
		switch_parameters['capacitance matrix'] = obj.capacitance_matrix 

		return _id, from_node, to_node, phases, phase_A_state, phase_B_state, phase_C_state, switch_parameters

	elif obj_type == 'regulator':
		tap_status = "Fixed" if tap_controls["Fixed"] else "Variable"
		band_center = getattr(obj, 'bandcenter', 120.0)
		bandwidth = getattr(obj, 'bandwidth', 2.0)
		pt_ratio = getattr(obj, 'pt_ratio')
		pt_phase = getattr(obj, 'pt_phase', None)
		ct_ratio = getattr(obj, 'ct_ratio', None)
		if ct_ratio:
			ct_prim = getattr(obj, 'ct_prim') if getattr(
				obj, 'ct_prim') is not None else 5 * ct_ratio
		else:
			ct_prim = getattr(obj, 'ct_prim') if getattr(
				obj, 'ct_prim') is not None else 5 * 140

		ct_phase = getattr(obj, 'ct_phase', None)
		highstep = getattr(obj, 'highstep', 16.0)
		lowstep = -getattr(obj, 'lowstep', 16.0)
		winding1 = obj.windings[0]
		winding2 = obj.windings[1]
		connected_transformer = getattr(obj, 'connected_transformer', None)
		connected_winding = getattr(obj, 'winding', None)

		try:
			voltage_setpoint = np.around(obj.per_unit, 3)
		except AttributeError:
			voltage_setpoint = 1.0

		if casetype == 0:
			phase_list = [phase.default_value for phase in obj.phases]
		else:
			phase_list = [phase.phase for phase in winding1.phase_windings]
		phases = _PHASE[convert_to_str(phase_list, '')]

		if connected_winding == 1:
			nom_voltage = winding1.nominal_voltage
			calculated_rated_power = np.around(
				ct_prim * nom_voltage * 3 * 10**-3, -2)
			rated_power = winding1.rated_power if winding1.rated_power is not None else calculated_rated_power
			compensator_resistance = getattr(winding1.phase_windings[0],
											 'compensator_r', 1.0)
			compensator_reactance = getattr(winding1.phase_windings[0],
											'compensator_x', 1.0)
		else:
			nom_voltage = winding2.nominal_voltage
			calculated_rated_power = np.around(
				ct_prim * nom_voltage * 3 * 10**-3, -2)
			rated_power = winding2.rated_power if winding2.rated_power is not None else calculated_rated_power
			compensator_resistance = getattr(winding2.phase_windings[0],
											 'compensator_r', 1.0)
			compensator_reactance = getattr(winding2.phase_windings[0],
											'compensator_x', 1.0)
		if casetype == 0:
			reg_type = winding1.reg_type if hasattr(winding1,
													'reg_type') else 'A'
		else:
			reg_type = 'B'

		if winding1.connection_type == "Y" and winding2.connection_type == "Y":
			connection_type = 1
		elif winding1.connection_type == "D" and winding2.connection_type == "D":
			connection_type = 2
		elif winding1.connection_type == "D" and winding2.connection_type == "Y":
			connection_type = 1
		else:
			connection_type = -1

		if pt_ratio is None:
			pt_ratio = np.around(winding1.nominal_voltage / 120, 1)

		if ct_prim is None:
			ct_prim = 700

		tapA = int(0)
		tapB = int(0)
		tapC = int(0)

		aR_max = getattr(obj, "maxtap")
		aR_min = getattr(obj, "mintap")
		num_tap_steps = getattr(obj, "numtaps")

		for winding in winding2.phase_windings:
			tapA = winding.tap_position if winding.phase == 'A' and winding.tap_position is not None else tapA
			tapB = winding.tap_position if winding.phase == 'B' and winding.tap_position is not None else tapB
			tapC = winding.tap_position if winding.phase == 'C' and winding.tap_position is not None else tapC

		if casetype != 0:
			if type(tapA) is float or type(tapB) is float or type(
					tapC) is float:
				# Change default taps if they are provided as per unit values
				# Use Type A regulator equations to calculate new tap starting points
				if phases & 0x1 == 1:
					tapA = int((tapA - 1) / 0.00625) if tapA > 1 else 0.0
				if phases & 0x2 == 2:
					tapB = int((tapB - 1) / 0.00625) if tapB > 1 else 0.0
				if phases & 0x4 == 4:
					tapC = int((tapC - 1) / 0.00625) if tapC > 1 else 0.0

		return _id, phases, from_node, to_node, connection_type, reg_type, \
			   tapA, tapB, tapC, pt_phase, ct_phase, connected_transformer, connected_winding, bustype, \
			   voltage_setpoint, nom_voltage, band_center, bandwidth, pt_ratio, ct_prim, highstep, lowstep, \
			   rated_power, compensator_resistance, compensator_reactance, tap_status, aR_max, aR_min, num_tap_steps

	elif obj_type == 'reactor':
		phase_list = [phase.phase for phase in obj.phase_reactors]
		phases = _PHASE[convert_to_str(phase_list, '')]
		connection_type = 0 if obj.connection_type == "Y" else 1
		for reactor in obj.phase_reactors:
			resistanceA = reactor.resistance if reactor.phase == 'A' else 0
			resistanceB = reactor.resistance if reactor.phase == 'B' else 0
			resistanceC = reactor.resistance if reactor.phase == 'C' else 0
			resistanceN = reactor.resistance if reactor.phase == 'N' else 0

			reactanceA = reactor.reactance if reactor.phase == 'A' else 0
			reactanceB = reactor.reactance if reactor.phase == 'B' else 0
			reactanceC = reactor.reactance if reactor.phase == 'C' else 0
			reactanceN = reactor.reactance if reactor.phase == 'N' else 0

		resistance = [resistanceA, resistanceB, resistanceC, resistanceN]
		reactance = [reactanceA, reactanceB, reactanceC, reactanceN]

		return _id, from_node, to_node, phases, connection_type, resistance, reactance

	elif obj_type == 'ibdg':
		phase_list = [phase.default_value for phase in obj.phases]
		phases = _PHASE[convert_to_str(phase_list, '')]
		is_triplex = getattr(obj, 'is_triplex', False)
		connection_type = getattr(obj, 'connection_type')
		rated_power = getattr(obj, 'rated_power')
		active_rating = getattr(obj, 'active_rating')
		reactive_rating = getattr(obj, 'reactive_rating')
		power_factor = getattr(obj, "power_factor")
		pv_max_power = getattr(obj, "pv_power_max")
		irradiance = getattr(obj, "irradiance")

		return _id, phases, nom_voltage, connection_type, rated_power, active_rating, reactive_rating, power_factor, \
			   pv_max_power, irradiance, is_triplex

	else:
		phase_list = [phase.default_value for phase in obj.phases]
		phases = _PHASE[convert_to_str(phase_list, '')]
		is_triplex = getattr(obj, 'is_triplex', None) if hasattr(
			obj, 'is_triplex') else False

		try:
			voltage_setpoint = np.around(obj.per_unit, 3)
		except AttributeError:
			voltage_setpoint = 1.0

		# Check the node connect type
		# 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
		if not is_triplex:
			if phases == _PHASE["ABCD"]:
				connection_type = 0
			else:
				connection_type = 1
		else:
			connection_type = 2

		return _id, phases, nom_voltage, bustype, is_triplex, connection_type, voltage_setpoint
