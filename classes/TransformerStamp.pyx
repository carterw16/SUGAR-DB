"""Implements transformer stamps.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 10-05-2020
  Updated Date: 10-10-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Development

"""
#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

cpdef stamp_Y(i, j, val, Y_val, Y_row, Y_col, idx):
	if val == 0.:
		return Y_val, Y_row, Y_col, idx

	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return Y_val, Y_row, Y_col, idx

cpdef stamp_primary(self,
                  Ylin_val,
                  Ylin_row,
                  Ylin_col,
                  idx_Y,
				  int node_Vr_pos_from,
				  int node_Vi_pos_from,
				  int node_Vr_neg_from,
				  int node_Vi_neg_from,
				  int node_Vr_neg_secondary,
				  int node_Vi_neg_secondary,
				  int node_Vr_pos_secondary,
				  int node_Vi_pos_secondary,
				  int node_Vr_primary,
				  int node_Vi_primary,
				  int phase,
				  dual_dict):
	if self.stamp_dual:
		node_Lr_pos_from = dual_dict['Lr_from_pos']
		node_Li_pos_from = dual_dict['Li_from_pos']
		node_Lr_neg_from = dual_dict['Lr_from_neg']
		node_Li_neg_from = dual_dict['Li_from_neg']
		node_Lr_neg_secondary = dual_dict['Lr_to_neg']
		node_Li_neg_secondary = dual_dict['Li_to_neg']
		node_Lr_pos_secondary = dual_dict['Lr_pos_secondary']
		node_Li_pos_secondary = dual_dict['Li_pos_secondary']
		node_Lr_primary = dual_dict['Lr_primary']
		node_Li_primary = dual_dict['Li_primary']
	# ==============================================================================
	#         Phase ANGLE != 0
	# ==============================================================================
	if self.ang:
		node_Vr_primary_link = self.nodes_primary_link[self.all_phases[phase]][0]
		node_Vi_primary_link = self.nodes_primary_link[self.all_phases[phase]][1]
		node_Vr_secondary_link = self.nodes_secondary_link[self.all_phases[phase]][0]
		node_Vi_secondary_link = self.nodes_secondary_link[self.all_phases[phase]][1]

		# Define phase shifts
		_tr_cos = self.tr * self.cosd(self.ang)
		_tr_sin = self.tr * self.sind(self.ang)

		# # PRIMARY VOLTAGES # #
		# Real #
		# Vr_prim = Vr_pos_from - tr * Vr_sec_pos *cos(theta)  + tr * Vr_sec_neg * cos(theta) - Vr_prim_link
		# Vr_prim_link = - Vr_prim + Vr_sec_link
		# Vr_sec_link = Vr_prim_link + tr * sin(theta) * Vi_sec_pos - tr * sin(theta) * Vi_sec_neg - Vr_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected
		# the link node is connecting secondary imaginary phase shift to the rest of the primary voltages

		# Primary
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_pos_from, 1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_pos_secondary, -_tr_cos,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_neg_secondary, _tr_cos,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_primary_link, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)
		# Primary Link
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary_link, node_Vr_primary, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary_link, node_Vr_secondary_link, 1,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		# Secondary Link
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_secondary_link, node_Vr_primary_link, 1,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_secondary_link, node_Vi_pos_secondary, _tr_sin,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_secondary_link, node_Vi_neg_secondary, -_tr_sin,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_secondary_link, node_Vr_neg_from, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		# Imaginary #
		# Vi_prim = Vi_pos_from - tr * Vi_sec_pos * cos(theta)  + tr * Vi_sec_neg * cos(theta) - Vi_prim_link
		# Vi_prim_link = - Vi_prim + Vi_sec_link
		# Vi_sec_link = Vi_prim - tr * sin(theta) * Vr_sec_pos + tr * sin(theta) * Vr_sec_neg - Vi_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected
		# the link node is connecting secondary imaginary phase shift to the rest of the primary voltages

		# Primary
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_pos_from, 1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_pos_secondary, -_tr_cos,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_neg_secondary, _tr_cos,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_primary_link, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		# Primary Link
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary_link, node_Vi_primary, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary_link, node_Vi_secondary_link, 1,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		# Secondary Link
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_secondary_link, node_Vi_primary_link, 1,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_secondary_link, node_Vr_pos_secondary, -_tr_sin,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_secondary_link, node_Vr_neg_secondary, _tr_sin,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_secondary_link, node_Vi_neg_from, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		# DUAL STAMP (TRANSPOSE OF REGULAR)
		if self.stamp_dual:
			node_Lr_primary_link = self.nodes_primary_link[self.all_phases[phase]][2]
			node_Li_primary_link = self.nodes_primary_link[self.all_phases[phase]][3]
			node_Lr_secondary_link = self.nodes_secondary_link[self.all_phases[phase]][2]
			node_Li_secondary_link = self.nodes_secondary_link[self.all_phases[phase]][3]
			
			# Primary - Real Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_pos_from, node_Lr_primary, 1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_pos_secondary, node_Lr_primary,  -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_neg_secondary, node_Lr_primary, _tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary_link, node_Lr_primary, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
			# Primary Link - Real Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_primary_link, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_secondary_link, node_Lr_primary_link, 1,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			
			# Secondary Link - Real Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary_link, node_Lr_secondary_link, 1,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_pos_secondary, node_Lr_secondary_link, _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_neg_secondary, node_Lr_secondary_link, -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_neg_from, node_Lr_secondary_link, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
		
			# Primary - Imaginary Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_pos_from, node_Li_primary, 1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_pos_secondary, node_Li_primary, -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_neg_secondary, node_Li_primary, _tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary_link, node_Li_primary, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
	
			# Primary Link - Imaginary Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_primary_link, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_secondary_link, node_Li_primary_link, 1,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
	
			# Secondary Link - Imaginary Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary_link, node_Li_secondary_link, 1,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_pos_secondary, node_Li_secondary_link, -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_neg_secondary, node_Li_secondary_link, _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_neg_from, node_Li_secondary_link, -1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)

	else:
		# ==============================================================================
		#         Phase ANGLE = 0
		# ==============================================================================

		# # PRIMARY VOLTAGES # #
		# Real #
		# Vr_primary = Vr_pos_from - tr * Vr_secondary_pos  + tr * Vr_secondary_neg - Vr_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_pos_from, 1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_pos_secondary, -self.tr,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_neg_secondary, self.tr,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_primary, node_Vr_neg_from, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		# Imaginary #
		# Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_pos_from, 1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_pos_secondary, -self.tr,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_neg_from, -1, Ylin_val,
													 Ylin_row, Ylin_col, idx_Y)

		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_primary, node_Vi_neg_secondary, self.tr,
													 Ylin_val, Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			# Real Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_pos_from, node_Lr_primary, 1, Ylin_val,
															  Ylin_row, Ylin_col, idx_Y)
	
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_pos_secondary, node_Lr_primary, -self.tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_neg_secondary, node_Lr_primary, self.tr,
														  Ylin_val, Ylin_row, Ylin_col, idx_Y)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_neg_from, node_Lr_primary, -1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)

			# Imaginary  Circuit Transpose
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_pos_from, node_Li_primary, 1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_pos_secondary, node_Li_primary, -self.tr,
														  Ylin_val, Ylin_row, Ylin_col, idx_Y)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_neg_from, node_Li_primary, -1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)
			
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_neg_secondary, node_Li_primary, self.tr,
														  Ylin_val, Ylin_row, Ylin_col, idx_Y)
	# # PRIMARY CURRENTS # #
	# Real #
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_from, node_Vr_primary, 1, Ylin_val,
												 Ylin_row, Ylin_col, idx_Y)

	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_from, node_Vr_primary, -1, Ylin_val,
												 Ylin_row, Ylin_col, idx_Y)

	# Imaginary #
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_from, node_Vi_primary, 1, Ylin_val,
												 Ylin_row, Ylin_col, idx_Y)

	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_from, node_Vi_primary, -1, Ylin_val,
												 Ylin_row, Ylin_col, idx_Y)
	
	if self.stamp_dual:
		# Real #
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_pos_from, 1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)
	
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_neg_from, -1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)
	
		# Imaginary #
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_pos_from, 1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)
	
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_neg_from,-1, Ylin_val,
														  Ylin_row, Ylin_col, idx_Y)  

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_secondary(self,
                    Ylin_val,
                    Ylin_row,
                    Ylin_col,
                    idx_Y,
					int node_Vr_pos_to,
					int node_Vi_pos_to,
					int node_Vr_neg_secondary,
					int node_Vi_neg_secondary,
					int node_Vr_pos_secondary,
					int node_Vi_pos_secondary,
					int node_Vr_primary,
					int node_Vi_primary,
					dual_dict):
	# Calculate the conductance, susceptance and transformer ration based on homotopy option
	if self.stamp_dual:
		node_Lr_pos_to = dual_dict['Lr_to_pos']
		node_Li_pos_to = dual_dict['Li_to_pos']
		node_Lr_neg_secondary = dual_dict['Lr_to_neg']
		node_Li_neg_secondary = dual_dict['Li_to_neg']
		node_Lr_pos_secondary = dual_dict['Lr_pos_secondary']
		node_Li_pos_secondary = dual_dict['Li_pos_secondary'] 
		node_Lr_primary = dual_dict['Lr_primary']
		node_Li_primary = dual_dict['Li_primary']
		
	_tr = self.tr
	_B = self.B
	_G = self.G
	if self.hasShunt:
		_Bshunt = self.Bshunt
		_Gshunt = self.Gshunt
	else:
		_Bshunt = 0
		_Gshunt = 0

	if _tr != 0:
		if self.ang:
			# ==============================================================================
			#         Phase ANGLE != 0
			# ==============================================================================
			_tr_sin = self.tr * self.sind(self.ang)
			_tr_cos = self.tr * self.cosd(self.ang)

			# # SECONDARY CURRENT SOURCES # #
			# Real #
			# Ir_pos_secondary = - tr * Ir_primary * cos(theta) - tr * Ii_primary * sin(theta)
			# Ir_neg_secondary = tr * Ir_primary * cos(theta) + tr * Ii_primary * sin(theta)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -_tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vi_primary, -_tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vr_primary, _tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vi_primary, _tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			# Imaginary #
			# Ii_pos_secondary = tr * Ir_primary * sin(theta) - tr * Ii_primary * cos(theta)
			# Ii_neg_secondary = -tr * Ir_primary * sin(theta) + tr * Ii_primary * cos(theta)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vr_primary, _tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -_tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vr_primary, -_tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vi_primary, _tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			if self.stamp_dual:
				# Real Circuit - Transpose
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_pos_secondary, -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Lr_pos_secondary, -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_neg_secondary, _tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Lr_neg_secondary,  _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				# Imaginary - Transpose

				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Li_pos_secondary,  _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_pos_secondary, -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Li_neg_secondary,  -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_neg_secondary,  _tr_cos,
																  Ylin_val, Ylin_row, Ylin_col, idx_Y) 
		else:
			# ==============================================================================
			#         Phase ANGLE = 0
			# ==============================================================================

			# # SECONDARY CURRENT SOURCES # #
			# Real #
			# Ir_pos_secondary = - tr * Ir_primary
			# Ir_neg_secondary = tr * Ir_primary
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -_tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vr_primary, _tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			# Imaginary #
			# Ii_pos_secondary = - tr * Ii_primary
			# Ii_neg_secondary = tr * Ii_primary
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -_tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vi_primary, _tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_pos_secondary, -_tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_neg_secondary, _tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_pos_secondary, -_tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_neg_secondary,  _tr,
																  Ylin_val, Ylin_row, Ylin_col, idx_Y) 
	
	# Branch stamping for the susceptances
	if _B != 0:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vi_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vi_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vr_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vr_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vi_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vi_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vr_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vr_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)

		if self.stamp_dual:
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Lr_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Lr_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Li_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Li_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Lr_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Lr_pos_secondary,  _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Li_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)     
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Li_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)

	# Adding resistance if existent
	if _G:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vr_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vr_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vi_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vi_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vr_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vr_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vi_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vi_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)

		if self.stamp_dual:
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Lr_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Lr_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Li_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Li_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Lr_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Lr_pos_secondary,  -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Li_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)     
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Li_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
	# ==============================================================================
	#         SHUNTS MISSING
	# ==============================================================================
	if self.hasShunt:
		# Stamp shunt on the secondary of the transformer
		if _Bshunt:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vr_pos_to, node_Vi_pos_to, -_Bshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vi_pos_to, node_Vr_pos_to, _Bshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			
			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Li_pos_to, node_Lr_pos_to, -_Bshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Lr_pos_to, node_Li_pos_to, _Bshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)

		if _Gshunt:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vr_pos_to, node_Vr_pos_to, _Gshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vi_pos_to, node_Vi_pos_to, _Gshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			
			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Lr_pos_to, node_Lr_pos_to, _Gshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Li_pos_to, node_Li_pos_to, _Gshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

# This function is only used for center-tap transformers right now
cpdef stamp_GB(double _G,
			  double _B,
			 Ylin_val,
			 Ylin_row,
			 Ylin_col,
			 idx_Y,
			  int pos_real,
			  int pos_imag,
			  int extra_real,
			  int extra_imag):
	# Branch stamping for the susceptances
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_real, extra_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_imag, extra_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_real, extra_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_real, pos_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_imag, extra_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Adding resistance if existent
	if _G != 0:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, extra_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, extra_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_real, extra_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_real, pos_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_imag, extra_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_imag, pos_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_GB_dual(double _G,
				  double _B, 
				  Ylin_val, 
				  Ylin_row, 
				  Ylin_col, 
				  idx_Y, 
				  int pos_real,
				  int pos_imag, 
				  int extra_real, 
				  int extra_imag):

	# Transpose of stampGB
	# Branch stamping for the susceptances
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_real, pos_imag,  _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_imag, extra_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_imag, extra_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		extra_real, extra_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	#
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
		pos_real, extra_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	

	# Adding resistance if existent
	if _G != 0:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_real, pos_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_imag, pos_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_real, extra_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, extra_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			extra_imag, extra_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, extra_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

# This function is only used for center-tap transformers right now
cpdef stamp_GB_shunt(double _G,
					double _B,
				   Ylin_val,
				   Ylin_row,
				   Ylin_col,
				   idx_Y,
					int pos_real,
					int pos_imag):
	if _B != 0:
		# Branch stamping for the susceptances
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Adding resistance if existent
	if _G != 0:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_GB_shunt_dual(double _G, 
						double _B, 
						Ylin_val, 
						Ylin_row, 
						Ylin_col, 
						idx_Y,
						int pos_real,
						int pos_imag):
		
	if _B != 0:
		# Branch stamping for the susceptances
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Adding resistance if existent
	if _G != 0:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_ground_neutral(self,
						   int neg_real_node,
						   int neg_imag_node,
						   int gnd_real_V,
						   int gnd_imag_V,
						 Ylin_val,
						 Ylin_row,
						 Ylin_col,
						 idx_Y,
						 stamped_ground,
						 grounded_node,
						 dual_dict):
	if grounded_node not in stamped_ground:
		# Stamp Vng = 0 + j0
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			gnd_real_V, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			gnd_imag_V, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		# Stamp Current Sources
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			neg_real_node, gnd_real_V, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			neg_imag_node, gnd_imag_V, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		stamped_ground.add(grounded_node)
		
		if self.stamp_dual:
			node_Lr_neg = dual_dict['Lr_neg']
			node_Li_neg = dual_dict['Li_neg']
			node_Lr_gnd = dual_dict['Lr_gnd']
			node_Li_gnd = dual_dict['Li_gnd']

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Lr_neg, node_Lr_gnd, 1, Ylin_val, Ylin_row, Ylin_col,
				idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Li_neg, node_Li_gnd, 1, Ylin_val, Ylin_row, Ylin_col,
				idx_Y)
			#
			# Stamp Current Sources
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Lr_gnd, node_Lr_neg, 1, Ylin_val, Ylin_row, Ylin_col,
				idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Li_gnd, node_Li_neg, 1, Ylin_val, Ylin_row, Ylin_col,
				idx_Y)
	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_linear(self,
				 Ylin_val,
				 Ylin_row,
				 Ylin_col,
				 idx_Y,
				 stamped_ground,
				 tx_stepping_enabled = False,
				 h_factor = 0,
				 g_homotopy = 0,
				 b_homotopy = 0):
	# Tx-stepping calls are now incorporates where ever primary/secondary calls are made
	# They are distinguished by the tx_stepping_enabled variable.
	if self.casetype == 1 and ((self.phases == 1 or self.phases == 2) or self.phases == 4):
		nodes_from, nodes_to, nodes_intermediate = self.get_nodes()
		# #########Stamp Primary Circuit##################
		# from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
		_pos = 0
		_neg = 1
		_winding = 0
		if self.stamp_dual:
			dualdict = {
				'Lr_from_pos': nodes_from['LR'][_pos], 
				'Li_from_pos': nodes_from['LI'][_pos],
				'Lr_from_neg': nodes_from['LR'][_neg], 
				'Li_from_neg': nodes_from['LI'][_neg],
				'Lr_to_pos': nodes_to['LR'][_pos], 
				'Li_to_pos': nodes_to['LI'][_pos],
				'Lr_to_neg': nodes_to['LR'][_neg], 
				'Li_to_neg': nodes_to['LI'][_neg],
				'Lr_pos_secondary': nodes_intermediate['LR POS SEC'][_winding],
				'Li_pos_secondary': nodes_intermediate['LI POS SEC'][_winding],
				'Lr_primary': nodes_intermediate['LR PRIM'][_winding],
				'Li_primary': nodes_intermediate['LI PRIM'][_winding]}
		else: 
			dualdict = None
				
		if tx_stepping_enabled:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_tx_stepping(
				Ylin_val, Ylin_row, Ylin_col, idx_Y,
				nodes_to['VR'][0], nodes_to['VI'][0],
				nodes_to['VR'][1], nodes_to['VI'][1],
				nodes_intermediate['VR POS SEC'][0],
				nodes_intermediate['VI POS SEC'][0],
				nodes_intermediate['VR PRIM'][0],
				nodes_intermediate['VI PRIM'][0],
				h_factor, g_homotopy, b_homotopy, dual_dict = dualdict)
		else:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
				Ylin_val, Ylin_row, Ylin_col, idx_Y,
				nodes_from['VR'][0], nodes_from['VI'][0],
				nodes_from['VR'][1], nodes_from['VI'][1],
				nodes_to['VR'][1], nodes_to['VI'][1],
				nodes_intermediate['VR POS SEC'][0],
				nodes_intermediate['VI POS SEC'][0],
				nodes_intermediate['VR PRIM'][0],
				nodes_intermediate['VI PRIM'][0], 0, dualdict)

			# #########Stamp Secondary Circuit##################
			# Without homotopy option OFF
			Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
				Ylin_val, Ylin_row, Ylin_col, idx_Y,
				nodes_to['VR'][0], nodes_to['VI'][0],
				nodes_to['VR'][1], nodes_to['VI'][1],
				nodes_intermediate['VR POS SEC'][0],
				nodes_intermediate['VI POS SEC'][0],
				nodes_intermediate['VR PRIM'][0],
				nodes_intermediate['VI PRIM'][0], dual_dict = dualdict)


	else:
		if not self.isSinglePhase:
			# Find the from (Primary) bus nodes to stamp
			nodeA_Vr_from = self.nodeA_Vr_from
			nodeA_Vi_from = self.nodeA_Vi_from
			nodeB_Vr_from = self.nodeB_Vr_from
			nodeB_Vi_from = self.nodeB_Vi_from
			nodeC_Vr_from = self.nodeC_Vr_from
			nodeC_Vi_from = self.nodeC_Vi_from

			wyePrimarySet = {
				self._WYE_WYE, self._GWYE_GWYE, self._WYE_DELTA,
				self._GWYE_DELTA
			}
			if self.connect_type in wyePrimarySet:
				nodeN_Vr_from = self.nodeN_Vr_from
				nodeN_Vi_from = self.nodeN_Vi_from

				if self.stamp_dual:
					nodeN_Lr_from = self.nodeN_Lr_from
					nodeN_Li_from = self.nodeN_Li_from

			# Find the to (Secondary) bus nodes to stamp
			nodeA_Vr_to = self.nodeA_Vr_to
			nodeA_Vi_to = self.nodeA_Vi_to
			nodeB_Vr_to = self.nodeB_Vr_to
			nodeB_Vi_to = self.nodeB_Vi_to
			nodeC_Vr_to = self.nodeC_Vr_to
			nodeC_Vi_to = self.nodeC_Vi_to

			# Collect the secondary positive terminal nodes
			nodeA_Vr_pos_secondary = self.nodeA_Vr_pos_secondary
			nodeB_Vr_pos_secondary = self.nodeB_Vr_pos_secondary
			nodeC_Vr_pos_secondary = self.nodeC_Vr_pos_secondary
			nodeA_Vi_pos_secondary = self.nodeA_Vi_pos_secondary
			nodeB_Vi_pos_secondary = self.nodeB_Vi_pos_secondary
			nodeC_Vi_pos_secondary = self.nodeC_Vi_pos_secondary

			nodeA_Vr_primary = self.nodeA_Vr_primary
			nodeB_Vr_primary = self.nodeB_Vr_primary
			nodeC_Vr_primary = self.nodeC_Vr_primary
			nodeA_Vi_primary = self.nodeA_Vi_primary
			nodeB_Vi_primary = self.nodeB_Vi_primary
			nodeC_Vi_primary = self.nodeC_Vi_primary

			wyeSecondarySet = {
				self._WYE_WYE, self._DELTA_GWYE, self._GWYE_GWYE,
				self._DELTA_WYE
			}
			if self.connect_type in wyeSecondarySet:
				nodeN_Vr_to = self.nodeN_Vr_to
				nodeN_Vi_to = self.nodeN_Vi_to

				if self.stamp_dual:
					nodeN_Lr_to = self.nodeN_Lr_to
					nodeN_Li_to = self.nodeN_Li_to
			
			if self.stamp_dual:
					nodeA_Lr_from = self.nodeA_Lr_from
					nodeA_Li_from = self.nodeA_Li_from
					nodeB_Lr_from = self.nodeB_Lr_from
					nodeB_Li_from = self.nodeB_Li_from
					nodeC_Lr_from = self.nodeC_Lr_from
					nodeC_Li_from = self.nodeC_Li_from
					
					nodeA_Lr_to = self.nodeA_Lr_to
					nodeA_Li_to = self.nodeA_Li_to
					nodeB_Lr_to = self.nodeB_Lr_to
					nodeB_Li_to = self.nodeB_Li_to
					nodeC_Lr_to = self.nodeC_Lr_to
					nodeC_Li_to = self.nodeC_Li_to
					
					nodeA_Lr_pos_secondary = self.nodeA_Lr_pos_secondary
					nodeA_Li_pos_secondary = self.nodeA_Li_pos_secondary 
					nodeB_Lr_pos_secondary = self.nodeB_Lr_pos_secondary
					nodeB_Li_pos_secondary = self.nodeB_Li_pos_secondary 
					nodeC_Lr_pos_secondary = self.nodeC_Lr_pos_secondary 
					nodeC_Li_pos_secondary = self.nodeC_Li_pos_secondary 

					nodeA_Lr_primary = self.nodeA_Lr_primary
					nodeA_Li_primary = self.nodeA_Li_primary 
					nodeB_Lr_primary = self.nodeB_Lr_primary 
					nodeB_Li_primary = self.nodeB_Li_primary 
					nodeC_Lr_primary = self.nodeC_Lr_primary 
					nodeC_Li_primary = self.nodeC_Li_primary 

			# Put together the set of nodes
			nodeVr_from = [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from]
			nodeVi_from = [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from]
			nodeVr_to = [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to]
			nodeVi_to = [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]
			node_Vr_pos_secondary = [nodeA_Vr_pos_secondary, nodeB_Vr_pos_secondary, nodeC_Vr_pos_secondary]
			node_Vi_pos_secondary = [nodeA_Vi_pos_secondary, nodeB_Vi_pos_secondary, nodeC_Vi_pos_secondary]
			node_Vr_primary = [nodeA_Vr_primary, nodeB_Vr_primary, nodeC_Vr_primary]
			node_Vi_primary = [nodeA_Vi_primary, nodeB_Vi_primary, nodeC_Vi_primary]
			
			if self.stamp_dual:
					nodeLr_from = [nodeA_Lr_from, nodeB_Lr_from, nodeC_Lr_from]
					nodeLi_from = [nodeA_Li_from, nodeB_Li_from, nodeC_Li_from]
					nodeLr_to = [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to]
					nodeLi_to = [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to]
					node_Lr_pos_secondary = [nodeA_Lr_pos_secondary, nodeB_Lr_pos_secondary, nodeC_Lr_pos_secondary]
					node_Li_pos_secondary = [nodeA_Li_pos_secondary, nodeB_Li_pos_secondary, nodeC_Li_pos_secondary]
					node_Lr_primary = [nodeA_Lr_primary, nodeB_Lr_primary, nodeC_Lr_primary]
					node_Li_primary = [nodeA_Li_primary, nodeB_Li_primary, nodeC_Li_primary]

			(_A, _B, _C) = (0, 1, 2)
			_N = 3

		elif self.isSinglePhase:
			pass

		try:
			# Individual phases of the transformer is implemented for this type
			# Option 1 - Stamp wye-wye
			if self.connect_type == self._WYE_WYE or self.connect_type == self._GWYE_GWYE:
				nodeVr_from.append(nodeN_Vr_from)
				nodeVi_from.append(nodeN_Vi_from)
				nodeVr_to.append(nodeN_Vr_to)
				nodeVi_to.append(nodeN_Vi_to)
				
				if self.stamp_dual:
					nodeLr_from.append(nodeN_Lr_from)
					nodeLi_from.append(nodeN_Li_from)
					nodeLr_to.append(nodeN_Lr_to)
					nodeLi_to.append(nodeN_Li_to)
					
				# Third number in the set is used to find what phases are on
				phases = [(_A, _N, self.phaseA), (_B, _N, self.phaseB),
						(_C, _N, self.phaseC)]
				for phaseSet in phases:
					if phaseSet[2] & self.phases == phaseSet[2]:
						# #########Stamp Primary Circuit##################
						# from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
						_pos = phaseSet[0]
						_neg = phaseSet[1]
						_winding = phaseSet[0]
						if self.stamp_dual:
							dualdict = {
								'Lr_from_pos': nodeLr_from[_pos], 
								'Li_from_pos': nodeLi_from[_pos],
								'Lr_from_neg': nodeLr_from[_neg], 
								'Li_from_neg': nodeLi_from[_neg],
								'Lr_to_pos': nodeLr_to[_pos], 
								'Li_to_pos': nodeLi_to[_pos],
								'Lr_to_neg': nodeLr_to[_neg], 
								'Li_to_neg': nodeLi_to[_neg],
								'Lr_pos_secondary': node_Lr_pos_secondary[_winding],
								'Li_pos_secondary': node_Li_pos_secondary[_winding],
								'Lr_primary': node_Lr_primary[_winding],
								'Li_primary': node_Li_primary[_winding]}
						else: 
							dualdict = None
							
							
						if tx_stepping_enabled:
							Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_tx_stepping(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodeVr_to[_pos], nodeVi_to[_pos],
								nodeVr_to[_neg], nodeVi_to[_neg],
								node_Vr_pos_secondary[_winding],
								node_Vi_pos_secondary[_winding],
								node_Vr_primary[_winding], node_Vi_primary[_winding],
								h_factor, g_homotopy, b_homotopy, dual_dict = dualdict)
						else:
							Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodeVr_from[_pos], nodeVi_from[_pos],
								nodeVr_from[_neg], nodeVi_from[_neg],
								nodeVr_to[_neg], nodeVi_to[_neg],
								node_Vr_pos_secondary[_winding],
								node_Vi_pos_secondary[_winding],
								node_Vr_primary[_winding], node_Vi_primary[_winding], _pos, dualdict)

							# #########Stamp Secondary Circuit##################
							# Without homotopy option OFF
							Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodeVr_to[_pos], nodeVi_to[_pos],
								nodeVr_to[_neg], nodeVi_to[_neg],
								node_Vr_pos_secondary[_winding],
								node_Vi_pos_secondary[_winding],
								node_Vr_primary[_winding], node_Vi_primary[_winding], dual_dict = dualdict)

			if self.connect_type == self._GWYE_GWYE:
				# #########Stamp Neutral to Ground##################
				_neg = _N
				if self.stamp_dual:
					dual_dict_prim = {
					'Lr_neg': nodeLr_from[_neg],
					'Li_neg': nodeLi_from[_neg],
					'Lr_gnd': self.nodeGnd_Lr_primary,
					'Li_gnd': self.nodeGnd_Li_primary,
					}
					dual_dict_sec = {
					'Lr_neg': nodeLr_to[_neg],
					'Li_neg': nodeLi_to[_neg],
					'Lr_gnd': self.nodeGnd_Lr_secondary,
					'Li_gnd': self.nodeGnd_Li_secondary,
					}
				else:
					dual_dict_prim = None
					dual_dict_sec = None
				# One for set of three phase
				# Stamp primary
				Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
					nodeVr_from[_neg], nodeVi_from[_neg], self.nodeGnd_Vr_primary,
					self.nodeGnd_Vi_primary, Ylin_val, Ylin_row, Ylin_col, idx_Y,
					stamped_ground, self.from_node, dual_dict = dual_dict_prim)
				# Stamp secondary
				Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
					nodeVr_to[_neg], nodeVi_to[_neg], self.nodeGnd_Vr_secondary,
					self.nodeGnd_Vi_secondary, Ylin_val, Ylin_row, Ylin_col, idx_Y,
					stamped_ground, self.to_node, dual_dict = dual_dict_sec)

			# Option 2 - Stamp Delta - Delta
			elif self.connect_type == self._DELTA_DELTA:
				phases = [(_A, _B), (_B, _C), (_C, _A)]
				for phaseSet in phases:
					# In delta - delta get a tuple of AB, BC, CA
					# #########Stamp Primary Circuit##################
					# from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
					_pos = phaseSet[0]
					_neg = phaseSet[1]
					_winding = phaseSet[0]
					if self.stamp_dual:
						dualdict = {
								'Lr_from_pos': nodeLr_from[_pos], 
								'Li_from_pos': nodeLi_from[_pos],
								'Lr_from_neg': nodeLr_from[_neg], 
								'Li_from_neg': nodeLi_from[_neg],
								'Lr_to_pos': nodeLr_to[_pos], 
								'Li_to_pos': nodeLi_to[_pos],
								'Lr_to_neg': nodeLr_to[_neg], 
								'Li_to_neg': nodeLi_to[_neg],
								'Lr_pos_secondary': node_Lr_pos_secondary[_winding],
								'Li_pos_secondary': node_Li_pos_secondary[_winding],
								'Lr_primary': node_Lr_primary[_winding],
								'Li_primary': node_Li_primary[_winding]}
					else: 
						dualdict = None
						
					if tx_stepping_enabled:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_tx_stepping(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos], nodeVi_to[_pos],
							nodeVr_to[_neg], nodeVi_to[_neg],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding],
							node_Vr_primary[_winding], node_Vi_primary[_winding],
							h_factor, g_homotopy, b_homotopy, dual_dict = dualdict)
					else:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_from[_pos], nodeVi_from[_pos],
							nodeVr_from[_neg], nodeVi_from[_neg],
							nodeVr_to[_neg],
							nodeVi_to[_neg], node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], _pos, dual_dict = dualdict)

						# #########Stamp Secondary Circuit##################
						# Without homotopy - option off
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos], nodeVi_to[_pos], nodeVr_to[_neg],
							nodeVi_to[_neg], node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], dual_dict = dualdict)


			# Option 3 or 7- Stamp Delta - GWye
			# ** This will have a phase shift
			elif self.connect_type == self._DELTA_GWYE or self.connect_type == self._DELTA_WYE:

				# The connection is as follows: (Step down)
				# from -> [(_C, _A), (_A, _B), (_B, _C)]
				# to ->   [(_N, _A), (_N, _B), (_N, _C)]

				# The connection is as follows: (Step up)
				# from -> [(_A, _B), (_B, _C), (_C, _A)]
				# to ->   [(_A, _N), (_B, _N), (_C, _N)]
				if self.is_stepdown:
					phases = [(_C, _A), (_A, _B), (_B, _C)]
				else:
					phases = [(_A, _B), (_B, _C), (_C, _A)]

				nodeVr_to.append(nodeN_Vr_to)
				nodeVi_to.append(nodeN_Vi_to)
				
				if self.stamp_dual:
					nodeLr_to.append(nodeN_Lr_to)
					nodeLi_to.append(nodeN_Li_to)
					
				for phaseSet in phases:
					# from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
					if self.is_stepdown:
						_pos_from = phaseSet[0]
						_pos_to = _N
						_neg_from = phaseSet[1]  # Primary is delta connected
						# *This unusual behavior is to get NA, NB, NC
						_neg_to = phaseSet[1]  # Secondary is wye connected
						_winding = phaseSet[1]
					else:
						_pos_from = phaseSet[0]
						_pos_to = phaseSet[0]
						_neg_from = phaseSet[1]  # Primary is delta connected
						_neg_to = _N  # Secondary is wye connected
						_winding = phaseSet[0]
					
					if self.stamp_dual:
						dualdict = {
								'Lr_from_pos': nodeLr_from[_pos_from], 
								'Li_from_pos': nodeLi_from[_pos_from],
								'Lr_from_neg': nodeLr_from[_neg_from], 
								'Li_from_neg': nodeLi_from[_neg_from],
								'Lr_to_pos': nodeLr_to[_pos_to], 
								'Li_to_pos': nodeLi_to[_pos_to],
								'Lr_to_neg': nodeLr_to[_neg_to], 
								'Li_to_neg': nodeLi_to[_neg_to],
								'Lr_pos_secondary': node_Lr_pos_secondary[_winding],
								'Li_pos_secondary': node_Li_pos_secondary[_winding],
								'Lr_primary': node_Lr_primary[_winding],
								'Li_primary': node_Li_primary[_winding]}
					else: 
						dualdict = None
						
					if tx_stepping_enabled:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_tx_stepping(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos_to], nodeVi_to[_pos_to],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding],
							node_Vr_primary[_winding], node_Vi_primary[_winding],
							h_factor, g_homotopy, b_homotopy, dual_dict = dualdict)
					else:
						# #########Stamp Primary Circuit - DELTA##################
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_from[_pos_from], nodeVi_from[_pos_from],
							nodeVr_from[_neg_from], nodeVi_from[_neg_from],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], _pos_from, dual_dict = dualdict)

						# #########Stamp Secondary Circuit - Gwye##################
						# Without homotopy - option off
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos_to], nodeVi_to[_pos_to],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], dual_dict = dualdict)

				if self.connect_type == self._DELTA_GWYE:
					# #########Stamp Neutral to Ground##################
					# One for set of three phase
					# Stamp secondary
					# Make sure to ground the neutral and not neg to in this case
					if self.stamp_dual:
						dual_dict_sec = {
						'Lr_neg': nodeLr_to[_N],
						'Li_neg': nodeLi_to[_N],
						'Lr_gnd': self.nodeGnd_Lr_secondary,
						'Li_gnd': self.nodeGnd_Li_secondary,
						}
					else:
						dual_dict_sec = None
						
					Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
						nodeVr_to[_N], nodeVi_to[_N], self.nodeGnd_Vr_secondary,
						self.nodeGnd_Vi_secondary, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						stamped_ground, self.to_node, dual_dict = dual_dict_sec)

			# Option 8 or 9- Stamp Wye- Delta or GWYE-Delta
			elif self.connect_type == self._WYE_DELTA or self.connect_type == self._GWYE_DELTA:
				# The connection is as follows: (Step down)
				# from -> [(_A, _N), (_B, _N), (_C, _N)]
				# to ->   [(_A, _B), (_B, _C), (_C, _A)]

				# The connection is as follows: (Step up)
				# from -> [(_A, _N), (_B, _N), (_C, _N)]
				# to ->   [(_A, _C), (_B, _A), (_C, _B)]
				if self.is_stepdown:
					phases = [(_A, _B), (_B, _C), (_C, _A)]
				else:
					phases = [(_A, _C), (_B, _A), (_C, _B)]
				nodeVr_from.append(nodeN_Vr_from)
				nodeVi_from.append(nodeN_Vi_from)
				
				if self.stamp_dual:
					nodeLr_from.append(nodeN_Lr_from)
					nodeLi_from.append(nodeN_Li_from)

				for phaseSet in phases:
					# from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
					_pos = phaseSet[0]
					_neg_from = _N  # Primary is wye connected
					_neg_to = phaseSet[1]  # Secondary is delta connected
					_winding = phaseSet[0]
					
					if self.stamp_dual:
						dualdict = {
								'Lr_from_pos': nodeLr_from[_pos], 
								'Li_from_pos': nodeLi_from[_pos],
								'Lr_from_neg': nodeLr_from[_neg_from], 
								'Li_from_neg': nodeLi_from[_neg_from],
								'Lr_to_pos': nodeLr_to[_pos], 
								'Li_to_pos': nodeLi_to[_pos],
								'Lr_to_neg': nodeLr_to[_neg_to], 
								'Li_to_neg': nodeLi_to[_neg_to],
								'Lr_pos_secondary': node_Lr_pos_secondary[_winding],
								'Li_pos_secondary': node_Li_pos_secondary[_winding],
								'Lr_primary': node_Lr_primary[_winding],
								'Li_primary': node_Li_primary[_winding]}
					else: 
						dualdict = None
						
					if tx_stepping_enabled:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_tx_stepping(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos], nodeVi_to[_pos],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding],
							node_Vr_primary[_winding], node_Vi_primary[_winding],
							h_factor, g_homotopy, b_homotopy, dual_dict = dualdict)
					else:
						# #########Stamp Primary Circuit - DELTA##################
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_from[_pos], nodeVi_from[_pos],
							nodeVr_from[_neg_from], nodeVi_from[_neg_from],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], _pos, dual_dict = dualdict)

						# #########Stamp Secondary Circuit - Gwye##################
						# Homotopy option - disabled
						Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
							Ylin_val, Ylin_row, Ylin_col, idx_Y,
							nodeVr_to[_pos], nodeVi_to[_pos],
							nodeVr_to[_neg_to], nodeVi_to[_neg_to],
							node_Vr_pos_secondary[_winding],
							node_Vi_pos_secondary[_winding], node_Vr_primary[_winding],
							node_Vi_primary[_winding], dual_dict = dualdict)

				if self.connect_type == self._GWYE_DELTA:
					# #########Stamp Neutral to Ground##################
					# One for set of three phase
					# Stamp secondary
					if self.stamp_dual:
						dual_dict_prim = {
							'Lr_neg': nodeLr_from[_neg_from],
							'Li_neg': nodeLi_from[_neg_from],
							'Lr_gnd': self.nodeGnd_Lr_primary,
							'Li_gnd': self.nodeGnd_Li_primary,
							}
					else:
						dual_dict_prim = None
						
					Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
						nodeVr_from[_neg_from], nodeVi_from[_neg_from],
						self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary, Ylin_val, Ylin_row,
						Ylin_col, idx_Y, stamped_ground, self.from_node, dual_dict = dual_dict_prim)
					
			elif self.connect_type == self._CENTER_TAP:
				if self.phases & 0x01 == 1:  # Phase A
					node1_Vr_from = self.node1_Vr_from
					node1_Vi_from = self.node1_Vi_from
					if self.stamp_dual:
						node1_Lr_from = self.node1_Lr_from 
						node1_Li_from = self.node1_Li_from 
				elif self.phases & 0x02 == 2:  # Phase B
					node1_Vr_from = self.node1_Vr_from
					node1_Vi_from = self.node1_Vi_from
					if self.stamp_dual:
						node1_Lr_from = self.node1_Lr_from 
						node1_Li_from = self.node1_Li_from 
				elif self.phases & 0x04 == 4:  # Phase C
					node1_Vr_from = self.node1_Vr_from
					node1_Vi_from = self.node1_Vi_from
					if self.stamp_dual:
						node1_Lr_from = self.node1_Lr_from 
						node1_Li_from = self.node1_Li_from 
						
				if not tx_stepping_enabled:
					# *This is the case when neutral is grounded (both on primary and secondary
					# Stamp the primary 4 current sources (between Xtra nodes and ground)
					# Primary Current Sources
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node_Vr_CT_primary, self.node1_Vr_secondary, -(1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node_Vr_CT_primary, self.node2_Vr_secondary, (1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node_Vi_CT_primary, self.node1_Vi_secondary, -(1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node_Vi_CT_primary, self.node2_Vi_secondary, (1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Lr_secondary, self.node_Lr_CT_primary, -(1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)
						#
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Lr_secondary, self.node_Lr_CT_primary, (1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)
						#
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Li_secondary, self.node_Li_CT_primary, -(1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)
						#
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Li_secondary, self.node_Li_CT_primary, (1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)
						
					# Stamp 4 voltage sources (with 8 stamps) (each has two elements)
					# Van_R - (VA_R/tr) = 0
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vr_secondary, self.node1_Vr_sending, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vr_secondary, self.nodeN_Vr_to, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vr_secondary, self.node_Vr_CT_primary, (-1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Lr_sending, self.node1_Lr_secondary, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.nodeN_Lr_to, self.node1_Lr_secondary, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node_Lr_CT_primary, self.node1_Lr_secondary,(-1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)
						
					# Vbn_R + (VA_R/tr) = 0
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vr_secondary, self.node2_Vr_sending, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vr_secondary, self.nodeN_Vr_to, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vr_secondary, self.node_Vr_CT_primary, (1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Lr_sending, self.node2_Lr_secondary, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.nodeN_Lr_to, self.node2_Lr_secondary, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node_Lr_CT_primary, self.node2_Lr_secondary, (1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)  
						
					# Van_I - (VA_I/tr) = 0
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vi_secondary, self.node1_Vi_sending, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vi_secondary, self.nodeN_Vi_to, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vi_secondary, self.node_Vi_CT_primary, (-1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Li_sending, self.node1_Li_secondary, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.nodeN_Li_to, self.node1_Li_secondary, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node_Li_CT_primary, self.node1_Li_secondary, (-1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)  
						
					# Van_I + (VA_I/tr) = 0
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vi_secondary, self.node2_Vi_sending, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vi_secondary, self.nodeN_Vi_to, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vi_secondary, self.node_Vi_CT_primary, (1 / self.tr),
						Ylin_val, Ylin_row, Ylin_col, idx_Y)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Li_sending, self.node2_Li_secondary, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.nodeN_Li_to, self.node2_Li_secondary, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node_Li_CT_primary, self.node2_Li_secondary, (1 / self.tr),
							Ylin_val, Ylin_row, Ylin_col, idx_Y)                        
					
					# Stamp the current sources for the voltage sources
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vr_sending, self.node1_Vr_secondary, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vr_sending, self.node2_Vr_secondary, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node1_Vi_sending, self.node1_Vi_secondary, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.node2_Vi_sending, self.node2_Vi_secondary, 1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)

					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.nodeN_Vr_to, self.node1_Vr_secondary, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.nodeN_Vi_to, self.node2_Vr_secondary, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.nodeN_Vr_to, self.node1_Vi_secondary, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
						self.nodeN_Vi_to, self.node2_Vi_secondary, -1, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Lr_secondary, self.node1_Lr_sending,  1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Lr_secondary, self.node2_Lr_sending, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Li_secondary, self.node1_Li_sending, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Li_secondary, self.node2_Li_sending, 1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)

						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Lr_secondary, self.nodeN_Lr_to, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Lr_secondary, self.nodeN_Li_to, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node1_Li_secondary, self.nodeN_Lr_to, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							self.node2_Li_secondary, self.nodeN_Li_to, -1, Ylin_val,
							Ylin_row, Ylin_col, idx_Y) 
						
					# Ground the neutrals if needed
					if self.stamp_dual:
						dual_dict_sec = {
							'Lr_neg': self.nodeN_Lr_to,
							'Li_neg': self.nodeN_Li_to,
							'Lr_gnd': self.nodeGnd_Lr_secondary,
							'Li_gnd': self.nodeGnd_Li_secondary,
							}
					else:
						dual_dict_sec = None
					Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
						self.nodeN_Vr_to, self.nodeN_Vi_to, self.nodeGnd_Vr_secondary,
						self.nodeGnd_Vi_secondary, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						stamped_ground, self.to_node, dual_dict = dual_dict_sec)
					
					
					# This is the impedance stamps shared by
					# Stamp Z0
					(self.G0, self.B0) = self.calc_G_B(self.r0, self.x0)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						self.G0, self.B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node_Vr_CT_primary, self.node_Vi_CT_primary, node1_Vr_from,
						node1_Vi_from)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							self.G0, self.B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node_Lr_CT_primary, self.node_Li_CT_primary, node1_Lr_from,
							node1_Li_from) 
					
					# Stamp shunt impedance
					if self.hasShunt:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_shunt(
							self.Gshunt, self.Bshunt, Ylin_val, Ylin_row,
							Ylin_col, idx_Y, node1_Vr_from, node1_Vi_from)
						if self.stamp_dual:
							Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_shunt_dual(
								self.Gshunt, self.Bshunt, Ylin_val, Ylin_row,
								Ylin_col, idx_Y, node1_Lr_from, node1_Li_from)                            

					# Stamp Z1 and Z2
					(self.G1, self.B1) = self.calc_G_B(self.r1, self.x1)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						self.G1, self.B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node1_Vr_sending, self.node1_Vi_sending, self.node1_Vr_to,
						self.node1_Vi_to)
					(self.G2, self.B2) = self.calc_G_B(self.r2, self.x2)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						self.G2, self.B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node2_Vr_sending, self.node2_Vi_sending, self.node2_Vr_to,
						self.node2_Vi_to)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							self.G1, self.B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node1_Lr_sending, self.node1_Li_sending, self.node1_Lr_to,
							self.node1_Li_to)
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							self.G2, self.B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node2_Lr_sending, self.node2_Li_sending, self.node2_Lr_to,
							self.node2_Li_to) 
							
				elif tx_stepping_enabled:
					_B0 = self.B0 * h_factor * b_homotopy
					_G0 = self.G0 * h_factor * g_homotopy
					_B1 = self.B1 * h_factor * b_homotopy
					_G1 = self.G1 * h_factor * g_homotopy
					_B2 = self.B2 * h_factor * b_homotopy
					_G2 = self.G2 * h_factor * g_homotopy

					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						_G0, _B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node_Vr_CT_primary, self.node_Vi_CT_primary, node1_Vr_from,
						node1_Vi_from)
					
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							_G0, _B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node_Lr_CT_primary, self.node_Li_CT_primary, node1_Lr_from,
							node1_Li_from)
					# Stamp Z1 and Z2
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						_G1, _B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node1_Vr_sending, self.node1_Vi_sending, self.node1_Vr_to,
						self.node1_Vi_to)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							_G1, _B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node_Lr_CT_primary, self.node_Li_CT_primary, node1_Lr_from,
							node1_Li_from) 
					
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB(
						_G2, _B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
						self.node2_Vr_sending, self.node2_Vi_sending, self.node2_Vr_to,
						self.node2_Vi_to)
					if self.stamp_dual:
						Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_GB_dual(
							_G2, _B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
							self.node_Lr_CT_primary, self.node_Li_CT_primary, node1_Lr_from,
							node1_Li_from) 
			else:
				raise ValueError
		except ValueError:
			print('Illegal connect_type for the transformer')

	return Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground

cpdef stamp_tx_stepping(self,
                    Ylin_val,
                    Ylin_row,
                    Ylin_col,
                    idx_Y,
					int node_Vr_pos_to,
					int node_Vi_pos_to,
					int node_Vr_neg_secondary,
					int node_Vi_neg_secondary,
					int node_Vr_pos_secondary,
					int node_Vi_pos_secondary,
					int node_Vr_primary,
					int node_Vi_primary,
					h_factor,
					g_homotopy,
					b_homotopy,
					dual_dict):
					
	if self.stamp_dual:
		node_Lr_pos_to = dual_dict['Lr_to_pos']
		node_Li_pos_to = dual_dict['Li_to_pos']
		node_Lr_neg_secondary = dual_dict['Lr_to_neg']
		node_Li_neg_secondary = dual_dict['Li_to_neg']
		node_Lr_pos_secondary = dual_dict['Lr_pos_secondary']
		node_Li_pos_secondary = dual_dict['Li_pos_secondary'] 
		node_Lr_primary = dual_dict['Lr_primary']
		node_Li_primary = dual_dict['Li_primary']
		
	# Copied stamp_secondary and made amendments
	# This function is called from stamp_linear in this class

	# B and G Homotopy values
	_B = self.B * h_factor * b_homotopy
	_G = self.G * h_factor * g_homotopy

	# B and G shunt values
	if self.hasShunt:
		_Bshunt = -h_factor * b_homotopy * self.Bshunt
		_Gshunt = -h_factor * g_homotopy * self.Gshunt

	# Tap changer homotopy value
	_tr = h_factor * (1 - self.tr)

	# Phase shifter homotopy value
	if self.ang:
		_ang = -h_factor * self.ang

	# The stamps for tx_stepping do not change from stamp_secondary
	if _tr != 0:
		if self.ang:
			# ==============================================================================
			#         Phase ANGLE != 0
			# ==============================================================================
			_tr_sin = _tr * self.sind(_ang)
			_tr_cos = _tr * self.cosd(_ang)

			# # SECONDARY CURRENT SOURCES # #
			# Real #
			# Ir_pos_secondary = - tr * Ir_primary * cos(theta) - tr * Ii_primary * sin(theta)
			# Ir_neg_secondary = tr * Ir_primary * cos(theta) + tr * Ii_primary * sin(theta)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -_tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vi_primary, -_tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vr_primary, _tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vi_primary, _tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			# Imaginary #
			# Ii_pos_secondary = tr * Ir_primary * sin(theta) - tr * Ii_primary * cos(theta)
			# Ii_neg_secondary = -tr * Ir_primary * sin(theta) + tr * Ii_primary * cos(theta)

			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vr_primary, _tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -_tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vr_primary, -_tr_sin,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vi_primary, _tr_cos,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			if self.stamp_dual:
				# Real Circuit - Transpose
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_pos_secondary, -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Lr_pos_secondary, -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_neg_secondary, _tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Lr_neg_secondary,  _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				# Imaginary - Transpose

				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Li_pos_secondary,  _tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_pos_secondary, -_tr_cos,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Li_neg_secondary,  -_tr_sin,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_neg_secondary,  _tr_cos,
																  Ylin_val, Ylin_row, Ylin_col, idx_Y)             
  
		else:
			# ==============================================================================
			#         Phase ANGLE = 0
			# ==============================================================================

			# # SECONDARY CURRENT SOURCES # #
			# Real #
			# Ir_pos_secondary = - tr * Ir_primary
			# Ir_neg_secondary = tr * Ir_primary
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -_tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vr_neg_secondary, node_Vr_primary, _tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			# Imaginary #
			# Ii_pos_secondary = - tr * Ii_primary
			# Ii_neg_secondary = tr * Ii_primary
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -_tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Vi_neg_secondary, node_Vi_primary, _tr,
														 Ylin_val, Ylin_row, Ylin_col, idx_Y)

			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_pos_secondary, -_tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Lr_primary, node_Lr_neg_secondary, _tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_pos_secondary, -_tr,
															  Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(node_Li_primary, node_Li_neg_secondary,  _tr,
																  Ylin_val, Ylin_row, Ylin_col, idx_Y)             
 

	# Branch stamping for the susceptances
	if _B != 0:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vi_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vi_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vr_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vr_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vi_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vi_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vr_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vr_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		if self.stamp_dual:
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Lr_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Lr_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Li_pos_to, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Li_pos_to, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Lr_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Lr_pos_secondary,  _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Li_pos_secondary, _B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)     
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Li_pos_secondary, -_B, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
	# Adding resistance if existent
	if _G != 0:
		# Branch stamping for the conductance
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vr_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_to, node_Vr_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vi_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_to, node_Vi_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vr_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vr_pos_secondary, node_Vr_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vi_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
			node_Vi_pos_secondary, node_Vi_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		
		if self.stamp_dual:
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Lr_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Lr_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Li_pos_to, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Li_pos_to, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_secondary, node_Lr_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Lr_pos_to, node_Lr_pos_secondary,  -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_secondary, node_Li_pos_secondary, _G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)     
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
							node_Li_pos_to, node_Li_pos_secondary, -_G, Ylin_val, Ylin_row, Ylin_col,
							idx_Y)
			
	if self.hasShunt:
		# Stamp shunt on the secondary of the transformer
		if _Bshunt:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vr_pos_to, node_Vi_pos_to, -_Bshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vi_pos_to, node_Vr_pos_to, _Bshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			
			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Li_pos_to, node_Lr_pos_to, -_Bshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Lr_pos_to, node_Li_pos_to, _Bshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)

		if _Gshunt:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vr_pos_to, node_Vr_pos_to, _Gshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
				node_Vi_pos_to, node_Vi_pos_to, _Gshunt, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			
			if self.stamp_dual:
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Lr_pos_to, node_Lr_pos_to, _Gshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stamp_Y(
					node_Li_pos_to, node_Li_pos_to, _Gshunt, Ylin_val, Ylin_row,
					Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y
