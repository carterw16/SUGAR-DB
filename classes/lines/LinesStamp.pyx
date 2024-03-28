"""Implements line stamps.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 10-05-2020
  Updated Date: 10-18-2020
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""
#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

cpdef stampY(i, j, val, Y_val, Y_row, Y_col, idx):
	if val == 0.:
		return Y_val, Y_row, Y_col, idx
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return Y_val, Y_row, Y_col, idx

cpdef stamp_phase_ckt(self, Ylin_val, Ylin_row, Ylin_col, idx_Y, nodeVr_from,
					  nodeVi_from, nodeVr_to, nodeVi_to, _self, _mutual1,
					  _mutual2, nodeLr_from, nodeLi_from, nodeLr_to, nodeLi_to):
	# FIRST ELEMENT ALWAYS REPRESENTS THE CURRENT INDEX FOR MNA
	# THEREFORE IT IS EITHER THE REAL OR IMAGINARY CURRENT INDEX
	# FOR THAT PHASE
	Bself = self.Ymatrix[_self, _self].imag
	if self.line_parameters['file type'] == 'DSS' and self.single_phase != -1:
		Bmutual1 = None
		Bmutual2 = None
	else:
		Bmutual1 = self.Ymatrix[_self, _mutual1].imag
		if not self.isTriplex:
			if self.line_parameters['file type'] == 'DSS' and self.phases != 'ACB':
				Bmutual2 = None
			else:
				Bmutual2 = self.Ymatrix[_self, _mutual2].imag
		else:
			Bmutual2 = None

	Bself_shunt = self.Yshunt[_self,
							  _self].imag * 0.5 if self.hasShunt else None
	
	if self.line_parameters['file type'] == 'DSS' and self.single_phase != -1:
		Bmutual1_shunt = None
		Bmutual2_shunt = None
	elif self.line_parameters['file type'] == 'DSS' and self.phases != 'ACB':
		Bmutual1_shunt = self.Yshunt[
						 _self, _mutual1].imag * 0.5 if self.hasShunt else None
		Bmutual2_shunt = None
	else:
		Bmutual1_shunt = self.Yshunt[
							_self, _mutual1].imag * 0.5 if self.hasShunt else None
		Bmutual2_shunt = self.Yshunt[
							_self, _mutual2].imag * 0.5 if self.hasShunt else None

	Gself = self.Ymatrix[_self, _self].real
	if self.line_parameters['file type'] == 'DSS' and self.single_phase != -1:
		Gmutual1 = None
		Gmutual2 = None
	else:
		Gmutual1 = self.Ymatrix[_self, _mutual1].real
		if not self.isTriplex:
			if self.line_parameters['file type'] == 'DSS' and self.phases != 'ACB':
				Gmutual2 = None
			else:
				Gmutual2 = self.Ymatrix[_self, _mutual2].real
		else:
			Gmutual2 = None

	# ======================================================================
	#              stamp_phase_ckt self susceptance for the phase
	# ======================================================================
	# Self Series

	# Real
	if Bself:
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeVr_from[_self], nodeVi_from[_self], -Bself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_to[_self], Bself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_self], -Bself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_from[_self], Bself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		# Imaginary Circuit
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_self], Bself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_to[_self], -Bself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_self], Bself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_from[_self], -Bself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual: #(should be the transpose of above)
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_self], -Bself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_from[_self], Bself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_self], -Bself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_to[_self], Bself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			# Imaginary Circuit
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_self], Bself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_from[_self], -Bself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_self], Bself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_to[_self], -Bself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)  
	# Self - Shunt
	if self.hasShunt:
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_from[_self], -Bself_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_self], Bself_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_self], -Bself_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_self], Bself_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_self], -Bself_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_self], Bself_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_self], -Bself_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_self], Bself_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)

	if Gself:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_from[_self], Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_from[_self], Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_to[_self], Gself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_to[_self], Gself, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_to[_self], -Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_to[_self], -Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_from[_self], -Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_from[_self], -Gself, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_from[_self], Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_from[_self], Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_to[_self], Gself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_to[_self], Gself, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_from[_self], -Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_from[_self], -Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_to[_self], -Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			  #
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_to[_self], -Gself, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
	# =====================================================================
	#         # stamp_phase_ckt mutual susceptance for the phases
	# =====================================================================
	# Mutual Series

	# Mutual 1 - Susceptance
	# Real Circuit
	if Bmutual1:
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_from[_mutual1], -Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_to[_mutual1], Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_mutual1], -Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_from[_mutual1], Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		# Imaginary Circuit
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_mutual1], Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_to[_mutual1], -Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_mutual1], Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_from[_mutual1], -Bmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_mutual1], -Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_from[_mutual1], Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_mutual1], -Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_to[_mutual1], Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			# Imaginary Circuit
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_mutual1], Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_from[_mutual1], -Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_mutual1], Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_to[_mutual1], -Bmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
	# Mutual 2 - Susceptance
	# Real Circuit
	if Bmutual2:
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_from[_mutual2], -Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_to[_mutual2], Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_mutual2], -Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_from[_mutual2], Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		# Imaginary Circuit
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_mutual2], Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_to[_mutual2], -Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_mutual2], Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_from[_mutual2], -Bmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_mutual2], -Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_from[_mutual2], Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_mutual2], -Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_to[_mutual2], Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			# Imaginary Circuit
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_mutual2], Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_from[_mutual2], -Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_mutual2], Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_to[_mutual2], -Bmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)

	# Mutual Shunts
	# Mutual 1
	if Bmutual1_shunt:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_from[_mutual1], -Bmutual1_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_mutual1], Bmutual1_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_mutual1], -Bmutual1_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_mutual1], Bmutual1_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		
		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_mutual1], -Bmutual1_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_mutual1], Bmutual1_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_mutual1], -Bmutual1_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_mutual1], Bmutual1_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
				
	if Bmutual2_shunt:
		# Mutual 2
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVi_from[_mutual2], -Bmutual2_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVr_from[_mutual2], Bmutual2_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVi_to[_mutual2], -Bmutual2_shunt,
			Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVr_to[_mutual2], Bmutual2_shunt, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			# Mutual 2
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLr_from[_mutual2], -Bmutual2_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLi_from[_mutual2], Bmutual2_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLr_to[_mutual2], -Bmutual2_shunt,
				Ylin_val, Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLi_to[_mutual2], Bmutual2_shunt, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
	
	if Gmutual1:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_from[_mutual1], Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_from[_mutual1], Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_to[_mutual1], Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_to[_mutual1], Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_to[_mutual1], -Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_to[_mutual1], -Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_from[_mutual1], -Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_from[_mutual1], -Gmutual1, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)

		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_from[_mutual1], Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_from[_mutual1], Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_to[_mutual1], Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_to[_mutual1], Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_from[_mutual1], -Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_from[_mutual1], -Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_to[_mutual1], -Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_to[_mutual1], -Gmutual1, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			
	if Gmutual2:
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_from[_mutual2], Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_from[_mutual2], Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_to[_mutual2], Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_to[_mutual2], Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_from[_self], nodeVr_to[_mutual2], -Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_from[_self], nodeVi_to[_mutual2], -Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVr_to[_self], nodeVr_from[_mutual2], -Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		#
		Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
			nodeVi_to[_self], nodeVi_from[_mutual2], -Gmutual2, Ylin_val,
			Ylin_row, Ylin_col, idx_Y)
		
		if self.stamp_dual:
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_from[_mutual2], Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_from[_mutual2], Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_to[_mutual2], Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_to[_mutual2], Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_to[_self], nodeLr_from[_mutual2], -Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_to[_self], nodeLi_from[_mutual2], -Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLr_from[_self], nodeLr_to[_mutual2], -Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			#
			Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
				nodeLi_from[_self], nodeLi_to[_mutual2], -Gmutual2, Ylin_val,
				Ylin_row, Ylin_col, idx_Y)
			
	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_linear(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, idx_Y):
	# Get the From/To bus nodes to stamp
	nodes_from, nodes_to = self.get_nodes(node, node_key)

	if self.line_parameters['file type'] == 'DSS':
		# In DSS we can have one, two, and three phase lines
		num_of_phases = len(nodes_from['VR']) - 1
		if num_of_phases == 3:
			phases = {0, 1, 2}
		elif num_of_phases == 2:
			phases = {0, 1}
		elif num_of_phases == 1:
			phases = {0}
		
		# Collect nodes	
		nodeVr_from = nodes_from['VR'][0:-1]
		nodeVi_from = nodes_from['VI'][0:-1]
		nodeVr_to = nodes_to['VR'][0:-1]
		nodeVi_to = nodes_to['VI'][0:-1]

		if self.stamp_dual:
			nodeLr_from = nodes_from['LR'][0:-1]
			nodeLi_from = nodes_from['LI'][0:-1]
			nodeLr_to = nodes_to['LR'][0:-1]
			nodeLi_to = nodes_to['LI'][0:-1]
		else:
			nodeLr_from = []
			nodeLi_from = []
			nodeLr_to = []
			nodeLi_to = []

		
	else:
		# Overhead or Underground line: 3 phases
		num_of_phases = 3
		phases = {0, 1, 2}

		# Collect nodes
		nodeVr_from = [nodes_from.VR[0], nodes_from.VR[1], nodes_from.VR[2]]
		nodeVi_from = [nodes_from.VI[0], nodes_from.VI[1], nodes_from.VI[2]]
		nodeVr_to = [nodes_to.VR[0], nodes_to.VR[1], nodes_to.VR[2]]
		nodeVi_to = [nodes_to.VI[0], nodes_to.VI[1], nodes_to.VI[2]]

		if self.stamp_dual:
			nodeLr_from = [nodes_from.LR[0], nodes_from.LR[1], nodes_from.LR[2]]
			nodeLi_from = [nodes_from.LI[0], nodes_from.LI[1], nodes_from.LI[2]]
			nodeLr_to = [nodes_to.LR[0], nodes_to.LR[1], nodes_to.LR[2]]
			nodeLi_to = [nodes_to.LI[0], nodes_to.LI[1], nodes_to.LI[2]]
		else:
			nodeLr_from = []
			nodeLi_from = []
			nodeLr_to = []
			nodeLi_to = []
		
	for _self in range(num_of_phases):
		phases = {0, 1, 2}
		phases.remove(_self)
		_mutual1 = phases.pop()
		_mutual2 = phases.pop()

		Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_phase_ckt(Ylin_val, Ylin_row, Ylin_col, idx_Y, nodeVr_from,
																   nodeVi_from, nodeVr_to, nodeVi_to, _self, _mutual1,
																   _mutual2, nodeLr_from, nodeLi_from, nodeLr_to, nodeLi_to)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_tx_stepping(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, idx_Y, h_factor, g_homotopy, b_homotopy):

	# USING BASIC CODE FROM STAMP_LINEAR + ADDITIONS TO IMPLEMENT TX_STEPPING
	# This function is called directly from Homotopy.stamp_tx_stepping()
	# This function can be used for overhead, underground, and triplex lines.

	nodes_from, nodes_to = self.get_nodes(node, node_key)


    # Determine whether or not there will be 2 or 3 phases
	if self.isTriplex == True:
		# Triplex line: 2 phases
		num_of_self = 2
	else:
		if self.line_parameters['file type'] == 'DSS':
			# In DSS we can have one, two, and three phase lines
			num_of_self = len(nodes_from) - 1

		else:
			# Overhead or Underground line: 3 phases
			num_of_self = 3

		
	
    # Although this should not be called if the homotopy factor, this is a check to prevent
    # it from being needlessly stamped if h_factor = 0
	if h_factor != 0:
		for _self in range(num_of_self):
			nodeVr_from = nodes_from.VR[_self]
			nodeVi_from = nodes_from.VI[_self]
			nodeVr_to = nodes_to.VR[_self]
			nodeVi_to = nodes_to.VI[_self]
			
			if self.stamp_dual:
				nodeLr_from = nodes_from.LR[_self]
				nodeLi_from = nodes_from.LI[_self]
				nodeLr_to = nodes_to.LR[_self]
				nodeLi_to = nodes_to.LI[_self]

			# The new G, B homotopy value are being created by multiplying the value of the homotopy factor by
			# G or B homotopy (can be user-defined input, otherwise the value is 400), and then the original
			# real or imaginary value from the admittance matrix, Y. This ensures that the new homotopy value is much
			# greater than the original when started from the trivial solution. The values change but the stamps otherwise stay the same.
			if self.Ymatrix[_self,_self].real:
				g_tx_stepping = h_factor*g_homotopy*self.Ymatrix[_self,_self].real
				
				# Real - From
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_from, nodeVr_from, g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_from, nodeVr_to, -g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				# Real - To
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_to, nodeVr_from, -g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_to, nodeVr_to, g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				# Imag - From
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_from, nodeVi_from, g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_from, nodeVi_to, -g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				# Imag - To
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_to, nodeVi_from, -g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_to, nodeVi_to, g_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				if self.stamp_dual:
					# Real - From
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_from, nodeLr_from, g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_from, nodeLr_to, -g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					
					# Imag - From
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_from, nodeLi_from, g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_from, nodeLi_to, -g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y) 
					# Real - To
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_to, nodeLr_to, g_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_to, nodeLr_from, -g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					# Imag - To
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_to, nodeLi_to, g_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_to, nodeLi_from, -g_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
			
			if self.Ymatrix[_self,_self].imag:
				b_tx_stepping = h_factor*b_homotopy*self.Ymatrix[_self,_self].imag
			
				# Real - From
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_from, nodeVi_from, -b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_from, nodeVi_to, b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				# Real - From
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_to, nodeVi_from, b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVr_to, nodeVi_to, -b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)       

				# Imag - From
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_from, nodeVr_from, b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_from, nodeVr_to, -b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				# Imag - To
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_to, nodeVr_from, -b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeVi_to, nodeVr_to, b_tx_stepping, Ylin_val, Ylin_row, Ylin_col, idx_Y)
				
				if self.stamp_dual:
					# Imaginary - From
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_from, nodeLr_from, -b_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)                    
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_from, nodeLr_to, b_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)
					
					# Imaginary - To
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_to, nodeLr_from, b_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_to, nodeLr_to, -b_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)

					# Real - From
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_from, nodeLi_from, b_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_from, nodeLi_to, -b_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y) 
					
					# Real - To
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_to, nodeLi_from, -b_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_to, nodeLi_to, b_tx_stepping, Ylin_val, Ylin_row,
						Ylin_col, idx_Y)
 				
			if self.hasShunt:
				# Shunts are decreased in homotopy to accommodate the virtual shorts
				b_shunt_tx_stepping = -h_factor*b_homotopy*self.Yshunt[_self,_self].imag*0.5
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
					nodeVr_from, nodeVi_from, -b_shunt_tx_stepping, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
					nodeVi_from, nodeVr_from, b_shunt_tx_stepping, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
					nodeVr_to, nodeVi_to, -b_shunt_tx_stepping, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
				#
				Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
					nodeVi_to, nodeVr_to, b_shunt_tx_stepping, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
				
				if self.stamp_dual:
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_from, nodeLr_from, -b_shunt_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_from, nodeLi_from, b_shunt_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLi_to, nodeLr_to, -b_shunt_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)
					#
					Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
						nodeLr_to, nodeLi_to, b_shunt_tx_stepping, Ylin_val,
						Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_neutral(self, int neutral_eqn1, int neutral_eqn2, int from_nodeVr_N,
					int to_nodeVr_N, int from_nodeVi_N, int to_nodeVi_N, Ylin_val,
					Ylin_row, Ylin_col, idx_Y):
	# Voltage sources
	# Real Neutral Connect
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		neutral_eqn1, from_nodeVr_N, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		neutral_eqn1, to_nodeVr_N, -1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	# Imag Neutral Connect
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		neutral_eqn2, from_nodeVi_N, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		neutral_eqn2, to_nodeVi_N, -1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	# Current Sources
	# Real Neutral Connect
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		from_nodeVr_N, neutral_eqn1, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		to_nodeVr_N, neutral_eqn1, -1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	# Imag Neutral Connect
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		from_nodeVi_N, neutral_eqn2, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
	Ylin_val, Ylin_row, Ylin_col, idx_Y = stampY(
		to_nodeVi_N, neutral_eqn2, -1, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y


cpdef stamp_linear_triplex(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, idx_Y):
	# Get the From/To bus nodes to stamp
	nodes_from, nodes_to = self.get_nodes(node, node_key)

	# Collect nodes
	nodeVr_from = [nodes_from.VR[0], nodes_from.VR[1]]
	nodeVi_from = [nodes_from.VI[0], nodes_from.VI[1]]
	nodeVr_to = [nodes_to.VR[0], nodes_to.VR[1]]
	nodeVi_to = [nodes_to.VI[0], nodes_to.VI[1]]
	
	if self.stamp_dual:
		nodeLr_from = [nodes_from.LR[0], nodes_from.LR[1]]
		nodeLi_from = [nodes_from.LI[0], nodes_from.LI[1]]
		nodeLr_to = [nodes_to.LR[0], nodes_to.LR[1]]
		nodeLi_to = [nodes_to.LI[0], nodes_to.LI[1]]
	else:
		nodeLr_from = []
		nodeLi_from = []
		nodeLr_to = []
		nodeLi_to = []

	for _self in range(2):
		phases = {0, 1}
		phases.remove(_self)
		_mutual1 = phases.pop()

		Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_phase_ckt(Ylin_val, Ylin_row, Ylin_col, idx_Y, nodeVr_from,
																   nodeVi_from, nodeVr_to, nodeVi_to, _self, _mutual1,
																   None, nodeLr_from, nodeLi_from, nodeLr_to, nodeLi_to)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

