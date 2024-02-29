#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

cpdef stampY(i, j, val, Ynlin_val, Ynlin_row, Ynlin_col, idx):
	Ynlin_val[idx] = val
	Ynlin_row[idx] = i
	Ynlin_col[idx] = j
	idx += 1
	return idx

cpdef stampJ(i, val, Jnlin_val, Jnlin_row, idx):
	Jnlin_val[idx] = val
	Jnlin_row[idx] = i
	idx += 1
	return idx

cpdef stamp_line_currents(self,
                        Ynlin_val,
                        Ynlin_row,
                        Ynlin_col,
                        Jnlin_val,
                        Jnlin_row,
                        idx_Y,
                        idx_J,
                        _self,
                        V_self_from,
                        V_self_to,
                        _mutual1,
                        V_mutual1_from,
                        V_mutual1_to,
                        Imag2,
                        Ir,
                        Ii,
                        V_mutual2_from=None,
                        V_mutual2_to=None,
                        _mutual2=None):
	# Get the kth from and to real and imaginary voltages
	Vr_self_from, Vi_self_from = V_self_from.real, V_self_from.imag
	Vr_self_to, Vi_self_to = V_self_to.real, V_self_to.imag

	Vr_mutual1_from, Vi_mutual1_from = V_mutual1_from.real, V_mutual1_from.imag
	Vr_mutual1_to, Vi_mutual1_to = V_mutual1_to.real, V_mutual1_to.imag

	if not self.isTriplex:
		Vr_mutual2_from, Vi_mutual2_from = V_mutual2_from.real, V_mutual2_from.imag
		Vr_mutual2_to, Vi_mutual2_to = V_mutual2_to.real, V_mutual2_to.imag
	else:
		Vr_mutual2_from, Vi_mutual2_from = 0, 0
		Vr_mutual2_to, Vi_mutual2_to = 0, 0

	# Compute the voltage drop across the lines
	Vr_self = Vr_self_from - Vr_self_to
	Vi_self = Vi_self_from - Vi_self_to

	Vr_mutual1 = Vr_mutual1_from - Vr_mutual1_to
	Vi_mutual1 = Vi_mutual1_from - Vi_mutual1_to

	Vr_mutual2 = Vr_mutual2_from - Vr_mutual2_to
	Vi_mutual2 = Vi_mutual2_from - Vi_mutual2_to

	# Collect the self and mutual series admittance values
	Bself = self.Y_series[_self, _self].imag
	Bmutual1 = self.Y_series[_self, _mutual1].imag
	Bmutual2 = self.Y_series[_self,
	                         _mutual2].imag if not self.isTriplex else 0.0

	Gself = self.Y_series[_self, _self].real
	Gmutual1 = self.Y_series[_self, _mutual1].real
	Gmutual2 = self.Y_series[_self,
	                         _mutual2].real if not self.isTriplex else 0.0

	# Collect the self and mutual shunt susceptance values
	Bself_shunt = self.Y_shunt[_self,
	                           _self].imag * 0.5 if self.hasShunt else 0.0
	Bmutual1_shunt = self.Y_shunt[
		                 _self, _mutual1].imag * 0.5 if self.hasShunt else 0.0
	Bmutual2_shunt = self.Y_shunt[
		                 _self, _mutual2].imag * 0.5 if self.hasShunt else 0.0

	# # # Compute shunt, series, line, and magnitude ^2 historical currents and partials # # #
	if self.type == 'XfmrCT' or self.type == 'Xfmr':
		Ir_shunt_hist = -(Bself_shunt * Vi_self_from)
		Ii_shunt_hist = Bself_shunt * Vr_self_from
		Ir_series_hist = Vr_self * Gself - Bself * Vi_self
		Ii_series_hist = Vi_self * Gself + Bself * Vr_self

	else:
		Ir_shunt_hist = -(Bself_shunt * Vi_self_from +
		                  Bmutual1_shunt * Vi_mutual1_from +
		                  Bmutual2_shunt * Vi_mutual2_from)
		Ii_shunt_hist = Bself_shunt * Vr_self_from + Bmutual1_shunt * Vr_mutual1_from + Bmutual2_shunt * \
		                Vr_mutual2_from

		Ir_series_hist = Vr_self * Gself - Bself * Vi_self + Vr_mutual1 * Gmutual1 - Vi_mutual1 * Bmutual1 + \
		                 Vr_mutual2 * Gmutual2 - Vi_mutual2 * Bmutual2

		Ii_series_hist = Vi_self * Gself + Bself * Vr_self + Gmutual1 * Vi_mutual1 + Bmutual1 * Vi_mutual1 + \
		                 Gmutual2 * Vi_mutual2 + Bmutual2 * Vr_mutual2

	# # Ir # #
	Ir_series_hist + Ir_shunt_hist - Ir

	dIr_hist_dVrsf = Gself
	dIr_hist_dVrst = -Gself
	dIr_hist_dVisf = -Bself_shunt - Bself
	dIr_hist_dVist = Bself

	dIr_hist_dVr1f = Gmutual1
	dIr_hist_dVr1t = -Gmutual1
	dIr_hist_dVi1f = -Bmutual1_shunt - Bmutual1
	dIr_hist_dVi1t = Bmutual1

	dIr_hist_dVr2f = Gmutual2
	dIr_hist_dVr2t = -Gmutual2
	dIr_hist_dVi2f = -Bmutual2_shunt - Bmutual2
	dIr_hist_dVi2t = Bmutual2

	dIr_hist_dIr = -1

	# Ir Stamps
	idx_Y = stampY(self.node_Ir, self.node_self_Vr_from,
	                    dIr_hist_dVrsf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    idx_Y)
	idx_Y = stampY(self.node_Ir, self.node_self_Vr_to, dIr_hist_dVrst,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(self.node_Ir, self.node_self_Vi_from,
	                    dIr_hist_dVisf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    idx_Y)
	idx_Y = stampY(self.node_Ir, self.node_self_Vi_to, dIr_hist_dVist,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	if self.type != 'Xfmr' or self.type != 'XfmrCT':
		idx_Y = stampY(self.node_Ir, self.node_mutual1_Vr_from,
		                    dIr_hist_dVr1f, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ir, self.node_mutual1_Vr_to,
		                    dIr_hist_dVr1t, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ir, self.node_mutual1_Vi_from,
		                    dIr_hist_dVi1f, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ir, self.node_mutual1_Vi_to,
		                    dIr_hist_dVi1t, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		if not self.isTriplex:
			idx_Y = stampY(self.node_Ir, self.node_mutual2_Vr_from,
			                    dIr_hist_dVr2f, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ir, self.node_mutual2_Vr_to,
			                    dIr_hist_dVr2t, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ir, self.node_mutual2_Vi_from,
			                    dIr_hist_dVi2f, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ir, self.node_mutual2_Vi_to,
			                    dIr_hist_dVi2t, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)

	idx_Y = stampY(self.node_Ir, self.node_Ir, dIr_hist_dIr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	# # Ii # #
	Ii_series_hist + Ii_shunt_hist - Ii

	dIi_hist_dVrsf = Bself + Bself_shunt
	dIi_hist_dVrst = -Bself
	dIi_hist_dVisf = Gself
	dIi_hist_dVist = -Gself

	dIi_hist_dVr1f = Bmutual1 + Bmutual1_shunt
	dIi_hist_dVr1t = -Bmutual1
	dIi_hist_dVi1f = Gmutual1
	dIi_hist_dVi1t = -Gmutual1

	dIi_hist_dVr2f = Bmutual2 + Bmutual2_shunt
	dIi_hist_dVr2t = -Bmutual2
	dIi_hist_dVi2f = Gmutual2
	dIi_hist_dVi2t = -Gmutual2

	dIi_hist_dIi = -1

	# Ii Stamps
	idx_Y = stampY(self.node_Ii, self.node_self_Vr_from,
	                    dIi_hist_dVrsf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    idx_Y)
	idx_Y = stampY(self.node_Ii, self.node_self_Vr_to, dIi_hist_dVrst,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(self.node_Ii, self.node_self_Vi_from,
	                    dIi_hist_dVisf, Ynlin_val, Ynlin_row, Ynlin_col,
	                    idx_Y)
	idx_Y = stampY(self.node_Ii, self.node_self_Vi_to, dIi_hist_dVist,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	if self.type != 'Xfmr' or self.type != 'XfmrCT':
		idx_Y = stampY(self.node_Ii, self.node_mutual1_Vr_from,
		                    dIi_hist_dVr1f, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ii, self.node_mutual1_Vr_to,
		                    dIi_hist_dVr1t, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ii, self.node_mutual1_Vi_from,
		                    dIi_hist_dVi1f, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		idx_Y = stampY(self.node_Ii, self.node_mutual1_Vi_to,
		                    dIi_hist_dVi1t, Ynlin_val, Ynlin_row, Ynlin_col,
		                    idx_Y)
		if not self.isTriplex:
			idx_Y = stampY(self.node_Ii, self.node_mutual2_Vr_from,
			                    dIi_hist_dVr2f, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ii, self.node_mutual2_Vr_to,
			                    dIi_hist_dVr2t, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ii, self.node_mutual2_Vi_from,
			                    dIi_hist_dVi2f, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)
			idx_Y = stampY(self.node_Ii, self.node_mutual2_Vi_to,
			                    dIi_hist_dVi2t, Ynlin_val, Ynlin_row,
			                    Ynlin_col, idx_Y)

	idx_Y = stampY(self.node_Ii, self.node_Ii, dIi_hist_dIi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	# # Imag2 # #
	Imag2_hist = Ir ** 2 + Ii ** 2 - Imag2

	dImag2_dIr = 2 * Ir
	dImag2_dIi = 2 * Ii
	dImag2_dImag2 = -1

	_Imag2 = Imag2_hist - dImag2_dIr * Ir - dImag2_dIi * Ii - dImag2_dImag2 * Imag2

	idx_Y = stampY(self.node_Imag2, self.node_Ir, dImag2_dIr,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(self.node_Imag2, self.node_Ii, dImag2_dIi,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(self.node_Imag2, self.node_Imag2, dImag2_dImag2,
	                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	idx_J = stampJ(self.node_Imag2, -_Imag2, Jnlin_val, Jnlin_row,
	                    idx_J)

	return idx_Y, idx_J
