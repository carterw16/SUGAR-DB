#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3
import numpy as np
cimport numpy as np

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


cpdef _phase_calculations(self,
                          double _cP,
                          double _cQ,
                          double _cG,
                          double _cB,
                          double _cIr,
                          double _cIi,
                          int node_Vr_from,
                          int node_Vi_from,
                          int node_Vr_to,
                          int node_Vi_to,
                          int Vr_node_from,
                          int Vi_node_from,
                          int Vr_node_to,
                          int Vi_node_to,
                          Ynlin_val,
                          Ynlin_row,
                          Ynlin_col,
                          Jnlin_val,
                          Jnlin_row,
                          idx_Y,
                          idx_J):
	"""
			Stamps the partial derivatives and their phase changes into the admittance (Y) matrix.

			:param _cP: Constant real power.
			:param _cQ: Constant reactive power.
			:param _cG: Constant conductance.
			:param _cB: Constant susceptance.
			:param _cIr: Constant real current.
			:param _cIi: Constant imaginary current.
			:param node_Vr_from: From node of the real circuit.
			:param node_Vi_from: From node of the imaginary circuit.
			:param node_Vr_to: To node of the real circuit.
			:param node_Vi_to: To node of the imaginary circuit.
			:param Vr_node_from: Voltage of the from node in the real circuit.
			:param Vi_node_from: Voltage of the from node in the imaginary circuit.
			:param Vr_node_to: Voltage of the to node in the real circuit.
			:param Vi_node_to: Voltage of the to node in the imaginary circuit.
			:param Ynlin_val: The value of the partial derivative to be stamped into the admittance matrix.
			:param Ynlin_row: The admittance matrix row.
			:param Ynlin_col: The admittance matrix column.
			:param Jnlin_val: The value of the partial derivative to be stamped into the excitation vector.
			:param Jnlin_row: The row of the excitation vector.
			:param _isGnd: Indicated whether load is grounded.
			:return: None
		"""

	if not self.isGnd:
		# Calculate partials for PQ loads
		Vr_across = Vr_node_from - Vr_node_to
		Vi_across = Vi_node_from - Vi_node_to

		# Calculate Partials from the voltage difference between from and to node
		(_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVrl, _dIildVil) = \
			self.calc_partials(_cP, _cQ, Vr_across, Vi_across)

		# 4 real and 4 imaginary stamps for the from circuit
		idx_Y = stampY(node_Vr_from, node_Vr_from, _dIrldVrl + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)  # G
		#
		idx_Y = stampY(node_Vr_from, node_Vr_to, -_dIrldVrl - _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vr_from, node_Vi_from, _dIrldVil - _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vr_from, node_Vi_to, -_dIrldVil + _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vi_from, _dIildVil + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vi_to, -_dIildVil - _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vr_from, _dIildVrl + _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vr_to, -_dIildVrl - _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		# 4 real and imaginary circuits for the to circuit
		idx_Y = stampY(node_Vr_to, node_Vr_from, -_dIrldVrl - _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vr_to, node_Vr_to, _dIrldVrl + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vr_to, node_Vi_from, -_dIrldVil + _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vr_to, node_Vi_to, _dIrldVil - _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_to, node_Vi_from, -_dIildVil - _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_to, node_Vi_to, _dIildVil + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_to, node_Vr_from, -_dIildVrl - _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_to, node_Vr_to, _dIildVrl + _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Independent current sources for from and to branches
		idx_J = stampJ(
			node_Vr_from,
			-((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
			  _cIr), Jnlin_val, Jnlin_row, idx_J)
		#
		idx_J = stampJ(
			node_Vi_from,
			-((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
			  _cIi), Jnlin_val, Jnlin_row, idx_J)
		#
		idx_J = stampJ(
			node_Vr_to,
			((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
			 _cIr), Jnlin_val, Jnlin_row, idx_J)
		#
		idx_J = stampJ(
			node_Vi_to,
			((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
			 _cIi), Jnlin_val, Jnlin_row, idx_J)
	else:
		# Calculate partials for PQ loads
		Vr_across = Vr_node_from
		Vi_across = Vi_node_from

		# Calculate Partials from the voltage difference between from and to node
		(_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVrl, _dIildVil) = \
			self.calc_partials(_cP, _cQ, Vr_across, Vi_across)

		# 4 real and 4 imaginary stamps for the from circuit
		idx_Y = stampY(node_Vr_from, node_Vr_from, _dIrldVrl + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)  # G
		#
		idx_Y = stampY(node_Vr_from, node_Vi_from, _dIrldVil - _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vi_from, _dIildVil + _cG,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		#
		idx_Y = stampY(node_Vi_from, node_Vr_from, _dIildVrl + _cB,
		                    Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# Independent current sources for from and to branches
		idx_J = stampJ(
			node_Vr_from,
			-((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
			  _cIr), Jnlin_val, Jnlin_row, idx_J)
		#
		idx_J = stampJ(
			node_Vi_from,
			-((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
			  _cIi), Jnlin_val, Jnlin_row, idx_J)

	return idx_Y, idx_J
