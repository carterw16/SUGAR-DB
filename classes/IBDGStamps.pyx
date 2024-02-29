#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3
import numpy as np
cimport numpy as np
ctypedef np.float64_t DTYPE_t

def stampY(int i, int j, double val, Ynlin_val, Ynlin_row, Ynlin_col, int idx):
	Ynlin_val[idx] = val
	Ynlin_row[idx] = i
	Ynlin_col[idx] = j
	idx += 1
	return idx

def stampJ(int i, double val, Jnlin_val, Jnlin_row, int idx):
	Jnlin_val[idx] = val
	Jnlin_row[idx] = i
	idx += 1
	return idx

def stamp_ibdg_partials(self,
                        dict Vseq,
                        phase,
                        dIdgrPos,
                        dIdgiPos,
                        dIdgrNeg,
                        dIdgiNeg,
                        IdgrPos_hist,
                        IdgiPos_hist,
                        IdgrNeg_hist,
                        IdgiNeg_hist,
                        np.ndarray[DTYPE_t, ndim=1] Ynlin_val,
                        np.ndarray[np.int64_t, ndim=1] Ynlin_row,
                        np.ndarray[np.int64_t, ndim=1] Ynlin_col,
                        np.ndarray[DTYPE_t, ndim=1] Jnlin_val,
                        np.ndarray[np.int64_t, ndim=1] Jnlin_row,
                        int idx_Y,
                        int idx_J):

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
	idx_Y = stampY(nodePos_Vr, nodePos_Vr, dIdgrPos.dVposr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vr, nodePos_Vi, dIdgrPos.dVposi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vr, nodeNeg_Vr, dIdgrPos.dVnegr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vr, nodeNeg_Vi, dIdgrPos.dVnegi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	# Idgi Pos Stamps
	idx_Y = stampY(nodePos_Vi, nodePos_Vr, dIdgiPos.dVposr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vi, nodePos_Vi, dIdgiPos.dVposi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vi, nodeNeg_Vr, dIdgiPos.dVnegr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodePos_Vi, nodeNeg_Vi, dIdgiPos.dVnegi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	# Idgr Neg Stamps
	idx_Y = stampY(nodeNeg_Vr, nodePos_Vr, dIdgrNeg.dVposr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vr, nodePos_Vi, dIdgrNeg.dVposi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vr, nodeNeg_Vr, dIdgrNeg.dVnegr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vr, nodeNeg_Vi, dIdgrNeg.dVnegi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	# Idgi Pos Stamps
	idx_Y = stampY(nodeNeg_Vi, nodePos_Vr, dIdgiNeg.dVposr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vi, nodePos_Vi, dIdgiNeg.dVposi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vi, nodeNeg_Vr, dIdgiNeg.dVnegr, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)
	idx_Y = stampY(nodeNeg_Vi, nodeNeg_Vi, dIdgiNeg.dVnegi, Ynlin_val,
	                    Ynlin_row, Ynlin_col, idx_Y)

	#  Q Stamps
	if not self.fixed_Q:
		idx_Y = stampY(nodePos_Vr, node_Q3, dIdgrPos.dQ3, Ynlin_val,
		                    Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = stampY(nodePos_Vi, node_Q3, dIdgiPos.dQ3, Ynlin_val,
		                    Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = stampY(nodeNeg_Vr, node_Q3, dIdgrNeg.dQ3, Ynlin_val,
		                    Ynlin_row, Ynlin_col, idx_Y)
		idx_Y = stampY(nodeNeg_Vi, node_Q3, dIdgiNeg.dQ3, Ynlin_val,
		                    Ynlin_row, Ynlin_col, idx_Y)
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

	idx_J = stampJ(nodePos_Vr, -IdgrPos, Jnlin_val, Jnlin_row, idx_J)
	idx_J = stampJ(nodePos_Vi, -IdgiPos, Jnlin_val, Jnlin_row, idx_J)
	idx_J = stampJ(nodeNeg_Vr, -IdgrNeg, Jnlin_val, Jnlin_row, idx_J)
	idx_J = stampJ(nodeNeg_Vi, -IdgiNeg, Jnlin_val, Jnlin_row, idx_J)

	return idx_Y, idx_J

def stamp_Q3_partials(self,
                        Vabc,
                        dFQ3,
                        FQ3_hist,
                        np.ndarray[DTYPE_t, ndim=1] Ynlin_val,
                        np.ndarray[np.int64_t, ndim=1] Ynlin_row,
                        np.ndarray[np.int64_t, ndim=1] Ynlin_col,
                        np.ndarray[DTYPE_t, ndim=1] Jnlin_val,
                        np.ndarray[np.int64_t, ndim=1] Jnlin_row,
                        int idx_Y,
                        int idx_J,
                        bint homotopy_enabled,
                        double h_factor):

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

	# Va_mag = Vabc.Vmag[0]
	# Vb_mag = Vabc.Vmag[1]
	# Vc_mag = Vabc.Vmag[2]

	Q3 = Vabc.Q[0]

	if homotopy_enabled and self.Qlim_type != 'Q-V_Switch':
		h_bound = (1 + (self.h_factor_init + 1) * h_factor)
	else:
		h_bound = 1
	# h_bound = 1

	if self.flag_Qmax and self.flag_Qlim:
		Qmax = self.Qmax * h_bound
		idx_Y = self.stamp_Y(node_Q3, node_Q3, 1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stamp_J(node_Q3, -Qmax, Jnlin_val, Jnlin_row, idx_J)
	elif self.flag_Qmin and self.flag_Qlim:
		Qmin = self.Qmin * h_bound
		idx_Y = self.stamp_Y(node_Q3, node_Q3, 1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		idx_J = self.stamp_J(node_Q3, -Qmin, Jnlin_val, Jnlin_row, idx_J)
	else:
		# if self.phases & 0x01 and self.phases & 0x02 and self.phases & 0x04:
		#     idx_Y = self.stampY(node_Q3, node_Var, dFQ3.dVr[0], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#     idx_Y = self.stampY(node_Q3, node_Vai, dFQ3.dVi[0], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#     idx_Y = self.stampY(node_Q3, node_Vbr, dFQ3.dVr[1], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#     idx_Y = self.stampY(node_Q3, node_Vbi, dFQ3.dVi[1], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#     idx_Y = self.stampY(node_Q3, node_Vcr, dFQ3.dVr[2], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#     idx_Y = self.stampY(node_Q3, node_Vci, dFQ3.dVi[2], Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)
		#
		#     idx_Y = self.stampY(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val,
		#                         Ynlin_row, Ynlin_col, idx_Y)

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

			# idx_Y = self.stampY(self.nodeA_Vmag2, nodeA_Vr, dVmag2_dVr,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeA_Vmag2, nodeA_Vi, dVmag2_dVi,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeA_Vmag2, self.nodeA_Vmag2, dVmag2_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeA_Vmag2, -_Vmag2, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# Vmag = Vabc.Vmag[0]
			#
			# Vmag_hist = np.sqrt(Vmag2) - Vmag
			#
			# dVmag_dVmag2 = (1 / (2 * np.sqrt(Vmag2)))
			# dVmag_dVmag = -1
			#
			# _Vmag = Vmag_hist - dVmag_dVmag * Vmag - dVmag_dVmag2 * Vmag2
			#
			# idx_Y = self.stampY(self.nodeA_Vmag, self.nodeA_Vmag2, dVmag_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeA_Vmag, self.nodeA_Vmag, dVmag_dVmag,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeA_Vmag, -_Vmag, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# # Stamp Q3 Va partials # #
			# idx_Y = self.stampY(node_Q3, self.nodeA_Vmag, dFQ3.dV[0], Ynlin_val,
			#                     Ynlin_row, Ynlin_col, idx_Y)

			idx_Y = self.stamp_Y(node_Q3, nodeA_Vr, dFQ3.dVr[0], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(node_Q3, nodeA_Vi, dFQ3.dVi[0], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

		# idx_Y = self.stampY(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val,
		#                     Ynlin_row, Ynlin_col, idx_Y)
		# FQ3 = FQ3_hist - dFQ3.dVr[0] * Var - dFQ3.dVi[0] * Vai - dFQ3.dQ * Q3
		# idx_J = self.stampJ(node_Q3, -FQ3, Jnlin_val, Jnlin_row, idx_J)

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

			# idx_Y = self.stampY(self.nodeB_Vmag2, nodeB_Vr, dVmag2_dVr,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeB_Vmag2, nodeB_Vi, dVmag2_dVi,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeB_Vmag2, self.nodeB_Vmag2, dVmag2_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeB_Vmag2, -_Vmag2, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# Vmag = Vabc.Vmag[1]
			#
			# Vmag_hist = np.sqrt(Vmag2) - Vmag
			#
			# dVmag_dVmag2 = (1 / (2 * np.sqrt(Vmag2)))
			# dVmag_dVmag = -1
			#
			# _Vmag = Vmag_hist - dVmag_dVmag * Vmag - dVmag_dVmag2 * Vmag2
			#
			# idx_Y = self.stampY(self.nodeB_Vmag, self.nodeB_Vmag2, dVmag_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeB_Vmag, self.nodeB_Vmag, dVmag_dVmag,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeB_Vmag, -_Vmag, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# # Stamp Q3 Vb partials # #
			# idx_Y = self.stampY(node_Q3, self.nodeB_Vmag, dFQ3.dV[1], Ynlin_val,
			#                     Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(node_Q3, nodeB_Vr, dFQ3.dVr[1], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(node_Q3, nodeB_Vi, dFQ3.dVi[1], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		# idx_Y = self.stampY(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val,
		#                     Ynlin_row, Ynlin_col, idx_Y)
		# FQ3 = FQ3_hist - dFQ3.dVr[1] * Vbr - dFQ3.dVi[1] * Vbi - dFQ3.dQ * Q3
		# idx_J = self.stampJ(node_Q3, -FQ3, Jnlin_val, Jnlin_row, idx_J)

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

			# idx_Y = self.stampY(self.nodeC_Vmag2, nodeC_Vr, dVmag2_dVr,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeC_Vmag2, nodeC_Vi, dVmag2_dVi,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeC_Vmag2, self.nodeC_Vmag2, dVmag2_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeC_Vmag2, -_Vmag2, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# Vmag = Vabc.Vmag[2]
			#
			# Vmag_hist = np.sqrt(Vmag2) - Vmag
			#
			# dVmag_dVmag2 = (1 / (2 * np.sqrt(Vmag2)))
			# dVmag_dVmag = -1
			#
			# _Vmag = Vmag_hist - dVmag_dVmag * Vmag - dVmag_dVmag2 * Vmag2
			#
			# idx_Y = self.stampY(self.nodeC_Vmag, self.nodeC_Vmag2, dVmag_dVmag2,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			# idx_Y = self.stampY(self.nodeC_Vmag, self.nodeC_Vmag, dVmag_dVmag,
			#                     Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			# idx_J = self.stampJ(self.nodeC_Vmag, -_Vmag, Jnlin_val, Jnlin_row,
			#                     idx_J)

			# # Stamp Q3 Vc partials # #
			# idx_Y = self.stampY(node_Q3, self.nodeC_Vmag, dFQ3.dV[2], Ynlin_val,
			#                     Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(node_Q3, nodeC_Vr, dFQ3.dVr[2], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			idx_Y = self.stamp_Y(node_Q3, nodeC_Vi, dFQ3.dVi[2], Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		# idx_Y = self.stampY(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val,
		#                     Ynlin_row, Ynlin_col, idx_Y)
		# FQ3 = FQ3_hist - dFQ3.dVr[2] * Vcr - dFQ3.dVi[2] * Vci - dFQ3.dQ * Q3
		# idx_J = self.stampJ(node_Q3, -FQ3, Jnlin_val, Jnlin_row, idx_J)

	idx_Y = self.stamp_Y(node_Q3, node_Q3, dFQ3.dQ, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

	FQ3 = FQ3_hist - dFQ3.dVr[0] * Var - dFQ3.dVi[0] * Vai - dFQ3.dVr[1] * Vbr - \
	      dFQ3.dVi[1] * Vbi - dFQ3.dVr[2] * Vcr - dFQ3.dVi[2] * Vci - dFQ3.dQ * Q3

	# FQ3 = FQ3_hist - dFQ3.dV[0] * Va_mag - dFQ3.dV[1] * Vb_mag - dFQ3.dV[2] * Vc_mag - dFQ3.dQ * Q3

	idx_J = self.stamp_J(node_Q3, -FQ3, Jnlin_val, Jnlin_row, idx_J)

	return idx_Y, idx_J