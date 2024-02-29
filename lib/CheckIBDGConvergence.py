import numpy as np
from itertools import count


class CheckIBDGConvergence:

	def __init__(self):
		self._converge_counter = count(0)
		self._check_count_limit = 20
		self.V_prev = None

	def check_ibdg_convergence(self,
	                           err_max,
	                           tol,
	                           node,
	                           V,
	                           V_prev_dict,
	                           homotopy_enabled,
	                           h_factor,
	                           past_h_factor,
	                           h_step,
	                           ibdgs):
		for ele in ibdgs:
			V = self.apply_limits(ele, node, V)
		# if not homotopy_enabled:
		# 	for ele in ibdgs:
		# 		V = self.apply_limits(ele, node, V)
		# 	if self._converge_counter.__next__() > self._check_count_limit and err_max > tol or err_max > 1e5:
		# 		self.V_prev = V_prev_dict[-1]
		# 		V = np.copy(self.V_prev)
		# 		self._converge_counter = count(0)
		#
		# else:
		# 	if err_max > tol and h_factor != 1:
		# 		converge_counter = self._converge_counter.__next__()
		# 		reset_IBDG_conv = False
		# 		num_steps = 1
		#
		# 		if converge_counter < 10:
		# 			reset_IBDG_conv = False
		# 		elif err_max > 5e-1 and converge_counter > 10:
		# 			reset_IBDG_conv = True
		# 		elif err_max > 5e-2 and converge_counter > 15:
		# 			reset_IBDG_conv = True
		# 		elif err_max > 5e-3 and converge_counter > 20:
		# 			reset_IBDG_conv = True
		# 		elif converge_counter > 25:
		# 			reset_IBDG_conv = True
		#
		# 		if reset_IBDG_conv:
		# 			h_factor, h_step, V = self.reset_IBDGs(ibdgs, node, past_h_factor, h_step, V_prev_dict, num_steps)

		return h_factor, h_step, V

	def reset_IBDGs(self, ibdgs, node, past_h_factor, h_step, V_prev_dict, num_steps):
		h_factor = past_h_factor[-1]
		h_step /= 5
		idx = len(past_h_factor) - num_steps
		Vsol = np.copy(V_prev_dict[idx])

		for ele in ibdgs:
			Vsol = self.apply_limits(ele, node, Vsol)

		return h_factor, h_step, Vsol

	def apply_limits(self, ibdg, node, Vsol):

		Vabc_k1 = ibdg.get_nodes_states(node, node_key, Vsol)
		Q3_k1 = Vabc_k1.Q[0]
		Q3_k = ibdg.Q3_prev

		Var = Vabc_k1.VR[0]
		Vbr = Vabc_k1.VR[1]
		Vcr = Vabc_k1.VR[2]

		Vai = Vabc_k1.VI[0]
		Vbi = Vabc_k1.VI[1]
		Vci = Vabc_k1.VI[2]

		Va = np.abs(Var + 1j * Vai)
		Vb = np.abs(Vbr + 1j * Vbi)
		Vc = np.abs(Vcr + 1j * Vci)

		Va_pu = Va / ibdg.Vnom_ph
		Vb_pu = Vb / ibdg.Vnom_ph
		Vc_pu = Vc / ibdg.Vnom_ph

		if ibdg.phases & 0x1 == 1:
			Vset = Va_pu
		elif ibdg.phases & 0x2 == 2:
			Vset = Vb_pu
		else:
			Vset = Vc_pu
		# V_pu = np.array([Va_pu, Vb_pu, Vc_pu])
		# Va_diff = abs(Va_pu - Vset_pu)
		# Vb_diff = abs(Vb_pu - Vset_pu)
		# Vc_diff = abs(Vc_pu - Vset_pu)
		# err = np.array([Va_diff, Vb_diff, Vc_diff])
		# err_max = np.amax(err)
		# tol = 1e-2

		if Q3_k1 == ibdg.Qmax:
			if Vset < ibdg.V4:
				ibdg.flag_Qlim = False
				ibdg.flag_Qmax = False

		elif Q3_k1 == ibdg.Qmin:
			if Vset > ibdg.V1:
				ibdg.flag_Qlim = False
				ibdg.flag_Qmin = False

		elif Q3_k1 < ibdg.Qmax:
			Q3_k1 = ibdg.Qmax
			Vsol[ibdg.node_Q3] = Q3_k1
			ibdg.flag_Qlim = True
			ibdg.flag_Qmax = True

		elif Q3_k1 > ibdg.Qmin:
			Q3_k1 = ibdg.Qmin
			Vsol[ibdg.node_Q3] = Q3_k1
			ibdg.flag_Qlim = True
			ibdg.flag_Qmin = True

		elif Q3_k1 == 0:
			if Vset < ibdg.V2:
				ibdg.flag_Qlim = True
				ibdg.flag_Qmin = True
			elif Vset > ibdg.V3:
				ibdg.flag_Qlim = True
				ibdg.flag_Qmax = True
			else:
				ibdg.flag_Qlim = False
				ibdg.flag_Qmin = False
				ibdg.flag_Qmax = False

		else:
			if Vset < ibdg.V1:
				ibdg.flag_Qlim = True
				ibdg.flag_Qmin = True
			elif Vset > ibdg.V4:
				ibdg.flag_Qlim = True
				ibdg.flag_Qmax = True
			else:
				ibdg.flag_Qlim = False
				ibdg.flag_Qmin = False
				ibdg.flag_Qmax = False
				if Q3_k1 > 0:
					if ibdg.V1 <= Vset <= ibdg.V2:
						pass
					else:
						V3 = ibdg.V3*ibdg.Vnom
						V4 = ibdg.V4*ibdg.Vnom
						V_k = Vset * ibdg.Vnom
						v3_4 = [V3, V4]
						Q3_4 = [0, ibdg.Qmax]
						m = (Q3_4[1] - Q3_4[0]) / (v3_4[1] - v3_4[0])
						b = (Q3_4[0] - (m * ibdg.V3))
						Vsol[ibdg.node_Q3] = ((m * V_k) + b)
				else:
					if ibdg.V3 <= Vset <= ibdg.V4:
						pass
					else:
						V1 = ibdg.V1*ibdg.Vnom
						V2 = ibdg.V2*ibdg.Vnom
						V_k = Vset * ibdg.Vnom
						v1_2 = [V1, V2]
						Q1_2 = [ibdg.Qmin, 0]
						m = (Q1_2[1] - Q1_2[0]) / (v1_2[1] - v1_2[0])
						b = (Q1_2[0] - (m * V1))
						Vsol[ibdg.node_Q3] = ((m * V_k) + b)


		return Vsol
