from classes.OperationalLimits.OperationalLimits_base import OperationalLimits

def stampY(i, j, val, Y_val, Y_row, Y_col, idx):
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return idx


def stampJ(i, val, J_val, J_row, idx):
	J_val[idx] = val
	J_row[idx] = i
	idx += 1

	return idx 


class VoltageLimits(OperationalLimits):
	pu_max2 = []
	pu_min2 = []

	def __init__(self, node, voltage_bound_settings):
		pu_max = voltage_bound_settings["Vmag Bound Max pu"]
		pu_min = voltage_bound_settings["Vmag Bound Min pu"]
		self.cs_eps = voltage_bound_settings['cs eps tol']
		OperationalLimits.__init__(self, pu_max, pu_min)
		self.v2max_pu = pu_max**2
		self.v2min_pu = pu_min**2
		self.connected_node_name = node.name
		self.node = node
		self.node_name = node.name
		self.vnom = node.Vnom
		self.vnom2 = self.vnom**2


	def init_bounding_variables(self, Vinit):
		if self.node.bustype == 3:
			return
		VoltageLimits.pu_max2.extend(([self.v2max_pu]*len(self.node.vmag2_inds)))
		VoltageLimits.pu_min2.extend([self.v2min_pu]*len(self.node.vmag2_inds))
		Vr = Vinit[self.node.vr_bounding_inds]
		Vi = Vinit[self.node.vi_bounding_inds]
		Vinit[self.node.vmag2_inds] = (Vr**2 + Vi**2)/self.vnom2
		Vinit[self.node.Lvmag2_inds] = 1
		Vinit[self.node.umax_vmag2_inds] = 1e-4
		Vinit[self.node.umin_vmag2_inds] = 1e-4
		# if not self.node.isTriplex:

		# 	nodeA_vmag2_norm = (Vinit[self.node.nodeA_Vr]**2 + Vinit[self.node.nodeA_Vi]**2)/self.vnom2
		# 	Vinit[self.node.nodeA_vmag2_index] = nodeA_vmag2_norm
		# 	nodeB_vmag2_norm = (Vinit[self.node.nodeB_Vr]**2 + Vinit[self.node.nodeB_Vi]**2)/self.vnom2
		# 	Vinit[self.node.nodeB_vmag2_index] = nodeB_vmag2_norm
		# 	nodeC_vmag2_norm = (Vinit[self.node.nodeC_Vr]**2 + Vinit[self.node.nodeC_Vi]**2)/self.vnom2
		# 	Vinit[self.node.nodeC_vmag2_index] = nodeC_vmag2_norm

		# 	Vinit[self.node.nodeA_Lvmag2_index] = 1
		# 	Vinit[self.node.nodeB_Lvmag2_index] = 1
		# 	Vinit[self.node.nodeC_Lvmag2_index] = 1

		# 	Vinit[self.node.nodeA_umax_vmag2_index] = 1e-4
		# 	Vinit[self.node.nodeB_umax_vmag2_index] = 1e-4
		# 	Vinit[self.node.nodeC_umax_vmag2_index] = 1e-4
		# 	Vinit[self.node.nodeA_umin_vmag2_index] = 1e-4
		# 	Vinit[self.node.nodeB_umin_vmag2_index] = 1e-4
		# 	Vinit[self.node.nodeC_umin_vmag2_index] = 1e-4
		# else:
		# 	VoltageLimits.pu_max2.extend([self.v2max_pu, self.v2max_pu])
		# 	VoltageLimits.pu_min2.extend([self.v2min_pu, self.v2min_pu])

		# 	node1_vmag2_norm = (Vinit[self.node.node1_Vr]**2 + Vinit[self.node.node1_Vi]**2)/self.vnom2
		# 	node2_vmag2_norm = (Vinit[self.node.node2_Vr]**2 + Vinit[self.node.node2_Vi]**2)/self.vnom2
		# 	Vinit[self.node.node1_vmag2_index] = node1_vmag2_norm
		# 	Vinit[self.node.node2_vmag2_index] = node2_vmag2_norm

		# 	Vinit[self.node.node1_Lvmag2_index] = 1
		# 	Vinit[self.node.node2_Lvmag2_index] = 1

		# 	Vinit[self.node.node1_umax_vmag2_index] = 1e-4
		# 	Vinit[self.node.node2_umax_vmag2_index] = 1e-4
		# 	Vinit[self.node.node1_umin_vmag2_index] = 1e-4
		# 	Vinit[self.node.node2_umin_vmag2_index] = 1e-4


	def stamp_nonlinear(self, node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
													 idx_Y, idx_J):
		if self.node.bustype == 3:
			# try not limiting voltage on slack bus for now?
			return (idx_Y, idx_J)
		# if len(Nodes.bounded_bus_list) > 0 and self.name not in Nodes.bounded_bus_list:
		# 	return
		# if self.node.isTriplex:
		# 	Vr_inds = [self.node.node1_Vr, self.node.node2_Vr]
		# 	Vi_inds = [self.node.node1_Vi, self.node.node2_Vi]
		# 	Lr_inds = [self.node.node1_dual_eq_var_r, self.node.node2_dual_eq_var_r]
		# 	Li_inds = [self.node.node1_dual_eq_var_i, self.node.node2_dual_eq_var_i]
		# 	vmag2_inds = [self.node.node1_vmag2_index, self.node.node2_vmag2_index]
		# 	Lvmag2_inds = [self.node.node1_Lvmag2_index, self.node.node2_Lvmag2_index]
		# 	umax_vmag2_inds = [self.node.node1_umax_vmag2_index, self.node.node2_umax_vmag2_index]
		# 	umin_vmag2_inds = [self.node.node1_umin_vmag2_index, self.node.node2_umin_vmag2_index]
		# else:
		# 	Vr_inds = [self.node.nodeA_Vr, self.node.nodeB_Vr, self.node.nodeC_Vr]
		# 	Vi_inds = [self.node.nodeA_Vi, self.node.nodeB_Vi, self.node.nodeC_Vi]
		# 	Lr_inds = [self.node.nodeA_dual_eq_var_r, self.node.nodeB_dual_eq_var_r, self.node.nodeC_dual_eq_var_r]
		# 	Li_inds = [self.node.nodeA_dual_eq_var_i, self.node.nodeB_dual_eq_var_i, self.node.nodeC_dual_eq_var_i]
		# 	vmag2_inds = [self.node.nodeA_vmag2_index, self.node.nodeB_vmag2_index, self.node.nodeC_vmag2_index]
		# 	Lvmag2_inds = [self.node.nodeA_Lvmag2_index, self.node.nodeB_Lvmag2_index, self.node.nodeC_Lvmag2_index]
		# 	umax_vmag2_inds = [self.node.nodeA_umax_vmag2_index, self.node.nodeB_umax_vmag2_index, self.node.nodeC_umax_vmag2_index]
		# 	umin_vmag2_inds = [self.node.nodeA_umin_vmag2_index, self.node.nodeB_umin_vmag2_index, self.node.nodeC_umin_vmag2_index]

		Vr_inds = self.node.vr_bounding_inds
		Vi_inds = self.node.vi_bounding_inds
		Lr_inds = self.node.Lr_bounding_inds
		Li_inds = self.node.Li_bounding_inds
		vmag2_inds = self.node.vmag2_inds
		Lvmag2_inds = self.node.Lvmag2_inds
		umax_vmag2_inds = self.node.umax_vmag2_inds
		umin_vmag2_inds = self.node.umin_vmag2_inds

		nPhases = len(Vr_inds)
		for i in range(nPhases):
			# if Vr_inds[i] not in unique_nonzero_index:
			# 	continue

			Vr = V[Vr_inds[i]]
			Vi = V[Vi_inds[i]]
			LR = V[Lr_inds[i]]
			LI = V[Li_inds[i]]
			vmag2 = V[vmag2_inds[i]]
			Lvmag2 = V[Lvmag2_inds[i]]
			umax_vmag2 = V[umax_vmag2_inds[i]]
			umin_vmag2 = V[umin_vmag2_inds[i]]

			#
			vnom2 = self.vnom2

			Fv2_hist = vmag2 - (Vr*Vr)/vnom2 - (Vi*Vi)/vnom2
			dFv2_dvmag2 = 1
			dFv2_dVr = -2/vnom2*Vr
			dFv2_dVi = -2/vnom2*Vi
			# idx_Y = stampY(i, j, val, Y_val, Y_row, Y_col, idx)
			idx_Y = stampY(vmag2_inds[i], vmag2_inds[i], dFv2_dvmag2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			
			#
			idx_Y = stampY(vmag2_inds[i], Vr_inds[i],
									dFv2_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

			#
			idx_Y = stampY(vmag2_inds[i], Vi_inds[i],
									dFv2_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		
			#
			Fv2_jval = -(Fv2_hist - dFv2_dvmag2*vmag2 - dFv2_dVr*Vr - dFv2_dVi*Vi)
			idx_J = stampJ(vmag2_inds[i], Fv2_jval, 
									Jnlin_val, Jnlin_row, idx_J)
		
			#### stamping dLagrange/dVr = -2/vnom2*Vr*Lvmag2
			Fvr_hist = -2/vnom2*Vr*Lvmag2
			dFvrdLv2 = -2/vnom2*Vr
			dFvrdVr = -2/vnom2*Lvmag2
			#
			idx_Y = stampY(Lr_inds[i], Vr_inds[i],
								dFvrdVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_Y = stampY(Lr_inds[i], Lvmag2_inds[i],
								dFvrdLv2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			Fvr_jval = -(Fvr_hist - dFvrdLv2*Lvmag2-dFvrdVr*Vr)
			idx_J = stampJ(Lr_inds[i], Fvr_jval, Jnlin_val, Jnlin_row, idx_J)
			
			#### stamping dLagrange/dVi = -2/vnom2*Vi*Lvmag2
			Fvi_hist = -2/vnom2*Vi*Lvmag2
			dFvidLv2 = -2/vnom2*Vi
			dFvidVi = -2/vnom2*Lvmag2
			#
			idx_Y = stampY(Li_inds[i], Vi_inds[i],
								dFvidVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_Y = stampY(Li_inds[i], Lvmag2_inds[i],
								dFvidLv2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_J = stampJ(Li_inds[i], -(Fvi_hist - dFvidLv2*Lvmag2-dFvidVi*Vi), 
									Jnlin_val, Jnlin_row, idx_J)
		
			# stamping dLagrange/dvmag2
			#
			idx_Y = stampY(Lvmag2_inds[i], Lvmag2_inds[i],
								1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_Y = stampY(Lvmag2_inds[i], umax_vmag2_inds[i],
								1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
		
			idx_Y = stampY(Lvmag2_inds[i], umin_vmag2_inds[i],
								-1, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

			# if abs(vmag2 - self.v2max_pu) < 1e-6:
			# 	print("%s, node %d: UB delta %.3e" % (self.node.name, i, (self.v2max_pu-vmag2)))
			umax_hist = umax_vmag2*(vmag2 - self.v2max_pu) + self.cs_eps
			dumax_dumax = vmag2 - self.v2max_pu
			dumax_dvmag2 = umax_vmag2
		
			#
			idx_Y = stampY(umax_vmag2_inds[i], umax_vmag2_inds[i],
								dumax_dumax, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_Y = stampY(umax_vmag2_inds[i], vmag2_inds[i],
								dumax_dvmag2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_J = stampJ(umax_vmag2_inds[i], -(umax_hist - dumax_dumax*umax_vmag2 - dumax_dvmag2*vmag2), 
									Jnlin_val, Jnlin_row, idx_J)
		
			### stamping lower bound
			# if abs(vmag2 - self.v2min_pu) < 1e-6:
				# print("%s, node %d: LB delta %.3e" % (self.node.name, i, (vmag2 - self.v2min_pu)))
			umin_hist = -umin_vmag2*(vmag2 - self.v2min_pu) + self.cs_eps
			dumin_dumin = -(vmag2 - self.v2min_pu)
			dumin_dvmag2 = -umin_vmag2
		
			#
			idx_Y = stampY(umin_vmag2_inds[i], umin_vmag2_inds[i],
									dumin_dumin, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_Y = stampY(umin_vmag2_inds[i], vmag2_inds[i],
									dumin_dvmag2, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
			#
			idx_J = stampJ(umin_vmag2_inds[i], -(umin_hist - dumin_dumin*umin_vmag2 - dumin_dvmag2*vmag2), 
									Jnlin_val, Jnlin_row, idx_J)
			
		return (idx_Y, idx_J)
	
	def calc_residuals(self, V, res_eqn):
		if self.node.bustype == 3:
			return
		# if len(Nodes.bounded_bus_list) > 0 and self.name not in Nodes.bounded_bus_list:
		# 	return
		# if self.node.isTriplex:
		# 	Vr_inds = [self.node.node1_Vr, self.node.node2_Vr]
		# 	Vi_inds = [self.node.node1_Vi, self.node.node2_Vi]
		# 	Lr_inds = [self.node.node1_dual_eq_var_r, self.node.node2_dual_eq_var_r]
		# 	Li_inds = [self.node.node1_dual_eq_var_i, self.node.node2_dual_eq_var_i]
		# 	vmag2_inds = [self.node.node1_vmag2_index, self.node.node2_vmag2_index]
		# 	Lvmag2_inds = [self.node.node1_Lvmag2_index, self.node.node2_Lvmag2_index]
		# 	umax_vmag2_inds = [self.node.node1_umax_vmag2_index, self.node.node2_umax_vmag2_index]
		# 	umin_vmag2_inds = [self.node.node1_umin_vmag2_index, self.node.node2_umin_vmag2_index]
		# else:
		# 	Vr_inds = [self.node.nodeA_Vr, self.node.nodeB_Vr, self.node.nodeC_Vr]
		# 	Vi_inds = [self.node.nodeA_Vi, self.node.nodeB_Vi, self.node.nodeC_Vi]
		# 	Lr_inds = [self.node.nodeA_dual_eq_var_r, self.node.nodeB_dual_eq_var_r, self.node.nodeC_dual_eq_var_r]
		# 	Li_inds = [self.node.nodeA_dual_eq_var_i, self.node.nodeB_dual_eq_var_i, self.node.nodeC_dual_eq_var_i]
		# 	vmag2_inds = [self.node.nodeA_vmag2_index, self.node.nodeB_vmag2_index, self.node.nodeC_vmag2_index]
		# 	Lvmag2_inds = [self.node.nodeA_Lvmag2_index, self.node.nodeB_Lvmag2_index, self.node.nodeC_Lvmag2_index]
		# 	umax_vmag2_inds = [self.node.nodeA_umax_vmag2_index, self.node.nodeB_umax_vmag2_index, self.node.nodeC_umax_vmag2_index]
		# 	umin_vmag2_inds = [self.node.nodeA_umin_vmag2_index, self.node.nodeB_umin_vmag2_index, self.node.nodeC_umin_vmag2_index]
		
		Vr_inds = self.node.vr_bounding_inds
		Vi_inds = self.node.vi_bounding_inds
		Lr_inds = self.node.Lr_bounding_inds
		Li_inds = self.node.Li_bounding_inds
		vmag2_inds = self.node.vmag2_inds
		Lvmag2_inds = self.node.Lvmag2_inds
		umax_vmag2_inds = self.node.umax_vmag2_inds
		umin_vmag2_inds = self.node.umin_vmag2_inds

		nPhases = len(Vr_inds)

		for i in range(nPhases):			
			Vr = V[Vr_inds[i]]
			Vi = V[Vi_inds[i]]
			LR = V[Lr_inds[i]]
			LI = V[Li_inds[i]]
			vmag2 = V[vmag2_inds[i]]
			Lvmag2 = V[Lvmag2_inds[i]]
			umax_vmag2 = V[umax_vmag2_inds[i]]
			umin_vmag2 = V[umin_vmag2_inds[i]]
			vnom2 = self.vnom2#[i]

			res_eqn[vmag2_inds[i]] += vmag2 - Vr*Vr/vnom2 - Vi*Vi/vnom2			
			res_eqn[umax_vmag2_inds[i]] += umax_vmag2*(vmag2 - self.v2max_pu) + self.cs_eps
			res_eqn[umin_vmag2_inds[i]] += -umin_vmag2*(vmag2 - self.v2min_pu) + self.cs_eps
			res_eqn[Lr_inds[i]] += -2/vnom2*Vr*Lvmag2
			res_eqn[Li_inds[i]] += -2/vnom2*Vi*Lvmag2
			res_eqn[Lvmag2_inds[i]] += Lvmag2 + umax_vmag2 - umin_vmag2
