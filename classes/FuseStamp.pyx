"""Implements fuse stamps.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 10-18-2020
  Updated Date: 10-18-2020
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""
#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3


cpdef stampY(i, j, val, Y_val, Y_row, Y_col, idx):
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return Y_val, Y_row, Y_col, idx

cpdef stamp_short_circuit(int node_fuse_Vr,
						  int node_Vr_from,
						  int node_Vr_to,
						  int node_fuse_Vi,
						  int node_Vi_from,
						  int node_Vi_to,
						  Ylin_val,
						  Ylin_row,
						  Ylin_col,
						  idx_Y):

	# Real Voltage Source Equation
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Vr, node_Vr_from, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Vr, node_Vr_to, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	# Imaginary Voltage Source Equation
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Vi, node_Vi_from, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Vi, node_Vi_to, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	# Currents for Real and Imag Nodes
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Vr_from, node_fuse_Vr, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Vr_to, node_fuse_Vr, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Vi_from, node_fuse_Vi, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Vi_to, node_fuse_Vi, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y

cpdef stamp_short_circuit_dual(int node_fuse_Lr, 
							int node_Lr_from, 
							int node_Lr_to, 
							int node_fuse_Li,
							int node_Li_from, 
							int node_Li_to, 
							Ylin_val,
							Ylin_row,
							Ylin_col, 
							idx_Y):
	# Real Voltage Source Equation
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Lr_from, node_fuse_Lr, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Lr_to, node_fuse_Lr, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	# Imaginary Voltage Source Equation
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Li_from, node_fuse_Li, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_Li_to, node_fuse_Li, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	# Currents for Real and Imag Nodes
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Lr, node_Lr_from, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Lr, node_Lr_to, -1, Ylin_val, Ylin_row,
						   Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Li, node_Li_from, 1, Ylin_val,
						   Ylin_row, Ylin_col, idx_Y)
	(Ylin_val, Ylin_row, Ylin_col,
	 idx_Y) = stampY(node_fuse_Li, node_Li_to,  -1, Ylin_val, Ylin_row,
							   Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, idx_Y
	
	
cpdef stamp_linear(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, idx_Y):
	# Find the from bus nodes to stamp
	nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr
	nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi
	nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr
	nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi
	nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
	nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi
	# Find the to bus nodes to stamp
	nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
	nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi
	nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
	nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi
	nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
	nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi
	# Stamp Phase A if exist
	if self.phases & 0x1 == 1 and self.phase_A_state == self._GOOD:
		(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit(self.nodeA_Ir, nodeA_Vr_from,
																	nodeA_Vr_to, self.nodeA_Ii,
																	nodeA_Vi_from, nodeA_Vi_to, Ylin_val,
																	Ylin_row, Ylin_col, idx_Y)
	# Stamp Phase B is exist
	if self.phases & 0x2 == 2 and self.phase_B_state == self._GOOD:
		(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit(self.nodeB_Ir, nodeB_Vr_from,
																	nodeB_Vr_to, self.nodeB_Ii,
																	nodeB_Vi_from, nodeB_Vi_to, Ylin_val,
																	Ylin_row, Ylin_col, idx_Y)
	# Stamp Phase C is exist
	if self.phases & 0x4 == 4 and self.phase_C_state == self._GOOD:
		(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit(self.nodeC_Ir, nodeC_Vr_from,
																	nodeC_Vr_to, self.nodeC_Ii,
																	nodeC_Vi_from, nodeC_Vi_to, Ylin_val,
																	Ylin_row, Ylin_col, idx_Y)
	if self.stamp_dual:
			nodeA_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
			nodeA_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i
			nodeB_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
			nodeB_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
			nodeC_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
			nodeC_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i
			# Find the to bus nodes to stamp
			nodeA_Lr_to = node[node_key[self.to_node]].nodeA_dual_eq_var_r
			nodeA_Li_to = node[node_key[self.to_node]].nodeA_dual_eq_var_i
			nodeB_Lr_to = node[node_key[self.to_node]].nodeB_dual_eq_var_r
			nodeB_Li_to = node[node_key[self.to_node]].nodeB_dual_eq_var_i
			nodeC_Lr_to = node[node_key[self.to_node]].nodeC_dual_eq_var_r
			nodeC_Li_to = node[node_key[self.to_node]].nodeC_dual_eq_var_i
			# Stamp Phase A if exist
			if self.phases & 0x1 == 1 and self.phase_A_state == self._GOOD:
				(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit_dual(
					self.nodeA_Lr_fuse, nodeA_Lr_from, nodeA_Lr_to,
					self.nodeA_Li_fuse, nodeA_Li_from, nodeA_Li_to, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
			# Stamp Phase B is exist
			if self.phases & 0x2 == 2 and self.phase_B_state == self._GOOD:
				(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit_dual(
					self.nodeB_Lr_fuse, nodeB_Lr_from, nodeB_Lr_to,
					self.nodeB_Li_fuse, nodeB_Li_from, nodeB_Li_to, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
			# Stamp Phase C is exist
			if self.phases & 0x4 == 4 and self.phase_C_state == self._GOOD:
				(Ylin_val, Ylin_row, Ylin_col, idx_Y) = stamp_short_circuit_dual(
					self.nodeC_Lr_fuse, nodeC_Lr_from, nodeC_Lr_to,
					self.nodeC_Li_fuse, nodeC_Li_from, nodeC_Li_to, Ylin_val,
					Ylin_row, Ylin_col, idx_Y)
					
	return Ylin_val, Ylin_row, Ylin_col, idx_Y
