#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

""" Nonlinear device stamps.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 06-05-2020
  Updated Date: 10-18-2020
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""
from classes.Regulators import VariableRegulator
import numpy as np
cimport numpy as np

ctypedef  np.float64_t DTYPE_t
DTYPE = np.int64
DTYPE_D = np.float64


cpdef stamp_nonlinear(node_key, nodes, regulator, load, triplexload, ibdg, voltage_bounds, current_meas, features, V, int idx_Y_init,
					  double lf, homotopy_enabled, homotopy, stamp_dual):
	cdef int n_load = len(load)
	cdef int n_tpload = len(triplexload)
	cdef int n_ibdg = len(ibdg)
	cdef int n_currmeas = len(current_meas)
	cdef int n_regulator = len(regulator)

	nnz_init = 2 * idx_Y_init
	cdef int nnz = nnz_init

	cdef np.ndarray[np.int64_t, ndim=1] Ynlin_row, Ynlin_col, Jnlin_row
	cdef np.ndarray[double, ndim=1] Ynlin_val, Jnlin_val

	Ynlin_row = np.zeros(nnz, dtype=DTYPE)
	Ynlin_col = np.zeros(nnz, dtype=DTYPE)
	Ynlin_val = np.zeros(nnz, dtype=DTYPE_D)
	Jnlin_row = np.zeros(nnz, dtype=DTYPE)
	Jnlin_val = np.zeros(nnz, dtype=DTYPE_D)

	cdef int idx_Y = 0
	cdef int idx_J = 0
	
	enable_CM = features['Current Meas']

	for ele in range(n_regulator):
		if isinstance(regulator[ele], VariableRegulator):
			(idx_Y, idx_J) = regulator[ele].stamp_nonlinear(node_key, nodes, V, Ynlin_val,
															Ynlin_row, Ynlin_col, Jnlin_val,
															Jnlin_row, idx_Y, idx_J,
															homotopy_enabled, homotopy.h_factor,
															homotopy.g_value, homotopy.b_value)

	# Stamp Loads
	for ele in range(n_load):
		(idx_Y, idx_J) = load[ele].stamp_nonlinear(node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
												   Jnlin_val, Jnlin_row, idx_Y, idx_J, enable_CM)

	for ele in range(n_tpload):
		(idx_Y, idx_J) = triplexload[ele].stamp_nonlinear(node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
														  Jnlin_val, Jnlin_row, idx_Y, idx_J)

	if ibdg:
		for ele in range(n_ibdg):
			(idx_Y, idx_J) = ibdg[ele].stamp_nonlinear(node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
													 idx_Y, idx_J)

	if voltage_bounds:
		for ele in range(len(voltage_bounds)):
			(idx_Y, idx_J) = voltage_bounds[ele].stamp_nonlinear(node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
													 idx_Y, idx_J)
				
	if enable_CM:
		for ele in range(n_currmeas):
			(idx_Y, idx_J) = current_meas[ele].stamp_nonlinear(V, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, 
																Jnlin_row, idx_Y, idx_J)
	
	return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J
