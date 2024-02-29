#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

""" Linear device stamps.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 06-05-2020
  Updated Date: 10-18-2020
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""

from classes.Regulators import FixedTapRegulator

from classes.Regulators import FixedTapRegulator, VariableRegulator

import numpy as np

cimport numpy as np



cpdef stamp_linear(node_key, slack, overheadline, undergroundline, triplexline, nodes, xfmr,
				   regulator, switch, fuse, capacitor, reactor, curr_meas, features, stamped_ground, 
				   int idx_Y_init, V, stamp_dual):

	cdef int n_ohline = len(overheadline)
	cdef int n_ugline = len(undergroundline)
	cdef int n_tpline = len(triplexline)
	cdef int n_xfmr = len(xfmr)
	cdef int n_reg = len(regulator)
	cdef int n_switch = len(switch)
	cdef int n_fuse = len(fuse)
	cdef int n_slack = len(slack)
	cdef int n_reactor = len(reactor)
	cdef int n_cap = len(capacitor)
	cdef int n_currmeas = len(curr_meas)

	nnz_init = 2 * idx_Y_init
	cdef int nnz = nnz_init

	cdef np.ndarray[np.int64_t, ndim=1] Ylin_row, Ylin_col, Jlin_row
	cdef np.ndarray[double, ndim=1] Ylin_val, Jlin_val

	Ylin_row = np.zeros(nnz, dtype=np.int64)
	Ylin_col = np.zeros(nnz, dtype=np.int64)
	Ylin_val = np.zeros(nnz, dtype=np.float64)
	Jlin_row = np.zeros(nnz, dtype=np.int64)
	Jlin_val = np.zeros(nnz, dtype=np.float64)

	cdef int idx_Y = 0
	cdef int idx_J = 0

	# Features
	enable_CM = features['Current Meas']
	fixed_taps = features['Tap Controls']['Fixed'] if 'Tap Controls' in features else True

	# Stamp Lines
	# Overhead Lines
	for ele in range(n_ohline):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = overheadline[ele].stamp_linear(node_key, nodes, Ylin_val, Ylin_row, Ylin_col,
																			 idx_Y)

	# Underground Lines
	for ele in range(n_ugline):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = undergroundline[ele].stamp_linear(node_key, nodes, Ylin_val, Ylin_row, Ylin_col,
																				idx_Y)

	# Triplex Lines
	for ele in range(n_tpline):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = triplexline[ele].stamp_linear_triplex(node_key, nodes, Ylin_val, Ylin_row,
																					Ylin_col, idx_Y)

	# Stamp Transformers
	for ele in range(n_xfmr):
		Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground = xfmr[ele].stamp_linear(Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground)
				

	# Stamp Regulators
	for ele in range(n_reg):
		if isinstance(regulator[ele], FixedTapRegulator):
			(Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J, stamped_ground) = regulator[ele].stamp_linear(Ylin_val, Ylin_row, 
																								Ylin_col, Jlin_row, Jlin_val, 
																								idx_Y, idx_J, stamped_ground)
		else:
			(idx_Y, stamped_ground) = regulator[ele].stamp_linear_ground(Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground)

	# Stamp Switch
	for ele in range(n_switch):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = switch[ele].stamp_linear(node_key, nodes, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp Fuse
	for ele in range(n_fuse):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = fuse[ele].stamp_linear(node_key, nodes, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp Slack
	for ele in range(n_slack):
		Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J = slack[ele].stamp_linear(node_key, nodes, Ylin_val,
																								 Ylin_row, Ylin_col,
																								 Jlin_val, Jlin_row,
																								 idx_Y, idx_J)
	# stamp reactors
	for ele in range(n_reactor):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = reactor[ele].stamp_linear(Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp capacitors
	for ele in range(n_cap):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = capacitor[ele].stamp_linear(node_key, nodes, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	return Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J, stamped_ground
