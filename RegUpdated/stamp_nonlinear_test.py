#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

"""
	stamp_nonlinear: stamps all the nonlinear three-phase powerflow models.
	__author__ = ["Amritanshu Pandey"]
	__copyright__ = "Copyright 2017-2020"
	__credits__ = ["Amritanshu Pandey"]
	__license__ = "GPL"
	__version__ = "1.0.1"
	__maintainer__ = ["Amritanshu Pandey", "Naeem Turner-Bandele"]
	__email__ = "nturnerb@andrew.cmu.edu"
	__status__ = "Development"
"""

import numpy as np


def stamp_nonlinear_test(nodes, load, triplexload, dg, current_meas, enable_CM, V, idx_Y_init, lf, \
                      pstep__,
                      pstepV__, stamp_dual):

	n_load = len(load)
	n_tpload = len(triplexload)
	n_dg = len(dg)
	n_currmeas = len(current_meas)

	nnz_init = 2 * idx_Y_init
	nnz = nnz_init

	Ynlin_row = np.zeros(nnz, dtype=np.int64)
	Ynlin_col = np.zeros(nnz, dtype=np.int64)
	Ynlin_val = np.zeros(nnz, dtype=np.float64)
	Jnlin_row = np.zeros(nnz, dtype=np.int64)
	Jnlin_val = np.zeros(nnz, dtype=np.float64)

	idx_Y = 0
	idx_J = 0

	# Stamp Loads
	for ele in range(len(load)):
		(idx_Y, idx_J) = load[ele].stamp_nonlinear(nodes, V, lf, pstep__, Ynlin_val, Ynlin_row, Ynlin_col,
												   Jnlin_val, Jnlin_row, idx_Y, idx_J, stamp_dual, enable_CM)

	for ele in range(len(triplexload)):
		(idx_Y, idx_J) = triplexload[ele].stamp_nonlinear(nodes, V, lf, pstep__, Ynlin_val, Ynlin_row, Ynlin_col,
														  Jnlin_val, Jnlin_row, idx_Y, idx_J)

	if dg:
		for ele in range(len(dg)):
			(idx_Y, idx_J) = dg[ele].stamp_nonlinear(nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
													 idx_Y, idx_J)

	if enable_CM:
		for ele in range(n_currmeas):
			Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J = current_meas[ele].stamp_nonlinear(
				V, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J)

	return Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J
