#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

"""
	stamp_linear: stamps all the linear three-phase powerflow models.
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


from classes.CurrentMeasurements import CurrentMeasurements, find_ammeter


def stamp_linear_test(slack, overheadline, undergroundline, triplexline, nodes, xfmr,
                   regulator, switch, fuse, capacitor, reactor, curr_meas, features, 
                   stamped_ground, manual_initialize, idx_Y_init, stamp_dual, node, V):

	n_ohline = len(overheadline)
	n_ugline = len(undergroundline)
	n_tpline = len(triplexline)
	n_xfmr = len(xfmr)
	n_reg = len(regulator)
	n_switch = len(switch)
	n_fuse = len(fuse)
	n_slack = len(slack)
	n_reactor = len(reactor)
	n_cap = len(capacitor)
	n_currmeas = len(curr_meas)

	nnz_init = 2 * idx_Y_init
	nnz = nnz_init

	Ylin_row = np.zeros(nnz, dtype=np.int64)
	Ylin_col = np.zeros(nnz, dtype=np.int64)
	Ylin_val = np.zeros(nnz, dtype=np.float64)
	Jlin_row = np.zeros(nnz, dtype=np.int64)
	Jlin_val = np.zeros(nnz, dtype=np.float64)

	idx_Y = 0
	idx_J = 0

	# Features
	enable_CM = features['Current Meas']
	enable_Cap = features['Capacitors']

	# Stamp Lines
	# Overhead Lines
	for ele in range(n_ohline):
		ammeter = find_ammeter(overheadline[ele], curr_meas, enable_CM)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = overheadline[ele].stamp_linear(nodes, Ylin_val, Ylin_row, Ylin_col,
		                                                                     idx_Y, enable_CM, ammeter, stamp_dual)

	# Underground Lines
	for ele in range(n_ugline):
		ammeter = find_ammeter(undergroundline[ele], curr_meas, enable_CM)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = undergroundline[ele].stamp_linear(nodes, Ylin_val, Ylin_row, Ylin_col,
		                                                                        idx_Y, enable_CM, ammeter, stamp_dual)

	# Triplex Lines
	for ele in range(n_tpline):
		ammeter = find_ammeter(triplexline[ele], curr_meas, enable_CM)
		Ylin_val, Ylin_row, Ylin_col, idx_Y = triplexline[ele].stamp_linear_triplex(nodes, Ylin_val, Ylin_row,
		                                                                            Ylin_col, idx_Y, enable_CM, ammeter, stamp_dual)

	# Stamp Transformers
	for ele in range(n_xfmr):
		ammeter = find_ammeter(xfmr[ele], curr_meas, enable_CM)
		Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground = xfmr[ele].stamp_linear(Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground, enable_CM,
		                                                             ammeter, stamp_dual = stamp_dual)

	# Stamp Regulators
	for ele in range(n_reg):
		Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground = regulator[ele].stamp_linear(node, V, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                                                    stamped_ground, manual_initialize)

	# Stamp Switch
	for ele in range(n_switch):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = switch[ele].stamp_linear(nodes, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp Fuse
	for ele in range(n_fuse):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = fuse[ele].stamp_linear(nodes, Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp Slack
	for ele in range(n_slack):
		Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J = slack[ele].stamp_linear(nodes, Ylin_val,
		                                                                                         Ylin_row, Ylin_col,
		                                                                                         Jlin_val, Jlin_row,
		                                                                                         idx_Y, idx_J, stamp_dual = stamp_dual)
	# stamp reactors
	for ele in range(n_reactor):
		Ylin_val, Ylin_row, Ylin_col, idx_Y = reactor[ele].stamp_linear(Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp Current Measurers
	if enable_CM:
		for ele in range(n_currmeas):
			Ylin_val, Ylin_row, Ylin_col, idx_Y = curr_meas[ele].stamp_linear(Ylin_val, Ylin_row, Ylin_col, idx_Y)

	# Stamp capacitors
	if enable_Cap:  # The purpose of the feature if turn on or off a given component
		for ele in range(n_cap):
			Ylin_val, Ylin_row, Ylin_col, idx_Y = capacitor[ele].stamp_linear(nodes, Ylin_val, Ylin_row, Ylin_col,
			                                                                  idx_Y, stamp_dual = stamp_dual)

	return Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J, stamped_ground
