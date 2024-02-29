# cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

"""
	stamp_linear_homotopy: homotopy stamp for the linear powerflow models.
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
from classes.Nodes import Nodes


def stamp_linear_homotopy_test(overheadline, undergroundline, triplexline, xfmr, regulator, features, node, V,
                               curr_meas, stamped_ground, idx_Y_init, h_factor, G_homotopy, B_homotopy, manual_init):
	nnz_init = 2 * idx_Y_init
	nnz = nnz_init

	Ylin_h_row = np.zeros(nnz, dtype=np.int64)
	Ylin_h_col = np.zeros(nnz, dtype=np.int64)
	Ylin_h_val = np.zeros(nnz, dtype=np.float64)
	Jlin_h_row = np.zeros(nnz, dtype=np.int64)
	Jlin_h_val = np.zeros(nnz, dtype=np.float64)

	idx_Y = 0
	idx_J = 0

	# Features
	enable_CM = features['Current Meas']

	# Stamp Lines
	# Overhead Lines
	for ele in range(len(overheadline)):
		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = overheadline[ele].stamp_linear_homotopy(node, Ylin_h_val,
		                                                                                    Ylin_h_row, Ylin_h_col,
		                                                                                    idx_Y, h_factor,
		                                                                                    G_homotopy, B_homotopy
		                                                                                    )

	# Underground Lines
	for ele in range(len(undergroundline)):
		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = undergroundline[ele].stamp_linear_homotopy(node, Ylin_h_val,
		                                                                                       Ylin_h_row, Ylin_h_col,
		                                                                                       idx_Y, h_factor,
		                                                                                       G_homotopy,
		                                                                                       B_homotopy)

	# Triplex Lines
	for ele in range(len(triplexline)):
		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = triplexline[ele].stamp_linear_triplex_homotopy(node, Ylin_h_val,
		                                                                                           Ylin_h_row,
		                                                                                           Ylin_h_col, idx_Y,
		                                                                                           h_factor,
		                                                                                           G_homotopy,
		                                                                                           B_homotopy)

	# Stamp Transformers
	for ele in range(len(xfmr)):
		homotopy_option = True
		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y, stamped_ground = xfmr[ele].stamp_linear(Ylin_h_val, Ylin_h_row,
		                                                                                   Ylin_h_col,
		                                                                                   idx_Y, stamped_ground, 
		                                                                                   homotopy_option, h_factor,
		                                                                                   G_homotopy, B_homotopy)


# 	for ele in range(len(node)):
#           if node[ele].bustype != 3:
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeA_Vr, node[ele].nodeA_Lr, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeB_Vr, node[ele].nodeB_Lr, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeC_Vr, node[ele].nodeC_Lr, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeN_Vr, node[ele].nodeN_Lr, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeA_Vi, node[ele].nodeA_Li, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeB_Vi, node[ele].nodeB_Li, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeC_Vi, node[ele].nodeC_Li, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)
#               Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = Nodes.stampY(node[ele].nodeN_Vi, node[ele].nodeN_Li, - 1/2*h_factor**2, Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y)


	return Ylin_h_val, Ylin_h_row, Ylin_h_col, Jlin_h_row, Jlin_h_val, idx_Y, idx_J, stamped_ground
