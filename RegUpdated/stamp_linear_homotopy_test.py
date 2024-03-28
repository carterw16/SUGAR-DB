#cython: boundscheck=False, wraparound=False, nonecheck=False, initializedcheck=False, language_level=3

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


def stamp_linear_homotopy_test(overheadline, undergroundline, triplexline, xfmr, regulator, features, node, V,
                            curr_meas, stamped_ground, idx_Y_init, h_factor, G_homotopy, B_homotopy, manual_init, stamp_dual):
    

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
   		ammeter = find_ammeter(overheadline[ele], curr_meas, enable_CM)
   		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = overheadline[ele].stamp_linear_homotopy(node, Ylin_h_val,
   		                                                                                    Ylin_h_row, Ylin_h_col,
   		                                                                                    idx_Y, h_factor,
   		                                                                                    G_homotopy, B_homotopy,
   		                                                                                    enable_CM, ammeter, stamp_dual)
   
    	# Underground Lines
    for ele in range(len(undergroundline)):
    		ammeter = find_ammeter(undergroundline[ele], curr_meas, enable_CM)
    		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = undergroundline[ele].stamp_linear_homotopy(node, Ylin_h_val,
    		                                                                                       Ylin_h_row, Ylin_h_col,
    		                                                                                       idx_Y, h_factor,
    		                                                                                       G_homotopy,
    		                                                                                       B_homotopy, enable_CM,
    		                                                                                       ammeter, stamp_dual)
    
    	# Triplex Lines
    for ele in range(len(triplexline)):
    		ammeter = find_ammeter(triplexline[ele], curr_meas, enable_CM)
    		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y = triplexline[ele].stamp_linear_triplex_homotopy(node, Ylin_h_val,
    		                                                                                           Ylin_h_row,
    		                                                                                           Ylin_h_col, idx_Y,
    		                                                                                           h_factor,
    		                                                                                           G_homotopy,
    		                                                                                           B_homotopy,
    		                                                                                           enable_CM, ammeter, stamp_dual)
    
    	# Stamp Transformers
    for ele in range(len(xfmr)):
    		homotopy_option = True
    		ammeter = find_ammeter(xfmr[ele], curr_meas, enable_CM)
    		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y, stamped_ground = xfmr[ele].stamp_linear(Ylin_h_val, Ylin_h_row, Ylin_h_col,
    		                                                                   idx_Y, stamped_ground, enable_CM, ammeter,
    		                                                                   homotopy_option, h_factor,
    		                                                                   G_homotopy, B_homotopy, stamp_dual = stamp_dual)
    
    	# Stamp Regulators
    for ele in range(len(regulator)):
    		homotopy_option = True

    		Ylin_h_val, Ylin_h_row, Ylin_h_col, idx_Y, stamped_ground = regulator[ele].stamp_linear(node, V, Ylin_h_val, Ylin_h_row, Ylin_h_col,
    		                                                                   idx_Y, stamped_ground, manual_init,
    		                                                                   homotopy_option, h_factor,
    		                                                                   G_homotopy, B_homotopy)
    

    return Ylin_h_val, Ylin_h_row, Ylin_h_col, Jlin_h_row, Jlin_h_val, idx_Y, idx_J
