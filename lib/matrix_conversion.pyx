#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False

import numpy as np
from scipy.sparse import csr_matrix
cimport numpy as np

DTYPE = np.float64
ctypedef np.float_t DTYPE_t

def matrix_conversion(np.ndarray[DTYPE_t, ndim=1] Y_val,
					  np.ndarray[np.int64_t, ndim=1] Y_row,
					  np.ndarray[np.int64_t, ndim=1] Y_col,
					  np.ndarray[np.int64_t, ndim=1] J_row,
					  np.ndarray[np.int64_t, ndim=1] J_col,
					  np.ndarray[DTYPE_t, ndim=1] J_val,
					  int sizeY):

	Y = csr_matrix((Y_val, (Y_row, Y_col)), shape=(sizeY, sizeY), dtype=np.float64)
	if len(J_row) > 0:
		J = csr_matrix((J_val, (J_row, J_col)), shape=(sizeY, 1), dtype=np.float64)
	else:
		J = None

	return Y, J
