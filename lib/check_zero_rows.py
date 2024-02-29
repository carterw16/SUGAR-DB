from __future__ import print_function

import numpy as np


def check_zero_rows(Y, J):
    #print ('Shape of Y without zero rows removed' + str(Y.shape))

    def remove_zero_rows(X):
        # X is a scipy sparse matrix. We want to remove all zero rows from it
        _nonzero_row_index, _ = X.nonzero()
        _unique_nonzero_index = np.unique(_nonzero_row_index)
        return X[_unique_nonzero_index], _unique_nonzero_index

    (Y_red, unique_nonzero_index) = remove_zero_rows(Y)
    J_red = J[unique_nonzero_index]
    #print ('Shape of Y with zero rows removed' + str(Y_red.shape))

    (Y_red, unique_nonzero_index) = remove_zero_rows(Y_red.transpose())
    # print ('Shape of Y with zero columns removed' + str(Y_red.shape))

    if Y_red.shape[0] != Y_red.shape[1]:
        print('Shape of Y is non-square' + str(Y_red.shape))
        # # Take CSR matrix
        num_nonzeros = np.diff(Y.indptr)
        rowZero = [i for i, x in enumerate(num_nonzeros) if not x]
        print('Rows with all zero' + str(rowZero))

        # # Convert csr to csc matrix
        Y = Y.tocsc()

        num_nonzeros = np.diff(Y.indptr)
        colZero = [i for i, x in enumerate(num_nonzeros) if not x]
        print('Columns with all zero' + str(colZero))

    Y_red = Y_red.transpose()

    # print('Arrays Equal', np.array_equal(colZero, rowZero))
    return Y_red, J_red, unique_nonzero_index
