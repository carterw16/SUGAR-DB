'''
    Elements: This class serves as a superclass for linear and nonlinear elements. Created to ensure all linear and
    nonlinear elements which need stamping have their corresponding stamping functions.
    __author__ = ["Amritanshu Pandey", "Naeem Turner-Bandele"]
	__copyright__ = "Copyright 2019"
	__credits__ = ["Amritanshu Pandey", "Naeem Turner-Bandele"]
	__license__ = "GPL"
	__version__ = "1.0.1"
	__maainer__ = ["Amritanshu Pandey", "Naeem Turner-Bandele"]
	__email__ = "nturnerb@andrew.cmu.edu"
	__status__ = "Development"
'''


class LinearElement:

    @staticmethod
    def stamp_Y(i, j, val, Y_val, Y_row, Y_col, idx):
        if val == 0.:
            return Y_val, Y_row, Y_col, idx
        Y_val[idx] = val
        Y_row[idx] = i
        Y_col[idx] = j
        idx += 1

        return Y_val, Y_row, Y_col, idx

    @staticmethod
    def stamp_J(i, val, J_val, J_row, idx):
        if val == 0.:
            return J_val, J_row, idx        
        J_val[idx] = val
        J_row[idx] = i
        idx += 1

        return J_val, J_row, idx

    def __init__(self):
        pass


class NonlinearElement:

    @staticmethod
    def stamp_Y(i, j, val, Ynlin_val, Ynlin_row, Ynlin_col, idx):
        if val == 0.:
            return idx
        Ynlin_val[idx] = val
        Ynlin_row[idx] = i
        Ynlin_col[idx] = j
        idx += 1
        return idx

    @staticmethod
    def stamp_J(i, val, Jnlin_val, Jnlin_row, idx):
        if val == 0.:
            return idx
        Jnlin_val[idx] = val
        Jnlin_row[idx] = i
        idx += 1
        return idx

    def __init__(self):
        pass
