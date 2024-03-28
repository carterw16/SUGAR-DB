from classes.OperationalLimits.OperationalLimits_base import OperationalLimits

def stampY(i, j, val, Y_val, Y_row, Y_col, idx):
	Y_val[idx] = val
	Y_row[idx] = i
	Y_col[idx] = j
	idx += 1

	return idx


def stampJ(i, val, J_val, J_row, idx):
	J_val[idx] = val
	J_row[idx] = i
	idx += 1

	return idx 


class BatteryLimits(OperationalLimits):

	def __init__(self, source, node, P_max = 3000000, P_min = 0, SOC_max = 1, SOC_min = 0):
		# these should be settings set by parser or SUGAR3 settings
		self.node = node
		self.P_max = P_max
		self.P_min = P_min
		self.SOC_max = 1
		self.SOC_min = 0
		# add SOC max
		self.cs_eps = 1e-6


	def init_bounding_variables(self, Vinit):
		pass


	def stamp_nonlinear(self, node_key, nodes, V, lf, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
													 idx_Y, idx_J):
		return (idx_Y, idx_J)
	
	def calc_residuals(self, V, res_eqn):
		pass
