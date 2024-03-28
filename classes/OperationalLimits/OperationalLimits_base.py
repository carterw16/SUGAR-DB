
class OperationalLimits:
	
	def __init__(self, upper_bound, lower_bound):
		'''

		Parameters
		----------
		stamp_dual : boolean variable that flags whether or not to incorporate dual variables
		obj : integer value of either 1 for the L1-norm or 2 for the L2-norm that denotes which sets of equations to stamp

		Returns
		-------
		None.

		'''
		self.upper_bound = upper_bound
		self.lower_bound = lower_bound
	
