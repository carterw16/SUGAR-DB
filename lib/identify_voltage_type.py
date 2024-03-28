"""Handles identification of LN or LL voltages.

Author(s): Naeem Turner-Bandele
Created Date: 03-26-2021
Updated Date: 03-26-2021
Email: nturnerb@cmu.edu
Status: Development

"""

import numpy as np


def identify_voltage_type(Vnom):
	"""Determines whether the voltage is line to neutral or line to line.

	Args:
		Vnom (float): nominal voltage

	Returns:
		voltage_type (str): specifies if voltage is Line to Line (LL) or Line to Neutral (LN)
		Vnom (float): nominal voltage as LN

	"""
	dist_voltages = [
		0.208, 0.240, 0.4, 0.416, 0.480, 4.16, 4.8, 12.47, 13.2, 13.8, 24.94, 34.5, 69.0, 115.0, 23, 
		138.0, 230.0
	]
	voltage_type = None

	for ele in dist_voltages:
		Vnom_kV = Vnom * 1E-3
		if np.isclose(Vnom_kV, ele, 1e-4):
			voltage_type = 'LL'
			Vnom = Vnom / np.sqrt(3)

	if voltage_type is None:
		voltage_type = 'LN'

	return voltage_type, Vnom
