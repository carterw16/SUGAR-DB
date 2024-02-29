import numpy as np
from classes.OperationalLimits.VoltageLimits import VoltageLimits
from classes.Nodes import Nodes

def initialize_voltage_bounding(casedata, node_key, node_index_, v_bound_settings):
	# TODO whitelist/blacklist which buses are bounded
	voltage_bounds = []
	for ele in casedata.node:
		ele.assign_voltage_limit_nodes(node_index_)
		voltage_bounds.append(VoltageLimits(ele, v_bound_settings))
	casedata.voltage_bounds = voltage_bounds
	Nodes.vmag2_index = np.array(Nodes.vmag2_index)
	Nodes.Lvmag2_index = np.array(Nodes.Lvmag2_index)
	Nodes.umax_vmag2_index = np.array(Nodes.umax_vmag2_index)
	Nodes.umin_vmag2_index = np.array(Nodes.umin_vmag2_index)

	return node_index_