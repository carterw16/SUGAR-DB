import numpy as np
from classes.Nodes import Nodes
from classes.OperationalLimits.VoltageLimits import VoltageLimits
import ipdb


def operational_limits_diode_limiting(Vsol, V, cs_eps = 1e-5, type='bus_voltages', d=0.95, normalize=False, P_max = 1800000, P_min = 0, bat = None):
	if type == 'bus_voltages':
		xmin = np.array(VoltageLimits.pu_min2).reshape(-1,1)
		xmax = np.array(VoltageLimits.pu_max2).reshape(-1,1)
		x_index = Nodes.vmag2_index
		umin_index = Nodes.umin_vmag2_index
		umax_index = Nodes.umax_vmag2_index
	elif type == 'battery':
		if (bat != None):
			x_power_index = bat.infeas_real_var_index
			x_soc_index = bat.Bt_nodes
			x_index = bat.infeas_real_var_index + bat.Bt_nodes
			xmax = np.vstack(( P_max * np.ones((np.size(x_power_index,0),1)), 1 * np.ones((np.size(x_soc_index,0),1)) ))
			xmin = np.vstack(( P_min * np.ones((np.size(x_power_index,0),1)), 0 * np.ones((np.size(x_soc_index,0),1)) ))
			umin_index = bat.mu_index + bat.mu_index_Bt
			umax_index = bat.mu_index_upper + bat.mu_index_Bt_upper
		else:
			raise Exception("battery limiting without providing battery object")
	else:
		raise NotImplementedError

	if len(xmin):
		normalize_vector = np.maximum(np.abs(xmax), np.abs(xmin))
	else:
		normalize_vector = np.abs(xmax)

	# x_index is the vector of input bounding variable 
	# Find index where V greater than max
	V_big_idx = np.where(V[x_index]>xmax)[0]
	V[x_index][V_big_idx] = xmax[V_big_idx] - 1e-6
	# Find index where V less than min
	if len(xmin):
		V_small_idx = np.where(V[x_index]<xmin)[0]
		V[x_index][V_small_idx] = xmin[V_small_idx] + 1e-6
	
	if normalize:
		x_k = np.divide(V[x_index], normalize_vector)
		x_k1 = np.divide(Vsol[x_index], normalize_vector)
		xmax_orig = np.copy(xmax)
		xmax = np.divide(xmax, normalize_vector)
		if len(xmin):
			xmin_orig = np.copy(xmin)
			xmin = np.divide(xmin, normalize_vector)
		delta_x = x_k1 -  x_k
	else:
		x_k = V[x_index]
		x_k1 = Vsol[x_index]
		delta_x = x_k1 -  x_k

	umax_k = V[umax_index]
	umax_k1 = Vsol[umax_index]
	delta_umax = umax_k1 - umax_k
	umin_k = V[umin_index]
	umin_k1 = Vsol[umin_index]
	delta_umin = umin_k1 - umin_k

	# Find all the mu's that are greater than zero
	pos_umax= np.where(delta_umax>0)[0]
	neg_umax= np.where(delta_umax<0)[0]
	pos_umin= np.where(delta_umin>0)[0]
	neg_umin= np.where(delta_umin<0)[0]

	neg_x = np.where(delta_x<0)[0]
	pos_x = np.where(delta_x>0)[0]

	zero_x = np.where(delta_x==0)[0]
	zero_umax = np.where(delta_umax ==0)[0]
	zero_umin = np.where(delta_umin ==0)[0]
	delta_x[zero_x] = 1e-9
	delta_umax[zero_umax] = 1e-9
	delta_umin[zero_umin] = 1e-9

	# Find all mu's that are greater than zero
	if normalize:
		k_fac = 1e-2
	else:
		k_fac = 1e-2
	e_vec = k_fac*cs_eps*(np.ones(len(delta_umax)).reshape(-1,1))
	#d = 0.95 # was 0.98
	amax_u = np.ones(len(delta_umax)).reshape(-1,1)
	amax_u[neg_umax] = d*np.divide((e_vec[neg_umax] - umax_k[neg_umax]), delta_umax[neg_umax])
	amin_u = np.ones(len(delta_umin)).reshape(-1,1)
	amin_u[neg_umin] = d*np.divide((e_vec[neg_umin] - umin_k[neg_umin]), delta_umin[neg_umin])

	if len(xmin) == len(x_k) and len(xmax) == len(x_k): 
		amax_x = d*np.divide((xmax - x_k), delta_x)
		amin_x = d*np.divide((xmin - x_k), delta_x)
		a_x = np.fmin(1, np.maximum(amin_x, amax_x))
	elif len(xmin) == len(x_k):
		amin_x = d*np.divide((xmin - x_k), delta_x)
		a_x = np.fmin(1, amin_x)
		a_x[pos_x] = 1
	elif len(xmax) == len(x_k):
		amax_x = d*np.divide((xmax - x_k), delta_x)
		a_x = np.fmin(1,amax_x)
		a_x[neg_x] = 1
	delta_x[zero_x] = 0.0
	#delta_x[zero_umax] = 0.0
	#delta_x[zero_umin] = 0.0
	x_lim = (x_k + np.multiply(delta_x, a_x)).reshape(-1,1)

	# delta_umax[zero_umax] = 0.0
	# delta_umin[zero_umin] = 0.0
	amax_u = np.fmin(1, amax_u)
	amin_u = np.fmin(1, amin_u)
	
	umax_lim = umax_k + np.multiply(delta_umax, amax_u) 
	umin_lim = umin_k + np.multiply(delta_umin, amin_u) 

	V_lim = np.copy(Vsol)
	if normalize:
		V_lim[x_index] = np.multiply(np.copy(x_lim), normalize_vector)
		xmax = np.copy(xmax_orig)
		if len(xmin):
			xmin = np.copy(xmin_orig)
	else:
		V_lim[x_index] = np.copy(x_lim)
	V_lim[umax_index] = np.copy(umax_lim)
	V_lim[umin_index] = np.copy(umin_lim)

	#print("mu upper values (after limiting)", V_lim[umin_index])
	#print("mu lower values (after limiting)", V_lim[umax_index])
	#print("P values after Pch/Pd", V[x_index])

	
	return (V_lim)

def calc_voltage_violations(V, node, vmax_pu, vmin_pu):
	upper_bound_violators = {}
	lower_bound_violators = {}
	for ele in node:
		ub = ele.Vnom*vmax_pu
		lb = ele.Vnom*vmin_pu
		if ele.isTriplex:
			V1mag = np.sqrt(V[ele.node1_Vr]**2 + V[ele.node1_Vi]**2)
			V2mag = np.sqrt(V[ele.node2_Vr]**2 + V[ele.node2_Vi]**2)
			Vmags = np.array([V1mag, V2mag]).reshape(-1)
		else:
			VAmag = np.sqrt(V[ele.nodeA_Vr]**2 + V[ele.nodeA_Vi]**2)
			VBmag = np.sqrt(V[ele.nodeB_Vr]**2 + V[ele.nodeB_Vi]**2)
			VCmag = np.sqrt(V[ele.nodeC_Vr]**2 + V[ele.nodeC_Vi]**2)
			Vmags = np.array([VAmag, VBmag, VCmag]).reshape(-1)
		# remove any 0s from mags to handle unconnected nodes
		Vmags_nonzero = Vmags[Vmags != 0.0]
		if len(Vmags_nonzero) == 0:
			# all phases must have been unconnected?
			continue
		# print(ele.name, ele.Vnom_ln, ele.Vnom, Vmags/ele.Vnom_ln)
		if np.amax(Vmags_nonzero) > ub:
			upper_bound_violators[ele.name] = np.amax(Vmags_nonzero)/ele.Vnom
		if np.amin(Vmags_nonzero) < lb:
			lower_bound_violators[ele.name] = np.amin(Vmags_nonzero)/ele.Vnom
	print("%d nodes had a voltage mag above %.2f pu of nominal val" % (len(upper_bound_violators), vmax_pu))
	if len(upper_bound_violators) > 0:
		worst_bus = max(upper_bound_violators, key=upper_bound_violators.get)
		worst_bus_val = upper_bound_violators[worst_bus]
		print("Worst upper violation: %.3f p.u. at bus %s" % (worst_bus_val, worst_bus))
	# pp.pprint(upper_bound_violators)
	print("%d nodes had a voltage mag below %.2f pu of nominal val" % (len(lower_bound_violators), vmin_pu))
	# pp.pprint(lower_bound_violators)
	if len(lower_bound_violators) > 0:
		worst_bus = min(lower_bound_violators, key=lower_bound_violators.get)
		worst_bus_val = lower_bound_violators[worst_bus]
		print("Worst lower violation: %.3f p.u. at bus %s" % (worst_bus_val, worst_bus))

