from classes.OperationalLimits.OperationalLimits_base import OperationalLimits
from classes.Nodes import Nodes
import math

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


class VoltageUnbalance(OperationalLimits):
	def __init__(self, voltage_unbalance_settings):
		# Upper bound on voltage unbalance should be either 2 or 3% depending on the definition
		# Lower bound should be 0 for voltage unbalance
		self.cs_tolerance = voltage_unbalance_settings['cs eps tol'] 
		OperationalLimits.__init__(self, voltage_unbalance_settings['upper bound'], voltage_unbalance_settings['lower bound'])
	
	def initialize_voltage_unbalance_variables(self, node, V_init):
		if node.phases == 7 or node.phases == 15:
			node_dict = self.create_node_dict(node)
			V_init[node_dict['node_Vunb_pos_real']] = math.sqrt(node.Vnom)
			V_init[node_dict['node_Vunb_pos_imag']] = math.sqrt(node.Vnom)
			V_init[node_dict['node_Vunb_neg_real']] = 1e-6
			V_init[node_dict['node_Vunb_neg_imag']] = 1e-6
			V_init[node_dict['node_Vunb_tracking_sqd']] = 1e-6 # math.sqrt((V_init[node_dict['node_Vunb_neg_real']]**2 + V_init[node_dict['node_Vunb_neg_imag']]**2)/
													   # (V_init[node_dict['node_Vunb_pos_real']]**2 + V_init[node_dict['node_Vunb_pos_imag']]**2))

			V_init[node_dict['node_unb_mu_upper']] = 1e-6

			V_init[node_dict['node_unb_Lpr']] = 1e-6
			V_init[node_dict['node_unb_Lpi']] = 1e-6
			V_init[node_dict['node_unb_Lnr']] = 1e-6			
			V_init[node_dict['node_unb_Lni']] = 1e-6
			V_init[node_dict['node_unb_L_tracking']] = 1e-6 # self.cs_tolerance/(2*V_init[node_dict['node_Vunb_tracking_sqd']]*(V_init[node_dict['node_Vunb_tracking_sqd']] - .02))

		return V_init
	
	def stamp_voltage_unbalance_constraints(self, node, V, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
		node_dict = self.create_node_dict(node)
		cs_eps = self.cs_tolerance# *node.Vnom

		partials = self.calc_partial_derivatives(node_dict, V, cs_eps)

		# Stamp the voltage unbalance lagrangian equations
		Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_stationarity_constrations(node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J)
		Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_equality_constraints(node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J)
		Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J = self.stamp_complementary_slackness_constraints(node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J)
		
		return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

	def calc_partial_derivatives(self, node_dict, V, cs_eps):
		partials = {}
		
		# TODO need to accommodate when A, B, C are not in 0, 120, 240 order
		Var = V[node_dict['nodeA_Vr']]
		Vai = V[node_dict['nodeA_Vi']]
		Vbr = V[node_dict['nodeB_Vr']]
		Vbi = V[node_dict['nodeB_Vi']]
		Vcr = V[node_dict['nodeC_Vr']]
		Vci = V[node_dict['nodeC_Vi']]

		if Vbi > 0:
			print("Vbr is greater than 0. B should lag A by 120.")

		Vpr = V[node_dict['node_Vunb_pos_real']]		
		Vpi = V[node_dict['node_Vunb_pos_imag']]
		Vnr = V[node_dict['node_Vunb_neg_real']]
		Vni = V[node_dict['node_Vunb_neg_imag']]
		Vunb2 = V[node_dict['node_Vunb_tracking_sqd']]

		mu_ub = V[node_dict['node_unb_mu_upper']]

		lambda_pr = V[node_dict['node_unb_Lpr']]
		lambda_pi = V[node_dict['node_unb_Lpi']]
		lambda_nr = V[node_dict['node_unb_Lnr']]
		lambda_ni = V[node_dict['node_unb_Lni']]
		lambda_unb = V[node_dict['node_unb_L_tracking']]

		# TAYLOR APPROXIMATION - VR, VI STATIONARITY CONSTRAINTs
		# f_dVar = lambda_pr + lambda_nr
		partials['d2L_dVar2'] = 0
		partials['d2L_dVardLpr'] = 1
		partials['d2L_dVardLnr'] = 1
		
		# f_dVai = lambda_pi + lambda_ni
		partials['d2L_dVai2'] = 0
		partials['d2L_dVaidLpi'] = 1
		partials['d2L_dVaidLni'] = 1

		# f_dVbr = (-0.5*lambda_pr + math.sqrt(3)/2*lambda_pi - 0.5*lambda_nr - math.sqrt(3)/2*lambda_ni)
		partials['d2L_dVbr2'] = 0 
		partials['d2L_dVbrdLpr'] = -0.5
		partials['d2L_dVbrdLpi'] = math.sqrt(3)/2
		partials['d2L_dVbrdLnr'] = -0.5
		partials['d2L_dVbrdLni'] = -math.sqrt(3)/2
		
		# f_dVbi = (-math.sqrt(3)/2*lambda_pr - 0.5*lambda_pi + math.sqrt(3)/2*lambda_nr - 0.5*lambda_ni)
		partials['d2L_dVbi2'] = 0 
		partials['d2L_dVbidLpr'] = -math.sqrt(3)/2
		partials['d2L_dVbidLpi'] =-0.5
		partials['d2L_dVbidLnr'] = math.sqrt(3)/2
		partials['d2L_dVbidLni'] = -0.5
		
		# f_dVcr = (-0.5*lambda_pr - math.sqrt(3)/2*lambda_pi - 0.5*lambda_nr + math.sqrt(3)/2*lambda_ni)
		partials['d2L_dVcr2'] = 0 
		partials['d2L_dVcrdLpr'] = -0.5
		partials['d2L_dVcrdLpi'] = -math.sqrt(3)/2
		partials['d2L_dVcrdLnr'] = -0.5
		partials['d2L_dVcrdLni'] = math.sqrt(3)/2
		
		# f_dVci = (math.sqrt(3)/2*lambda_pr - 0.5*lambda_pi - math.sqrt(3)/2*lambda_nr - 0.5*lambda_ni)
		partials['d2L_dVci2'] = 0 
		partials['d2L_dVcidLpr'] = math.sqrt(3)/2
		partials['d2L_dVcidLpi'] = -0.5
		partials['d2L_dVcidLnr'] = -math.sqrt(3)/2
		partials['d2L_dVcidLni'] = -0.5
		
		# TAYLOR APPROXIMATION - VP, VN STATIONARITY CONSTRAINTS
		Vp_mag = Vpr**2 + Vpi**2
		Vp_mag2 = Vp_mag**2
		Vp_mag3 = Vp_mag2*Vp_mag

		Vn_mag = Vnr**2 + Vni**2

		# dL/dVpr = -3*lambda_pr - lambda_unb*(2*Vpr*(Vn_mag))/Vp_mag2
		f_dVpr =  -3*lambda_pr - lambda_unb*(2*Vpr*(Vn_mag))/Vp_mag2
		partials['d2L_dVpr2'] = lambda_unb*Vn_mag*(-2*Vp_mag + 8*Vpr**2)/Vp_mag3
		partials['d2L_dVprdVpi'] = lambda_unb*Vn_mag*(8*Vpr*Vpi)/Vp_mag3
		partials['d2L_dVprdVnr'] = lambda_unb*(-4*Vpr*Vnr)/Vp_mag2
		partials['d2L_dVprdVni'] = lambda_unb*(-4*Vpr*Vni)/Vp_mag2
		partials['d2L_dVprdLunb'] = (-2*Vpr*Vn_mag)/Vp_mag2
		partials['d2L_dVprdLpr'] = -3
		partials['dL_dVpr_hist'] = (-f_dVpr + partials['d2L_dVpr2']*Vpr + partials['d2L_dVprdVpi']*Vpi + 
							  partials['d2L_dVprdVnr']*Vnr + partials['d2L_dVprdVni']*Vni + 
							  partials['d2L_dVprdLunb']*lambda_unb + partials['d2L_dVprdLpr']*lambda_pr) 

		# dL/dVpi = -3*lambda_pi - lambda_unb*(2*Vpi*(Vn_mag))/Vp_mag2	     
		f_dVpi =  -3*lambda_pi  - lambda_unb*(2*Vpi*(Vn_mag))/Vp_mag2
		partials['d2L_dVpidVpr'] = lambda_unb*Vn_mag*(8*Vpr*Vpi)/Vp_mag3
		partials['d2L_dVpi2'] = lambda_unb*Vn_mag*(-2*Vp_mag + 8*Vpi**2)/Vp_mag3
		partials['d2L_dVpidVnr'] = lambda_unb*(-4*Vpi*Vnr)/Vp_mag2
		partials['d2L_dVpidVni'] = lambda_unb*(-4*Vpi*Vni)/Vp_mag2
		partials['d2L_dVpidLunb'] = (-2*Vpi*Vn_mag)/Vp_mag2
		partials['d2L_dVpidLpi'] = -3
		partials['dL_dVpi_hist'] = (-f_dVpi + partials['d2L_dVpidVpr']*Vpr + partials['d2L_dVpi2']*Vpi + 
							  partials['d2L_dVpidVnr']*Vnr + partials['d2L_dVpidVni']*Vni + 
							  partials['d2L_dVpidLunb']*lambda_unb  + partials['d2L_dVpidLpi']*lambda_pi) 
		
		# dL/dVnr = -3*lambda_nr + 2*Vnr/Vp_mag*lambda_unb
		f_dVnr =  -3*lambda_nr + (2*lambda_unb*Vnr)/Vp_mag
		partials['d2L_dVnrdVpr'] = lambda_unb*(-4*Vpr*Vnr)/Vp_mag2
		partials['d2L_dVnrdVpi'] = lambda_unb*(-4*Vpi*Vnr)/Vp_mag2
		partials['d2L_dVnr2'] = (2*lambda_unb)/Vp_mag2
		partials['d2L_dVnrdLunb'] = 2*Vnr/Vp_mag
		partials['d2L_dVnrdLnr'] = -3
		partials['dL_dVnr_hist'] = (-f_dVnr + partials['d2L_dVnrdVpr']*Vpr + partials['d2L_dVnrdVpi']*Vpi + 
							  partials['d2L_dVnr2']*Vnr + partials['d2L_dVnrdLunb']*lambda_unb + partials['d2L_dVnrdLnr']*lambda_nr) 

		# dL/dVni = -3*lambda_ni + 2*Vni/Vp_mag*lambda_unb	     
		f_dVni =  -3*lambda_ni + (2*lambda_unb*Vni)/Vp_mag
		partials['d2L_dVnidVpr'] = lambda_unb*(-4*Vpr*Vni)/Vp_mag2
		partials['d2L_dVnidVpi'] = lambda_unb*(-4*Vpi*Vni)/Vp_mag2
		partials['d2L_dVni2'] = (2*lambda_unb)/Vp_mag2
		partials['d2L_dVnidLunb'] = 2*Vnr/Vp_mag
		partials['d2L_dVnidLni'] = -3
		partials['dL_dVni_hist'] = (-f_dVni + partials['d2L_dVnidVpr']*Vpr + partials['d2L_dVnidVpi']*Vpi + 
							  partials['d2L_dVni2']*Vni + partials['d2L_dVnidLunb']*lambda_unb + partials['d2L_dVnidLni']*lambda_ni) 
		
		# TAYLOR APPROXIMATION - V_UNB STATIONARITY CONSTRAINT
		# f_dVunb = -lambda_unb + mu
		f_dVunb = -lambda_unb + mu_ub*(100)**2
		partials['dL_dVunb2'] = 0
		partials['dL_dVunbdLunb'] = -1
		partials['dL_dVunbdmu'] = (100)**2
		partials['dL_dVunb_hist'] = (-f_dVunb + partials['dL_dVunb2']*Vunb2
								+ partials['dL_dVunbdLunb']*lambda_unb + partials['dL_dVunbdmu']*mu_ub)

		# TAYLOR APPROXIMATION - V_UNB EQUALITY CONSTRAINT
		# f_unb: -Vunb2 + Vn_mag/Vp_mag = 0
		f_Vunb_hist = -Vunb2 + Vn_mag/Vp_mag
		partials['fVunb_dVunb'] = -1
		partials['fVunb_dVpr'] = (-2*Vpr*Vn_mag)/Vp_mag2
		partials['fVunb_dVpi'] = (-2*Vpi*Vn_mag)/Vp_mag2
		partials['fVunb_dVnr'] = 2*Vnr/Vp_mag
		partials['fVunb_dVni'] = 2*Vni/Vp_mag
		partials['fVunb_hist'] = (-f_Vunb_hist + partials['fVunb_dVunb']*Vunb2 + 
							partials['fVunb_dVpr']*Vpr + partials['fVunb_dVpi']*Vpi + 
							partials['fVunb_dVnr']*Vnr + partials['fVunb_dVni']*Vni)

		
		# TAYLOR APPROXIMATION - COMPLEMENTARY SLACKNESS
		# f_unb_ub = mu_ub*(Vunb - .02**2) - eps = 0
		f_unb_ub = mu_ub*(Vunb2 - (self.upper_bound + cs_eps)**2/100**2) 
		partials['dcs_ub_dVunb'] = mu_ub
		partials['dcs_ub_dmu_ub'] = (Vunb2 - (self.upper_bound + cs_eps)**2/100**2)
		partials['dcs_ub_hist'] = (-f_unb_ub + partials['dcs_ub_dVunb']*Vunb2 + 
							+ partials['dcs_ub_dmu_ub']*mu_ub)

		return partials
	
	def calc_residuals(self, V, node, res_eqn):
		if node.phases == 7 or node.phases == 15:
			Vpr = V[node.extra_unbalance_node_pos_real]
			Vpi = V[node.extra_unbalance_node_pos_imag]
			Vnr = V[node.extra_unbalance_node_neg_real]
			Vni = V[node.extra_unbalance_node_neg_imag]

			Lpr = V[node.dual_unb_eq_var_pr]
			Lpi = V[node.dual_unb_eq_var_pi]
			Lnr = V[node.dual_unb_eq_var_nr]
			Lni = V[node.dual_unb_eq_var_ni]

			Vunb2 = V[node.extra_voltage_unbalance_node_sqd]
			L_unb = V[node.dual_unb_eq_var_total]
			mu_unb = V[node.dual_unb_ineq_var_ub]

			Var = V[node.nodeA_Vr]
			Vai = V[node.nodeA_Vi]
			Vbr = V[node.nodeB_Vr]
			Vbi = V[node.nodeB_Vi]
			Vcr = V[node.nodeC_Vr]
			Vci = V[node.nodeC_Vi]

			# Equality Constraints
			res_eqn[node.extra_unbalance_node_pos_real] += -Vpr + Var - 0.5*Vbr - math.sqrt(3)/2*Vbi - 0.5*Vcr + math.sqrt(3)/2*Vci
			res_eqn[node.extra_unbalance_node_pos_imag] += -Vpi + Vai + math.sqrt(3)/2*Vbr - 0.5*Vbi - math.sqrt(3)/2*Vcr - 0.5*Vci
			res_eqn[node.extra_unbalance_node_neg_real] += -Vnr + Var - 0.5*Vbr + math.sqrt(3)/2*Vbi - 0.5*Vcr - math.sqrt(3)/2*Vci
			res_eqn[node.extra_unbalance_node_neg_imag] += -Vni + Vai - math.sqrt(3)/2*Vbr - 0.5*Vbi + math.sqrt(3)/2*Vcr - 0.5*Vci

			res_eqn[node.extra_voltage_unbalance_node_sqd] += -Vunb2 + (Vnr**2 + Vni**2)/(Vpr**2 + Vpi**2)

			# Inequality Constraints [Complementary Slackness]
			res_eqn[node.dual_unb_ineq_var_ub] += mu_unb*(Vunb2 - self.upper_bound**2) + self.cs_tolerance

			# Stationarity Constraints
			# f_dVar = lambda_pr + lambda_nr
			res_eqn[node.nodeA_dual_eq_var_r] = Lpr + Lnr
			# f_dVai = lambda_pi + lambda_ni
			res_eqn[node.nodeA_dual_eq_var_i] = Lpi + Lni
			# f_dVbr = (-0.5*lambda_pr + math.sqrt(3)/2*lambda_pi - 0.5*lambda_nr - math.sqrt(3)/2*lambda_ni)
			res_eqn[node.nodeB_dual_eq_var_r] = -0.5*(Lpr + Lnr) + math.sqrt(3)/2*(Lpi - Lni)
			# f_dVbi = (-math.sqrt(3)/2*lambda_pr - 0.5*lambda_pi + math.sqrt(3)/2*lambda_nr - 0.5*lambda_ni)
			res_eqn[node.nodeB_dual_eq_var_i] = -0.5*(Lpi + Lni) - math.sqrt(3)/2*(Lpr - Lnr)
			# f_dVcr = (-0.5*lambda_pr - math.sqrt(3)/2*lambda_pi - 0.5*lambda_nr + math.sqrt(3)/2*lambda_ni)
			res_eqn[node.nodeC_dual_eq_var_r] = -0.5*(Lpr + Lnr) - math.sqrt(3)/2*(Lpi - Lni)
			# f_dVci = (math.sqrt(3)/2*lambda_pr - 0.5*lambda_pi - math.sqrt(3)/2*lambda_nr - 0.5*lambda_ni)
			res_eqn[node.nodeC_dual_eq_var_i] = -0.5*(Lpi + Lni) + math.sqrt(3)/2*(Lpr - Lnr)
			
			# dL/dVpr = -3*lambda_pr - lambda_unb*(2*Vpr*(Vn_mag))/Vp_mag2
			res_eqn[node.dual_unb_eq_var_pr] = -3*Lpr - L_unb*(2*Vpr*(Vnr**2 + Vni**2)/(Vpr**2 + Vpi**2))
			# dL/dVpi = -3*lambda_pi - lambda_unb*(2*Vpi*(Vn_mag))/Vp_mag2
			res_eqn[node.dual_unb_eq_var_pi] = -3*Lpi - L_unb*(2*Vpi*(Vnr**2 + Vni**2)/(Vpr**2 + Vpi**2))
			# dL/dVnr = -3*lambda_nr + 2*Vnr/Vp_mag*lambda_unb
			res_eqn[node.dual_unb_eq_var_nr] = -3*Lnr + 2*L_unb*Vnr
			# dL/dVni = -3*lambda_ni + 2*Vni/Vp_mag*lambda_unb	
			res_eqn[node.dual_unb_eq_var_ni] = -3*Lni + 2*L_unb*Vni
			# f_dVunb = -lambda_unb + mu
			res_eqn[node.dual_unb_eq_var_pr] = -L_unb + mu_unb
		
		return res_eqn
	
	def create_node_dict(self, node):
		# Note: in future, could amend dict for node_dict['node_Vr] = [A, B, C] for easier looping
		node_dict = {}

		# Primal nodes for real/imag Vp, Vn (unbalance nodes)
		node_dict['node_Vunb_pos_real'] = node.extra_unbalance_node_pos_real
		node_dict['node_Vunb_pos_imag'] = node.extra_unbalance_node_pos_imag
		node_dict['node_Vunb_neg_real'] = node.extra_unbalance_node_neg_real
		node_dict['node_Vunb_neg_imag'] = node.extra_unbalance_node_neg_imag

		# Primal node for tracking unbalance
		node_dict['node_Vunb_tracking_sqd'] = node.extra_voltage_unbalance_node_sqd
		node_dict['node_unb_L_tracking'] = node.dual_unb_eq_var_total
		

		# Primal nodes for Vr, Vi
		node_dict['nodeA_Vr'] = node.nodeA_Vr
		node_dict['nodeA_Vi'] = node.nodeA_Vi
		node_dict['nodeB_Vr'] = node.nodeB_Vr
		node_dict['nodeB_Vi'] = node.nodeB_Vi
		node_dict['nodeC_Vr'] = node.nodeC_Vr
		node_dict['nodeC_Vi'] = node.nodeC_Vi

		# Dual nodes for dL/dVp, dL/dVn
		node_dict['node_unb_Lpr'] = node.dual_unb_eq_var_pr
		node_dict['node_unb_Lpi'] = node.dual_unb_eq_var_pi
		node_dict['node_unb_Lnr'] = node.dual_unb_eq_var_nr
		node_dict['node_unb_Lni'] = node.dual_unb_eq_var_ni

		# Dual nodes for dL/dVr, dL/dVi
		node_dict['nodeA_Lr'] = node.nodeA_dual_eq_var_r
		node_dict['nodeA_Li'] = node.nodeA_dual_eq_var_i
		node_dict['nodeB_Lr'] = node.nodeB_dual_eq_var_r
		node_dict['nodeB_Li'] = node.nodeB_dual_eq_var_i
		node_dict['nodeC_Lr'] = node.nodeC_dual_eq_var_r
		node_dict['nodeC_Li'] = node.nodeC_dual_eq_var_i

		# Dual inquality constraint variables
		# Nodes associated with complementary slackness equations
		node_dict['node_unb_mu_upper'] = node.dual_unb_ineq_var_ub
	    
		return node_dict

	def stamp_equality_constraints(self, node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
		# PRIMAL FEASIBILITY - EQUALITY CONSTRAINTS 
		# Stamping the equality constraint where Vp = f(Va, Vb, Vc) to Vp nodes
		# f_unb_pos_real: (-3*Vpr + Var - 0.5*Vbr - math.sqrt(3)/2*Vbr - 0.5*Vcr + math.sqrt(3)/2*Vcr) = 0
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['node_Vunb_pos_real'], -1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['nodeA_Vr'], 1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['nodeB_Vr'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['nodeB_Vi'], -math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['nodeC_Vr'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_real'], node_dict['nodeC_Vi'], math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)

		# f_unb_pos_imag: (-3*Vpi + Vai + math.sqrt(3)/2*Vbr - 0.5*Vbr - math.sqrt(3)/2*Vcr -0.5*Vcr) = 0
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['node_Vunb_pos_imag'], -1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['nodeA_Vi'], 1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['nodeB_Vr'], math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['nodeB_Vi'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['nodeC_Vr'], -math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_pos_imag'], node_dict['nodeC_Vi'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		
		# f_unb_neg_real: (-3*Vpr + Var - 0.5*Vbr + math.sqrt(3)/2*Vbr - 0.5*Vcr - math.sqrt(3)/2*Vcr) = 0
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['node_Vunb_neg_real'], -1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['nodeA_Vr'], 1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['nodeB_Vr'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['nodeB_Vi'], math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['nodeC_Vr'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_real'], node_dict['nodeC_Vi'], -math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)

		# f_unb_neg_imag: (-3*Vpi + Vai - math.sqrt(3)/2*Vbr - 0.5*Vbr + math.sqrt(3)/2*Vcr -0.5*Vcr) = 0
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['node_Vunb_neg_imag'], -1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['nodeA_Vi'], 1, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['nodeB_Vr'], -math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['nodeB_Vi'], -0.5, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['nodeC_Vr'], math.sqrt(3)/2, Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_neg_imag'], node_dict['nodeC_Vi'], -0.5, Y_val, Y_row, Y_col, idx_Y)

		# f_unb: -Vunb + Vn_mag/Vp_mag = 0 (NONLINEAR)
		idx_Y = stampY(node_dict['node_Vunb_tracking_sqd'], node_dict['node_Vunb_tracking_sqd'], partials['fVunb_dVunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_tracking_sqd'], node_dict['node_Vunb_pos_real'], partials['fVunb_dVpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_tracking_sqd'], node_dict['node_Vunb_pos_imag'], partials['fVunb_dVpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_tracking_sqd'], node_dict['node_Vunb_neg_real'], partials['fVunb_dVnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_Vunb_tracking_sqd'], node_dict['node_Vunb_neg_imag'], partials['fVunb_dVni'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_Vunb_tracking_sqd'], partials['fVunb_hist'], J_val, J_row, idx_J)

		return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J
	
	def stamp_complementary_slackness_constraints(self, node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
		# COMPLEMENTARY SLACKNESS (NONLINEAR)
		# Stamping the upper bound equations: mu_upper_bound * (Vunb - upper_bound) + eps = 0
		idx_Y = stampY(node_dict['node_unb_mu_upper'], node_dict['node_unb_mu_upper'], partials['dcs_ub_dmu_ub'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_mu_upper'], node_dict['node_Vunb_tracking_sqd'], partials['dcs_ub_dVunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_unb_mu_upper'], partials['dcs_ub_hist'], J_val, J_row, idx_J)
		
		return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J

	def stamp_stationarity_constrations(self, node_dict, partials, Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J):
		# STATIONARITY CONSTRAINTS (dL/dVar, dL/dVai, dL/dVbr, dL/dVbi, dL/dVcr, dL/dVci, dL/dVp, dL/dVn)
		# TODO - write with matrix format? loop?

		# Amending the dL/dVr equations to dual_real nodes (LINEAR EQUATIONS)
		idx_Y = stampY(node_dict['nodeA_Lr'], node_dict['node_unb_Lpr'], partials['d2L_dVardLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeA_Lr'], node_dict['node_unb_Lnr'], partials['d2L_dVardLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeA_Lr'], node_dict['nodeA_Vr'], partials['d2L_dVar2'], Y_val, Y_row, Y_col, idx_Y)

		idx_Y = stampY(node_dict['nodeB_Lr'], node_dict['node_unb_Lpr'], partials['d2L_dVbrdLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Lr'], node_dict['node_unb_Lpi'], partials['d2L_dVbrdLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Lr'], node_dict['node_unb_Lnr'], partials['d2L_dVbrdLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Lr'], node_dict['node_unb_Lni'], partials['d2L_dVbrdLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Lr'], node_dict['nodeB_Vr'], partials['d2L_dVbr2'], Y_val, Y_row, Y_col, idx_Y)
		
		idx_Y = stampY(node_dict['nodeC_Lr'], node_dict['node_unb_Lpr'], partials['d2L_dVcrdLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Lr'], node_dict['node_unb_Lpi'], partials['d2L_dVcrdLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Lr'], node_dict['node_unb_Lnr'], partials['d2L_dVcrdLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Lr'], node_dict['node_unb_Lni'], partials['d2L_dVcrdLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Lr'], node_dict['nodeC_Vr'], partials['d2L_dVcr2'], Y_val, Y_row, Y_col, idx_Y)
		
		# Amending the dL/dVi equations to dual_imag nodes (LINEAR EQUATIONS)
		idx_Y = stampY(node_dict['nodeA_Li'], node_dict['node_unb_Lpi'], partials['d2L_dVaidLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeA_Li'], node_dict['node_unb_Lni'], partials['d2L_dVaidLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeA_Li'], node_dict['nodeA_Vi'], partials['d2L_dVai2'], Y_val, Y_row, Y_col, idx_Y)
		
		idx_Y = stampY(node_dict['nodeB_Li'], node_dict['node_unb_Lpr'], partials['d2L_dVbidLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Li'], node_dict['node_unb_Lpi'], partials['d2L_dVbidLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Li'], node_dict['node_unb_Lnr'], partials['d2L_dVbidLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Li'], node_dict['node_unb_Lni'], partials['d2L_dVbidLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeB_Li'], node_dict['nodeB_Vi'], partials['d2L_dVbi2'], Y_val, Y_row, Y_col, idx_Y)
		
		idx_Y = stampY(node_dict['nodeC_Li'], node_dict['node_unb_Lpr'], partials['d2L_dVcidLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Li'], node_dict['node_unb_Lpi'], partials['d2L_dVcidLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Li'], node_dict['node_unb_Lnr'], partials['d2L_dVcidLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Li'], node_dict['node_unb_Lni'], partials['d2L_dVcidLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['nodeC_Li'], node_dict['nodeC_Vi'], partials['d2L_dVci2'], Y_val, Y_row, Y_col, idx_Y)
		
		# Stamping the dL/dVp equations to lambda_p nodes (NONLINEAR EQUATIONS)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_Vunb_pos_real'], partials['d2L_dVpr2'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_Vunb_pos_imag'], partials['d2L_dVprdVpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_Vunb_neg_real'], partials['d2L_dVprdVnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_Vunb_neg_imag'], partials['d2L_dVprdVni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_unb_L_tracking'], partials['d2L_dVprdLunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpr'], node_dict['node_unb_Lpr'], partials['d2L_dVprdLpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_unb_Lpr'], partials['dL_dVpr_hist'], J_val, J_row, idx_J)
		
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_Vunb_pos_real'], partials['d2L_dVpidVpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_Vunb_pos_imag'], partials['d2L_dVpi2'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_Vunb_neg_real'], partials['d2L_dVpidVnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_Vunb_neg_imag'], partials['d2L_dVpidVni'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_unb_L_tracking'], partials['d2L_dVpidLunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lpi'], node_dict['node_unb_Lpi'], partials['d2L_dVpidLpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_unb_Lpi'], partials['dL_dVpi_hist'], J_val, J_row, idx_J)
		
		idx_Y = stampY(node_dict['node_unb_Lnr'], node_dict['node_Vunb_pos_real'], partials['d2L_dVnrdVpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lnr'], node_dict['node_Vunb_pos_imag'], partials['d2L_dVnrdVpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lnr'], node_dict['node_Vunb_neg_real'], partials['d2L_dVnr2'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lnr'], node_dict['node_unb_L_tracking'], partials['d2L_dVnrdLunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lnr'], node_dict['node_unb_Lnr'], partials['d2L_dVnrdLnr'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_unb_Lnr'], partials['dL_dVnr_hist'], J_val, J_row, idx_J)
		
		idx_Y = stampY(node_dict['node_unb_Lni'], node_dict['node_Vunb_pos_real'], partials['d2L_dVnidVpr'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lni'], node_dict['node_Vunb_pos_imag'], partials['d2L_dVnidVpi'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lni'], node_dict['node_Vunb_neg_imag'], partials['d2L_dVni2'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lni'], node_dict['node_unb_L_tracking'], partials['d2L_dVnidLunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_Lni'], node_dict['node_unb_Lni'], partials['d2L_dVnidLni'], Y_val, Y_row, Y_col, idx_Y)
		idx_J = stampJ(node_dict['node_unb_Lni'], partials['dL_dVni_hist'], J_val, J_row, idx_J)
		
		# dL_dVunb
		idx_Y = stampY(node_dict['node_unb_L_tracking'], node_dict['node_unb_L_tracking'], partials['dL_dVunbdLunb'], Y_val, Y_row, Y_col, idx_Y)
		idx_Y = stampY(node_dict['node_unb_L_tracking'], node_dict['node_unb_mu_upper'], partials['dL_dVunbdmu'], Y_val, Y_row, Y_col, idx_Y)
		
		return Y_val, Y_row, Y_col, J_val, J_row, idx_Y, idx_J