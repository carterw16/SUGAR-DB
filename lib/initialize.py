"""Initialize distribution grid device models.

Author(s): Amrit Pandey, Naeem Turner-Bandele
Created Date: 12-19-2016
Updated Date: 10-16-2020
Email: nturnerb@cmu.edu
Status: Development

"""
import pdb
from typing import Optional
from classes.Regulators import VariableRegulator
import numpy as np
import ipdb

def calc_infeas_i(Vr, Vi, P, Q):
    infeas_i_real = (P*Vr + Q*Vi)/(Vr**2 + Vi**2) if (Vr**2 + Vi**2) != 0 else 0
    infeas_i_imag = (P*Vi - Q*Vr)/(Vr**2 + Vi**2) if (Vr**2 + Vi**2) != 0 else 0
    return infeas_i_real, infeas_i_imag

def calc_infeas_i_lambdas(Vr, Vi, Ir, Ii):
    # this is assuming the objective function is sum of square powers
    infeas_i_lambda_real = -2*Ir*(Vr**2 + Vi**2) 
    infeas_i_lambda_imag = -2*Ii*(Vr**2 + Vi**2)
    return infeas_i_lambda_real, infeas_i_lambda_imag

def calc_infeas_GB(Vr, Vi, P, Q, real):
    V_mag = (Vr**2 + Vi**2)
    if real:
        if np.any(V_mag == 0):
            infeas_GB = np.zeros_like(V_mag)
            for index in range(0, len(V_mag)):
                infeas_GB[index] = P[index]/V_mag[index]
        else:
            infeas_GB = np.divide(P,(Vr**2 + Vi**2))
    else:
        if np.any(V_mag == 0):
            infeas_GB = np.zeros_like(V_mag)
            for index in range(0, len(V_mag)):
                infeas_GB[index] = -Q[index]/V_mag[index]
        else:
            infeas_GB = -np.divide(Q,(Vr**2 + Vi**2))

    if np.any(abs(infeas_GB) > 1e20):
        infeas_GB[np.where(abs(infeas_GB) > 1e20)[0]] = 1e-15
    
    return infeas_GB

def calc_infeas_GB_lambdas(Vr, Vi, G, B):
    ipdb.set_trace()
    infeas_GB_lambda_real = 2*(Vr**2 + Vi**2)*(B*Vi - G*Vr)
    infeas_GB_lambda_imag = -2*(Vr**2 + Vi**2)*(B*Vr + G*Vi)

def initialize(node, node_key, regulator, ibdg, battery_sources, op_bounds, Vinit, lf,
               stamp_dual=False, obj: Optional[int] = 2, source_type = 'current', stamp_slack_bus = False, stamp_tplx_bus_flag = True, stamp_neutral_infeas_source = False):
    default_dual_eq_var_value = 1e-4
    default_dual_ineq_var_value = 1e-4

    for ele in range(len(node)):
        infeas_init_complex_power = complex(1e-8*node[ele].Vnom, 1e-8*node[ele].Vnom)
        if not node[ele].isTriplex:
            node_type = node[ele].bustype
            # Initial Conditions for i) Node ABCN voltages ii) Generator Vars iii) Xtra xfmr nodes can have any value
            if node[ele].phases & 0x1 == 1:
                Vinit[node[ele].nodeA_Vr] = node[ele].Va.real
                Vinit[node[ele].nodeA_Vi] = node[ele].Va.imag if node[ele].Va.imag != 0 else 1e-10
            else:
                try:
                    Vinit[node[ele].nodeA_Vr] = 0
                    Vinit[node[ele].nodeA_Vi] = 0
                except:
                    print("No Phase A on bus ", node[ele].name) 
                    continue
            
            if node[ele].phases & 0x2 == 2:
                Vinit[node[ele].nodeB_Vr] = node[ele].Vb.real
                Vinit[node[ele].nodeB_Vi] = node[ele].Vb.imag
            else:
                try:
                    Vinit[node[ele].nodeB_Vr] = 0
                    Vinit[node[ele].nodeB_Vi] = 0
                except:
                    print("No Phase B on bus ", node[ele].name)
                    continue 
            
            if node[ele].phases & 0x4 == 4:
                Vinit[node[ele].nodeC_Vr] = node[ele].Vc.real
                Vinit[node[ele].nodeC_Vi] = node[ele].Vc.imag
            else:
                try:
                    Vinit[node[ele].nodeC_Vr] = 0
                    Vinit[node[ele].nodeC_Vi] = 0
                except:
                    print("No Phase C on bus ", node[ele].name)
                    continue
            
            try:
                Vinit[node[ele].nodeN_Vr] = 0.0
                Vinit[node[ele].nodeN_Vi] = 0.0
            except:
                print("No neutral node on ", node[ele].name)  
                continue           
            
            if stamp_dual == True and (node[ele].bustype != 3 or stamp_slack_bus == True):
                if node[ele].phases & 0x1 == 1:
                    Vinit[node[ele].nodeA_dual_eq_var_r] = default_dual_eq_var_value
                    Vinit[node[ele].nodeA_dual_eq_var_i] = default_dual_eq_var_value
                else:
                    try:
                        Vinit[node[ele].nodeA_dual_eq_var_r] = 0
                        Vinit[node[ele].nodeA_dual_eq_var_i] = 0
                    except:
                        continue
                
                if node[ele].phases & 0x2 == 2:
                    Vinit[node[ele].nodeB_dual_eq_var_r] = default_dual_eq_var_value
                    Vinit[node[ele].nodeB_dual_eq_var_i] = default_dual_eq_var_value
                else:
                    try:
                        Vinit[node[ele].nodeB_dual_eq_var_r] = 0
                        Vinit[node[ele].nodeB_dual_eq_var_i] = 0
                    except:
                        continue
                
                if node[ele].phases & 0x4 == 4:
                    Vinit[node[ele].nodeC_dual_eq_var_r] = default_dual_eq_var_value
                    Vinit[node[ele].nodeC_dual_eq_var_i] = default_dual_eq_var_value
                else:
                    try:
                        Vinit[node[ele].nodeC_dual_eq_var_r] = 0
                        Vinit[node[ele].nodeC_dual_eq_var_i] = 0
                    except:
                        continue
                
                try:                       
                    Vinit[node[ele].nodeN_dual_eq_var_r] = default_dual_eq_var_value
                    Vinit[node[ele].nodeN_dual_eq_var_i] = default_dual_eq_var_value
                except:
                    continue

                if obj == 'L2': 
                    if source_type == 'current':
                        P = infeas_init_complex_power.real
                        Q = infeas_init_complex_power.imag

                        if node[ele].phases & 0x1 == 1:
                            Vinit[node[ele].nodeA_if_r], Vinit[node[ele].nodeA_if_i] = calc_infeas_i(Vinit[node[ele].nodeA_Vr], Vinit[node[ele].nodeA_Vi], P, Q)
                        else:
                            try:
                                Vinit[node[ele].nodeA_if_r], Vinit[node[ele].nodeA_if_i] = [0, 0]
                            except:
                                continue
                        if node[ele].phases & 0x2 == 2:
                            Vinit[node[ele].nodeB_if_r], Vinit[node[ele].nodeB_if_i] = calc_infeas_i(Vinit[node[ele].nodeB_Vr], Vinit[node[ele].nodeB_Vi], P, Q)
                        else:
                            try:
                                Vinit[node[ele].nodeB_if_r], Vinit[node[ele].nodeB_if_i] = [0, 0]
                            except:
                                continue
                        if node[ele].phases & 0x4 == 4:
                            Vinit[node[ele].nodeC_if_r], Vinit[node[ele].nodeC_if_i] = calc_infeas_i(Vinit[node[ele].nodeC_Vr], Vinit[node[ele].nodeC_Vi], P, Q)
                        else:
                            try:
                                Vinit[node[ele].nodeC_if_r], Vinit[node[ele].nodeC_if_i] = [0,0]
                            except:
                                continue
                        if stamp_neutral_infeas_source:
                            Vinit[node[ele].nodeN_if_r], Vinit[node[ele].nodeN_if_i] = [0, 0]
                        else:
                            continue
                    
                    elif source_type == 'PQ' or source_type == 'Q':
                        try:
                            if node[ele].phases & 0x1 == 1:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeA_P] = infeas_init_complex_power.real
                                Vinit[node[ele].nodeA_Q] = infeas_init_complex_power.imag
                            else:
                                try:
                                    if source_type == 'PQ':
                                        Vinit[node[ele].nodeA_P] = 0
                                    Vinit[node[ele].nodeA_Q] = 0
                                except:
                                    continue
                                
                            if node[ele].phases & 0x2 == 2:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeB_P] = infeas_init_complex_power.real
                                Vinit[node[ele].nodeB_Q] = infeas_init_complex_power.imag
                            else:
                                try:
                                    if source_type == 'PQ':
                                        Vinit[node[ele].nodeB_P] = 0
                                    Vinit[node[ele].nodeB_Q] = 0
                                except:
                                    continue
                                
                            if node[ele].phases & 0x4 == 4:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeC_P] = infeas_init_complex_power.real
                                Vinit[node[ele].nodeC_Q] = infeas_init_complex_power.imag
                            else:
                                try:
                                    if source_type == 'PQ':
                                        Vinit[node[ele].nodeC_P] = 0
                                    Vinit[node[ele].nodeC_Q] = 0
                                except:
                                    continue
                            
                            if stamp_neutral_infeas_source:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeN_P] = infeas_init_complex_power.real
                                Vinit[node[ele].nodeN_Q] = infeas_init_complex_power.imag
    
                        except:
                            if source_type == 'PQ':
                                Vinit[node[ele].nodeA_P] = infeas_init_complex_power.real
                            Vinit[node[ele].Q_nodes] = infeas_init_complex_power.imag
                    
                    elif source_type == 'GB' or source_type == 'B':
                        P = infeas_init_complex_power.real
                        Q = infeas_init_complex_power.imag

                        Vr = []
                        Vi = []
                        if node[ele].phases & 0x1 == 1:
                            Vr.append(Vinit[node[ele].nodeA_Vr])
                            Vi.append(Vinit[node[ele].nodeA_Vi])
                        if node[ele].phases & 0x2 == 2:
                            Vr.append(Vinit[node[ele].nodeB_Vr])
                            Vi.append(Vinit[node[ele].nodeB_Vi])
                        if node[ele].phases & 0x4 == 4:
                            Vr.append(Vinit[node[ele].nodeC_Vr])
                            Vi.append(Vinit[node[ele].nodeC_Vi])
                        if stamp_neutral_infeas_source:
                            Vr.append(Vinit[node[ele].nodeN_Vr])
                            Vi.append(Vinit[node[ele].nodeN_Vi])
                        Vr = np.array(Vr)
                        Vi = np.array(Vi)
                        
                        if source_type == 'GB':
                            try:
                                Vinit[node[ele].G_nodes] = calc_infeas_GB(Vr, Vi, P, Q, True)
                            except:
                                "Error intializing G infeas source nodes (intialize.py)"
                        try:        
                            Vinit[node[ele].B_nodes] = calc_infeas_GB(Vr, Vi, P, Q, False)
                        except:
                            "Error initializing B infeas source nodes (intialize.py)"
                if obj == 'L1' and source_type == 'current':
                    P = infeas_init_complex_power.real
                    Q = infeas_init_complex_power.imag

                    if node[ele].phases & 0x1 == 1:
                        Ir_temp, Ii_temp = calc_infeas_i(Vinit[node[ele].nodeA_Vr], Vinit[node[ele].nodeA_Vi], P, Q)
                        Ir = Ir_temp if Ir_temp >= 0 else 1e-2
                        Ii = Ii_temp if Ii_temp >= 0 else 1e-2

                        Vinit[node[ele].nodeA_if_r_plus] = Ir
                        Vinit[node[ele].nodeA_if_i_plus] = Ii
                        Vinit[node[ele].nodeA_if_r_minus] = Ir
                        Vinit[node[ele].nodeA_if_i_minus] = Ii

                        Vinit[node[ele].nodeA_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeA_dual_ineq_i_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeA_dual_ineq_r_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeA_dual_ineq_i_minus] = default_dual_ineq_var_value

                    else:
                        # In development - initializing phases where there are no connections to 0 (or not stamping at all to those phases)
                        try:
                            Vinit[node[ele].nodeA_if_r_plus] = 1e-2
                            Vinit[node[ele].nodeA_if_i_plus] = 1e-2
                            Vinit[node[ele].nodeA_if_r_minus] = 1e-2
                            Vinit[node[ele].nodeA_if_i_minus] = 1e-2

                            Vinit[node[ele].nodeA_dual_ineq_r_plus] =  1e-4
                            Vinit[node[ele].nodeA_dual_ineq_i_plus] =  1e-4
                            Vinit[node[ele].nodeA_dual_ineq_r_minus] =  1e-4
                            Vinit[node[ele].nodeA_dual_ineq_i_minus] =  1e-4
                        except:
                            continue

                    
                    if node[ele].phases & 0x2 == 2:
                        Ir_temp, Ii_temp = calc_infeas_i(Vinit[node[ele].nodeB_Vr], Vinit[node[ele].nodeB_Vi], P, Q)
                        Ir = Ir_temp if Ir_temp >= 0 else 1e-2
                        Ii = Ii_temp if Ii_temp >= 0 else 1e-2

                        Vinit[node[ele].nodeB_if_r_plus] = Ir
                        Vinit[node[ele].nodeB_if_i_plus] = Ii
                        Vinit[node[ele].nodeB_if_r_minus] = Ir
                        Vinit[node[ele].nodeB_if_i_minus] = Ii

                        Vinit[node[ele].nodeB_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeB_dual_ineq_i_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeB_dual_ineq_r_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeB_dual_ineq_i_minus] = default_dual_ineq_var_value

                    else:
                        # In development - initializing phases where there are no connections to 0 (or not stamping at all to those phases)
                        try:
                            Vinit[node[ele].nodeB_if_r_plus] = 1e-2
                            Vinit[node[ele].nodeB_if_i_plus] = 1e-2
                            Vinit[node[ele].nodeB_if_r_minus] = 1e-2
                            Vinit[node[ele].nodeB_if_i_minus] = 1e-2

                            Vinit[node[ele].nodeB_dual_ineq_r_plus] =  1e-4
                            Vinit[node[ele].nodeB_dual_ineq_i_plus] =  1e-4
                            Vinit[node[ele].nodeB_dual_ineq_r_minus] =  1e-4
                            Vinit[node[ele].nodeB_dual_ineq_i_minus] =  1e-4
                        except:
                            continue
                    
                    if node[ele].phases & 0x4 == 4:
                        Ir_temp, Ii_temp = calc_infeas_i(Vinit[node[ele].nodeB_Vr], Vinit[node[ele].nodeB_Vi], P, Q)
                        Ir = Ir_temp if Ir_temp >= 0 else 1e-2
                        Ii = Ii_temp if Ii_temp >= 0 else 1e-2

                        Vinit[node[ele].nodeC_if_r_plus] = Ir
                        Vinit[node[ele].nodeC_if_i_plus] = Ii
                        Vinit[node[ele].nodeC_if_r_minus] = Ir
                        Vinit[node[ele].nodeC_if_i_minus] = Ii

                        Vinit[node[ele].nodeC_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeC_dual_ineq_i_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeC_dual_ineq_r_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeC_dual_ineq_i_minus] = default_dual_ineq_var_value 

                    else:
                        # In development - initializing phases where there are no connections to 0 (or not stamping at all to those phases)
                        try:
                            Vinit[node[ele].nodeC_if_r_plus] = 1e-2
                            Vinit[node[ele].nodeC_if_i_plus] = 1e-2
                            Vinit[node[ele].nodeC_if_r_minus] = 1e-2
                            Vinit[node[ele].nodeC_if_i_minus] = 1e-2

                            Vinit[node[ele].nodeC_dual_ineq_r_plus] =  1e-4
                            Vinit[node[ele].nodeC_dual_ineq_i_plus] =  1e-4
                            Vinit[node[ele].nodeC_dual_ineq_r_minus] =  1e-4
                            Vinit[node[ele].nodeC_dual_ineq_i_minus] =  1e-4 
                        except:
                            continue

                    if stamp_neutral_infeas_source:
                        Vinit[node[ele].nodeN_if_r_plus] = 1e-2
                        Vinit[node[ele].nodeN_if_i_plus] = 1e-2
                        Vinit[node[ele].nodeN_if_r_minus] = 1e-2
                        Vinit[node[ele].nodeN_if_i_minus] = 1e-2

                        Vinit[node[ele].nodeN_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeN_dual_ineq_i_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeN_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].nodeN_dual_ineq_i_minus] = default_dual_ineq_var_value


                if obj == 'L1' and (source_type == 'PQ' or source_type == 'Q'):
                    try:
                        if node[ele].phases & 0x1 == 1:
                            if source_type == 'PQ':
                                Vinit[node[ele].nodeA_p_plus] = 0.001 #(debug) #infeas_init_complex_power.real
                                Vinit[node[ele].nodeA_p_minus] = 1500000 # (debug) #infeas_init_complex_power.real/2
                                Vinit[node[ele].nodeA_dual_ineq_var_p_plus] = default_dual_ineq_var_value
                                Vinit[node[ele].nodeA_dual_ineq_var_p_minus] = default_dual_ineq_var_value/2

                                # (debug) upper bounds
                                #Vinit[node[ele].nodeA_dual_ineq_r_plus_upper] = default_dual_ineq_var_value
                                #Vinit[node[ele].nodeA_dual_ineq_r_minus_upper] = default_dual_ineq_var_value


                            Vinit[node[ele].nodeA_q_plus] = infeas_init_complex_power.imag
                            Vinit[node[ele].nodeA_q_minus] = infeas_init_complex_power.imag/2
                            Vinit[node[ele].nodeA_dual_ineq_var_q_plus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeA_dual_ineq_var_q_minus] = default_dual_ineq_var_value/2
                        else:
                            try:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeA_p_plus] = 0
                                    Vinit[node[ele].nodeA_p_minus] = 0
                                    Vinit[node[ele].nodeA_dual_ineq_var_p_plus] = 0
                                    Vinit[node[ele].nodeA_dual_ineq_var_p_minus] = 0
                                Vinit[node[ele].nodeA_q_plus] = 0
                                Vinit[node[ele].nodeA_q_minus] = 0
                                Vinit[node[ele].nodeA_dual_ineq_var_q_plus] = 0
                                Vinit[node[ele].nodeA_dual_ineq_var_q_minus] = 0
                            except:
                                continue
                                
                        if node[ele].phases & 0x2 == 2:
                            if source_type == 'PQ':
                                Vinit[node[ele].nodeB_p_plus] = 0.001 #(debug) #infeas_init_complex_power.real
                                Vinit[node[ele].nodeB_p_minus] = 1500000 #(debug) #infeas_init_complex_power.real/2
                                Vinit[node[ele].nodeB_dual_ineq_var_p_plus] = default_dual_ineq_var_value
                                Vinit[node[ele].nodeB_dual_ineq_var_p_minus] = default_dual_ineq_var_value/2

                                # (debug) upper bounds
                                #Vinit[node[ele].nodeB_dual_ineq_r_plus_upper] = default_dual_ineq_var_value
                                #Vinit[node[ele].nodeB_dual_ineq_r_minus_upper] = default_dual_ineq_var_value

                            Vinit[node[ele].nodeB_q_plus] = infeas_init_complex_power.imag
                            Vinit[node[ele].nodeB_q_minus] = infeas_init_complex_power.imag/2
                            Vinit[node[ele].nodeB_dual_ineq_var_q_plus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeB_dual_ineq_var_q_minus] = default_dual_ineq_var_value/2
                        else:
                            try:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeB_p_plus] = 0
                                    Vinit[node[ele].nodeB_p_minus] = 0
                                    Vinit[node[ele].nodeB_dual_ineq_var_p_plus] = 0
                                    Vinit[node[ele].nodeB_dual_ineq_var_p_minus] = 0
                                Vinit[node[ele].nodeB_q_plus] = 0
                                Vinit[node[ele].nodeB_q_minus] = 0
                                Vinit[node[ele].nodeB_dual_ineq_var_q_plus] = 0
                                Vinit[node[ele].nodeB_dual_ineq_var_q_minus] = 0
                            except:
                                continue
                                
                        if node[ele].phases & 0x4 == 4:
                            if source_type == 'PQ':
                                Vinit[node[ele].nodeC_p_plus] = 0.001 #(debug) #infeas_init_complex_power.real
                                Vinit[node[ele].nodeC_p_minus] = 1500000 #(debug) #infeas_init_complex_power.real/2
                                Vinit[node[ele].nodeC_dual_ineq_var_p_plus] = default_dual_ineq_var_value
                                Vinit[node[ele].nodeC_dual_ineq_var_p_minus] = default_dual_ineq_var_value/2

                                # (debug) upper bounds
                                #Vinit[node[ele].nodeC_dual_ineq_r_plus_upper] = default_dual_ineq_var_value
                                #Vinit[node[ele].nodeC_dual_ineq_r_minus_upper] = default_dual_ineq_var_value

                            Vinit[node[ele].nodeC_q_plus] = infeas_init_complex_power.imag
                            Vinit[node[ele].nodeC_q_minus] = infeas_init_complex_power.imag/2
                            Vinit[node[ele].nodeC_dual_ineq_var_q_plus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeC_dual_ineq_var_q_minus] = default_dual_ineq_var_value/2
                        else:
                            try:
                                if source_type == 'PQ':
                                    Vinit[node[ele].nodeC_p_plus] = 0
                                    Vinit[node[ele].nodeC_p_minus] = 0
                                    Vinit[node[ele].nodeC_dual_ineq_var_p_plus] = 0
                                    Vinit[node[ele].nodeC_dual_ineq_var_p_minus] = 0
                                Vinit[node[ele].nodeC_q_plus] = 0
                                Vinit[node[ele].nodeC_q_minus] = 0
                                Vinit[node[ele].nodeC_dual_ineq_var_q_plus] = 0
                                Vinit[node[ele].nodeC_dual_ineq_var_q_minus] = 0
                            except:
                                continue
                                
                    except:               
                        Vinit[node[ele].Q_plus_nodes] = infeas_init_complex_power.imag
                        Vinit[node[ele].Q_minus_nodes] = infeas_init_complex_power.imag
                        Vinit[node[ele].dual_ineq_i_plus_nodes] = default_dual_ineq_var_value
                        Vinit[node[ele].dual_ineq_i_minus_nodes] = default_dual_ineq_var_value

                        if source_type == 'PQ':
                            Vinit[node[ele].P_plus_nodes] = infeas_init_complex_power.real
                            Vinit[node[ele].P_minus_nodes] = infeas_init_complex_power.real
                            Vinit[node[ele].dual_ineq_r_plus_nodes] = default_dual_ineq_var_value
                            Vinit[node[ele].dual_ineq_r_minus_nodes] = default_dual_ineq_var_value                     
                    
            if node_type == 2:
                pass
        else:
            node_type = node[ele].bustype
            Vinit[node[ele].node1_Vr] = node[ele].V1.real
            Vinit[node[ele].node1_Vi] = -node[ele].V1.imag
            Vinit[node[ele].node2_Vr] = node[ele].V2.real
            Vinit[node[ele].node2_Vi] = -node[ele].V2.imag
            Vinit[node[ele].nodeN_Vr] = 0.0
            Vinit[node[ele].nodeN_Vi] = 0.0

            if stamp_dual == True and obj == 'L2' and (node[ele].bustype != 3 or stamp_slack_bus == True):
                Vinit[node[ele].node1_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].node1_dual_eq_var_i] = default_dual_eq_var_value
                Vinit[node[ele].node2_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].node2_dual_eq_var_i] = default_dual_eq_var_value
                Vinit[node[ele].nodeN_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].nodeN_dual_eq_var_i] = default_dual_eq_var_value
                
                if stamp_tplx_bus_flag:
                    if source_type == 'current':
                        P = infeas_init_complex_power.real
                        Q = infeas_init_complex_power.imag

                        Ir1_temp, Ii1_temp = calc_infeas_i(Vinit[node[ele].node1_Vr], Vinit[node[ele].node1_Vi], P, Q)
                        Ir2_temp, Ii2_temp = calc_infeas_i(Vinit[node[ele].node2_Vr], Vinit[node[ele].node2_Vi], P, Q)

                        Vinit[node[ele].node1_if_r], Vinit[node[ele].node1_if_i] = calc_infeas_i(Vinit[node[ele].node1_Vr], Vinit[node[ele].node1_Vi], P, Q)
                        Vinit[node[ele].node2_if_r], Vinit[node[ele].node2_if_i] = calc_infeas_i(Vinit[node[ele].node2_Vr], Vinit[node[ele].node2_Vi], P, Q)
                        if stamp_neutral_infeas_source:
                            Vinit[node[ele].nodeN_if_r],Vinit[node[ele].nodeN_if_i] = calc_infeas_i(Vinit[node[ele].nodeN_Vr], Vinit[node[ele].nodeN_Vi], P, Q)
                    
                    elif source_type == 'PQ' or source_type == 'Q':
                        if source_type == 'PQ':
                            Vinit[node[ele].P_nodes] = infeas_init_complex_power.real
                        Vinit[node[ele].Q_nodes] = infeas_init_complex_power.imag
                    
                    elif source_type == 'GB' or source_type == 'B':
                        P = infeas_init_complex_power.real
                        Q = infeas_init_complex_power.imag
                        if stamp_neutral_infeas_source:
                            Vr = np.array([Vinit[node[ele].node1_Vr], Vinit[node[ele].node2_Vr], Vinit[node[ele].nodeN_Vr]])
                            Vi = np.array([Vinit[node[ele].node1_Vi], Vinit[node[ele].node2_Vi], Vinit[node[ele].nodeN_Vi]])
                        else:
                            Vr = np.array([Vinit[node[ele].node1_Vr], Vinit[node[ele].node2_Vr]])
                            Vi = np.array([Vinit[node[ele].node1_Vi], Vinit[node[ele].node2_Vi]])
                           
                        if source_type == 'GB':
                            Vinit[node[ele].G_nodes] = calc_infeas_GB(Vr, Vi, P, Q, True)
                        Vinit[node[ele].B_nodes] = calc_infeas_GB(Vr, Vi, P, Q, False)

            if stamp_dual == True and obj == 'L1' and (node[ele].bustype != 3 or stamp_slack_bus == True):
                Vinit[node[ele].node1_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].node1_dual_eq_var_i] = default_dual_eq_var_value
                Vinit[node[ele].node2_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].node2_dual_eq_var_i] = default_dual_eq_var_value
                Vinit[node[ele].nodeN_dual_eq_var_r] = default_dual_eq_var_value
                Vinit[node[ele].nodeN_dual_eq_var_i] = default_dual_eq_var_value

                if stamp_tplx_bus_flag:
                    P = infeas_init_complex_power.real
                    Q = infeas_init_complex_power.imag
                    if source_type == 'current':
                        Ir1_temp, Ii1_temp = calc_infeas_i(Vinit[node[ele].node1_Vr], Vinit[node[ele].node1_Vi], P, Q)
                        Ir2_temp, Ii2_temp = calc_infeas_i(Vinit[node[ele].node2_Vr], Vinit[node[ele].node2_Vi], P, Q)

                        Ir1 = Ir1_temp if Ir1_temp >= 0 else 1e-2
                        Ii1 = Ii1_temp if Ii1_temp >= 0 else 1e-2
                        Ir2 = Ir2_temp if Ir2_temp >= 0 else 1e-2
                        Ii2 = Ii2_temp if Ii2_temp >= 0 else 1e-2
                        
                        Vinit[node[ele].node1_if_r_plus] = Ir1
                        Vinit[node[ele].node1_if_i_plus] = Ii1
                        Vinit[node[ele].node2_if_r_plus] = Ir2
                        Vinit[node[ele].node2_if_i_plus] = Ii2

                        Vinit[node[ele].node1_if_r_minus] = Ir1
                        Vinit[node[ele].node1_if_i_minus] = Ii1
                        Vinit[node[ele].node2_if_r_minus] = Ir2
                        Vinit[node[ele].node2_if_i_minus] = Ii2

                        Vinit[node[ele].node1_dual_ineq_r_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].node1_dual_ineq_i_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].node2_dual_ineq_r_minus] = default_dual_ineq_var_value
                        Vinit[node[ele].node2_dual_ineq_i_minus] = default_dual_ineq_var_value

                        Vinit[node[ele].node1_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].node1_dual_ineq_i_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].node2_dual_ineq_r_plus] = default_dual_ineq_var_value
                        Vinit[node[ele].node2_dual_ineq_i_plus] = default_dual_ineq_var_value

                        if stamp_neutral_infeas_source:
                            Irn_temp, Iin_temp = calc_infeas_i(Vinit[node[ele].nodeN_Vr], Vinit[node[ele].nodeN_Vi], P, Q)
                            Irn = Irn_temp if Irn_temp >= 0 else 1e-2
                            Iin = Ii1_temp if Iin_temp >= 0 else 1e-2
                            Vinit[node[ele].nodeN_if_r_plus] = Irn
                            Vinit[node[ele].nodeN_if_i_plus] = Iin
                            Vinit[node[ele].nodeN_if_r_minus] = Irn
                            Vinit[node[ele].nodeN_if_i_minus] = Iin

                            Vinit[node[ele].nodeN_dual_ineq_r_minus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeN_dual_ineq_i_minus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeN_dual_ineq_r_plus] = default_dual_ineq_var_value
                            Vinit[node[ele].nodeN_dual_ineq_i_plus] = default_dual_ineq_var_value
                    
                    if source_type == 'PQ' or source_type == 'Q':
                        if source_type == 'PQ':
                            Vinit[node[ele].P_plus_nodes] = infeas_init_complex_power.real
                            Vinit[node[ele].P_minus_nodes] =infeas_init_complex_power.real/2
                            Vinit[node[ele].dual_ineq_r_plus_nodes] = default_dual_ineq_var_value
                            Vinit[node[ele].dual_ineq_r_minus_nodes] = default_dual_ineq_var_value/2

                        Vinit[node[ele].Q_plus_nodes] = infeas_init_complex_power.imag
                        Vinit[node[ele].Q_minus_nodes] =infeas_init_complex_power.imag/2
                        Vinit[node[ele].dual_ineq_i_plus_nodes] = default_dual_ineq_var_value
                        Vinit[node[ele].dual_ineq_i_minus_nodes] = default_dual_ineq_var_value/2


    if battery_sources:
        for ele in battery_sources:
            Vinit = ele.initialize(Vinit)

    if regulator:
        for ele in regulator:
            if isinstance(ele, VariableRegulator):
                Vinit = ele.initialize_regulator(Vinit)

    if ibdg:
        for ele in ibdg:
            Vinit = ele.initialize_ibdg(node, node_key, Vinit, lf)
    
    if op_bounds:
        # in the future op_bounds could contain other limits like line/transformer
        # flow limits
        from classes.OperationalLimits.VoltageLimits import VoltageLimits
        from classes.Nodes import Nodes
        for ele in op_bounds:
            ele.init_bounding_variables(Vinit)
            if len(Nodes.vmag2_index) != len(VoltageLimits.pu_max2):
                test = 2


    return Vinit
