#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 14:24:32 2021

@author: emfoster
"""
import pandas as pd
import numpy as np
from classes.Nodes import Nodes


def calc_PQ_current(P, Q, Vr, Vi):
    # Purpose: this function is used to find the current from the PQ slack source
    V_mag = (np.square(Vr) + np.square(Vi))
    if np.any(V_mag == 0):
        Ir = np.zeros_like(V_mag)
        Ii = np.zeros_like(V_mag)
        for index in range(0, len(V_mag)):
            Ir[index] = (P[index]*Vr[index] + Q[index]*Vi[index])/V_mag[index] if V_mag[index] != 0 else 0
            Ii[index] = (P[index]*Vi[index] - Q[index]*Vr[index])/V_mag[index] if V_mag[index] != 0 else 0
    else:
        # Ir = (P*Vr + Q*Vi)/(Vr**2 + Vi**2)
        Ir = np.divide((np.multiply(P,Vr) + np.multiply(Q,Vi)), V_mag) 
        # Ii = (P*Vi - Q*Vr)/(Vr**2 + Vi**2)
        Ii = np.divide((np.multiply(P,Vi) - np.multiply(Q,Vr)), V_mag) 

    return Ir, Ii

def calc_GB_current(G, B, Vr, Vi):
    # Purpose: this function is used to find the current from the GB slack source
    # Ir = G*Vr - B*Vi
    Ir = np.multiply(G,Vr) - np.multiply(B,Vi)
    # Ii = G*Vi + B*Vr
    Ii = np.multiply(G,Vi) + np.multiply(B,Vr)

    return Ir, Ii

def dual_analysis(V_out, node, obj, stamp_slack, stamp_tplx, source_type):
    find_magnitude = True
   
    num_of_tplx = 0
    tplx_index = []
    num_of_ABC = 0
    ABC_index = []
    for i in range(0, len(node)):
        if node[i].isTriplex == False:
            num_of_ABC += 1
            ABC_index.append(i)
        elif node[i].isTriplex and stamp_tplx:
            num_of_tplx += 1
            tplx_index.append(i)
            

    ## IDENTIFYING INFEASIBILITY CURRENTS
    ABC_column_names = ['Node Name', 
                            'Infeas Ir (A)', 'Infeas Ii (A)', 
                            'Infeas Ir (B)', 'Infeas Ii (B)', 
                            'Infeas Ir (C)', 'Infeas Ii (C)']
    if source_type == 'Q':
        ABC_column_names.extend(['Infeas Q (A)', 'Infeas Q (B)', 'Infeas Q (C)'])
    elif source_type == 'PQ':
        ABC_column_names.extend(['Infeas P (A)', 'Infeas P (B)', 'Infeas P (C)'])
        ABC_column_names.extend(['Infeas Q (A)', 'Infeas Q (B)', 'Infeas Q (C)'])
    elif source_type == 'B':
        ABC_column_names.extend(['Infeas B (A)', 'Infeas B (B)', 'Infeas B (C)'])
    elif source_type == 'GB':
        ABC_column_names.extend(['Infeas G (A)', 'Infeas G (B)', 'Infeas G (C)'])
        ABC_column_names.extend(['Infeas B (A)', 'Infeas B (B)', 'Infeas B (C)'])

    ABC_dual_info = pd.DataFrame(columns=ABC_column_names, index = range(0, num_of_ABC))
    
    if stamp_tplx:
        tplx_column_names = ['Node Name','if1,r Value', 'if1,i Value', 'if2,r Value', 'if2,i Value', 'ifN,r Value', 'ifN,i Value']
        tplx_dual_info = pd.DataFrame(columns=tplx_column_names, index = range(0, num_of_tplx))

    max_per_node = []
    max_power = []

    print('Assuming sources are on all three phases (dual_analysis)')
    for i in range(0, num_of_ABC):
        temp_node_info = node[ABC_index[i]]
        if (temp_node_info.bustype!= 3) or (temp_node_info.bustype== 3 and stamp_slack):
            if obj == 'L2':
                if source_type == 'current':
                    Ir = [V_out[temp_node_info.nodeA_if_r], V_out[temp_node_info.nodeB_if_r], V_out[temp_node_info.nodeC_if_r]]
                    Ii = [V_out[temp_node_info.nodeA_if_i], V_out[temp_node_info.nodeB_if_i], V_out[temp_node_info.nodeC_if_i]]
                    
                elif source_type == 'PQ' or source_type == 'Q':
                    Q = V_out[temp_node_info.Q_nodes]
                    if source_type == 'PQ':
                        P = V_out[temp_node_info.P_nodes]
                    else:
                        P = np.zeros_like(Q)
                    Vr = V_out[temp_node_info.Vr_nodes]
                    if len(Vr) > len(Q):
                        Vr = V_out[temp_node_info.Vr_nodes][0:3]
                    Vi = V_out[temp_node_info.Vi_nodes][0:3]
                    if len(Vi) > len(Q):
                        Vi = V_out[temp_node_info.Vi_nodes][0:3]

                    Ir, Ii = calc_PQ_current(P, Q, Vr, Vi)
                
                elif source_type == 'GB' or source_type == 'B':
                    B = V_out[temp_node_info.B_nodes]
                    if source_type == 'GB':
                        G = V_out[temp_node_info.G_nodes]
                    else:
                        G = np.zeros_like(B)
                    Vr = V_out[temp_node_info.Vr_nodes]
                    if len(Vr) > len(B):
                        Vr = V_out[temp_node_info.Vr_nodes][0:3]
                    Vi = V_out[temp_node_info.Vi_nodes][0:3]
                    if len(Vi) > len(B):
                        Vi = V_out[temp_node_info.Vi_nodes][0:3]

                    Ir, Ii = calc_PQ_current(G, B, Vr, Vi)

            elif obj == 'L1' or obj == 1:
                if source_type == 'current':
                    Ir = [V_out[temp_node_info.nodeA_if_r_plus] - V_out[temp_node_info.nodeA_if_r_minus],
                          V_out[temp_node_info.nodeB_if_r_plus] - V_out[temp_node_info.nodeB_if_r_minus],
                          V_out[temp_node_info.nodeC_if_r_plus] - V_out[temp_node_info.nodeC_if_r_minus]]
                    Ii = [V_out[temp_node_info.nodeA_if_i_plus] - V_out[temp_node_info.nodeA_if_i_minus],
                          V_out[temp_node_info.nodeB_if_i_plus] - V_out[temp_node_info.nodeB_if_i_minus],
                          V_out[temp_node_info.nodeC_if_i_plus] - V_out[temp_node_info.nodeC_if_i_minus]]
                elif source_type == 'PQ' or source_type == 'Q':
                    Q = V_out[temp_node_info.Q_plus_nodes] - V_out[temp_node_info.Q_minus_nodes]
                    if source_type == 'PQ':
                        P = V_out[temp_node_info.P_plus_nodes] - V_out[temp_node_info.P_minus_nodes]
                    else:
                        P = np.zeros_like(Q)
                    Vr = V_out[temp_node_info.Vr_nodes]
                    if len(Vr) > len(Q):
                        Vr = V_out[temp_node_info.Vr_nodes][0:3]
                    Vi = V_out[temp_node_info.Vi_nodes][0:3]
                    if len(Vi) > len(Q):
                        Vi = V_out[temp_node_info.Vi_nodes][0:3]

                    Ir, Ii = calc_PQ_current(P, Q, Vr, Vi)
                   
                elif source_type == 'GB' or source_type == 'B':
                    print("Dual analysis not built out for G/B sources with L1")
                
            ABC_dual_info.loc[i, 'Infeas Ir (A)'] = Ir[0]
            ABC_dual_info.loc[i, 'Infeas Ii (A)'] = Ii[0]
            ABC_dual_info.loc[i, 'Infeas Ir (B)'] = Ir[1]
            ABC_dual_info.loc[i, 'Infeas Ii (B)'] = Ii[1]
            ABC_dual_info.loc[i, 'Infeas Ir (C)'] = Ir[2]
            ABC_dual_info.loc[i, 'Infeas Ii (C)'] = Ii[2]
            ABC_dual_info.loc[i, 'Node Name'] = temp_node_info.name

            if find_magnitude:
                Imag_A = np.sqrt(Ir[0]** 2 + Ii[0]** 2)
                Imag_B = np.sqrt(Ir[1]** 2 + Ii[1]** 2)
                Imag_C = np.sqrt(Ir[2]** 2 + Ii[2]** 2)

                max_if = max(Imag_A, Imag_B, Imag_C)
                min_if = min(Imag_A, Imag_B, Imag_C)

                if abs(min_if) > max_if:
                    max_per_node.append(min_if)
                else:
                    max_per_node.append(max_if)
                
                Vr = V_out[temp_node_info.Vr_nodes]
                if len(Vr) > len(Ir):
                    Vr = V_out[temp_node_info.Vr_nodes][0:3]
                Vi = V_out[temp_node_info.Vi_nodes][0:3]
                if len(Vi) > len(Ir):
                    Vi = V_out[temp_node_info.Vi_nodes][0:3]
                
                Vmag_A = np.sqrt(Vr[0]**2 + Vi[0]**2)
                Vmag_B = np.sqrt(Vr[1]**2 + Vi[1]**2)
                Vmag_C = np.sqrt(Vr[2]**2 + Vi[2]**2)
                relative_power = [Vmag_A*Imag_A, Vmag_B*Imag_B, Vmag_C*Imag_C]
                max_power.append(max(relative_power))
                
    for i in range(0, num_of_tplx):
        temp_node_info = node[tplx_index[i]]

        if temp_node_info.isTriplex == 1 and obj == 2:
            tplx_dual_info['if1,r Value'][i] = V_out[temp_node_info.node1_if_r]
            tplx_dual_info['if1,i Value'][i] = V_out[temp_node_info.node1_if_i]
            tplx_dual_info['if2,r Value'][i] = V_out[temp_node_info.node2_if_r]
            tplx_dual_info['if2,i Value'][i] = V_out[temp_node_info.node2_if_i]
            tplx_dual_info['ifN,r Value'][i] = V_out[temp_node_info.nodeN_if_r]
            tplx_dual_info['ifN,i Value'][i] = V_out[temp_node_info.nodeN_if_i]
            tplx_dual_info['Node Name'][i] = temp_node_info.name
            
            if find_magnitude:
                mag_1 = np.sqrt(V_out[temp_node_info.node1_if_r] ** 2 + V_out[temp_node_info.node1_if_i] ** 2)
                mag_2 = np.sqrt(V_out[temp_node_info.node2_if_r] ** 2 + V_out[temp_node_info.node2_if_i] ** 2)
                mag_N = np.sqrt(V_out[temp_node_info.nodeN_if_r] ** 2 + V_out[temp_node_info.nodeN_if_i] ** 2)

                max_if = max(mag_1, mag_2, mag_N)
                min_if = min(mag_1, mag_2, mag_N)

                if abs(min_if) > max_if:
                    max_per_node.append(min_if)
                else:
                    max_per_node.append(max_if)

        elif temp_node_info.isTriplex == 1 and obj == 1:
            tplx_dual_info['if1,r Value'][i] = V_out[temp_node_info.node1_if_r_plus] + V_out[temp_node_info.node1_if_r_plus]
            tplx_dual_info['if1,i Value'][i] = V_out[temp_node_info.node1_if_i_minus] + V_out[temp_node_info.node1_if_i_minus]
            tplx_dual_info['if2,r Value'][i] = V_out[temp_node_info.node2_if_r_plus] + V_out[temp_node_info.node2_if_r_plus]
            tplx_dual_info['if2,i Value'][i] = V_out[temp_node_info.node2_if_i_minus] + V_out[temp_node_info.node1_if_i_minus]
            tplx_dual_info['ifN,r Value'][i] = V_out[temp_node_info.nodeN_if_r_plus] + V_out[temp_node_info.nodeN_if_r_plus]
            tplx_dual_info['ifN,i Value'][i] = V_out[temp_node_info.nodeN_if_i_minus] + V_out[temp_node_info.node1_if_i_minus]
            tplx_dual_info['Node Name'][i] = temp_node_info.name

            if find_magnitude:
                mag_1 = np.sqrt((V_out[temp_node_info.node1_if_r_plus] - V_out[temp_node_info.node1_if_r_minus]) ** 2 +
                                (V_out[temp_node_info.node1_if_i_plus] - V_out[temp_node_info.node1_if_i_minus]) ** 2)
                mag_2 = np.sqrt((V_out[temp_node_info.node2_if_r_plus] - V_out[temp_node_info.node2_if_r_minus]) ** 2 +
                                (V_out[temp_node_info.node2_if_i_plus] - V_out[temp_node_info.node2_if_i_minus]) ** 2)
                mag_N = np.sqrt((V_out[temp_node_info.nodeN_if_r_plus] - V_out[temp_node_info.nodeN_if_r_minus]) ** 2 +
                                (V_out[temp_node_info.nodeN_if_i_plus] - V_out[temp_node_info.nodeN_if_i_minus]) ** 2)

                max_if = max(mag_1, mag_2, mag_N)
                min_if = min(mag_1, mag_2, mag_N)

                if abs(min_if) > max_if:
                    max_per_node.append(min_if)
                else:
                    max_per_node.append(max_if)

    if stamp_tplx == False:
        tplx_dual_info = None

    max_per_node = np.array(max_per_node)
    print_currents = True
    if print_currents:
        print()
        print('Maximum absolute infeasibility current per node')
        print('Number of infeasibility currents > 1e3:', np.size(np.where(abs(max_per_node) > 1e3)[0]))
        print('Number of infeasibility currents > 1e1:', np.size(np.where(abs(max_per_node) > 1e1)[0]))
        print('Number of infeasibility currents > 1e-1:', np.size(np.where(abs(max_per_node) > 1e-1)[0]))
        print('Number of infeasibility currents > 1e-3:', np.size(np.where(abs(max_per_node) > 1e-3)[0]))
        print('Number of infeasibility currents > 1e-6:', np.size(np.where(abs(max_per_node) > 1e-6)[0]))
        print('Number of infeasibility currents > 1e-9:', np.size(np.where(abs(max_per_node) > 1e-9)[0]))

    return ABC_dual_info, tplx_dual_info
