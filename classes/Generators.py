"""Conventional grid power sources.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 04-10-2017
  Updated Date: 10-18-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Development

"""

from __future__ import division

import numpy as np

from classes.Elements import LinearElement
from classes.Nodes import Nodes
from lib.identify_voltage_type import identify_voltage_type


class Slack(LinearElement):

    @staticmethod
    def almostEqual(d1, d2, epsilon=10**-7):
        return abs(d2 - d1) < epsilon

    @staticmethod
    def _totuple(a):
        try:
            return tuple(i for i in a)
        except TypeError:
            return

    @staticmethod
    def _convertline2phase(Vab, Vbc, Vca):
        _Veq = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
        _VLL = np.array([Vab, Vbc, Vca])
        _Vabc = np.linalg.solve(_Veq, _VLL)
        return Slack._totuple(_Vabc)

    @staticmethod
    def sind(x):
        return np.sin(float(x) * np.pi / 180)

    @staticmethod
    def cosd(x):
        return np.cos(float(x) * np.pi / 180)

    def __init__(self,
                 ID,
                 name,
                 phases,
                 nominal_voltage,
                 phase_angle,
                 voltage_setpoint,
                 connect_type,
                 slack_parameters,
                 voltage_A=None,
                 voltage_B=None,
                 voltage_C=None,
                 stamp_dual=False):

        super(Slack, self).__init__()

        self.ID = ID
        self.name = name
        self.phases = int(phases)
        # check if voltage is LL or LN and adjust nom voltage if needed
        voltage_type, nominal_voltage = identify_voltage_type(nominal_voltage)

        self.Vnom = nominal_voltage
        self.ang = phase_angle
        self.Vset = voltage_setpoint
        # Check the node connect type
        # 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
        self.connect_type = connect_type
        self.stamp_dual = stamp_dual
        self.slack_parameters = slack_parameters

        if slack_parameters['file type'] == 'DSS':
            # Calculate slack generator voltages depending on if line or phase.
            # self.Va = self.Vnom * self.Vset * (self.cosd(-0.4) + 1j * self.sind(-0.4)) if voltage_A is None else voltage_A#(self.Vnom *
                    #self.Vset) + 0 * 1j if voltage_A is None else voltage_A
            # self.Vb = (self.Vnom * self.Vset) * (self.cosd(-120.5) +
            #    self.sind(-120.5) * 1j) if voltage_B is None else voltage_B
                #self.cosd(-120) +
                #self.sind(-120) * 1j) if voltage_B is None else voltage_B
            #self.Vc = (self.Vnom * self.Vset) * (self.cosd(119.6) +
            #    self.sind(119.6) * 1j) if voltage_C is None else voltage_C
                #self.cosd(120) +
                #self.sind(120) * 1j) if voltage_C is None else voltage_C
            self.impedance_matrix = slack_parameters['impedance_matrix']
            self.admittance_matrix = np.linalg.inv(self.impedance_matrix)
            self.G_self = self.admittance_matrix[0,0].real
            self.G_mutual_1 = self.admittance_matrix[0,1].real
            self.G_mutual_2 = self.admittance_matrix[0,2].real
            self.B_self = -self.admittance_matrix[0,0].imag
            self.B_mutual_1 = -self.admittance_matrix[0,1].imag
            self.B_mutual_2 = -self.admittance_matrix[0,2].imag
        
        # Calculate slack generator voltages depending on if line or phase.
        self.Va = (self.Vnom * self.Vset) + 0 * 1j if voltage_A is None else voltage_A
        self.Vb = (self.Vnom * self.Vset) * (self.cosd(-120) + 
            self.sind(-120) * 1j) if voltage_B is None else voltage_B
        self.Vc = (self.Vnom * self.Vset) * (self.cosd(120) +
            self.sind(120) * 1j) if voltage_C is None else voltage_C

        # if delta connected
        if self.connect_type == 0:
            self.Vab = self.Va - self.Vb
            self.Vbc = self.Vb - self.Vc
            self.Vca = self.Vc - self.Va

        self.Vn = complex(0, 0)
        self.Ia = complex(0, 0)
        self.Ib = complex(0, 0)
        self.Ic = complex(0, 0)
        self.Sa = complex(0, 0)
        self.Sb = complex(0, 0)
        self.Sc = complex(0, 0)

        self.S = self.Sa + self.Sb + self.Sc

    def stamp_linear(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, Jlin_val,
                     Jlin_row, idx_Y, idx_J):

        slack_nodeA_Vr = node[node_key[self.ID]].slack_nodeA_Vr
        slack_nodeB_Vr = node[node_key[self.ID]].slack_nodeB_Vr
        slack_nodeC_Vr = node[node_key[self.ID]].slack_nodeC_Vr

        slack_nodeA_Vi = node[node_key[self.ID]].slack_nodeA_Vi
        slack_nodeB_Vi = node[node_key[self.ID]].slack_nodeB_Vi
        slack_nodeC_Vi = node[node_key[self.ID]].slack_nodeC_Vi

        nodeA_Vr = node[node_key[self.ID]].nodeA_Vr
        nodeB_Vr = node[node_key[self.ID]].nodeB_Vr
        nodeC_Vr = node[node_key[self.ID]].nodeC_Vr

        nodeA_Vi = node[node_key[self.ID]].nodeA_Vi
        nodeB_Vi = node[node_key[self.ID]].nodeB_Vi
        nodeC_Vi = node[node_key[self.ID]].nodeC_Vi

        Va_real = self.Va.real
        Vb_real = self.Vb.real
        Vc_real = self.Vc.real

        Va_imag = self.Va.imag
        Vb_imag = self.Vb.imag
        Vc_imag = self.Vc.imag

        if self.slack_parameters['file type'] == 'DSS':
            real_slack_current_nodes = [slack_nodeA_Vr, slack_nodeB_Vr, slack_nodeC_Vr]
            imag_slack_current_nodes = [slack_nodeA_Vi, slack_nodeB_Vi, slack_nodeC_Vi]

            Vr_nodes = [nodeA_Vr, nodeB_Vr, nodeC_Vr]
            Vi_nodes = [nodeA_Vi, nodeB_Vi, nodeC_Vi]

            V_set_real = [Va_real, Vb_real, Vc_real]
            V_set_imag = [Va_imag, Vb_imag, Vc_imag]

            for phase in range(0, len(Vr_nodes)):
                if phase == 0:
                    mutual_phase_index = [1, 2]
                elif phase == 1:
                    mutual_phase_index = [0, 2]
                else:
                    mutual_phase_index = [0, 1]
                
                r_self = self.impedance_matrix[phase, phase].real
                x_self = self.impedance_matrix[phase, phase].imag
                G_self = self.admittance_matrix[phase,phase].real
                B_self = self.admittance_matrix[phase,phase].imag
                # # STAMPING THE VOLTAGE EQUATION IN THE SLACK CURRENT NODE ROW
                # # Vx - Vset - Is*Zs
                # # REAL
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], Vr_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], real_slack_current_nodes[phase], -r_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], imag_slack_current_nodes[phase], x_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)

                # Jlin_val, Jlin_row, idx_J = Slack.stamp_J(real_slack_current_nodes[phase], V_set_real[phase], Jlin_val, Jlin_row, idx_J)
                
                # # IMAGINARY
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], Vi_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], real_slack_current_nodes[phase], -x_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], imag_slack_current_nodes[phase], -r_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                
                # Jlin_val, Jlin_row, idx_J = Slack.stamp_J(imag_slack_current_nodes[phase], V_set_imag[phase], Jlin_val, Jlin_row, idx_J)
                
                # # MUTUAL STAMPS
                # # STAMPING THE VOLTAGE NODE
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], real_slack_current_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], imag_slack_current_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
               
                # for mutual_phase in range(0, len(mutual_phase_index)):
                #     B_mutual = 0.5*self.impedance_matrix[phase, mutual_phase_index[mutual_phase]].imag
                #     G_mutual = 0.5*self.impedance_matrix[phase, mutual_phase_index[mutual_phase]].real
                #     # REAL
                #     # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], Vr_nodes[phase], 0.5*G_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     #Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], Vi_nodes[phase], -B_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], Vr_nodes[mutual_phase_index[mutual_phase]], G_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], Vi_nodes[mutual_phase_index[mutual_phase]], -B_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    
                #     # IMAGINARY
                #     # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], Vr_nodes[mutual_phase_index[mutual_phase]], B_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     # Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], Vi_nodes[mutual_phase_index[mutual_phase]], G_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], real_slack_current_nodes[mutual_phase_index[mutual_phase]], B_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                #     Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], imag_slack_current_nodes[mutual_phase_index[mutual_phase]], G_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)

                # STAMPING THE CURRENTS IN THE VOLTAGE NODE
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vr_nodes[phase], real_slack_current_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(Vi_nodes[phase], imag_slack_current_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
               
                # STAMPING THE VOLTAGE EQUATION IN THE SLACK CURRENT NODE ROW
                # REAL
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], Vr_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], real_slack_current_nodes[phase], -r_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], imag_slack_current_nodes[phase], x_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)

                # IMAGINARY
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], Vi_nodes[phase], 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], real_slack_current_nodes[phase], -x_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], imag_slack_current_nodes[phase], -r_self, Ylin_val, Ylin_row, Ylin_col, idx_Y)

                # MUTUAL STAMPS
                for mutual_phase in range(0, len(mutual_phase_index)):
                    r_mutual = self.impedance_matrix[phase, mutual_phase_index[mutual_phase]].real
                    x_mutual = self.impedance_matrix[phase, mutual_phase_index[mutual_phase]].imag
                    
                    # REAL
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], real_slack_current_nodes[mutual_phase_index[mutual_phase]], -r_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(real_slack_current_nodes[phase], imag_slack_current_nodes[mutual_phase_index[mutual_phase]], x_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    
                    # IMAGINARY
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], real_slack_current_nodes[mutual_phase_index[mutual_phase]], -x_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(imag_slack_current_nodes[phase], imag_slack_current_nodes[mutual_phase_index[mutual_phase]], -r_mutual, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(real_slack_current_nodes[phase], V_set_real[phase], Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(imag_slack_current_nodes[phase], V_set_imag[phase], Jlin_val, Jlin_row, idx_J)

        else:

            # Current Values
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeA_Vr, nodeA_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeA_Vi, nodeA_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeB_Vr, nodeB_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeB_Vi, nodeB_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeC_Vr, nodeC_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                slack_nodeC_Vi, nodeC_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)

            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeA_Vr, slack_nodeA_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeA_Vi, slack_nodeA_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeB_Vr, slack_nodeB_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeB_Vi, slack_nodeB_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeC_Vr, slack_nodeC_Vr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                nodeC_Vi, slack_nodeC_Vi, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)

            # Voltage set value equations
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeA_Vr, Va_real,
                                                    Jlin_val, Jlin_row, idx_J)
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeA_Vi, Va_imag,
                                                    Jlin_val, Jlin_row, idx_J)
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeB_Vr, Vb_real,
                                                    Jlin_val, Jlin_row, idx_J)
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeB_Vi, Vb_imag,
                                                    Jlin_val, Jlin_row, idx_J)
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeC_Vr, Vc_real,
                                                    Jlin_val, Jlin_row, idx_J)
            Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeC_Vi, Vc_imag,
                                                    Jlin_val, Jlin_row, idx_J)
            
            # Currently not stamping the dual variable at the slack nodes
            if self.stamp_dual:
                slack_nodeA_Lr = node[node_key[self.ID]].slack_nodeA_Lr
                slack_nodeB_Lr = node[node_key[self.ID]].slack_nodeB_Lr
                slack_nodeC_Lr = node[node_key[self.ID]].slack_nodeC_Lr

                slack_nodeA_Li = node[node_key[self.ID]].slack_nodeA_Li
                slack_nodeB_Li = node[node_key[self.ID]].slack_nodeB_Li
                slack_nodeC_Li = node[node_key[self.ID]].slack_nodeC_Li

                nodeA_Lr = node[node_key[self.ID]].nodeA_dual_eq_var_r
                nodeB_Lr = node[node_key[self.ID]].nodeB_dual_eq_var_r
                nodeC_Lr = node[node_key[self.ID]].nodeC_dual_eq_var_r

                nodeA_Li = node[node_key[self.ID]].nodeA_dual_eq_var_i
                nodeB_Li = node[node_key[self.ID]].nodeB_dual_eq_var_i
                nodeC_Li = node[node_key[self.ID]].nodeC_dual_eq_var_i

                # Current Values
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeA_Lr, slack_nodeA_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeA_Li, slack_nodeA_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeB_Lr, slack_nodeB_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeB_Li, slack_nodeB_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeC_Lr, slack_nodeC_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    nodeC_Li, slack_nodeC_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)

                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeA_Lr, nodeA_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeA_Li, nodeA_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeB_Lr, nodeB_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeB_Li, nodeB_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeC_Lr, nodeC_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                Ylin_val, Ylin_row, Ylin_col, idx_Y = Slack.stamp_Y(
                    slack_nodeC_Li, nodeC_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeA_Lr, 0, Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeA_Li, 0, Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeB_Lr, 0, Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeB_Li, 0, Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeC_Lr, 0, Jlin_val, Jlin_row, idx_J)
                Jlin_val, Jlin_row, idx_J = Slack.stamp_J(slack_nodeC_Li, 0, Jlin_val, Jlin_row, idx_J)


        return Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J
