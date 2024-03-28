"""Provides Triplex Load class for distribution system analysis.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-11-2017
  Updated Date: 10-11-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""
from __future__ import division

from classes.Elements import NonlinearElement
from classes.GlobalVars import connected_nodes


class TriplexLoads(NonlinearElement):

    @staticmethod
    def calcPQpartials(P, Q, Vr, Vi, Lr, Li, stamp_dual):
        QVrVi = Q * Vr * Vi
        PVrVi = P * Vr * Vi
        Vr2Vi2 = Vr * Vr - Vi * Vi
        Vmag2 = Vr * Vr + Vi * Vi
        Vmag4 = Vmag2 * Vmag2
        Irl_hist = (P * Vr + Q * Vi) / Vmag2
        Iil_hist = (P * Vi - Q * Vr) / Vmag2
        dIrldVrl = -(P * Vr2Vi2 + 2.0 * QVrVi) / Vmag4  # dIr/dVr
        dIrldVil = (Q * Vr2Vi2 - 2.0 * PVrVi) / Vmag4  # dIr/dVi
        dIildVil = (P * Vr2Vi2 + 2.0 * QVrVi) / Vmag4  # dIi/dVi
        dIildVrl = dIrldVil  # dIi/dVi
        
        if stamp_dual:
            Vmag6 = Vmag4*Vmag2
            Lr_hist = Lr*dIrldVrl + Li*dIildVrl
            Li_hist = Lr*dIrldVil + Li*dIildVil
            # HESSIAN TERMS FROM TAYLOR APPROXIMATION
            # dIr/dVr, dIr/dVr*dVi, dIr/dVi*dVr, dIr/Vi (all second derivatives)
            dIr2_dVr2 =    2*(P*Vr**3 + 3*Q*Vr**2*Vi + -3*P*Vr*Vi**2 - Q*Vi**3)/Vmag6
            dIr2_dVrdVi = -2*(Q*Vr**3 - 3*P*Vr**2*Vi - 3*Q*Vr*Vi**2 + P*Vi**3)/Vmag6
            dIr2_dVidVr = dIr2_dVrdVi
            dIr2_dVi2 =  -dIr2_dVr2  
            
            # dIi/dVr, dIi/dVr*dVi, dIi/dVi*dVr, dIi/dVi (all second derivatives)
            dIi2_dVr2 =    dIr2_dVrdVi
            dIi2_dVrdVi = -dIr2_dVr2
            dIi2_dVidVr = -dIr2_dVr2
            dIi2_dVi2 =   -dIr2_dVrdVi
            
            # TAYLOR APPROXIMATION [LINEAR TERMS]
            dLr_dVr = Lr*dIr2_dVr2 + Li*dIi2_dVr2
            dLr_dVi = Lr*dIr2_dVrdVi + Li*dIi2_dVrdVi
            dLr_dLr = dIrldVrl
            dLr_dLi = dIildVrl
            
            dLi_dVr = Lr*dIr2_dVidVr + Li*dIi2_dVidVr
            dLi_dVi = Li*dIr2_dVi2 + Li*dIi2_dVi2
            dLi_dLr = dIrldVil
            dLi_dLi = dIildVil
        else:
            Lr_hist = 0
            Li_hist = 0
            dLr_dVr = 0
            dLr_dVi = 0
            dLr_dLr = 0
            dLr_dLi = 0
            dLi_dVr = 0
            dLi_dVi = 0
            dLi_dLr = 0
            dLi_dLi = 0
        return (Irl_hist, Iil_hist, dIrldVrl, dIrldVil, dIildVrl, dIildVil,
                Lr_hist, Li_hist, 
                dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi,
                dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi)

    @staticmethod
    def calc_G_B(r, x):
        if r or x:
            G = r / (r**2 + x**2)
            B = -x / (r**2 + x**2)
            return G, B
        else:
            return 0, 0

    def __init__(self,
                 name,
                 ID,
                 parent,
                 phases,
                 current1,
                 current12,
                 current2,
                 currentN,
                 impedance1,
                 impedance12,
                 impedance2,
                 nominal_voltage,
                 power1=None,
                 power12=None,
                 power2=None,
                 shunt1=None,
                 shunt12=None,
                 shunt2=None,
                 stamp_dual = False):
        super(TriplexLoads, self).__init__()
        self.name = name
        self.ID = ID
        self.parent = parent
        self.phases = phases
        self.nominal_voltage = int(nominal_voltage)
        self.stamp_dual = stamp_dual
        # Triplex Node Load Parameters
        self.current1 = current1 if current1 else 0 + 1j * 0
        self.current2 = current2 if current2 else 0 + 1j * 0
        self.current12 = current12 if current12 else 0 + 1j * 0
        self.currentN = currentN if currentN else 0 + 1j * 0
        self.impedance1 = impedance1 if impedance1 else 0 + 1j * 0
        self.impedance2 = impedance2 if impedance2 else 0 + 1j * 0
        self.impedance12 = impedance12 if impedance12 else 0 + 1j * 0
        self.power1 = power1 if power1 else 0 + 1j * 0
        self.power2 = power2 if power2 else 0 + 1j * 0
        self.power12 = power12 if power12 else 0 + 1j * 0
        self.shunt1 = shunt1 if shunt1 else 0 + 1j * 0
        self.shunt2 = shunt2 if shunt2 else 0 + 1j * 0
        self.shunt12 = shunt12 if shunt12 else 0 + 1j * 0
        # Calculate conductance/susceptance if any
        (self.cG_12,
         self.cB_12) = TriplexLoads.calc_G_B(self.impedance12.real,
                                             self.impedance12.imag)
        (self.cG_1, self.cB_1) = TriplexLoads.calc_G_B(self.impedance1.real,
                                                       self.impedance1.imag)
        (self.cG_2, self.cB_2) = TriplexLoads.calc_G_B(self.impedance2.real,
                                                       self.impedance2.imag)
        # Check if any load between 1 and 2
        if self.power12 or self.shunt12 or self.current12 or self.impedance12:
            self.type12 = True
            self.isGnd = False
            # Getting additional nodes to ground neutral
        else:
            self.type12 = False
            self.isGnd = True
        # Check if any load between 1 and Gnd
        if self.power1 or self.shunt1 or self.current1 or self.impedance1:
            self.type1 = True
            self.isGnd = True
        else:
            self.type1 = False
            self.isGnd = False

        # Check if any load between 2 and Gnd
        if self.power2 or self.shunt2 or self.current2 or self.impedance2:
            self.type2 = True
            self.isGnd = True
        else:
            self.type2 = False
            self.isGnd = False

        # Add connected nodes
        connected_nodes.add(self.ID)

    def _phase_calculations(self, _cP, _cQ, _cG, _cB, _cIr, _cIi, _cGshunt,
                            _cBshunt, node_Vr_from, node_Vi_from, node_Vr_to,
                            node_Vi_to, Vr_node_from, Vi_node_from, Vr_node_to,
                            Vi_node_to, Ynlin_val, Ynlin_row, Ynlin_col,
                            Jnlin_val, Jnlin_row, idx_Y, idx_J, 
                            node_Lr_from, node_Li_from, node_Lr_to, node_Li_to,
                            Lr_from, Li_from, Lr_to, Li_to):

        if not self.isGnd:
            # Calculate partials for PQ loads
            Vr_across = Vr_node_from - Vr_node_to
            Vi_across = Vi_node_from - Vi_node_to
            
            if self.stamp_dual:
                Lr_across = Lr_from - Lr_to
                Li_across = Lr_to - Li_to
            else:
                Lr_across = []
                Li_across = []

            # Calculate Partials from the voltage difference between from and to node
            (_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVrl, _dIildVil, 
                 Lr_hist, Li_hist, 
                 dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi,
                 dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi) = \
             TriplexLoads.calcPQpartials(_cP, _cQ, Vr_across, Vi_across, Lr_across, Li_across, self.stamp_dual)
             

            # 4 real and 4 imaginary stamps for the from circuit
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from,
                                _dIrldVrl + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)  # G
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_to,
                                -_dIrldVrl - (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from,
                                _dIrldVil - (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_to,
                                -_dIrldVil + (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from,
                                _dIildVil + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_to,
                                -_dIildVil - (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from,
                                _dIildVrl + (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_to,
                                -_dIildVrl - (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            # 4 real and imaginary circuits for the to circuit
            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_from,
                                -_dIrldVrl - (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_to,
                                _dIrldVrl + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_from,
                                -_dIrldVil + (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_to,
                                _dIrldVil - (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_from,
                                -_dIildVil - (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_to,
                                _dIildVil + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_from,
                                -_dIildVrl - (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_to,
                                _dIildVrl + (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)

            # Independent current sources for from and to branches
            idx_J = self.stamp_J(
                node_Vr_from,
                -((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
                  _cIr), Jnlin_val, Jnlin_row, idx_J)
            #
            idx_J = self.stamp_J(
                node_Vi_from,
                -((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
                  _cIi), Jnlin_val, Jnlin_row, idx_J)
            #
            idx_J = self.stamp_J(
                node_Vr_to,
                ((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
                 _cIr), Jnlin_val, Jnlin_row, idx_J)
            #
            idx_J = self.stamp_J(
                node_Vi_to,
                ((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
                 _cIi), Jnlin_val, Jnlin_row, idx_J)
            
            if self.stamp_dual:
                idx_Y = self.stamp_Y(node_Lr_from, node_Vr_from, dLr_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Vr_to, -dLr_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Vi_from, dLr_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Vi_to, -dLr_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)               
                
                
                idx_Y = self.stamp_Y(node_Lr_from, node_Lr_from, dLr_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Lr_to, -dLr_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)               
                idx_Y = self.stamp_Y(node_Lr_from, node_Li_from, dLr_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Li_to, -dLr_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                idx_J = self.stamp_J(node_Lr_from, ( \
                            - Lr_hist  \
                            + Vr_across*dLr_dVr + Vi_across*dLr_dVi  \
                            + Lr_across*dLr_dLr + Li_across*dLr_dLi\
                                  ), Jnlin_val, Jnlin_row, idx_J)
                # idx_J = self.stamp_J(node_Lr_from, ( \
                #             - Lr_across*_dIrldVrl - Li_across*_dIildVrl  \
                #             + Vr_across*dLr_dVr + Vi_across*dLr_dVi  \
                #             + Lr_across*dLr_dLr + Li_across*dLr_dLi\
                #                   ), Jnlin_val, Jnlin_row, idx_J)
                # LR_to
                idx_Y = self.stamp_Y(node_Lr_to, node_Vr_from, -dLr_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_to, node_Vr_to, dLr_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_to, node_Vi_from, -dLr_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_to, node_Vi_to, dLr_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)               
                
                
                idx_Y = self.stamp_Y(node_Lr_to, node_Lr_from, -dLr_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_to, node_Lr_to, dLr_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)               
                idx_Y = self.stamp_Y(node_Lr_to, node_Li_from, -dLr_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_to, node_Li_to, dLr_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                idx_J = self.stamp_J(node_Lr_to, ( \
                            Lr_hist  \
                            - Vr_across*dLr_dVr - Vi_across*dLr_dVi  \
                            - Lr_across*dLr_dLr - Li_across*dLr_dLi \
                                  ), Jnlin_val, Jnlin_row, idx_J)
                
                    # idx_J = self.stamp_J(node_Lr_to, ( \
                #             - (-(Lr_across)*_dIrldVrl) - (-(Li_across)*_dIildVrl)  \
                #             + (-(Vr_across)*dLr_dVr) + (-(Vi_across)*dLr_dVi)  \
                #             + (-(Lr_across)*dLr_dLr) + (-(Li_across)*dLr_dLi) \
                #                   ), Jnlin_val, Jnlin_row, idx_J)
                
                # LI_from
                idx_Y = self.stamp_Y(node_Li_from, node_Vr_from, dLi_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Vr_to, -dLi_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Vi_from, dLi_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Vi_to, -dLi_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                idx_Y = self.stamp_Y(node_Li_from, node_Lr_from, dLi_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Lr_to, -dLi_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Li_from, dLi_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Li_to, -dLi_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                 
                idx_J = self.stamp_J(node_Li_from, ( \
                                    - Li_hist \
                                    + Vr_across*dLi_dVr + Vi_across*dLi_dVi \
                                    + Lr_across*dLi_dLr + Li_across*dLi_dLi \
                                        ), Jnlin_val, Jnlin_row, idx_J)  
                # idx_J = self.stamp_J(node_Li_from, ( \
                #                     - Lr_across*_dIrldVil - Li_across*_dIildVil \
                #                     + Vr_across*dLi_dVr + Vi_across*dLi_dVi \
                #                     + Lr_across*dLi_dLr + Li_across*dLi_dLi \
                #                         ), Jnlin_val, Jnlin_row, idx_J)                   
                #LI_to
                idx_Y = self.stamp_Y(node_Li_to, node_Vr_from, -dLi_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Vr_to, dLi_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Vi_from, -dLi_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Vi_to, dLi_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                idx_Y = self.stamp_Y(node_Li_to, node_Lr_from, -dLi_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Lr_to, dLi_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Li_from, -dLi_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_to, node_Li_to, dLi_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                idx_J = self.stamp_J(node_Li_to, ( \
                                    Li_hist \
                                    - Vr_across*dLi_dVr - Vi_across*dLi_dVi \
                                    - Lr_across*dLi_dLr - Li_across*dLi_dLi \
                                        ), Jnlin_val, Jnlin_row, idx_J)                
                # idx_J = self.stamp_J(node_Li_to, ( \
                #                     - (-(Lr_across)*_dIrldVil) - (-(Li_across)*_dIildVil) \
                #                     + (-(Vr_across)*dLi_dVr) + (-(Vi_across)*dLi_dVi) \
                #                     + (-(Lr_across)*dLi_dLr) + (-(Li_across)*dLi_dLi) \
                #                         ), Jnlin_val, Jnlin_row, idx_J)
 
        else:
            # Calculate partials for PQ TriplexLoads
            Vr_across = Vr_node_from
            Vi_across = Vi_node_from
            if self.stamp_dual:
                Lr_across = Lr_from
                Li_across = Li_from
            else:
                Lr_across = 0
                Li_across = 0

            # Calculate Partials from the voltage difference between from and to node
            (_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVrl, _dIildVil, 
                         Lr_hist, Li_hist, dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi,
                            dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi) = \
                         TriplexLoads.calcPQpartials(_cP, _cQ, Vr_across, Vi_across, Lr_across, Li_across, self.stamp_dual)
             
            # 4 real and 4 imaginary stamps for the from circuit
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from,
                                _dIrldVrl + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)  # G
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from,
                                _dIrldVil - (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from,
                                _dIildVil + (_cG + _cGshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from,
                                _dIildVrl + (_cB + _cBshunt), Ynlin_val,
                                Ynlin_row, Ynlin_col, idx_Y)

            # Independent current sources for from and to branches
            idx_J = self.stamp_J(
                node_Vr_from,
                -((_Irl_hist - _dIrldVrl * Vr_across - _dIrldVil * Vi_across) +
                  _cIr), Jnlin_val, Jnlin_row, idx_J)
            #
            idx_J = self.stamp_J(
                node_Vi_from,
                -((_Iil_hist - _dIildVrl * Vr_across - _dIildVil * Vi_across) +
                  _cIi), Jnlin_val, Jnlin_row, idx_J)
            
            if self.stamp_dual:
                # Real
                idx_Y = self.stamp_Y(node_Lr_from, node_Vr_from, dLr_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Vi_from, dLr_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Lr_from, dLr_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Lr_from, node_Li_from, dLr_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                # Imaginary
                idx_Y = self.stamp_Y(node_Li_from, node_Vr_from, dLi_dVr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Vi_from, dLi_dVi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Lr_from, dLi_dLr, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                idx_Y = self.stamp_Y(node_Li_from, node_Li_from, dLi_dLi, Ynlin_val,
                                    Ynlin_row, Ynlin_col, idx_Y)
                
                # J
                idx_J = self.stamp_J(node_Lr_from, (-Lr_from*_dIrldVrl - Li_from*_dIildVrl + \
                                    Vr_across*dLr_dVr + Vi_across*dLr_dVi + \
                                    Lr_from*dLr_dLr + Li_from*dLr_dLi), Jnlin_val, Jnlin_row, idx_J)
                idx_J = self.stamp_J(node_Li_from, (-Lr_from*_dIrldVil - Li_from*_dIildVil + \
                                    Vr_across*dLi_dVr + Vi_across*dLi_dVi + \
                                    Lr_from*dLi_dLr + Li_from*dLi_dLi), Jnlin_val, Jnlin_row, idx_J)
                                               
        return idx_Y, idx_J

    def stamp_nonlinear(self, node_key, node, V, lf, Ynlin_val, Ynlin_row,
                        Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J):

        # Assuming neutral is grounded
        # Find the nodes to stamp
        node1_Vr = node[node_key[self.ID]].node1_Vr
        node1_Vi = node[node_key[self.ID]].node1_Vi
        node2_Vr = node[node_key[self.ID]].node2_Vr
        node2_Vi = node[node_key[self.ID]].node2_Vi
        nodeN_Vr = node[node_key[self.ID]].nodeN_Vr
        nodeN_Vi = node[node_key[self.ID]].nodeN_Vi

       
        # Find the node voltages
        (Vr_node1, Vi_node1, Vr_node2, Vi_node2, Vr_nodeN, Vi_nodeN) = \
         (V[node1_Vr], V[node1_Vi], V[node2_Vr], V[node2_Vi], V[nodeN_Vr],
          V[nodeN_Vi])
         
        if self.stamp_dual:
            node1_Lr = node[node_key[self.ID]].node1_dual_eq_var_r
            node1_Li = node[node_key[self.ID]].node1_dual_eq_var_i
            node2_Lr = node[node_key[self.ID]].node2_dual_eq_var_r
            node2_Li = node[node_key[self.ID]].node2_dual_eq_var_i
            nodeN_Lr = node[node_key[self.ID]].nodeN_dual_eq_var_r
            nodeN_Li = node[node_key[self.ID]].nodeN_dual_eq_var_i
            
            (Lr_node1, Li_node1, Lr_node2, Li_node2, Lr_nodeN, \
                 Li_nodeN) = (V[node1_Lr], V[node1_Li], V[node2_Lr], \
                                  V[node2_Li], V[nodeN_Lr], V[nodeN_Li])
        
        else:
            node1_Lr = []
            node1_Li = []
            node2_Lr = []
            node2_Li = []
            nodeN_Lr = []
            nodeN_Li = []
           
            Lr_node1 = []
            Lr_node2 = []
            Lr_nodeN = []
            
            Li_node1 = []
            Li_node2 = []
            Li_nodeN = []
         

        # If the load is connect between phases 1 and 2
        if self.type12:
            _cP = self.power12.real * lf
            _cQ = self.power12.imag * lf
            _cG = self.cG_12 * lf
            _cB = self.cB_12 * lf
            _cIr = self.current12.real * lf
            _cIi = self.current12.imag * lf
            _cGshunt = self.shunt12.real
            _cBshunt = self.shunt12.imag
            node_Vr_from = node1_Vr
            node_Vi_from = node1_Vi
            node_Vr_to = node2_Vr
            node_Vi_to = node2_Vi
            Vr_node_from = Vr_node1
            Vi_node_from = Vi_node1
            Vr_node_to = Vr_node2
            Vi_node_to = Vi_node2

            if _cP or _cQ or _cG or _cB or _cIr or _cIi or _cGshunt or _cBshunt:
                idx_Y, idx_J = self._phase_calculations(
                    _cP, _cQ, _cG, _cB, _cIr, _cIi, _cGshunt, _cBshunt,
                    node_Vr_from, node_Vi_from, node_Vr_to, node_Vi_to,
                    Vr_node_from, Vi_node_from, Vr_node_to, Vi_node_to,
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J, 
                    node1_Lr, node1_Li, node2_Lr, node2_Li,
                    Lr_node1, Li_node1, Lr_node2, Li_node2) 
                

        # If the load is connect between phases 1
        if self.type1:
            _cP = self.power1.real * lf
            _cQ = self.power1.imag * lf
            _cG = self.cG_1 * lf
            _cB = self.cB_1 * lf
            _cIr = self.current1.real * lf
            _cIi = self.current1.imag * lf
            _cGshunt = self.shunt1.real
            _cBshunt = self.shunt1.imag
            node_Vr_from = node1_Vr
            node_Vi_from = node1_Vi
            node_Vr_to = nodeN_Vr
            node_Vi_to = nodeN_Vi
            Vr_node_from = Vr_node1
            Vi_node_from = Vi_node1
            Vr_node_to = Vr_nodeN
            Vi_node_to = Vi_nodeN

            if _cP or _cQ or _cG or _cB or _cIr or _cIi or _cGshunt or _cBshunt:
                idx_Y, idx_J = self._phase_calculations(
                    _cP, _cQ, _cG, _cB, _cIr, _cIi, _cGshunt, _cBshunt,
                    node_Vr_from, node_Vi_from, node_Vr_to, node_Vi_to,
                    Vr_node_from, Vi_node_from, Vr_node_to, Vi_node_to,
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J,
                    node1_Lr, node1_Li, nodeN_Lr, nodeN_Li,
                    Lr_node1, Li_node1, Lr_nodeN, Li_nodeN) 

        # If the load is connect between phases 2
        if self.type2:
            _cP = self.power2.real * lf
            _cQ = self.power2.imag * lf
            _cG = self.cG_2 * lf
            _cB = self.cB_2 * lf
            _cIr = self.current2.real * lf
            _cIi = self.current2.imag * lf
            _cGshunt = self.shunt2.real
            _cBshunt = self.shunt2.imag
            node_Vr_from = node2_Vr
            node_Vi_from = node2_Vi
            node_Vr_to = nodeN_Vr
            node_Vi_to = nodeN_Vi
            Vr_node_from = Vr_node2
            Vi_node_from = Vi_node2
            Vr_node_to = Vr_nodeN
            Vi_node_to = Vi_nodeN

            if _cP or _cQ or _cG or _cB or _cIr or _cIi or _cGshunt or _cBshunt:
                idx_Y, idx_J = self._phase_calculations(
                    _cP, _cQ, _cG, _cB, _cIr, _cIi, _cGshunt, _cBshunt,
                    node_Vr_from, node_Vi_from, node_Vr_to, node_Vi_to,
                    Vr_node_from, Vi_node_from, Vr_node_to, Vi_node_to,
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J, 
                    node2_Lr, node2_Li, nodeN_Lr, nodeN_Li,
                    Lr_node2, Li_node2, Lr_nodeN, Li_nodeN) 

        return idx_Y, idx_J

    def calc_residual(self, V, node_key, node, res_eqn, lf):
        node1_Vr = node[node_key[self.ID]].node1_Vr
        node1_Vi = node[node_key[self.ID]].node1_Vi
        node2_Vr = node[node_key[self.ID]].node2_Vr
        node2_Vi = node[node_key[self.ID]].node2_Vi
        nodeN_Vr = node[node_key[self.ID]].nodeN_Vr
        nodeN_Vi = node[node_key[self.ID]].nodeN_Vi

       
        # Find the node voltages
        (Vr_node1, Vi_node1, Vr_node2, Vi_node2, Vr_nodeN, Vi_nodeN) = \
         (V[node1_Vr], V[node1_Vi], V[node2_Vr], V[node2_Vi], V[nodeN_Vr],
          V[nodeN_Vi])        
        
        if self.stamp_dual:
            node1_Lr = node[node_key[self.ID]].node1_dual_eq_var_r
            node1_Li = node[node_key[self.ID]].node1_dual_eq_var_i
            node2_Lr = node[node_key[self.ID]].node2_dual_eq_var_r
            node2_Li = node[node_key[self.ID]].node2_dual_eq_var_i
            nodeN_Lr = node[node_key[self.ID]].nodeN_dual_eq_var_r
            nodeN_Li = node[node_key[self.ID]].nodeN_dual_eq_var_i
            
            (Lr_node1, Li_node1, Lr_node2, Li_node2, Lr_nodeN, Li_nodeN) = \
             (V[node1_Lr], V[node1_Li], V[node2_Lr], V[node2_Li], V[nodeN_Lr],
              V[nodeN_Li]) 
             
        if self.type12:
            _cP = self.power12.real * lf 
            _cQ = self.power12.imag * lf
            _cG = self.cG_12 * lf 
            _cB = self.cB_12 * lf 
            _cIr = self.current12.real * lf
            _cIi = self.current12.imag * lf 
            _cGshunt = self.shunt12.real
            _cBshunt = self.shunt12.imag
            node_Vr_from = node1_Vr
            node_Vi_from = node1_Vi
            node_Vr_to = node2_Vr
            node_Vi_to = node2_Vi
            Vr_from = Vr_node1
            Vi_from = Vi_node1
            Vr_to = Vr_node2
            Vi_to = Vi_node2   
            
            if self.stamp_dual:
                node_Lr_from = node1_Lr
                node_Li_from = node1_Li
                node_Lr_to = node2_Lr
                node_Li_to = node2_Li
                Lr_from = Lr_node1
                Li_from = Li_node1
                Lr_to = Lr_node2
                Li_to = Li_node2
            
            if self.isGnd:  
                # Equality
                ipq_from_denom = (Vr_from)**2 + (Vi_from)**2
                ipq_from_to_r = _cP*(Vr_from) + _cQ*(Vi_from)
                ipq_from_to_i = _cP*(Vi_from) - _cQ*(Vr_from) 
                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * (Vr_from * Vr_from - Vi_from * Vi_from) - 2.0 * _cP * Vr_from * Vi_from)  # dIr/dVi
                    dipq_dvi_from_i = (_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_from_denom = ipq_from_denom**2       
                
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_from_denom + (_cG + _cGshunt) * Vr_from - (_cB + _cBshunt)* Vi_from
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_from_denom + (_cB + _cBshunt) * Vr_from + (_cG + _cGshunt)* Vi_from
                    
                    if self.stamp_dual:
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_from_denom * Lr_from + dipq_dvi_from_r/dipq_from_denom * Li_from
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_from_denom * Lr_from + dipq_dvi_from_i/dipq_from_denom * Li_from             

            else:   
                Vr_across = Vr_from - Vr_to
                Vi_across = Vi_from - Vi_to 
                
                ipq_denom = Vr_across**2 + Vi_across**2
                ipq_from_to_r = _cP*Vr_across + _cQ*Vi_across
                ipq_from_to_i = _cP*Vi_across - _cQ*Vr_across
                
                ipq_to_from_r = -ipq_from_to_r
                ipq_to_from_i = -ipq_from_to_i
                                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * (Vr_across**2 - Vi_across**2) + 2.0 * _cQ * Vr_across * Vi_across)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * (Vr_across**2 - Vi_across**2) - 2.0 * _cP * Vr_across * Vi_across)  # dIr/dVi
                    dipq_dvi_from_i = -dipq_dvr_from_r # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_denom = ipq_denom**2                          
                
                if abs(ipq_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_denom + (_cG + _cGshunt) * Vr_across - (_cB + _cBshunt) * Vi_across
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_denom + (_cB + _cBshunt) * Vr_across + (_cG + _cGshunt) * Vi_across
                    res_eqn[node_Vr_to] += -_cIr - ipq_from_to_r/ipq_denom - (_cG + _cGshunt) * Vr_across + (_cB + _cBshunt) * Vi_across
                    res_eqn[node_Vi_to] += -_cIi - ipq_from_to_i/ipq_denom - (_cB + _cBshunt) * Vr_across - (_cG + _cGshunt) * Vi_across
                    
                    if self.stamp_dual:
                        Lr_across = Lr_from - Lr_to 
                        Li_across = Li_from - Li_to
                        
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_denom * (Lr_across) + dipq_dvi_from_r/dipq_denom * (Li_across)
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_denom * (Lr_across) + dipq_dvi_from_i/dipq_denom * (Li_across)
                
                        res_eqn[node_Lr_to] += -(dipq_dvr_from_r/dipq_denom * (Lr_across) + dipq_dvi_from_r/dipq_denom * (Li_across))#dipq_dvr_to_r/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_r/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase])
                        res_eqn[node_Li_to] += -(dipq_dvr_from_i/dipq_denom * (Lr_across) + dipq_dvi_from_i/dipq_denom * (Li_across))#dipq_dvr_to_i/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_i/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase]) 
                              
        if self.type1:
            _cP = self.power1.real * lf 
            _cQ = self.power1.imag * lf 
            _cG = self.cG_1 * lf 
            _cB = self.cB_1 * lf 
            _cIr = self.current1.real * lf 
            _cIi = self.current1.imag * lf 
            _cGshunt = self.shunt1.real
            _cBshunt = self.shunt1.imag
            node_Vr_from = node1_Vr
            node_Vi_from = node1_Vi
            node_Vr_to = nodeN_Vr
            node_Vi_to = nodeN_Vi
            Vr_from = Vr_node1
            Vi_from = Vi_node1
            Vr_to = Vr_nodeN
            Vi_to = Vi_nodeN
            
            if self.stamp_dual:
                node_Lr_from = node1_Lr
                node_Li_from = node1_Li
                node_Lr_to = nodeN_Lr
                node_Li_to = nodeN_Li
                Lr_from = Lr_node1
                Li_from = Li_node1
                Lr_to = Lr_nodeN
                Li_to = Li_nodeN
                
            if not self.isGnd: 
                ipq_from_denom = (Vr_from - Vr_to)**2 + (Vi_from - Vi_to)**2
                ipq_from_to_r = _cP*(Vr_from-Vr_to) + _cQ*(Vi_from - Vi_to)
                ipq_from_to_i = _cP*(Vi_from-Vi_to) - _cQ*(Vr_from - Vr_to)
                
                ipq_to_denom = (Vr_to - Vr_from)**2 + (Vi_to - Vi_from)**2
                ipq_to_from_r = _cP*(Vr_to - Vr_from) + _cQ*(Vi_to - Vi_from)
                ipq_to_from_i = _cP*(Vi_to - Vi_from) - _cQ*(Vr_to - Vr_from)
                                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) + 2.0 * _cQ * Vr_from * Vi_from)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) - 2.0 * _cP * Vr_from * Vi_from)  # dIr/dVi
                    dipq_dvi_from_i = (_cP * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) + 2.0 * _cQ * Vr_from * Vi_from)  # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_from_denom = ipq_from_denom**2                          
                
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_from_denom + (_cG + _cGshunt) * (Vr_from - Vr_to) - (_cB) * (Vi_from - Vi_to)
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_from_denom + (_cB + _cBshunt) * (Vr_from - Vr_to) + (_cG) * (Vi_from - Vi_to)
                    res_eqn[node_Vr_to] += -_cIr + ipq_to_from_r/ipq_to_denom - (_cG + _cGshunt) * (Vr_from - Vr_to) + (_cB) * (Vi_from - Vi_to)
                    res_eqn[node_Vi_to] += -_cIi + ipq_to_from_i/ipq_to_denom - (_cB + _cBshunt) * (Vr_from - Vr_to) - (_cG) * (Vi_from - Vi_to)
                    
                    if self.stamp_dual:
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_r/dipq_from_denom * (Li_from - Li_to)
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_i/dipq_from_denom * (Li_from - Li_to)
                
                        res_eqn[node_Lr_to] += -(dipq_dvr_from_r/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_r/dipq_from_denom * (Li_from - Li_to))
                        res_eqn[node_Li_to] += -(dipq_dvr_from_i/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_i/dipq_from_denom * (Li_from - Li_to))                  
            else: 
                # Equality
                ipq_from_denom = (Vr_from)**2 + (Vi_from)**2
                ipq_from_to_r = _cP*(Vr_from) + _cQ*(Vi_from)
                ipq_from_to_i = _cP*(Vi_from) - _cQ*(Vr_from) 
                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * (Vr_from * Vr_from - Vi_from * Vi_from) - 2.0 * _cP * Vr_from * Vi_from)  # dIr/dVi
                    dipq_dvi_from_i = (_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_from_denom = ipq_from_denom**2       
                
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_from_denom + (_cG + _cGshunt) * Vr_from - (_cB + _cBshunt) * Vi_from
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_from_denom + (_cB + _cBshunt) * Vr_from + (_cG + _cGshunt) * Vi_from
                    
                    if self.stamp_dual:
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_from_denom * Lr_from + dipq_dvi_from_r/dipq_from_denom * Li_from
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_from_denom * Lr_from + dipq_dvi_from_i/dipq_from_denom * Li_from    
             
        if self.type2:
            _cP = self.power2.real * lf 
            _cQ = self.power2.imag * lf 
            _cG = self.cG_2 * lf 
            _cB = self.cB_2 * lf 
            _cIr = self.current2.real * lf 
            _cIi = self.current2.imag * lf 
            _cGshunt = self.shunt2.real
            _cBshunt = self.shunt2.imag
            node_Vr_from = node2_Vr
            node_Vi_from = node2_Vi
            node_Vr_to = nodeN_Vr
            node_Vi_to = nodeN_Vi
            Vr_from = Vr_node2
            Vi_from = Vi_node2
            Vr_to = Vr_nodeN
            Vi_to = Vi_nodeN
            
            if self.stamp_dual:
                node_Lr_from = node2_Lr
                node_Li_from = node2_Li
                node_Lr_to = nodeN_Lr
                node_Li_to = nodeN_Li
                Lr_from = Lr_node2
                Li_from = Li_node2
                Lr_to = Lr_nodeN
                Li_to = Li_nodeN
            
            if self.isGnd:                   
                # Equality
                ipq_from_denom = (Vr_from)**2 + (Vi_from)**2
                ipq_from_to_r = _cP*(Vr_from) + _cQ*(Vi_from)
                ipq_from_to_i = _cP*(Vi_from) - _cQ*(Vr_from) 
                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * (Vr_from * Vr_from - Vi_from * Vi_from) - 2.0 * _cP * Vr_from * Vi_from)  # dIr/dVi
                    dipq_dvi_from_i = (_cP * (Vr_from * Vr_from - Vi_from * Vi_from) + 2.0 * _cQ * Vr_from * Vi_from)  # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_from_denom = ipq_from_denom**2       
                
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_from_denom + (_cG + _cGshunt) * Vr_from - (_cB + _cBshunt) * Vi_from
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_from_denom + (_cB + _cBshunt) * Vr_from + (_cG + _cGshunt) * Vi_from
                    
                    if self.stamp_dual:
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_from_denom * Lr_from + dipq_dvi_from_r/dipq_from_denom * Li_from
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_from_denom * Lr_from + dipq_dvi_from_i/dipq_from_denom * Li_from 
               
            else:                                    
                ipq_from_denom = (Vr_from - Vr_to)**2 + (Vi_from - Vi_to)**2
                ipq_from_to_r = _cP*(Vr_from-Vr_to) + _cQ*(Vi_from - Vi_to)
                ipq_from_to_i = _cP*(Vi_from-Vi_to) - _cQ*(Vr_from - Vr_to)
                
                ipq_to_denom = (Vr_to - Vr_from)**2 + (Vi_to - Vi_from)**2
                ipq_to_from_r = _cP*(Vr_to - Vr_from) + _cQ*(Vi_to - Vi_from)
                ipq_to_from_i = _cP*(Vi_to - Vi_from) - _cQ*(Vr_to - Vr_from)
                                
                if self.stamp_dual:
                    # Stationarity
                    dipq_dvr_from_r = -(_cP * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) + 2.0 * _cQ * Vr_from * Vi_from)  # dIr/dVr
                    dipq_dvi_from_r = (_cQ * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) - 2.0 * _cP * Vr_from * Vi_from)  # dIr/dVi
                    dipq_dvi_from_i = (_cP * ((Vr_from - Vr_to)**2 - (Vi_from - Vi_to)**2) + 2.0 * _cQ * Vr_from * Vi_from)  # dIi/dVi
                    dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                    dipq_from_denom = ipq_from_denom**2                          
                
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from] += _cIr + ipq_from_to_r/ipq_from_denom + (_cG + _cGshunt) * (Vr_from - Vr_to) - (_cB + _cBshunt) * (Vi_from - Vi_to)
                    res_eqn[node_Vi_from] += _cIi + ipq_from_to_i/ipq_from_denom + (_cB + _cBshunt) * (Vr_from - Vr_to) + (_cG + _cGshunt) * (Vi_from - Vi_to)
                    res_eqn[node_Vr_to] += -_cIr + ipq_to_from_r/ipq_to_denom - (_cG + _cGshunt) * (Vr_from - Vr_to) + (_cB + _cBshunt) * (Vi_from - Vi_to)
                    res_eqn[node_Vi_to] += -_cIi + ipq_to_from_i/ipq_to_denom - (_cB + _cBshunt) * (Vr_from - Vr_to) - (_cG + _cGshunt) * (Vi_from - Vi_to)
                    
                    if self.stamp_dual:
                        res_eqn[node_Lr_from] += dipq_dvr_from_r/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_r/dipq_from_denom * (Li_from - Li_to)
                        res_eqn[node_Li_from] += dipq_dvr_from_i/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_i/dipq_from_denom * (Li_from - Li_to)
                
                        res_eqn[node_Lr_to] += -(dipq_dvr_from_r/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_r/dipq_from_denom * (Li_from - Li_to))#dipq_dvr_to_r/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_r/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase])
                        res_eqn[node_Li_to] += -(dipq_dvr_from_i/dipq_from_denom * (Lr_from - Lr_to) + dipq_dvi_from_i/dipq_from_denom * (Li_from - Li_to))#dipq_dvr_to_i/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_i/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase]) 
  
           
        return res_eqn
