"""Reactors are constant impedance elements.

  Author(s): Naeem Turner-Bandele
  Created Date: 07-10-2020
  Updated Date: 07-20-2020
  Email: nturnerb@cmu.edu
  Status: Validated

  The Reactors class creates a constant impedance series element that can be used to control the impedance of a
  particular path and control current flow.

  Typical usage example:

  reactor = Reactors()
"""

from __future__ import division

from classes.Elements import LinearElement
from classes.Nodes import Nodes


class Reactors(LinearElement):

    def __init__(self, node, ID, name, from_node, to_node, connection_type,
                 phases, resistance, reactance, stamp_dual = False):
        super(Reactors, self).__init__()
        self.name = name
        self.ID = ID
        self.from_node = from_node
        self.to_node = to_node
        self.type = connection_type
        self.phases = phases
        self.G = []
        self.B = []
        self.stamp_dual = stamp_dual
        # calculate impedances
        if self.phases & 0x1 == 1:
            self.impedance_A = complex(resistance[0], reactance[0])
            Ga, Ba = self.calc_G_B(self.impedance_A.real, self.impedance_A.imag)
            self.G.append(Ga)
            self.B.append(Ba)
        else:
            self.G.append(0)
            self.B.append(0)
        if self.phases & 0x2 == 2:
            self.impedance_B = complex(resistance[1], reactance[1])
            Gb, Bb = self.calc_G_B(self.impedance_B.real, self.impedance_B.imag)
            self.G.append(Gb)
            self.B.append(Bb)
        else:
            self.G.append(0)
            self.B.append(0)
        if self.phases & 0x4 == 4:
            self.impedance_C = complex(resistance[2], reactance[2])
            Gc, Bc = self.calc_G_B(self.impedance_C.real, self.impedance_C.imag)
            self.G.append(Gc)
            self.B.append(Bc)
        else:
            self.G.append(0)
            self.B.append(0)
        if self.phases & 0x8 == 8:
            self.impedance_N = complex(resistance[3], reactance[3])
            Gn, Bn = self.calc_G_B(self.impedance_N.real, self.impedance_N.imag)
            self.G.append(Gn)
            self.B.append(Bn)
        else:
            self.G.append(0)
            self.B.append(0)

    @staticmethod
    def calc_G_B(r, x):
        if r == 0:
            G = None
        else:
            G = r / (r**2 + x**2)
        if x == 0:
            B = None
        else:
            B = -x / (r**2 + x**2)
        return G, B

    def assign_nodes(self, node_key, node):
        if self.phases & 0x1 == 1:  # Check for phase A
            self.from_nodeA_Vr = node[node_key[self.from_node]].nodeA_Vr
            self.from_nodeA_Vi = node[node_key[self.from_node]].nodeA_Vi
            self.to_nodeA_Vr = node[node_key[self.to_node]].nodeA_Vr
            self.to_nodeA_Vi = node[node_key[self.to_node]].nodeA_Vi
        if self.phases & 0x2 == 2:  # Check for phase B
            self.from_nodeB_Vr = node[node_key[self.from_node]].nodeB_Vr
            self.from_nodeB_Vi = node[node_key[self.from_node]].nodeB_Vi
            self.to_nodeB_Vr = node[node_key[self.to_node]].nodeB_Vr
            self.to_nodeB_Vi = node[node_key[self.to_node]].nodeB_Vi
        if self.phases & 0x4 == 4:  # Check for phase C
            self.from_nodeC_Vr = node[node_key[self.from_node]].nodeC_Vr
            self.from_nodeC_Vi = node[node_key[self.from_node]].nodeC_Vi
            self.to_nodeC_Vr = node[node_key[self.to_node]].nodeC_Vr
            self.to_nodeC_Vi = node[node_key[self.to_node]].nodeC_Vi

        self.from_nodeN_Vr = node[node_key[self.from_node]].nodeN_Vr
        self.from_nodeN_Vi = node[node_key[self.from_node]].nodeN_Vi
        self.to_nodeN_Vr = node[node_key[self.to_node]].nodeN_Vr
        self.to_nodeN_Vi = node[node_key[self.to_node]].nodeN_Vi
        
        if self.stamp_dual:
            if self.phases & 0x1 == 1:  # Check for phase A
                self.from_nodeA_Lr = node[node_key[self.from_node]].nodeA_dual_eq_var_r
                self.from_nodeA_Li = node[node_key[self.from_node]].nodeA_dual_eq_var_i
                self.to_nodeA_Lr = node[node_key[self.to_node]].nodeA_dual_eq_var_r
                self.to_nodeA_Li = node[node_key[self.to_node]].nodeA_dual_eq_var_i
            if self.phases & 0x2 == 2:  # Check for phase B
                self.from_nodeB_Lr = node[node_key[self.from_node]].nodeB_dual_eq_var_r
                self.from_nodeB_Li = node[node_key[self.from_node]].nodeB_dual_eq_var_i
                self.to_nodeB_Lr = node[node_key[self.to_node]].nodeB_dual_eq_var_r
                self.to_nodeB_Li = node[node_key[self.to_node]].nodeB_dual_eq_var_i
            if self.phases & 0x4 == 4:  # Check for phase C
                self.from_nodeC_Lr = node[node_key[self.from_node]].nodeC_dual_eq_var_r
                self.from_nodeC_Li = node[node_key[self.from_node]].nodeC_dual_eq_var_i
                self.to_nodeC_Lr = node[node_key[self.to_node]].nodeC_dual_eq_var_r
                self.to_nodeC_Li = node[node_key[self.to_node]].nodeC_dual_eq_var_i
    
            self.from_nodeN_Lr = node[node_key[self.from_node]].nodeN_dual_eq_var_r
            self.from_nodeN_Li = node[node_key[self.from_node]].nodeN_dual_eq_var_i
            self.to_nodeN_Lr = node[node_key[self.to_node]].nodeN_dual_eq_var_r
            self.to_nodeN_Li = node[node_key[self.to_node]].nodeN_dual_eq_var_i
    
    def stamp_impedance(self, node_Vr_from, node_Vi_from, node_Vr_to,
                        node_Vi_to, Ylin_val, Ylin_row, Ylin_col, idx_Y, phase):
        G = self.G[phase]
        B = self.B[phase]
        if G:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from, G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_from, node_Vr_to, -G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_to, node_Vr_to, G, Ylin_val, Ylin_row, Ylin_col,
                                                               idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_to, node_Vr_from, -G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from, G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_from, node_Vi_to, -G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_to, node_Vi_to, G, Ylin_val, Ylin_row, Ylin_col,
                                                               idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_to, node_Vi_from, -G, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
        if B:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from, -B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_from, node_Vi_to, B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from, B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_from, node_Vr_to, -B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_to, node_Vi_to, -B, Ylin_val, Ylin_row, Ylin_col,
                                                               idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_to, node_Vi_from, B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_to, node_Vr_to, B, Ylin_val, Ylin_row, Ylin_col,
                                                               idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_to, node_Vr_from, -B, Ylin_val, Ylin_row,
                                                               Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_impedance_dual(self, node_Lr_from, node_Li_from, node_Lr_to,
                        node_Li_to, Ylin_val, Ylin_row, Ylin_col, idx_Y, phase):
        G = self.G[phase]
        B = self.B[phase]
        if G:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_from, node_Lr_from, G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_to, node_Lr_from, -G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_to, node_Lr_to, G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                 node_Lr_from, node_Lr_to,-G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_from, node_Li_from, G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_to, node_Li_from, -G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_to, node_Li_to, G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_from, node_Li_to, -G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
        if B:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_from,node_Lr_from,  -B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_to, node_Lr_from, B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_from, node_Li_from,B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_to, node_Li_from,-B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_to, node_Lr_to, -B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Li_from, node_Lr_to, B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_to, node_Li_to, B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Lr_from, node_Li_to, -B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_linear(self, Ylin_val, Ylin_row, Ylin_col, idx_Y):
        if self.phases & 0x1 == 1:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance(
                self.from_nodeA_Vr, self.from_nodeA_Vi, self.to_nodeA_Vr,
                self.to_nodeA_Vi, Ylin_val, Ylin_row, Ylin_col, idx_Y, 0)
        if self.phases & 0x2 == 2:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance(
                self.from_nodeB_Vr, self.from_nodeB_Vi, self.to_nodeB_Vr,
                self.to_nodeB_Vi, Ylin_val, Ylin_row, Ylin_col, idx_Y, 1)
        if self.phases & 0x4 == 4:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance(
                self.from_nodeC_Vr, self.from_nodeC_Vi, self.to_nodeC_Vr,
                self.to_nodeC_Vi, Ylin_val, Ylin_row, Ylin_col, idx_Y, 2)

        if self.phases & 0x8 == 8:
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance(
                self.from_nodeN_Vr, self.from_nodeN_Vi, self.to_nodeN_Vr,
                self.to_nodeN_Vi, Ylin_val, Ylin_row, Ylin_col, idx_Y, 3)
            
        if self.stamp_dual:
            if self.phases & 0x1 == 1:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance_dual(
                    self.from_nodeA_Lr, self.from_nodeA_Li, self.to_nodeA_Lr,
                    self.to_nodeA_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y, 0)
            if self.phases & 0x2 == 2:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance_dual(
                    self.from_nodeB_Lr, self.from_nodeB_Li, self.to_nodeB_Lr,
                    self.to_nodeB_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y, 1)
            if self.phases & 0x4 == 4:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance_dual(
                    self.from_nodeC_Lr, self.from_nodeC_Li, self.to_nodeC_Lr,
                    self.to_nodeC_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y, 2)
    
            if self.phases & 0x8 == 8:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_impedance_dual(
                    self.from_nodeN_Lr, self.from_nodeN_Li, self.to_nodeN_Lr,
                    self.to_nodeN_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y, 3)
        return Ylin_val, Ylin_row, Ylin_col, idx_Y
