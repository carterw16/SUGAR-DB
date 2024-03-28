"""Capacitors are objects used for reactive power compensation and voltage support.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 04-11-2017
  Updated Date: 10-19-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Development

  The capacitors class creates and implements a set of switchable capacitors. The phase value of each capacitor is
  in VARs.

  Typical usage example:

  capacitor = Capacitors()
"""

from __future__ import division

from classes.Elements import LinearElement
from classes.GlobalVars import connected_nodes
from classes.Nodes import Nodes


class Capacitors(LinearElement):

    def __init__(self, name, parent, ID, phases, nominal_voltage, connect_type,
                 capacitorA, capacitorB, capacitorC, pt_phase, delay,
                 upper_limit, lower_limit, stamp_dual = False):
        super(Capacitors, self).__init__()
        self.name = name
        if parent:
            self.ID = parent  # Load bus number
        else:
            self.ID = ID  # Load bus number

        self.phases = phases
        self.capacitor_phases = [0x01, 0x02, 0x04]
        self.connect_type = connect_type  # WYE (0) or DELTA (1)
        self.capacitorA = float(capacitorA) if capacitorA else 0.0  # VARs
        self.capacitorB = float(capacitorB) if capacitorB else 0.0  # VARs
        self.capacitorC = float(capacitorC) if capacitorC else 0.0  # VARs
        # Other parameters
        self.nominal_voltage = nominal_voltage
        self.pt_phase = pt_phase
        # Upper and Lower Limits could be V or VAr. Not specifying units right now.
        self.upper_limit = upper_limit
        self.lower_limit = lower_limit
        self.delay = delay
        (self.Ba, self.Bb,
         self.Bc) = map(self.calcB,
                        [self.capacitorA, self.capacitorB, self.capacitorC],
                        3 * [self.nominal_voltage])
        # Add to set connected nodes
        connected_nodes.add(self.ID)
        self.stamp_dual = stamp_dual

    @staticmethod
    def calcB(VARs, Vll):
        B = VARs / Vll**2
        return B

    def stamp_linear(self, node_key, node, Ylin_val, Ylin_row, Ylin_col, idx_Y):
        nodeA_Vr = node[node_key[self.ID]].nodeA_Vr
        nodeA_Vi = node[node_key[self.ID]].nodeA_Vi
        nodeB_Vr = node[node_key[self.ID]].nodeB_Vr
        nodeB_Vi = node[node_key[self.ID]].nodeB_Vi
        nodeC_Vr = node[node_key[self.ID]].nodeC_Vr
        nodeC_Vi = node[node_key[self.ID]].nodeC_Vi
        
        if self.stamp_dual:
            nodeA_Lr = node[node_key[self.ID]].nodeA_dual_eq_var_r
            nodeA_Li = node[node_key[self.ID]].nodeA_dual_eq_var_i
            nodeB_Lr = node[node_key[self.ID]].nodeB_dual_eq_var_r
            nodeB_Li = node[node_key[self.ID]].nodeB_dual_eq_var_i
            nodeC_Lr = node[node_key[self.ID]].nodeC_dual_eq_var_r
            nodeC_Li = node[node_key[self.ID]].nodeC_dual_eq_var_i

        # For WYE-Connection -- Neutral Grounded
        if self.connect_type == 0:  # Wye Connected
            _B = [self.Ba, self.Bb, self.Bc]
            _nodeVr = [nodeA_Vr, nodeB_Vr, nodeC_Vr]
            _nodeVi = [nodeA_Vi, nodeB_Vi, nodeC_Vi]
           
            if self.stamp_dual:
                _nodeLr = [nodeA_Lr, nodeB_Lr, nodeC_Lr]
                _nodeLi = [nodeA_Li, nodeB_Li, nodeC_Li]
            
            for _phase_idx in range(3):
                if self.capacitor_phases[_phase_idx] & self.phases == int(
                        self.capacitor_phases[_phase_idx]):
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVr[_phase_idx], _nodeVi[_phase_idx],
                                                                          -_B[_phase_idx], Ylin_val, Ylin_row, 
                                                                          Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVi[_phase_idx], _nodeVr[_phase_idx],
                                                                          _B[_phase_idx], Ylin_val, Ylin_row, 
                                                                          Ylin_col, idx_Y)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLi[_phase_idx], _nodeLr[_phase_idx], 
                                                                              -_B[_phase_idx], Ylin_val, Ylin_row, 
                                                                              Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLr[_phase_idx], _nodeLi[_phase_idx],
                                                                              _B[_phase_idx], Ylin_val, Ylin_row, 
                                                                              Ylin_col, idx_Y)   

        # Delta connected
        elif self.connect_type == 1:
            _B = [self.Ba, self.Bb, self.Bc]
            _nodeVr_from = [nodeA_Vr, nodeB_Vr, nodeC_Vr]
            _nodeVi_from = [nodeA_Vi, nodeB_Vi, nodeC_Vi]
            _nodeVr_to = [nodeB_Vr, nodeC_Vr, nodeA_Vr]
            _nodeVi_to = [nodeB_Vi, nodeC_Vi, nodeA_Vi]
            
            if self.stamp_dual:
                _nodeLr_from = [nodeA_Lr, nodeB_Lr, nodeC_Lr]
                _nodeLi_from = [nodeA_Li, nodeB_Li, nodeC_Li]
                _nodeLr_to = [nodeB_Lr, nodeC_Lr, nodeA_Lr]
                _nodeLi_to = [nodeB_Li, nodeC_Li, nodeA_Li]
                
            for _phase_idx in range(3):
                if self.capacitor_phases[_phase_idx] & self.phases == int(
                        self.capacitor_phases[_phase_idx]):
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVr_from[_phase_idx], 
                                                                          _nodeVi_from[_phase_idx], -_B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVr_from[_phase_idx], 
                                                                          _nodeVi_to[_phase_idx], _B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVi_from[_phase_idx], 
                                                                          _nodeVr_from[_phase_idx], _B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVi_from[_phase_idx],
                                                                          _nodeVr_to[_phase_idx], -_B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVr_to[_phase_idx], 
                                                                          _nodeVi_to[_phase_idx], -_B[_phase_idx],
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVr_to[_phase_idx], 
                                                                          _nodeVi_from[_phase_idx], _B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeVi_to[_phase_idx], 
                                                                          _nodeVr_to[_phase_idx], _B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y( _nodeVi_to[_phase_idx], 
                                                                          _nodeVr_from[_phase_idx], -_B[_phase_idx], 
                                                                          Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLi_from[_phase_idx], _nodeLr_from[_phase_idx],
                                                                              -_B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLi_to[_phase_idx], _nodeLr_from[_phase_idx],
                                                                              _B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLr_from[_phase_idx], _nodeLi_from[_phase_idx],
                                                                              _B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLr_to[_phase_idx], _nodeLi_from[_phase_idx],
                                                                              -_B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLi_to[_phase_idx], _nodeLr_to[_phase_idx],
                                                                              -_B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLi_from[_phase_idx], _nodeLr_to[_phase_idx],
                                                                              _B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLr_to[_phase_idx], _nodeLi_to[_phase_idx],
                                                                              _B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = Capacitors.stamp_Y(_nodeLr_from[_phase_idx], _nodeLi_to[_phase_idx],
                                                                              -_B[_phase_idx], Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        
        return Ylin_val, Ylin_row, Ylin_col, idx_Y
