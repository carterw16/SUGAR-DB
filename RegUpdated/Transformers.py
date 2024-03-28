"""Power delivery element that steps-up or steps-down voltage.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-10-2017
  Updated Date: 07-20-2020
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

  A power deliver element with two or more windings that steps-up or steps-down voltage in power system. Types of
  transformers currently implemented include Wye-Wye, Gwye-Gwye, Delta-Delta, Delta-Gwye, Single Phase, Center-Tap,
  Delta-Wye, Wye-Delta, and Gwye-Delta. Phase shifting transformers have not yet been implemented.

"""

from __future__ import division

from itertools import count
from types import SimpleNamespace

import numpy as np
from termcolor import colored

from classes.Elements import LinearElement
from classes.GlobalVars import _XFMR
from classes.Nodes import Nodes


class Transformer(LinearElement):
    _ids = count(0)

    # Create variable names
    _WYE_WYE = -1
    _GWYE_GWYE = 1
    _DELTA_DELTA = 2
    _DELTA_GWYE = 3
    _SINGLE_PHASE = 4
    _CENTER_TAP = 5
    # _GWYE_GWYE = 6
    _DELTA_WYE = 7
    _WYE_DELTA = 8
    _GWYE_DELTA = 9
    # Phase name
    phaseA = 0x01
    phaseB = 0x02
    phaseC = 0x04

    @staticmethod
    def sind(x):
        return np.sin(float(x) * np.pi / 180)

    @staticmethod
    def cosd(x):
        return np.cos(float(x) * np.pi / 180)

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

    def __init__(self,
                 node_index_,
                 node,
                 name,
                 ID,
                 phases,
                 from_node,
                 to_node,
                 nominal_voltage,
                 connect_type,
                 install_type,
                 power_rating,
                 primary_voltage,
                 secondary_voltage,
                 resistance,
                 reactance,
                 shunt_impedance,
                 is_center_tap=False,
                 is_grounded=False,
                 full_load_loss=None,
                 no_load_loss=None,
                 phase_shift=False,
                 stamp_dual = False):

        super(Transformer, self).__init__()

        self.name = name
        self.id = self._ids.__next__()
        self.ID = ID
        self.phases = int(phases)
        self.from_node = from_node
        self.to_node = to_node
        self.nominal_voltage = float(
            nominal_voltage) if nominal_voltage else None
        self.connect_type = connect_type
        self.install_type = install_type
        self.primary_voltage = primary_voltage
        self.secondary_voltage = secondary_voltage
        self.phase_shift = float(phase_shift) if phase_shift else None
        self.is_center_tap = is_center_tap
        self.is_grounded = is_grounded
        self.full_load_loss = float(full_load_loss) if full_load_loss else None
        self.no_load_loss = float(no_load_loss) if no_load_loss else None
        self.shunt_impedance = shunt_impedance if shunt_impedance else None
        self.stamp_dual = stamp_dual

        # Calculated Parameters
        if self.connect_type == _XFMR['WYE_WYE'] or self.connect_type == _XFMR[
                'GWYE_GWYE']:  # Wye-Wye or Gwye-Gwye
            self.Vpri = self.primary_voltage
            self.Vsec = self.secondary_voltage
            self.ang = None
        elif self.connect_type == _XFMR[
                'DELTA_DELTA']:  # If primary voltage is given line to line
            self.Vpri = self.primary_voltage * np.sqrt(3)
            self.Vsec = self.secondary_voltage * np.sqrt(3)
            self.ang = None
        elif self.connect_type == _XFMR[
                'DELTA_GWYE']:  # Delta - wye or Delta - Gwye
            self.Vpri = self.primary_voltage * np.sqrt(3)
            self.Vsec = self.secondary_voltage
            self.ang = None
        elif self.connect_type == _XFMR['WYE_DELTA']:  # Wye-Delta or Gwye-Delta
            self.Vpri = self.primary_voltage
            self.Vsec = self.secondary_voltage * np.sqrt(3)
            self.ang = None
        elif self.connect_type == _XFMR[
                'SINGLE_PHASE_CENTER_TAPPED']:  # center-tap transformer
            self.Vpri = self.primary_voltage
            self.Vsec = self.secondary_voltage
            self.ang = None
        else:
            print('Invalid Transformer Type')

        # Assign Nodes
        self.nodeVr_meas = {'A': 0, 'B': 0, 'C': 0}
        self.nodeVi_meas = {'A': 0, 'B': 0, 'C': 0}
        if self.connect_type not in [self._SINGLE_PHASE, self._CENTER_TAP]:
            self.xfmr_type = 'Three-Phase'
            self.isSinglePhase = False
            self.assign_nodes_3phase(node_index_, node, stamp_dual)
        else:
            self.xfmr_type = 'Center Tap'
            self.isSinglePhase = True
            self.assign_nodes_centertap(node_index_, node)

        # Calculate turns ratio r and x (For single phase do not use sqrt(3))
        self.tr = self.Vpri / self.Vsec
        if self.tr > 1:
            self.is_stepdown = True
        else:
            self.is_stepdown = False

        if self.is_center_tap and len(power_rating) > 1:
            self._r = resistance[0] / 0.5
            self._x = reactance[0] / 0.8
            if self.phases & 0x1 == 1:  # Check for phase A
                self.power_rating = power_rating[0]
            elif self.phases & 0x2 == 2:
                self.power_rating = power_rating[1]
            else:
                self.power_rating = power_rating[2]
        elif self.is_center_tap and len(power_rating) == 1:
            self._r = resistance[0] / 0.5
            self._x = reactance[0] / 0.8
            self.power_rating = power_rating[0]
        else:
            self._r = resistance
            self._x = reactance
            self.power_rating = power_rating

        self.r = (self._r * self.Vsec**2) / self.power_rating
        self.x = (self._x * self.Vsec**2) / self.power_rating

        if self.connect_type == 5:  # Reference Kersting (2007)
            _ZbaseH = self.Vpri**2 / self.power_rating
            _ZbaseL = self.Vsec**2 / self.power_rating
            # Correct this if impedance, impedance1 and impedance2 is given instead
            self._r0 = resistance[0]
            self._r1 = resistance[1]
            self._r2 = resistance[2]
            self._x0 = reactance[0]
            self._x1 = reactance[1]
            self._x2 = reactance[2]
            # Convert the impedance to ohms
            self.r0 = self._r0 * _ZbaseH
            self.x0 = self._x0 * _ZbaseH
            self.r1 = self._r1 * _ZbaseL
            self.x1 = self._x1 * _ZbaseL
            self.r2 = self._r2 * _ZbaseL
            self.x2 = self._x2 * _ZbaseL

            (self.G0, self.B0) = self.calc_G_B(self.r0, self.x0)
            (self.G1, self.B1) = self.calc_G_B(self.r1, self.x1)
            (self.G2, self.B2) = self.calc_G_B(self.r2, self.x2)

        (self.G, self.B) = self.calc_G_B(self.r, self.x)

        if self.shunt_impedance:
            if self.connect_type != 5:
                self.xShunt = (self.shunt_impedance *
                               self.Vsec**2) / self.power_rating
            else:
                self.xShunt = (self.shunt_impedance *
                               self.Vpri**2) / self.power_rating
            self.hasShunt = True
        else:
            self.hasShunt = False
        if self.hasShunt:
            (self.Gshunt, self.Bshunt) = self.calc_G_B(self.xShunt.real,
                                                       self.xShunt.imag)
            self.hasShunt = True
        else:
            self.hasShunt = False

    def assign_nodes_3phase(self, node_index_, node, stamp_dual):
        # There are likely to be 12 additional nodes for the Transformer
        # 6 additional rows for the voltage source eqns. (No Phase Shifters
        # 12 additional rows for the voltage sources eqn (w Phase Shifters
        self.nodeA_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vr
        self.nodeA_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vi
        self.nodeB_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vr
        self.nodeB_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vi
        self.nodeC_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vr
        self.nodeC_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vi
        self.nodeN_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vr
        self.nodeN_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vi
        self.nodeA_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vr
        self.nodeA_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vi
        self.nodeB_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vr
        self.nodeB_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vi
        self.nodeC_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vr
        self.nodeC_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vi
        self.nodeN_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vr
        self.nodeN_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vi
        
        if stamp_dual:
            self.nodeA_LR_from = node[Nodes.nodeKey[self.from_node]].nodeA_Lr
            self.nodeA_LI_from = node[Nodes.nodeKey[self.from_node]].nodeA_Li
            self.nodeB_LR_from = node[Nodes.nodeKey[self.from_node]].nodeB_Lr
            self.nodeB_LI_from = node[Nodes.nodeKey[self.from_node]].nodeB_Li
            self.nodeC_LR_from = node[Nodes.nodeKey[self.from_node]].nodeC_Lr
            self.nodeC_LI_from = node[Nodes.nodeKey[self.from_node]].nodeC_Li
            self.nodeN_LR_from = node[Nodes.nodeKey[self.from_node]].nodeN_Lr
            self.nodeN_LI_from = node[Nodes.nodeKey[self.from_node]].nodeN_Li
            self.nodeA_LR_to = node[Nodes.nodeKey[self.to_node]].nodeA_Lr
            self.nodeA_LI_to = node[Nodes.nodeKey[self.to_node]].nodeA_Li
            self.nodeB_LR_to = node[Nodes.nodeKey[self.to_node]].nodeB_Lr
            self.nodeB_LI_to = node[Nodes.nodeKey[self.to_node]].nodeB_Li
            self.nodeC_LR_to = node[Nodes.nodeKey[self.to_node]].nodeC_Lr
            self.nodeC_LI_to = node[Nodes.nodeKey[self.to_node]].nodeC_Li
            self.nodeN_LR_to = node[Nodes.nodeKey[self.to_node]].nodeN_Lr
            self.nodeN_LI_to = node[Nodes.nodeKey[self.to_node]].nodeN_Li
        #############################################################
        # Add ground voltage source if needed
        if self.connect_type in [
                self._DELTA_DELTA, self._WYE_WYE, self._DELTA_GWYE,
                self._GWYE_GWYE, self._DELTA_WYE, self._WYE_DELTA,
                self._GWYE_DELTA
        ]:
            # Add for the three phase extra node in the secondary
            # circuit (for both real and imaginary)
            self.nodeXtraA_rs = node_index_.__next__()
            self.nodeXtraA_is = node_index_.__next__()
            self.nodeXtraB_rs = node_index_.__next__()
            self.nodeXtraB_is = node_index_.__next__()
            self.nodeXtraC_rs = node_index_.__next__()
            self.nodeXtraC_is = node_index_.__next__()

            # Add for the three phases extra node for the voltage sources
            self.nodeVsA_rp = node_index_.__next__()
            self.nodeVsA_ip = node_index_.__next__()
            self.nodeVsB_rp = node_index_.__next__()
            self.nodeVsB_ip = node_index_.__next__()
            self.nodeVsC_rp = node_index_.__next__()
            self.nodeVsC_ip = node_index_.__next__()
            
            if stamp_dual:
                self.nodeXtraA_Lrs = node_index_.__next__()
                self.nodeXtraA_Lis = node_index_.__next__()
                self.nodeXtraB_Lrs = node_index_.__next__()
                self.nodeXtraB_Lis = node_index_.__next__()
                self.nodeXtraC_Lrs = node_index_.__next__()
                self.nodeXtraC_Lis = node_index_.__next__()

                self.nodeA_Lrp = node_index_.__next__()
                self.nodeA_Lip = node_index_.__next__()
                self.nodeB_Lrp = node_index_.__next__()
                self.nodeB_Lip = node_index_.__next__()
                self.nodeC_Lrp = node_index_.__next__()
                self.nodeC_Lip = node_index_.__next__()
                
                self.LRnode_meas = {
                    'A': self.nodeA_Lrp,
                    'B': self.nodeB_Lrp,
                    'C': self.nodeC_Lrp
                }
                self.LInode_meas = {
                    'A': self.nodeA_Lip,
                    'B': self.nodeB_Lip,
                    'C': self.nodeC_Lip
                }
                

            # Add voltage source nodes to the meas dictionary
            self.nodeVr_meas = {
                'A': self.nodeVsA_rp,
                'B': self.nodeVsB_rp,
                'C': self.nodeVsC_rp
            }
            self.nodeVi_meas = {
                'A': self.nodeVsA_ip,
                'B': self.nodeVsB_ip,
                'C': self.nodeVsC_ip
            }

            # Add additional voltage sources for DELTA-WYE and WYE-DELTA
            if self.ang:
                self.nodeVsA_rp_2 = node_index_.__next__()
                self.nodeVsA_ip_2 = node_index_.__next__()
                self.nodeVsB_rp_2 = node_index_.__next__()
                self.nodeVsB_ip_2 = node_index_.__next__()
                self.nodeVsC_rp_2 = node_index_.__next__()
                self.nodeVsC_ip_2 = node_index_.__next__()
                
                if stamp_dual:
                    self.nodeA_Lrp_2 = node_index_.__next__()
                    self.nodeA_Lip_2 = node_index_.__next__()
                    self.nodeB_Lrp_2 = node_index_.__next__()
                    self.nodeB_Lip_2 = node_index_.__next__()
                    self.nodeC_Lrp_2 = node_index_.__next__()
                    self.nodeC_Lip_2 = node_index_.__next__()

            # Check for ground and add further nodes
            if self.connect_type in [self._GWYE_DELTA, self._GWYE_GWYE]:
                self.nodeVGnd_rp = node_index_.__next__()
                self.nodeVGnd_ip = node_index_.__next__()
                
                if stamp_dual:
                    self.nodeGnd_Lrp = node_index_.__next__()
                    self.nodeGnd_Lip = node_index_.__next__()
            if self.connect_type in [self._DELTA_GWYE, self._GWYE_GWYE]:
                self.nodeVGnd_rs = node_index_.__next__()
                self.nodeVGnd_is = node_index_.__next__()
                
                if stamp_dual:
                    self.nodeGnd_Lrs = node_index_.__next__()
                    self.nodeGnd_Lis = node_index_.__next__()
        else:
            print('Incorrect transformer type for assigning nodes')
        return None

    # A separate function to assign node for center-tap transformer
    # from node should be AN, BN (regular node)
    # to node should be triplex node (1, 2, N)
    def assign_nodes_centertap(self, node_index_, node):
        # First find the Vr and Vi of the single phase that is connected
        if self.phases & 0x01 == 1:  # Phase A
            self.node1_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vr
            self.node1_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vi
        elif self.phases & 0x02 == 2:  # Phase B
            self.node1_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vr
            self.node1_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vi
        elif self.phases & 0x04 == 4:  # Phase C
            self.node1_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vr
            self.node1_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vi
        self.nodeN_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vr
        self.nodeN_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vi
        # Secondary Hot Wires
        self.to_node1_Vr = node[Nodes.nodeKey[self.to_node]].node1_Vr
        self.to_node1_Vi = node[Nodes.nodeKey[self.to_node]].node1_Vi
        self.to_node2_Vr = node[Nodes.nodeKey[self.to_node]].node2_Vr
        self.to_node2_Vi = node[Nodes.nodeKey[self.to_node]].node2_Vi
        self.nodeN_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vr
        self.nodeN_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vi
        
        
        if self.stamp_dual:
            if self.phases & 0x01 == 1:  # Phase A
                self.node1_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Lr
                self.node1_Li_from = node[Nodes.nodeKey[self.from_node]].nodeA_Li
            elif self.phases & 0x02 == 2:  # Phase B
                self.node1_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Lr
                self.node1_Li_from = node[Nodes.nodeKey[self.from_node]].nodeB_Li
            elif self.phases & 0x04 == 4:  # Phase C
                self.node1_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Lr
                self.node1_Li_from = node[Nodes.nodeKey[self.from_node]].nodeC_Li
            self.nodeN_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeN_Lr
            self.nodeN_Li_from = node[Nodes.nodeKey[self.from_node]].nodeN_Li
            # Secondary Hot Wires
            self.to_node1_Lr = node[Nodes.nodeKey[self.to_node]].node1_Lr
            self.to_node1_Li = node[Nodes.nodeKey[self.to_node]].node1_Li
            self.to_node2_Lr = node[Nodes.nodeKey[self.to_node]].node2_Lr
            self.to_node2_Li = node[Nodes.nodeKey[self.to_node]].node2_Li
            self.nodeN_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeN_Lr
            self.nodeN_Li_to = node[Nodes.nodeKey[self.to_node]].nodeN_Li           
        # Assign new nodes
        try:
            if self.connect_type == self._CENTER_TAP:
                if not self.phase_shift:
                    # 2 additional nodes on the primary side
                    self.nodeXtraP_rp = node_index_.__next__()
                    self.nodeXtraP_ip = node_index_.__next__()
                    # 4 additional extra nodes on sending side
                    self.nodeXtra1_rs = node_index_.__next__()
                    self.nodeXtra1_is = node_index_.__next__()
                    self.nodeXtra2_rs = node_index_.__next__()
                    self.nodeXtra2_is = node_index_.__next__()

                    # 4 additional rows for voltage sources
                    self.nodeVs1_rs = node_index_.__next__()
                    self.nodeVs1_is = node_index_.__next__()
                    self.nodeVs2_rs = node_index_.__next__()
                    self.nodeVs2_is = node_index_.__next__()

                    # 2 additional variables from grounding neutrals
                    self.gnd_real_Vsec = node_index_.__next__()
                    self.gnd_imag_Vsec = node_index_.__next__()
                    self.gnd_real_Vpri = node_index_.__next__()
                    self.gnd_imag_Vpri = node_index_.__next__()
                    
                    if self.stamp_dual:
                        # 2 additional nodes on the primary side
                        self.nodeLrXtraP_rp = node_index_.__next__()
                        self.nodeLiXtraP_ip = node_index_.__next__()
                        # 4 additional extra nodes on sending side
                        self.nodeLrXtra1_rs = node_index_.__next__()
                        self.nodeLiXtra1_is = node_index_.__next__()
                        self.nodeLrXtra2_rs = node_index_.__next__()
                        self.nodeLiXtra2_is = node_index_.__next__()
    
                        # 4 additional rows for voltage sources
                        self.nodeLr1_rs = node_index_.__next__()
                        self.nodeLi1_is = node_index_.__next__()
                        self.nodeLr2_rs = node_index_.__next__()
                        self.nodeLi2_is = node_index_.__next__()
    
                        # 2 additional variables from grounding neutrals
                        self.gnd_real_Lrsec = node_index_.__next__()
                        self.gnd_imag_Lisec = node_index_.__next__()
                        self.gnd_real_Lrpri = node_index_.__next__()
                        self.gnd_imag_Lipri = node_index_.__next__()   
                else:
                    pass
            else:
                raise Exception
        except Exception:
            print('Undefined Transformer Type')
        return None

    def getNodes(self, node, V, stamp_dual):
        """
		Find the indices of the solution vector and the values at those indices.
		:param node: The vector of all system nodes.
		:param V: The solution vector.
		:return: Three dictionaries which hold the phase and sequence nodes.
		"""
        # Find the from bus nodes to stamp
        nodeA_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vr
        nodeA_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vi
        nodeB_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vr
        nodeB_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vi
        nodeC_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vr
        nodeC_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vi

        # Find the to bus nodes to stamp
        if not self.isSinglePhase:
            nodeA_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vr
            nodeA_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vi
            nodeB_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vr
            nodeB_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vi
            nodeC_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vr
            nodeC_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vi

            (Var_to, Vai_to, Vbr_to, Vbi_to, Vcr_to, Vci_to) = \
             (V[nodeA_Vr_to], V[nodeA_Vi_to], V[nodeB_Vr_to], V[nodeB_Vi_to], V[nodeC_Vr_to],
              V[nodeC_Vi_to])
             
            if stamp_dual:
                nodeA_LR_from = node[Nodes.nodeKey[self.from_node]].nodeA_Lr
                nodeA_LI_from = node[Nodes.nodeKey[self.from_node]].nodeA_Li
                nodeB_LR_from = node[Nodes.nodeKey[self.from_node]].nodeB_Lr
                nodeB_LI_from = node[Nodes.nodeKey[self.from_node]].nodeB_Li
                nodeC_LR_from = node[Nodes.nodeKey[self.from_node]].nodeC_Lr
                nodeC_LI_from = node[Nodes.nodeKey[self.from_node]].nodeC_Li
                
                nodeA_LR_to = node[Nodes.nodeKey[self.to_node]].nodeA_Lr
                nodeA_LI_to = node[Nodes.nodeKey[self.to_node]].nodeA_Li
                nodeB_LR_to = node[Nodes.nodeKey[self.to_node]].nodeB_Lr
                nodeB_LI_to = node[Nodes.nodeKey[self.to_node]].nodeB_Li
                nodeC_LR_to = node[Nodes.nodeKey[self.to_node]].nodeC_Lr
                nodeC_LI_to = node[Nodes.nodeKey[self.to_node]].nodeC_Li
              
               
                lam_from = {'LR': [V[nodeA_LR_from], V[nodeB_LR_from], V[nodeC_LR_from]],
                            'LRnode': [nodeA_LR_from, nodeB_LR_from, nodeC_LR_from],
                            'LI': [V[nodeA_LI_from], V[nodeB_LI_from], V[nodeC_LI_from]],
                            'LInode': [nodeA_LI_from, nodeB_LI_from, nodeC_LI_from]}
                
                lam_to = {'LR': [V[nodeA_LR_to], V[nodeB_LR_to], V[nodeC_LR_to]],
                            'LRnode': [nodeA_LR_to, nodeB_LR_to, nodeC_LR_to],
                            'LI': [V[nodeA_LI_to], V[nodeB_LI_to], V[nodeC_LI_to]],
                            'LInode': [nodeA_LI_to, nodeB_LI_to, nodeC_LI_to]}

            else:
                lam_from = []
                lam_to = []

            Vabc_to = {
                'VR': [Var_to, Vbr_to, Vcr_to],
                'VI': [Vai_to, Vbi_to, Vci_to],
                'nodeVR': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to],
                'nodeVI': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]
            }

        else:
            node1_Vr = node[Nodes.nodeKey[self.to_node]].node1_Vr
            node1_Vi = node[Nodes.nodeKey[self.to_node]].node1_Vi
            node2_Vr = node[Nodes.nodeKey[self.to_node]].node2_Vr
            node2_Vi = node[Nodes.nodeKey[self.to_node]].node2_Vi
            nodeN_Vr = node[Nodes.nodeKey[self.to_node]].nodeN_Vr
            nodeN_Vi = node[Nodes.nodeKey[self.to_node]].nodeN_Vi

            (V1r, V1i, V2r, V2i, VNr, VNi) = \
             (V[node1_Vr], V[node1_Vi], V[node2_Vr], V[node2_Vi], V[nodeN_Vr],
              V[nodeN_Vi])

            Vabc_to = {
                'VR': [V1r, V2r, VNr],
                'VI': [V1i, V2i, VNi],
                'nodeVR': [node1_Vr, node2_Vr, nodeN_Vr],
                'nodeVI': [node1_Vi, node2_Vi, nodeN_Vi]
            }
            
            lam_to = {'LR': None, 'LRnode': None, 'LI': None, 'LInode': None}
            lam_from = {'LR': None, 'LRnode': None, 'LI': None, 'LInode': None}

        # Find the node voltages and reactive power

        (Var_from, Vai_from, Vbr_from, Vbi_from, Vcr_from, Vci_from) = \
         (V[nodeA_Vr_from], V[nodeA_Vi_from], V[nodeB_Vr_from], V[nodeB_Vi_from], V[nodeC_Vr_from],
          V[nodeC_Vi_from])

        Vabc_from = {
            'VR': [Var_from, Vbr_from, Vcr_from],
            'VI': [Vai_from, Vbi_from, Vci_from],
            'nodeVR': [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from],
            'nodeVI': [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from]
        }

        Vabc_from = SimpleNamespace(**Vabc_from)
        Vabc_to = SimpleNamespace(**Vabc_to)
        lam_to = SimpleNamespace(**lam_to)
        lam_from = SimpleNamespace(**lam_from)

        return Vabc_from, Vabc_to, lam_from, lam_to

    def stampPrimary(self,
                     Ylin_val,
                     Ylin_row,
                     Ylin_col,
                     idx_Y,
                     from_pos_real,
                     from_pos_imag,
                     from_neg_real,
                     from_neg_imag,
                     to_pos_real,
                     to_pos_imag,
                     to_neg_real,
                     to_neg_imag,
                     to_extra_real,
                     to_extra_imag,
                     Vsource_real,
                     Vsource_imag,
                     from_mid_real=None,
                     from_mid_imag=None,
                     Vsource_real_2=None,
                     Vsource_imag_2=None,
                     dual_dict = None,
                     stamp_dual = 0):
        # ==============================================================================
        #         Phase ANGLE = 0
        # ==============================================================================

        if stamp_dual:
            LR_from_p = dual_dict['LR_from_p']
            LI_from_p = dual_dict['LI_from_p']
            LR_from_n = dual_dict['LR_from_n']
            LI_from_n = dual_dict['LI_from_n']
            LR_to_n = dual_dict['LR_to_n']
            LI_to_n = dual_dict['LI_to_n']
            LR_extra_node = dual_dict['LR_extra_node']
            LI_extra_node = dual_dict['LI_extra_node']
            LR_source = dual_dict['LR_source']
            LI_source = dual_dict['LI_source']
            
            
        # PRIMARY VOLTAGES
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_real, from_pos_real, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_imag, from_pos_imag, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_real, from_neg_real, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_imag, from_neg_imag, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_real, to_extra_real, -self.tr, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_imag, to_extra_imag, -self.tr, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_real, to_neg_real, self.tr, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            Vsource_imag, to_neg_imag, self.tr, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        # PRIMARY CURRENTS
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            from_pos_real, Vsource_real, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            from_pos_imag, Vsource_imag, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            from_neg_real, Vsource_real, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            from_neg_imag, Vsource_imag, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        # DUAL STAMP (TRANSPOSE OF REGULAR)
        if stamp_dual:
            # PRIMARY VOLTAGES
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_from_p, LR_source, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_from_p, LI_source, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_from_n, LR_source, -1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_from_n, LI_source, -1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_extra_node, LR_source, -self.tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_extra_node, LI_source, -self.tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_to_n, LR_source, self.tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_to_n, LI_source, self.tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            # PRIMARY CURRENTS
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_source, LR_from_p, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_source, LI_from_p, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LR_source, LR_from_n, -1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                LI_source, LI_from_n, -1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
        # ==============================================================================
        #         Phase ANGLE = 30 (WYE-DELTA and DELTA-WYE)
        # ==============================================================================
        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stampSecondary(self,
                       Ylin_val,
                       Ylin_row,
                       Ylin_col,
                       idx_Y,
                       from_pos_real,
                       from_pos_imag,
                       from_neg_real,
                       from_neg_imag,
                       to_pos_real,
                       to_pos_imag,
                       to_neg_real,
                       to_neg_imag,
                       to_extra_real,
                       to_extra_imag,
                       Vsource_real,
                       Vsource_imag,
                       homotopy_option,
                       B_factor=None,
                       G_factor=None,
                       dual_dict = None,
                       stamp_dual = 0):
        # Calculate the conductance, susceptance and transformer ration based on homotopy option
        if stamp_dual:
            LR_to_p = dual_dict['LR_to_p']
            LI_to_p = dual_dict['LI_to_p']
            LR_to_n = dual_dict['LR_to_n']
            LI_to_n = dual_dict['LI_to_n']
            LR_extra_node = dual_dict['LR_extra_node']
            LI_extra_node = dual_dict['LI_extra_node'] 
            LR_source = dual_dict['LR_source']
            LI_source = dual_dict['LI_source']
            
        if not homotopy_option:
            _tr = self.tr
            _B = self.B
            _G = self.G
            if self.hasShunt:
                _Bshunt = self.Bshunt
                _Gshunt = self.Gshunt
        else:
            _tr = 0
            _B = B_factor * self.B
            _G = G_factor * self.G
            _Bshunt = 0
            _Gshunt = 0

        if _tr != 0:
            # ==============================================================================
            #         Phase ANGLE = 0
            # ==============================================================================
            # Secondary current sources
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_real, Vsource_real, -_tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_imag, Vsource_imag, -_tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_neg_real, Vsource_real, _tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_neg_imag, Vsource_imag, _tr, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            
            if stamp_dual:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_source, LR_extra_node, -_tr, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_source, LI_extra_node, -_tr, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_source, LR_to_n, _tr, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_source, LI_to_n, _tr, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)

        # Branch stamping for the susceptances
        if _B != 0:
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_real, to_pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_real, to_extra_imag, _B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_imag, to_pos_real, _B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_imag, to_extra_real, -_B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_real, to_extra_imag, -_B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_real, to_pos_imag, _B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_imag, to_extra_real, _B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_imag, to_pos_real, -_B, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)

            if stamp_dual:
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_to_p, LR_to_p, -_B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_extra_node, LR_to_p, _B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_to_p, LI_to_p, _B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_extra_node, LI_to_p, -_B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_extra_node, LR_extra_node, -_B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_to_p, LR_extra_node, _B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_extra_node, LI_extra_node, _B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_to_p, LI_extra_node, -_B, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
        # Adding resistance if existent
        if _G != 0:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_real, to_pos_real, _G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_real, to_extra_real, -_G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_imag, to_pos_imag, _G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_pos_imag, to_extra_imag, -_G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_real, to_extra_real, _G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_real, to_pos_real, -_G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_imag, to_extra_imag, _G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                to_extra_imag, to_pos_imag, -_G, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            if stamp_dual:
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_to_p, LR_to_p, _G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_extra_node, LR_to_p, -_G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_to_p, LI_to_p, _G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_extra_node, LI_to_p, -_G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_extra_node, LR_extra_node, _G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_to_p, LR_extra_node, -_G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_extra_node, LI_extra_node, _G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_to_p, LI_extra_node, -_G, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
        # ==============================================================================
        #         SHUNTS MISSING
        # ==============================================================================
        if self.hasShunt:
            # Stamp shunt on the secondary of the transformer
            if _Bshunt:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    to_pos_real, to_pos_imag, -_Bshunt, Ylin_val, Ylin_row,
                    Ylin_col, idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    to_pos_imag, to_pos_real, _Bshunt, Ylin_val, Ylin_row,
                    Ylin_col, idx_Y)
                
                if stamp_dual:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        LR_to_p, LI_to_p, _Bshunt, Ylin_val, Ylin_row,
                        Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        LI_to_p, LR_to_p, -_Bshunt, Ylin_val, Ylin_row,
                        Ylin_col, idx_Y)

            if _Gshunt:
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    to_pos_real, to_pos_real, _Gshunt, Ylin_val, Ylin_row,
                    Ylin_col, idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    to_pos_imag, to_pos_imag, _Gshunt, Ylin_val, Ylin_row,
                    Ylin_col, idx_Y)
                
                if stamp_dual:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        LR_to_p, LR_to_p, _Gshunt, Ylin_val, Ylin_row,
                        Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        LI_to_p, LI_to_p, _Gshunt, Ylin_val, Ylin_row,
                        Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    # This function is only used for center-tap transformers right now
    def stampGB(self, _G, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y, pos_real,
                pos_imag, extra_real, extra_imag, dual_dict = None, stamp_dual = 0):
        
        # Branch stamping for the susceptances
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_real, extra_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_imag, extra_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_real, extra_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_real, pos_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_imag, extra_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        

        # Adding resistance if existent
        if _G != 0:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, extra_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, extra_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_real, extra_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_real, pos_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_imag, extra_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_imag, pos_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stampGB_dual(self, _G, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y, pos_real,
                pos_imag, extra_real, extra_imag, dual_dict = None, stamp_dual = 0):

        # Branch stamping for the susceptances
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_real, pos_imag,  _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_imag, extra_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_imag, extra_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            extra_real, extra_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            pos_real, extra_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        

        # Adding resistance if existent
        if _G != 0:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_real, pos_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_imag, pos_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_real, extra_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, extra_real, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                extra_imag, extra_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, extra_imag, -_G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    # This function is only used for center-tap transformers right now
    def stampGB_shunt(self, _G, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                      pos_real, pos_imag, dual_dict = None, stamp_dual = 0):
            
        if _B != 0:
            # Branch stamping for the susceptances
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_imag, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_real, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        # Adding resistance if existent
        if _G != 0:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stampGB_shunt_dual(self, _G, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                      pos_real, pos_imag, dual_dict = None, stamp_dual = 0):
            
        if _B != 0:
            # Branch stamping for the susceptances
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_real, -_B, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_imag, _B, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        # Adding resistance if existent
        if _G != 0:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_real, pos_real, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                pos_imag, pos_imag, _G, Ylin_val, Ylin_row, Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_ground_neutral(self, neg_real_node, neg_imag_node, gnd_real_V,
                             gnd_imag_V, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                             grounded_node, dual_dict = None, stamp_dual = 0):
        if grounded_node not in stamped_ground:
            # Stamp Vng = 0 + j0
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                gnd_real_V, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                gnd_imag_V, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            # Stamp Current Sources
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                neg_real_node, gnd_real_V, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                neg_imag_node, gnd_imag_V, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            stamped_ground.add(grounded_node)
            
            if stamp_dual:
                LR_n = dual_dict['LR_n']
                LI_n = dual_dict['LI_n']
                LR_gnd = dual_dict['LR_gnd']
                LI_gnd = dual_dict['LI_gnd']

                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_n, LR_gnd, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_n, LI_gnd, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                # Stamp Current Sources
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LR_gnd, LR_n, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                #
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                    LI_gnd, LI_n, 1, Ylin_val, Ylin_row, Ylin_col,
                    idx_Y)
                
        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_ground_neutral_dual(self, neg_real_node, neg_imag_node, gnd_real_Lr,
                              gnd_imag_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                              grounded_node):
        if grounded_node not in stamped_ground:
            # Stamp Vng = 0 + j0
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                neg_real_node, gnd_real_Lr, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                neg_imag_node, gnd_imag_Li, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            # Stamp Current Sources
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                gnd_real_Lr, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                gnd_imag_Li, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col,
                idx_Y)
            #
            stamped_ground.add(grounded_node)
                
        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_linear(self,
                     Ylin_val,
                     Ylin_row,
                     Ylin_col,
                     idx_Y,
                     stamped_ground,
                     enable_CM,
                     curr_meas,
                     homotopy_option=False,
                     h_factor=None,
                     G_homotopy=None,
                     B_homotopy=None,
                     stamp_dual = 0):

        if not self.isSinglePhase:
            # Find the from (Primary) bus nodes to stamp
            nodeA_Vr_from = self.nodeA_Vr_from
            nodeA_Vi_from = self.nodeA_Vi_from
            nodeB_Vr_from = self.nodeB_Vr_from
            nodeB_Vi_from = self.nodeB_Vi_from
            nodeC_Vr_from = self.nodeC_Vr_from
            nodeC_Vi_from = self.nodeC_Vi_from

            wyePrimarySet = {
                self._WYE_WYE, self._GWYE_GWYE, self._WYE_DELTA,
                self._GWYE_DELTA
            }
            if self.connect_type in wyePrimarySet:
                nodeN_Vr_from = self.nodeN_Vr_from
                nodeN_Vi_from = self.nodeN_Vi_from
                
                if stamp_dual:
                    nodeN_LR_from = self.nodeN_LR_from
                    nodeN_LI_from = self.nodeN_LI_from

            # Find the to (Secondary) bus nodes to stamp
            nodeA_Vr_to = self.nodeA_Vr_to
            nodeA_Vi_to = self.nodeA_Vi_to
            nodeB_Vr_to = self.nodeB_Vr_to
            nodeB_Vi_to = self.nodeB_Vi_to
            nodeC_Vr_to = self.nodeC_Vr_to
            nodeC_Vi_to = self.nodeC_Vi_to

            # Find the Xtra Secondary Nodes
            nodeA_Vr_extra = self.nodeXtraA_rs
            nodeB_Vr_extra = self.nodeXtraB_rs
            nodeC_Vr_extra = self.nodeXtraC_rs
            nodeA_Vi_extra = self.nodeXtraA_is
            nodeB_Vi_extra = self.nodeXtraB_is
            nodeC_Vi_extra = self.nodeXtraC_is

            nodeVsA_rp = self.nodeVsA_rp
            nodeVsB_rp = self.nodeVsB_rp
            nodeVsC_rp = self.nodeVsC_rp
            nodeVsA_ip = self.nodeVsA_ip
            nodeVsB_ip = self.nodeVsB_ip
            nodeVsC_ip = self.nodeVsC_ip

            wyeSecondarySet = {
                self._WYE_WYE, self._DELTA_GWYE, self._GWYE_GWYE,
                self._DELTA_WYE
            }
            if self.connect_type in wyeSecondarySet:
                nodeN_Vr_to = self.nodeN_Vr_to
                nodeN_Vi_to = self.nodeN_Vi_to
                
                if stamp_dual:
                    nodeN_LR_to = self.nodeN_LR_to
                    nodeN_LI_to = self.nodeN_LI_to
            
            if stamp_dual:
                nodeA_LR_from = self.nodeA_LR_from
                nodeA_LI_from = self.nodeA_LI_from
                nodeB_LR_from = self.nodeB_LR_from
                nodeB_LI_from = self.nodeB_LI_from
                nodeC_LR_from = self.nodeC_LR_from
                nodeC_LI_from = self.nodeC_LI_from
                
                nodeA_LR_to = self.nodeA_LR_to
                nodeA_LI_to = self.nodeA_LI_to
                nodeB_LR_to = self.nodeB_LR_to
                nodeB_LI_to = self.nodeB_LI_to
                nodeC_LR_to = self.nodeC_LR_to
                nodeC_LI_to = self.nodeC_LI_to
                
                nodeA_LR_extra = self.nodeXtraA_Lrs 
                nodeA_LI_extra = self.nodeXtraA_Lis 
                nodeB_LR_extra = self.nodeXtraB_Lrs 
                nodeB_LI_extra = self.nodeXtraB_Lis 
                nodeC_LR_extra = self.nodeXtraC_Lrs 
                nodeC_LI_extra = self.nodeXtraC_Lis 

                nodeA_LR_rp = self.nodeA_Lrp 
                nodeA_LI_ip = self.nodeA_Lip 
                nodeB_LR_rp = self.nodeB_Lrp 
                nodeB_LI_ip = self.nodeB_Lip 
                nodeC_LR_rp = self.nodeC_Lrp 
                nodeC_LI_ip = self.nodeC_Lip 

            # Put together the set of nodes
            nodeVr_from = [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from]
            nodeVi_from = [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from]
            nodeVr_to = [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to]
            nodeVi_to = [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]
            nodeVr_extra_to = [nodeA_Vr_extra, nodeB_Vr_extra, nodeC_Vr_extra]
            nodeVi_extra_to = [nodeA_Vi_extra, nodeB_Vi_extra, nodeC_Vi_extra]
            Vsource_Real = [nodeVsA_rp, nodeVsB_rp, nodeVsC_rp]
            Vsource_Imag = [nodeVsA_ip, nodeVsB_ip, nodeVsC_ip]
            
            if stamp_dual:
                nodeLR_from = [nodeA_LR_from, nodeB_LR_from, nodeC_LR_from]
                nodeLI_from = [nodeA_LI_from, nodeB_LI_from, nodeC_LI_from]
                nodeLR_to = [nodeA_LR_to, nodeB_LR_to, nodeC_LR_to]
                nodeLI_to = [nodeA_LI_to, nodeB_LI_to, nodeC_LI_to]
                nodeLR_extra_to = [nodeA_LR_extra, nodeB_LR_extra, nodeC_LR_extra]
                nodeLI_extra_to = [nodeA_LI_extra, nodeB_LI_extra, nodeC_LI_extra]
                nodeLR_source = [nodeA_LR_rp, nodeB_LR_rp, nodeC_LR_rp]
                nodeLI_source = [nodeA_LI_ip, nodeB_LI_ip, nodeC_LI_ip]

            # In future if there is phase angle
            if self.ang:
                Vsource2_Real = [
                    self.nodeVsA_rp_2, self.nodeVsB_rp_2, self.nodeVsC_rp_2
                ]
                Vsource2_Imag = [
                    self.nodeVsA_ip_2, self.nodeVsB_ip_2, self.nodeVsC_ip_2
                ]

            (_A, _B, _C) = (0, 1, 2)
            _N = 3

        elif self.isSinglePhase:
            pass

        try:
            # Individual phases of the transformer is implemented for this type
            # Option 1 - Stamp wye-wye
            if self.connect_type == self._WYE_WYE or self.connect_type == self._GWYE_GWYE:
                nodeVr_from.append(nodeN_Vr_from)
                nodeVi_from.append(nodeN_Vi_from)
                nodeVr_to.append(nodeN_Vr_to)
                nodeVi_to.append(nodeN_Vi_to)
                
                if stamp_dual:
                    nodeLR_from.append(nodeN_LR_from)
                    nodeLI_from.append(nodeN_LI_from)
                    nodeLR_to.append(nodeN_LR_to)
                    nodeLI_to.append(nodeN_LI_to)
                    
                # Third number in the set is used to find what phases are on
                phases = [(_A, _N, self.phaseA), (_B, _N, self.phaseB),
                          (_C, _N, self.phaseC)]
                for phaseSet in phases:
                    if phaseSet[2] & self.phases == phaseSet[2]:
                        # #########Stamp Primary Circuit##################
                        # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                        _pos = phaseSet[0]
                        _neg = phaseSet[1]
                        _winding = phaseSet[0]
                        if stamp_dual:
                            dual_dict = {
                                'LR_from_p': nodeLR_from[_pos], 
                                'LI_from_p': nodeLI_from[_pos],
                                'LR_from_n': nodeLR_from[_neg], 
                                'LI_from_n': nodeLI_from[_neg],
                                'LR_to_p': nodeLR_to[_pos], 
                                'LI_to_p': nodeLI_to[_pos],
                                'LR_to_n': nodeLR_to[_neg], 
                                'LI_to_n': nodeLI_to[_neg],
                                'LR_extra_node': nodeLR_extra_to[_winding],
                                'LI_extra_node': nodeLI_extra_to[_winding],
                                'LR_source': nodeLR_source[_winding],
                                'LI_source': nodeLI_source[_winding]}
                        else: 
                            dual_dict = None
                            
                        if not homotopy_option:
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampPrimary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodeVr_from[_pos], nodeVi_from[_pos],
                                nodeVr_from[_neg], nodeVi_from[_neg],
                                nodeVr_to[_pos], nodeVi_to[_pos],
                                nodeVr_to[_neg], nodeVi_to[_neg],
                                nodeVr_extra_to[_winding],
                                nodeVi_extra_to[_winding],
                                Vsource_Real[_winding], Vsource_Imag[_winding],
                                dual_dict = dual_dict, stamp_dual = stamp_dual)

                            # #########Stamp Secondary Circuit##################
                            # Without homotopy option OFF
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodeVr_from[_pos], nodeVi_from[_pos],
                                nodeVr_from[_neg], nodeVi_from[_neg],
                                nodeVr_to[_pos], nodeVi_to[_pos],
                                nodeVr_to[_neg], nodeVi_to[_neg],
                                nodeVr_extra_to[_winding],
                                nodeVi_extra_to[_winding],
                                Vsource_Real[_winding], Vsource_Imag[_winding],
                                homotopy_option, 
                                dual_dict = dual_dict, stamp_dual = stamp_dual)
                        else:
                            G_factor = h_factor * G_homotopy
                            B_factor = h_factor * B_homotopy
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodeVr_from[_pos], nodeVi_from[_pos],
                                nodeVr_from[_neg], nodeVi_from[_neg],
                                nodeVr_to[_pos], nodeVi_to[_pos],
                                nodeVr_to[_neg], nodeVi_to[_neg],
                                nodeVr_extra_to[_winding],
                                nodeVi_extra_to[_winding],
                                Vsource_Real[_winding], Vsource_Imag[_winding],
                                homotopy_option, G_factor, B_factor, dual_dict = dual_dict,
                                stamp_dual = stamp_dual)

            if self.connect_type == self._GWYE_GWYE:
                # #########Stamp Neutral to Ground##################
                _neg = _N
                # One for set of three phase
                if stamp_dual:
                    dual_dict_prim = {
                    'LR_n': nodeLR_from[_neg],
                    'LI_n': nodeLI_from[_neg],
                    'LR_gnd': self.nodeGnd_Lrp,
                    'LI_gnd': self.nodeGnd_Lip,
                    }
                    dual_dict_sec = {
                    'LR_n': nodeLR_to[_neg],
                    'LI_n': nodeLI_to[_neg],
                    'LR_gnd': self.nodeGnd_Lrs,
                    'LI_gnd': self.nodeGnd_Lis,
                    }
                else:
                    dual_dict_prim = None
                    dual_dict_sec = None
                # Stamp primary
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                    nodeVr_from[_neg], nodeVi_from[_neg], self.nodeVGnd_rp,
                    self.nodeVGnd_ip, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                    self.from_node, dual_dict = dual_dict_prim, stamp_dual = stamp_dual)
                # Stamp secondary
                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                    nodeVr_to[_neg], nodeVi_to[_neg], self.nodeVGnd_rs,
                    self.nodeVGnd_is, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                    self.to_node, dual_dict = dual_dict_sec, stamp_dual = stamp_dual)

            # Option 2 - Stamp Delta - Delta
            elif self.connect_type == self._DELTA_DELTA:
                phases = [(_A, _B), (_B, _C), (_C, _A)]
                for phaseSet in phases:
                    # In delta - delta get a tuple of AB, BC, CA
                    # #########Stamp Primary Circuit##################
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]
                    
                    if stamp_dual:
                            dualdict = {
                                'LR_from_p': nodeLR_from[_pos], 
                                'LI_from_p': nodeLI_from[_pos],
                                'LR_from_n': nodeLR_from[_neg], 
                                'LI_from_n': nodeLI_from[_neg],
                                'LR_to_p': nodeLR_to[_pos], 
                                'LI_to_p': nodeLI_to[_pos],
                                'LR_to_n': nodeLR_to[_neg], 
                                'LI_to_n': nodeLI_to[_neg],
                                'LR_extra_node': nodeLR_extra_to[_winding],
                                'LI_extra_node': nodeLI_extra_to[_winding],
                                'LR_source': nodeLR_source[_winding],
                                'LI_source': nodeLI_source[_winding]}
                    else: 
                        dualdict = None
                            
                    if not homotopy_option:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampPrimary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg], nodeVi_from[_neg],
                            nodeVr_to[_pos], nodeVi_to[_pos], nodeVr_to[_neg],
                            nodeVi_to[_neg], nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding],
                            dual_dict = dualdict, stamp_dual = stamp_dual)

                        # #########Stamp Secondary Circuit##################
                        # Without homotopy - option off
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg], nodeVi_from[_neg],
                            nodeVr_to[_pos], nodeVi_to[_pos], nodeVr_to[_neg],
                            nodeVi_to[_neg], nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, 
                            dual_dict = dualdict, stamp_dual = stamp_dual)
                    else:
                        G_factor = h_factor * G_homotopy
                        B_factor = h_factor * B_homotopy
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg], nodeVi_from[_neg],
                            nodeVr_to[_pos], nodeVi_to[_pos], nodeVr_to[_neg],
                            nodeVi_to[_neg], nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, G_factor,
                            B_factor, dual_dict = dualdict, stamp_dual = stamp_dual)

            # Option 3 or 7- Stamp Delta - GWye
            # ** This will have a phase shift
            elif self.connect_type == self._DELTA_GWYE or self.connect_type == self._DELTA_WYE:

                # The connection is as follows: (Step down)
                # from -> [(_C, _A), (_A, _B), (_B, _C)]
                # to ->   [(_N, _A), (_N, _B), (_N, _C)]

                # The connection is as follows: (Step up)
                # from -> [(_A, _B), (_B, _C), (_C, _A)]
                # to ->   [(_A, _N), (_B, _N), (_C, _N)]
                if self.is_stepdown:
                    phases = [(_C, _A), (_A, _B), (_B, _C)]
                else:
                    phases = [(_A, _B), (_B, _C), (_C, _A)]

                nodeVr_to.append(nodeN_Vr_to)
                nodeVi_to.append(nodeN_Vi_to)
                
                if stamp_dual:
                    nodeLR_to.append(nodeN_LR_to)
                    nodeLI_to.append(nodeN_LI_to)

                for phaseSet in phases:
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    if self.is_stepdown:
                        _pos_from = phaseSet[0]
                        _pos_to = _N
                        _neg_from = phaseSet[1]  # Primary is delta connected
                        # *This unusual behavior is to get NA, NB, NC
                        _neg_to = phaseSet[1]  # Secondary is wye connected
                        _winding = phaseSet[1]
                    else:
                        _pos_from = phaseSet[0]
                        _pos_to = phaseSet[0]
                        _neg_from = phaseSet[1]  # Primary is delta connected
                        _neg_to = _N  # Secondary is wye connected
                        _winding = phaseSet[0]
                        
                    if stamp_dual:
                        dual_dict = {
                            'LR_from_p': nodeLR_from[_pos_from], 
                            'LI_from_p': nodeLI_from[_pos_from],
                            'LR_from_n': nodeLR_from[_neg_from], 
                            'LI_from_n': nodeLI_from[_neg_from],
                            'LR_to_p': nodeLR_to[_pos_to], 
                            'LI_to_p': nodeLI_to[_pos_to],
                            'LR_to_n': nodeLR_to[_neg_to], 
                            'LI_to_n': nodeLI_to[_neg_to],
                            'LR_extra_node': nodeLR_extra_to[_winding],
                            'LI_extra_node': nodeLI_extra_to[_winding],
                            'LR_source': nodeLR_source[_winding],
                            'LI_source': nodeLI_source[_winding]}
                    else: 
                        dual_dict = None
                    
                    if not homotopy_option:
                        # #########Stamp Primary Circuit - DELTA##################
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampPrimary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos_from], nodeVi_from[_pos_from],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos_to], nodeVi_to[_pos_to],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], dual_dict = dual_dict, stamp_dual = stamp_dual)

                        # #########Stamp Secondary Circuit - Gwye##################
                        # Without homotopy - option off
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos_from], nodeVi_from[_pos_from],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos_to], nodeVi_to[_pos_to],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, dual_dict = dual_dict, stamp_dual = stamp_dual)
                    else:
                        G_factor = h_factor * G_homotopy
                        B_factor = h_factor * B_homotopy
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos_from], nodeVi_from[_pos_from],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos_to], nodeVi_to[_pos_to],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, G_factor,
                            B_factor, dual_dict = dual_dict, stamp_dual = stamp_dual)

                if self.connect_type == self._DELTA_GWYE:
                    # #########Stamp Neutral to Ground##################
                    # One for set of three phase
                    # Stamp secondary
                    # Make sure to ground the neutral and not neg to in this case
                    if stamp_dual:
                        dual_dict_sec = {
                        'LR_n': nodeLR_to[_N],
                        'LI_n': nodeLI_to[_N],
                        'LR_gnd': self.nodeGnd_Lrs,
                        'LI_gnd': self.nodeGnd_Lis,
                        }
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                        nodeVr_to[_N], nodeVi_to[_N], self.nodeVGnd_rs,
                        self.nodeVGnd_is, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                        self.to_node, dual_dict = dual_dict_sec, stamp_dual = stamp_dual)

            # Option 8 or 9- Stamp Wye- Delta or GWYE-Delta
            elif self.connect_type == self._WYE_DELTA or self.connect_type == self._GWYE_DELTA:
                # The connection is as follows: (Step down)
                # from -> [(_A, _N), (_B, _N), (_C, _N)]
                # to ->   [(_A, _B), (_B, _C), (_C, _A)]

                # The connection is as follows: (Step up)
                # from -> [(_A, _N), (_B, _N), (_C, _N)]
                # to ->   [(_A, _C), (_B, _A), (_C, _B)]
                if self.is_stepdown:
                    phases = [(_A, _B), (_B, _C), (_C, _A)]
                else:
                    phases = [(_A, _C), (_B, _A), (_C, _B)]
                nodeVr_from.append(nodeN_Vr_from)
                nodeVi_from.append(nodeN_Vi_from)
                
                if stamp_dual:
                    nodeLR_from.append(nodeN_LR_from)
                    nodeLI_from.append(nodeN_LI_from)

                for phaseSet in phases:
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    _pos = phaseSet[0]
                    _neg_from = _N  # Primary is wye connected
                    _neg_to = phaseSet[1]  # Secondary is delta connected
                    _winding = phaseSet[0]
                    
                    if stamp_dual:
                        dual_dict = {
                            'LR_from_p': nodeLR_from[_pos], 
                            'LI_from_p': nodeLI_from[_pos],
                            'LR_from_n': nodeLR_from[_neg_from], 
                            'LI_from_n': nodeLI_from[_neg_from],
                            'LR_to_p': nodeLR_to[_pos], 
                            'LI_to_p': nodeLI_to[_pos],
                            'LR_to_n': nodeLR_to[_neg_to], 
                            'LI_to_n': nodeLI_to[_neg_to],
                            'LR_extra_node': nodeLR_extra_to[_winding],
                            'LI_extra_node': nodeLI_extra_to[_winding],
                            'LR_source': nodeLR_source[_winding],
                            'LI_source': nodeLI_source[_winding]}
                    else: 
                        dual_dict = None
                        
                    if not homotopy_option:
                        # #########Stamp Primary Circuit - DELTA##################
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampPrimary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos], nodeVi_to[_pos],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], dual_dict = dual_dict, stamp_dual = stamp_dual)

                        # #########Stamp Secondary Circuit - Gwye##################
                        # Homotopy option - disabled
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos], nodeVi_to[_pos],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, 
                            dual_dict = dual_dict, stamp_dual = stamp_dual)
                    else:
                        G_factor = h_factor * G_homotopy
                        B_factor = h_factor * B_homotopy
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampSecondary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodeVr_from[_pos], nodeVi_from[_pos],
                            nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                            nodeVr_to[_pos], nodeVi_to[_pos],
                            nodeVr_to[_neg_to], nodeVi_to[_neg_to],
                            nodeVr_extra_to[_winding],
                            nodeVi_extra_to[_winding], Vsource_Real[_winding],
                            Vsource_Imag[_winding], homotopy_option, G_factor,
                            B_factor, dual_dict = dual_dict, stamp_dual = stamp_dual)

                if self.connect_type == self._GWYE_DELTA:
                    # #########Stamp Neutral to Ground##################
                    # One for set of three phase
                    # Stamp secondary
                    if stamp_dual:
                        dual_dict_prim = {
                            'LR_n': nodeLR_from[_neg_from],
                            'LI_n': nodeLI_from[_neg_from],
                            'LR_gnd': self.nodeGnd_Lrp,
                            'LI_gnd': self.nodeGnd_Lip,
                            }
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                        nodeVr_from[_neg_from], nodeVi_from[_neg_from],
                        self.nodeVGnd_rp, self.nodeVGnd_ip, Ylin_val, Ylin_row,
                        Ylin_col, idx_Y, stamped_ground, self.from_node, dual_dict = dual_dict_prim,
                        stamp_dual = stamp_dual)
            elif self.connect_type == self._CENTER_TAP:
                if self.phases & 0x01 == 1:  # Phase A
                    node1_Vr_from = self.node1_Vr_from if not enable_CM else curr_meas.nodeA_Vr_to
                    node1_Vi_from = self.node1_Vi_from if not enable_CM else curr_meas.nodeA_Vi_to
                    
                    if self.stamp_dual:
                        node1_Lr_from = self.node1_Lr_from if not enable_CM else curr_meas.nodeA_Lr_to
                        node1_Li_from = self.node1_Li_from if not enable_CM else curr_meas.nodeA_Li_to                        
                
                elif self.phases & 0x02 == 2:  # Phase B
                    node1_Vr_from = self.node1_Vr_from if not enable_CM else curr_meas.nodeB_Vr_to
                    node1_Vi_from = self.node1_Vi_from if not enable_CM else curr_meas.nodeB_Vi_to
                    
                    if self.stamp_dual:
                        node1_Lr_from = self.node1_Lr_from if not enable_CM else curr_meas.nodeB_Lr_to
                        node1_Li_from = self.node1_Li_from if not enable_CM else curr_meas.nodeB_Li_to                        
                
                elif self.phases & 0x04 == 4:  # Phase C
                    node1_Vr_from = self.node1_Vr_from if not enable_CM else curr_meas.nodeC_Vr_to
                    node1_Vi_from = self.node1_Vi_from if not enable_CM else curr_meas.nodeC_Vi_to
                    
                    if self.stamp_dual:
                        node1_Lr_from = self.node1_Lr_from if not enable_CM else curr_meas.nodeC_Lr_to
                        node1_Li_from = self.node1_Li_from if not enable_CM else curr_meas.nodeC_Li_to 
                        
                if not homotopy_option:
                    # *This is the case when neutral is grounded (both on primary and secondary
                    # Stamp the primary 4 current sources (between Xtra nodes and ground)
                    # Primary Current Sources
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtraP_rp, self.nodeVs1_rs, -(1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtraP_rp, self.nodeVs2_rs, (1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtraP_ip, self.nodeVs1_is, -(1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    #
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtraP_ip, self.nodeVs2_is, (1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr1_rs, self.nodeLrXtraP_rp, -(1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr2_rs, self.nodeLrXtraP_rp, (1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi1_is, self.nodeLiXtraP_ip, -(1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        #
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi2_is, self.nodeLiXtraP_ip, (1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    # Stamp 4 voltage sources (with 8 stamps) (each has two elements)
                    # Van_R - (VA_R/tr) = 0
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_rs, self.nodeXtra1_rs, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_rs, self.nodeN_Vr_to, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_rs, self.nodeXtraP_rp, (-1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLrXtra1_rs, self.nodeLr1_rs, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeN_Lr_to, self.nodeLr1_rs, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLrXtraP_rp, self.nodeLr1_rs,(-1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)                        
                    # Vbn_R + (VA_R/tr) = 0
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_rs, self.nodeXtra2_rs, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_rs, self.nodeN_Vr_to, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_rs, self.nodeXtraP_rp, (1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLrXtra2_rs, self.nodeLr2_rs, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeN_Lr_to, self.nodeLr2_rs, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLrXtraP_rp, self.nodeLr2_rs, (1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)                        
                    # Van_I - (VA_I/tr) = 0
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_is, self.nodeXtra1_is, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_is, self.nodeN_Vi_to, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs1_is, self.nodeXtraP_ip, (-1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLiXtra1_is, self.nodeLi1_is, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeN_Li_to, self.nodeLi1_is, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLiXtraP_ip, self.nodeLi1_is, (-1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)                        
                    # Van_I + (VA_I/tr) = 0
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_is, self.nodeXtra2_is, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_is, self.nodeN_Vi_to, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeVs2_is, self.nodeXtraP_ip, (1 / self.tr),
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLiXtra2_is, self.nodeLi2_is, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeN_Li_to, self.nodeLi2_is, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLiXtraP_ip, self.nodeLi2_is, (1 / self.tr),
                            Ylin_val, Ylin_row, Ylin_col, idx_Y)                        
                    # Stamp the current sources for the voltage sources
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtra1_rs, self.nodeVs1_rs, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtra2_rs, self.nodeVs2_rs, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtra1_is, self.nodeVs1_is, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeXtra2_is, self.nodeVs2_is, 1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)

                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeN_Vr_to, self.nodeVs1_rs, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeN_Vi_to, self.nodeVs2_rs, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeN_Vr_to, self.nodeVs1_is, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                        self.nodeN_Vi_to, self.nodeVs2_is, -1, Ylin_val,
                        Ylin_row, Ylin_col, idx_Y)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr1_rs, self.nodeLrXtra1_rs,  1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr2_rs, self.nodeLrXtra2_rs, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi1_is, self.nodeLiXtra1_is, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi2_is, self.nodeLiXtra2_is, 1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
    
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr1_rs, self.nodeN_Lr_to, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLr2_rs, self.nodeN_Li_to, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi1_is, self.nodeN_Lr_to, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                            self.nodeLi2_is, self.nodeN_Li_to, -1, Ylin_val,
                            Ylin_row, Ylin_col, idx_Y)                        
                    # Ground the neutrals if needed
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                        self.nodeN_Vr_to, self.nodeN_Vi_to, self.gnd_real_Vsec,
                        self.gnd_imag_Vsec, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                        self.to_node)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral_dual(
                            self.nodeN_Lr_to, self.nodeN_Li_to, self.gnd_real_Lrsec,
                            self.gnd_imag_Lisec, Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground,
                            self.to_node)    
                        
                    # This is the impedance stamps shared by
                    # Stamp Z0
                    (self.G0, self.B0) = self.calc_G_B(self.r0, self.x0)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        self.G0, self.B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtraP_rp, self.nodeXtraP_ip, node1_Vr_from,
                        node1_Vi_from)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_dual(
                            self.G0, self.B0, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtraP_rp, self.nodeLiXtraP_ip, node1_Lr_from,
                            node1_Li_from)                    
                    # Stamp shunt impedance
                    if self.hasShunt:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_shunt(
                            self.Gshunt, self.Bshunt, Ylin_val, Ylin_row,
                            Ylin_col, idx_Y, node1_Vr_from, node1_Vi_from)
                        
                        if self.stamp_dual:
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_shunt_dual(
                                self.Gshunt, self.Bshunt, Ylin_val, Ylin_row,
                                Ylin_col, idx_Y, node1_Lr_from, node1_Li_from)                            
                    # Stamp Z1 and Z2
                    (self.G1, self.B1) = self.calc_G_B(self.r1, self.x1)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        self.G1, self.B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtra1_rs, self.nodeXtra1_is, self.to_node1_Vr,
                        self.to_node1_Vi)
                    (self.G2, self.B2) = self.calc_G_B(self.r2, self.x2)
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        self.G2, self.B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtra2_rs, self.nodeXtra2_is, self.to_node2_Vr,
                        self.to_node2_Vi)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_dual(
                            self.G1, self.B1, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtra1_rs, self.nodeLiXtra1_is, self.to_node1_Lr,
                            self.to_node1_Li)
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_dual(
                            self.G2, self.B2, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtra2_rs, self.nodeLiXtra2_is, self.to_node2_Lr,
                            self.to_node2_Li)                       
                    
                else:  # Stamp all the homotopy functions
                    G_factor = h_factor * G_homotopy
                    B_factor = h_factor * B_homotopy
                    # Stamp homotopy Z0
                    G0_cont = self.G0 * G_factor
                    B0_cont = self.B0 * B_factor
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        G0_cont, B0_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtraP_rp, self.nodeXtraP_ip, node1_Vr_from,
                        node1_Vi_from)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_dual(
                            G0_cont, B0_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtraP_rp, self.nodeLiXtraP_ip, node1_Lr_from,
                            node1_Li_from)
                    # Stamp homotopy Z1
                    G1_cont = self.G1 * G_factor
                    B1_cont = self.B1 * B_factor
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        G1_cont, B1_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtraP_rp, self.nodeXtraP_ip, node1_Vr_from,
                        node1_Vi_from)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB_dual(
                            G1_cont, B1_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtraP_rp, self.nodeLiXtraP_ip, node1_Lr_from,
                            node1_Li_from)                    
                    # Stamp homotopy Z2
                    G2_cont = self.G2 * G_factor
                    B2_cont = self.B2 * B_factor
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                        G2_cont, B2_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        self.nodeXtraP_rp, self.nodeXtraP_ip, node1_Vr_from,
                        node1_Vi_from)
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampGB(
                            G2_cont, B2_cont, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            self.nodeLrXtraP_rp, self.nodeLiXtraP_ip, node1_Lr_from,
                            node1_Li_from)                       
            # TODO Call Homotopy Method using same function
            else:
                raise ValueError
        except ValueError:
            print(colored('Illegal connect_type for the transformer', 'red'))

        return Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground
