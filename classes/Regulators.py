"""Regulator implementation for distribution systems analysis.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 04-11-2017
  Updated Date: 10-28-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Fixed Regulators Validated, Variable WIP

  Regulators consist of an autotransformer and a load tap changing (LTC) mechanism. They are used to regulate the
  voltage on a particular bus. This class specifies the core properties of regulator objects and the stamps for
  regulators. Regulators are split between fixed regulators where the tap does not change and variable regulators
  where the taps change with voltage. See W.H. Kersting  Distribution System Analysis for details on
  how regulators work.
"""

from __future__ import division

from itertools import count
from types import SimpleNamespace

import numpy as np

from .Elements import LinearElement, NonlinearElement
from .Nodes import Nodes
from .RegulatorControl import RegulatorRatioControl


class FixedTapRegulator(LinearElement):
    # Create global class variables
    _GWYE_GWYE = 1
    _DELTA_DELTA = 2
    _ids = count(0)
    _regulator_locations = dict()

    def __init__(self,
                 name,
                 ID,
                 phases,
                 from_node,
                 to_node,
                 connect_type,
                 reg_type,
                 band_center,
                 bandwidth,
                 rated_power,
                 lowstep,
                 highstep,
                 num_tap_steps,
                 pt_ratio,
                 pt_phase,
                 ct_phase,
                 connected_transformer,
                 connected_winding,
                 tapA=1.0,
                 tapB=1.0,
                 tapC=1.0, 
                 stamp_dual = False):
        """

        Args:
            node_index_ (list): index of nodes for the network
            node (list): collection of node objects in the network
            name (str): Name of the regulator.
            ID  (str): Regulator ID.
            phases (int): The phases present in the regulator.
            from_node (str): From node that the element is connected to.
            to_node (str): To node that the element is connected to.
            connect_type (int): Type 1 or 2. Indicates Wye connection or Delta Connection.
            reg_type (str): Indicates if the regulator is Type A or B.
            band_center (float): the band center voltage of the regulator
            bandwidth (float): specifies how far the voltage can deviate from the band center
            rated_power (float): rated power of the regulator
            lowstep (int): minimum number of tap steps to lower by from neutral position
            highstep (int): maximum number of tap steps to raise by from neutral position
            num_tap_steps (int): total number of possible tap steps
            pt_ratio (float): The turns ratio used for the power transducer with a line-drop compensator.
            pt_phase (str): The phase of the power transducer used to monitor the voltage
            ct_phase (str): The phase of the current transformer
            connected_transformer (object): The transformer the voltage regulator is connected to
            connected_winding (int): the winding of the connected transformer being monitored by the regulator
            tapA (float): Phase A tap value of the regulator.
            tapB (float): Phase B tap value of the regulator.
            tapC (float): Phase C tap value of the regulator.
        """

        super(FixedTapRegulator, self).__init__()

        # Basic Properties
        self.name = name
        self.ID = ID
        self.from_node = from_node
        self.to_node = to_node
        self.connected_transformer = connected_transformer
        self.connected_winding = connected_winding
        self.reg_type = reg_type
        self.connect_type = connect_type
        self.phases = phases if phases else (pt_phase | ct_phase)
        self.stamp_dual = stamp_dual
        # grounded or not
        if self.connect_type == 1:
            self.isGnd = True
        else:
            self.isGnd = False

        # Regulator Settings
        self.pt_ratio = pt_ratio

        # Tap Values
        self.tapA = tapA
        self.tapB = tapB
        self.tapC = tapC

        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        aR_max = 1.1
        aR_min = 0.9
        self.aR = [1.0, 1.0, 1.0]
        self.aR_step = (aR_max - aR_min) / num_tap_steps
        self.band_center = band_center / self.pt_ratio
        self.bandwidth = bandwidth / self.pt_ratio
        self.Vreg_step = 120.0 * self.aR_step
        self.Vmax = (self.band_center + (self.bandwidth / 2))
        self.Vmin = (self.band_center - (self.bandwidth / 2))
        Va_reg = self.band_center
        Vb_reg = self.band_center
        Vc_reg = self.band_center
        if self.phases & 0x01 == int(0x01):
            Va_reg = self.band_center
        if self.phases & 0x02 == int(0x02):
            Vb_reg = self.band_center
        if self.phases & 0x04 == int(0x04):
            Vc_reg = self.band_center
        Vmag_reg = [Va_reg, Vb_reg, Vc_reg]
        self.raise_taps, self.lower_taps = self.check_regulator(Vmag_reg)

        self.rated_power = rated_power
        self.highstep = highstep
        self.lowstep = lowstep

        self.tap_positions = [self.tapA, self.tapB, self.tapC]
        self.initial_taps = [self.tapA, self.tapB, self.tapC]
        self.compute_effective_regulator_ratios()

        # Create and assign nodes
        self.nodeA_Vr_from = None
        self.nodeA_Vr_to = None
        self.nodeA_Vi_from = None
        self.nodeA_Vi_to = None

        self.nodeB_Vr_from = None
        self.nodeB_Vr_to = None
        self.nodeB_Vi_from = None
        self.nodeB_Vi_to = None

        self.nodeC_Vr_from = None
        self.nodeC_Vr_to = None
        self.nodeC_Vi_from = None
        self.nodeC_Vi_to = None

        self.nodeN_Vr_from = None
        self.nodeN_Vi_from = None
        self.nodeN_Vr_to = None
        self.nodeN_Vi_to = None

        self.nodeA_Vr_secondary = None
        self.nodeA_Vi_secondary = None
        self.nodeA_Vr_primary = None
        self.nodeA_Vi_primary = None

        self.nodeB_Vr_secondary = None
        self.nodeB_Vi_secondary = None
        self.nodeB_Vr_primary = None
        self.nodeB_Vi_primary = None

        self.nodeC_Vr_secondary = None
        self.nodeC_Vi_secondary = None
        self.nodeC_Vr_primary = None
        self.nodeC_Vi_primary = None

        self.nodeGnd_Vr_primary = None
        self.nodeGnd_Vi_primary = None
        self.nodeGnd_Vr_secondary = None
        self.nodeGnd_Vi_secondary = None

        self.node_set = []
        self.node_Vr_meas = {'A': 0, 'B': 0, 'C': 0}
        self.node_Vi_meas = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_diff = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_loadCenter = {'A': [], 'B': [], 'C': []}
        self.nodes_V_to = {'A': [], 'B': [], 'C': []}
        self.nodes_V_source = {'A': [], 'B': [], 'C': []}

        # Default loss terms
        self.G = 1e5
        self.B = 0

        # Count the nos. of changes
        self.taps_changed = False

        # Initialize Voltages
        self.Va_mag = 0
        self.Vb_mag = 0
        self.Vc_mag = 0
        self.Va_norm = 0
        self.Vb_norm = 0
        self.Vc_norm = 0
        if self.connect_type == self._DELTA_DELTA:
            self.Vab_mag = 0
            self.Vbc_mag = 0
            self.Vca_mag = 0

        # Initialize Currents
        self.Ia_mag = 0.0
        self.Ia_ang = 0.0
        self.Ib_mag = 0.0
        self.Ib_ang = 0.0
        self.Ic_mag = 0.0
        self.Ic_ang = 0.0

        # Initialiaze Effective Regulator Ratios
        self.aR_A = self.aR[0]
        self.aR_B = self.aR[1]
        self.aR_C = self.aR[2]

        self.isTriplex = False
        self.all_phases = {0: "A", 1: "B", 2: "C", 3: "N"}

    def assign_nodes(self,
                     node_key,
                     node_index_,
                     node):

        # There are likely to be 12 additional nodes for the XFMRs
        # 6 additional rows for the voltage source eqns. (No Phase Shifters
        # 12 additional rows for the voltage sources eqn (w Phase Shifters
        if self.phases & 0x1 == 1:  # Check for phase A
            self.nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr
            self.nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi
            self.nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
            self.nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi

            if self.stamp_dual:
                self.nodeA_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
                self.nodeA_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i
                self.nodeA_Lr_to = node[node_key[self.to_node]].nodeA_dual_eq_var_r
                self.nodeA_Li_to = node[node_key[self.to_node]].nodeA_dual_eq_var_i

        if self.phases & 0x2 == 2:  # Check for phase B
            self.nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr
            self.nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi
            self.nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
            self.nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi

            if self.stamp_dual:
                self.nodeB_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
                self.nodeB_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
                self.nodeB_Lr_to = node[node_key[self.to_node]].nodeB_dual_eq_var_r
                self.nodeB_Li_to = node[node_key[self.to_node]].nodeB_dual_eq_var_i
                
        if self.phases & 0x4 == 4:  # Check for phase C
            self.nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
            self.nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi
            self.nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
            self.nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi

            if self.stamp_dual:
                self.nodeC_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
                self.nodeC_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i
                self.nodeC_Lr_to = node[node_key[self.to_node]].nodeC_dual_eq_var_r
                self.nodeC_Li_to = node[node_key[self.to_node]].nodeC_dual_eq_var_i
                
        self.nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
        self.nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi
        self.nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
        self.nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi

        if self.stamp_dual:
            self.nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
            self.nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i
            self.nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
            self.nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i
       
        #############################################################
        # Add ground voltage source if needed
        if self.connect_type == self._GWYE_GWYE or self.connect_type == self._DELTA_DELTA:  # Wye-Wye
            # Add for the three phase extra node in the secondary
            # Add for the three phases extra node for the voltage sources (primarySource)
            # circuit (for both real and imaginary)
            if self.phases & 0x1 == 1:  # Check for phase A
                self.nodeA_Vr_secondary = node_index_.__next__()
                self.nodeA_Vi_secondary = node_index_.__next__()
                self.nodeA_Vr_primary = node_index_.__next__()
                self.nodeA_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["A"] = self.nodeA_Vr_primary
                self.node_Vi_meas["A"] = self.nodeA_Vi_primary
                Nodes.I_index.extend([self.nodeA_Vr_primary, self.nodeA_Vi_primary])
                self.node_set.extend([self.nodeA_Vr_secondary, self.nodeA_Vi_secondary, self.nodeA_Vr_primary,
                                      self.nodeA_Vi_primary])
                # self.node_set.extend([self.nodeA_Vr_primary, self.nodeA_Vi_primary])
                if self.stamp_dual:
                    self.nodeA_Lr_secondary = node_index_.__next__()
                    self.nodeA_Li_secondary = node_index_.__next__()
                    self.nodeA_Lr_primary = node_index_.__next__()
                    self.nodeA_Li_primary = node_index_.__next__()
                    Nodes.Lr_index.extend([self.nodeA_Lr_primary])
                    Nodes.Li_index.extend([self.nodeA_Li_primary])
                    
            if self.phases & 0x2 == 2:  # Check for phase B
                self.nodeB_Vr_secondary = node_index_.__next__()
                self.nodeB_Vi_secondary = node_index_.__next__()
                self.nodeB_Vr_primary = node_index_.__next__()
                self.nodeB_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["B"] = self.nodeB_Vr_primary
                self.node_Vi_meas["B"] = self.nodeB_Vi_primary
                Nodes.I_index.extend([self.nodeB_Vr_primary, self.nodeB_Vi_primary])
                self.node_set.extend([self.nodeB_Vr_secondary, self.nodeB_Vi_secondary, self.nodeB_Vr_primary,
                                      self.nodeB_Vi_primary])
                # self.node_set.extend([self.nodeB_Vr_primary, self.nodeB_Vi_primary])
                
                if self.stamp_dual:
                    self.nodeB_Lr_secondary = node_index_.__next__()
                    self.nodeB_Li_secondary = node_index_.__next__()
                    self.nodeB_Lr_primary = node_index_.__next__()
                    self.nodeB_Li_primary = node_index_.__next__()
                    Nodes.Lr_index.extend([self.nodeB_Lr_primary])
                    Nodes.Li_index.extend([self.nodeB_Li_primary])
                    
            if self.phases & 0x4 == 4:  # Check for phase C
                self.nodeC_Vr_secondary = node_index_.__next__()
                self.nodeC_Vi_secondary = node_index_.__next__()
                self.nodeC_Vr_primary = node_index_.__next__()
                self.nodeC_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["C"] = self.nodeC_Vr_primary
                self.node_Vi_meas["C"] = self.nodeC_Vi_primary
                Nodes.I_index.extend([self.nodeC_Vr_primary, self.nodeC_Vi_primary])
                self.node_set.extend([self.nodeC_Vr_secondary, self.nodeC_Vi_secondary, self.nodeC_Vr_primary,
                                      self.nodeC_Vi_primary])
                self.node_set.extend([self.nodeC_Vr_primary, self.nodeC_Vi_primary])
                
                if self.stamp_dual:
                    self.nodeC_Lr_secondary = node_index_.__next__()
                    self.nodeC_Li_secondary = node_index_.__next__()
                    self.nodeC_Lr_primary = node_index_.__next__()
                    self.nodeC_Li_primary = node_index_.__next__()
                    Nodes.Lr_index.extend([self.nodeC_Lr_primary])
                    Nodes.Li_index.extend([self.nodeC_Li_primary])
            # Check for ground and add further nodes
            if self.isGnd:
                if (self.from_node,
                    self.to_node) not in self._regulator_locations:
                    self.nodeGnd_Vr_primary = node_index_.__next__()
                    self.nodeGnd_Vi_primary = node_index_.__next__()
                    self.nodeGnd_Vr_secondary = node_index_.__next__()
                    self.nodeGnd_Vi_secondary = node_index_.__next__()
                    
                    if self.stamp_dual:
                        self.nodeGnd_Lr_primary = node_index_.__next__()
                        self.nodeGnd_Li_primary = node_index_.__next__()
                        self.nodeGnd_Lr_secondary = node_index_.__next__()
                        self.nodeGnd_Li_secondary = node_index_.__next__()
                else:
                    self.nodeGnd_Vr_primary = self._regulator_locations[(
                        self.from_node, self.to_node)][0]
                    self.nodeGnd_Vi_primary = self._regulator_locations[(
                        self.from_node, self.to_node)][1]
                    self.nodeGnd_Vr_secondary = self._regulator_locations[(
                        self.from_node, self.to_node)][2]
                    self.nodeGnd_Vi_secondary = self._regulator_locations[(
                        self.from_node, self.to_node)][3]
                    
                    if self.stamp_dual:
                        self.nodeGnd_Lr_primary = self._regulator_locations[(
							self.from_node, self.to_node)][0]
                        self.nodeGnd_Li_primary = self._regulator_locations[(
							self.from_node, self.to_node)][1]
                        self.nodeGnd_Lr_secondary = self._regulator_locations[(
							self.from_node, self.to_node)][2]
                        self.nodeGnd_Li_secondary = self._regulator_locations[(
							self.from_node, self.to_node)][3]

        return node_index_

    def check_connections(self):
        # check if regulator connections already exist
        if (self.from_node, self.to_node) not in self._regulator_locations:
            if self.connect_type == 1:
                self._regulator_locations[(self.from_node, self.to_node)] = [
                    self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary, self.nodeGnd_Vr_secondary,
                    self.nodeGnd_Vi_secondary
                ]
            else:
                pass

    def compute_effective_regulator_ratios(self):
        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        if self.reg_type == 'A':  # Type A
            self.aR = [(1 + (self.aR_step * tap)) ** -1 for tap in self.tap_positions]
        else:  # Type B
            self.aR = [1 - (self.aR_step * tap) for tap in self.tap_positions]

    def check_regulator(self,
                        Vabc):

        self.Va_mag = np.abs(Vabc[0])
        self.Vb_mag = np.abs(Vabc[1])
        self.Vc_mag = np.abs(Vabc[2])

        # Set raise or lower tap flags depending on connection type
        raise_tapA = False
        lower_tapA = False
        raise_tapB = False
        lower_tapB = False
        raise_tapC = False
        lower_tapC = False
        if self.phases & 0x01 == int(0x01):
            raise_tapA = True if self.Va_mag < self.Vmin else False
            lower_tapA = True if self.Va_mag > self.Vmax else False
        if self.phases & 0x02 == int(0x02):
            raise_tapB = True if self.Vb_mag < self.Vmin else False
            lower_tapB = True if self.Vb_mag > self.Vmax else False
        if self.phases & 0x04 == int(0x04):
            raise_tapC = True if self.Vc_mag < self.Vmin else False
            lower_tapC = True if self.Vc_mag > self.Vmax else False
        raise_taps = [raise_tapA, raise_tapB, raise_tapC]
        lower_taps = [lower_tapA, lower_tapB, lower_tapC]

        return raise_taps, lower_taps

    @staticmethod
    def calc_G_B(r,
                 x):
        if r == 0:
            G = None
        else:
            G = r / (r ** 2 + x ** 2)
        if x == 0:
            B = None
        else:
            B = -x / (r ** 2 + x ** 2)
        return G, B

    def append_nodes(self):
        # Put together the set of nodes
        nodes_Vr_from = []
        nodes_Vi_from = []
        nodes_Vr_to = []
        nodes_Vi_to = []
        nodes_Vr_secondary = []
        nodes_Vi_secondary = []
        Ir_sec = []
        Ii_sec = []
        
        if self.stamp_dual:
            nodes_Lr_from = []
            nodes_Li_from = []
            nodes_Lr_to = []
            nodes_Li_to = []
            nodes_Lr_secondary = []
            nodes_Li_secondary = []
            nodes_Lr_primary = []
            nodes_Li_primary = []            

        # Find the from (Primary) bus nodes to stamp
        if self.phases & 0x01 == int(0x01):
            nodes_Vr_from.append(self.nodeA_Vr_from)
            nodes_Vi_from.append(self.nodeA_Vi_from)
            nodes_Vr_to.append(self.nodeA_Vr_to)
            nodes_Vi_to.append(self.nodeA_Vi_to)
            nodes_Vr_secondary.append(self.nodeA_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeA_Vi_secondary)
            Ir_sec.append(self.nodeA_Vr_primary)
            Ii_sec.append(self.nodeA_Vi_primary)
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeA_Lr_from)
                nodes_Li_from.append(self.nodeA_Li_from)
                nodes_Lr_to.append(self.nodeA_Lr_to)
                nodes_Li_to.append(self.nodeA_Li_to)
                nodes_Lr_secondary.append(self.nodeA_Lr_secondary)
                nodes_Li_secondary.append(self.nodeA_Li_secondary)
                nodes_Lr_primary.append(self.nodeA_Lr_primary)
                nodes_Li_primary.append(self.nodeA_Li_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
			
            if self.stamp_dual:
                nodes_Lr_from.append(-1)
                nodes_Li_from.append(-1)
                nodes_Lr_to.append(-1)
                nodes_Li_to.append(-1)
                nodes_Lr_secondary(-1)
                nodes_Li_secondary(-1)
                nodes_Lr_primary(-1)
                nodes_Li_primary(-1)
                
        if self.phases & 0x02 == int(0x02):
            nodes_Vr_from.append(self.nodeB_Vr_from)
            nodes_Vi_from.append(self.nodeB_Vi_from)
            nodes_Vr_to.append(self.nodeB_Vr_to)
            nodes_Vi_to.append(self.nodeB_Vi_to)
            nodes_Vr_secondary.append(self.nodeB_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeB_Vi_secondary)
            Ir_sec.append(self.nodeB_Vr_primary)
            Ii_sec.append(self.nodeB_Vi_primary)
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeB_Lr_from)
                nodes_Li_from.append(self.nodeB_Li_from)
                nodes_Lr_to.append(self.nodeB_Lr_to)
                nodes_Li_to.append(self.nodeB_Li_to)
                nodes_Lr_secondary.append(self.nodeB_Lr_secondary)
                nodes_Li_secondary.append(self.nodeB_Li_secondary)
                nodes_Lr_primary.append(self.nodeB_Lr_primary)
                nodes_Li_primary.append(self.nodeB_Li_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
            
            if self.stamp_dual:
                nodes_Lr_from.append(-1)
                nodes_Li_from.append(-1)
                nodes_Lr_to.append(-1)
                nodes_Li_to.append(-1)
                nodes_Lr_secondary(-1)
                nodes_Li_secondary(-1)
                nodes_Lr_primary(-1)
                nodes_Li_primary(-1)
                
        if self.phases & 0x04 == int(0x04):
            nodes_Vr_from.append(self.nodeC_Vr_from)
            nodes_Vi_from.append(self.nodeC_Vi_from)
            nodes_Vr_to.append(self.nodeC_Vr_to)
            nodes_Vi_to.append(self.nodeC_Vi_to)
            nodes_Vr_secondary.append(self.nodeC_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeC_Vi_secondary)
            Ir_sec.append(self.nodeC_Vr_primary)
            Ii_sec.append(self.nodeC_Vi_primary)
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeC_Lr_from)
                nodes_Li_from.append(self.nodeC_Li_from)
                nodes_Lr_to.append(self.nodeC_Lr_to)
                nodes_Li_to.append(self.nodeC_Li_to)
                nodes_Lr_secondary.append(self.nodeC_Lr_secondary)
                nodes_Li_secondary.append(self.nodeC_Li_secondary)
                nodes_Lr_primary.append(self.nodeC_Lr_primary)
                nodes_Li_primary.append(self.nodeC_Li_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
            
            if self.stamp_dual:
                nodes_Lr_from.append(-1)
                nodes_Li_from.append(-1)
                nodes_Lr_to.append(-1)
                nodes_Li_to.append(-1)
                nodes_Lr_secondary(-1)
                nodes_Li_secondary(-1)
                nodes_Lr_primary(-1)
                nodes_Li_primary(-1)
                
        if self.connect_type == self._GWYE_GWYE:
            nodes_Vr_from.append(self.nodeN_Vr_from)
            nodes_Vi_from.append(self.nodeN_Vi_from)
            nodes_Vr_to.append(self.nodeN_Vr_to)
            nodes_Vi_to.append(self.nodeN_Vi_to)
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeN_Lr_from)
                nodes_Li_from.append(self.nodeN_Li_from)
                nodes_Lr_to.append(self.nodeN_Lr_to)
                nodes_Li_to.append(self.nodeN_Li_to)

        if self.stamp_dual:
            L_dict = {'nodes_Lr_from': nodes_Lr_from,
					  'nodes_Li_from': nodes_Li_from,
					  'nodes_Lr_to': nodes_Lr_to,
					  'nodes_Li_to': nodes_Li_to,
                      'nodes Lr_secondary': nodes_Lr_secondary,
                      'nodes Li_secondary': nodes_Li_secondary,
					  'nodes_Lr_primary': nodes_Lr_primary,
					  'nodes_Li_primary': nodes_Li_primary}
        else:
            L_dict = None
        
        return nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
               nodes_Vr_secondary, nodes_Vi_secondary, Ir_sec, Ii_sec, L_dict

    def update_values(self, V):
        # # Update Regulator Currents # #
        deg = True
        Ia = complex(V[self.nodeA_Vr_primary], V[self.nodeA_Vi_primary]) if hasattr(self, 'nodeA_Vr_primary') else 0
        self.Ia_mag = np.around(np.abs(Ia), 2)
        self.Ia_ang = np.around(np.angle(Ia, deg=deg), 2)
        Ib = complex(V[self.nodeB_Vr_primary], V[self.nodeB_Vi_primary]) if hasattr(self, 'nodeB_Vr_primary') else 0
        self.Ib_mag = np.around(np.abs(Ib), 2)
        self.Ib_ang = np.around(np.angle(Ib, deg=deg), 2)
        Ic = complex(V[self.nodeC_Vr_primary], V[self.nodeC_Vi_primary]) if hasattr(self, 'nodeC_Vr_primary') else 0
        self.Ic_mag = np.around(np.abs(Ic), 2)
        self.Ic_ang = np.around(np.angle(Ic, deg=deg), 2)

    def stamp_primary(self,
                      Ylin_val,
                      Ylin_row,
                      Ylin_col,
                      idx_Y,
                      node_Vr_pos_from,
                      node_Vi_pos_from,
                      node_Vr_neg_from,
                      node_Vi_neg_from,
                      node_Vr_neg_secondary,
                      node_Vi_neg_secondary,
                      node_Vr_pos_secondary,
                      node_Vi_pos_secondary,
                      node_Vr_primary,
                      node_Vi_primary,
                      phase):

        # # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]

        # # PRIMARY VOLTAGES # #
        # Real #
        # Vr_primary = Vr_pos_from - tr * Vr_secondary_pos  + tr * Vr_secondary_neg - Vr_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected

        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_primary, node_Vr_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_primary, node_Vr_pos_secondary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_primary, node_Vr_neg_secondary, aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_primary, node_Vr_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        # Imaginary #
        # Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_primary, node_Vi_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_primary, node_Vi_pos_secondary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_primary, node_Vi_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_primary, node_Vi_neg_secondary, aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)

        # PRIMARY CURRENTS
        # Real #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_pos_from, node_Vr_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vr_neg_from, node_Vr_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        # Imaginary #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_pos_from, node_Vi_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
            node_Vi_neg_from, node_Vi_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_primary_dual(self, Ylin_val, Ylin_row, Ylin_col, idx_Y,
					 node_Lr_pos_from, node_Li_pos_from, node_Lr_neg_from,
					 node_Li_neg_from,  node_Lr_neg_secondary,
					 node_Li_neg_secondary, node_Lr_pos_to,
					 node_Li_pos_to, node_Lr_primary, node_Li_primary,
					 phase):

		# # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]
		
        # DUAL OF PRIMARY VOLTAGES - TRANSPOSE OF LINEAR STAMPS
		# Real #
		# Vr_primary = Vr_pos_from - tr * Vr_secondary_pos  + tr * Vr_secondary_neg - Vr_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_pos_from, node_Lr_primary, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_pos_to, node_Lr_primary, -aR, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_neg_secondary, node_Lr_primary, aR, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_neg_from, node_Lr_primary, -1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)

		# Imaginary #
		# Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
		# neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_pos_from, node_Li_primary, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_pos_to, node_Li_primary, -aR, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_neg_from, node_Li_primary, -1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_neg_secondary, node_Li_primary,  aR, Ylin_val, Ylin_row,
			Ylin_col, idx_Y)

		# PRIMARY CURRENTS
		# Real #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_primary, node_Lr_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Lr_primary, node_Lr_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
		# Imaginary #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_primary, node_Li_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			node_Li_primary, node_Li_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
			idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_secondary(self,
                        Ylin_val,
                        Ylin_row,
                        Ylin_col,
                        idx_Y,
                        node_Vr_pos_to,
                        node_Vi_pos_to,
                        node_Vr_neg_secondary,
                        node_Vi_neg_secondary,
                        node_Vr_pos_secondary,
                        node_Vi_pos_secondary,
                        node_Vr_primary,
                        node_Vi_primary,
                        phase):

        # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]
        # Calculate conductance and susceptance of regulator if homotopy
        _G = self.G

        if aR:
            # # SECONDARY CURRENT SOURCES # #
            # Real #
            # Ir_pos_secondary = - aR * Ir_primary
            # Ir_neg_secondary = aR * Ir_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -aR, Ylin_val,
                                                               Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vr_neg_secondary, node_Vr_primary, aR, Ylin_val,
                                                               Ylin_row, Ylin_col, idx_Y)

            # Imaginary #
            # Ii_pos_secondary = - aR * Ii_primary
            # Ii_neg_secondary = aR * Ii_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -aR, Ylin_val,
                                                               Ylin_row, Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(node_Vi_neg_secondary, node_Vi_primary, aR, Ylin_val,
                                                               Ylin_row, Ylin_col, idx_Y)

        if _G:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vr_pos_to, node_Vr_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vr_pos_to, node_Vr_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vi_pos_to, node_Vi_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vi_pos_to, node_Vi_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vr_pos_secondary, node_Vr_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vr_pos_secondary, node_Vr_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vi_pos_secondary, node_Vi_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
                node_Vi_pos_secondary, node_Vi_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_secondary_dual(self, Ylin_val, Ylin_row, Ylin_col, idx_Y,
					   node_Lr_pos_to, node_Li_pos_to,
					   node_Lr_neg_secondary, node_Li_neg_secondary,
					   node_Lr_primary, node_Li_primary, phase,
					   homotopy = None, B_factor = None, G_factor = None):

        # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]
        # Calculate conductance and susceptance of regulator if homotopy
        _G = self.G
        if aR:
			# # SECONDARY CURRENT SOURCES # #
			# Real #
			# Ir_pos_secondary = - aR * Ir_primary
			# Ir_neg_secondary = aR * Ir_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
				node_Lr_primary, node_Lr_pos_to, -aR, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			   node_Lr_primary,  node_Lr_neg_secondary, aR, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)

			# Imaginary #
			# Ii_pos_secondary = - aR * Ii_primary
			# Ii_neg_secondary = aR * Ii_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
				node_Li_primary, node_Li_pos_to, -aR, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
				node_Li_primary, node_Li_neg_secondary, aR, Ylin_val, Ylin_row,
				Ylin_col, idx_Y)

		# Adding resistance if existent
#         if _G:
# 			# Branch stamping for the conductance
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Lr_pos_to, node_Lr_pos_to, _G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Lr_pos_secondary,node_Lr_pos_to,  -_G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Li_pos_to, node_Li_pos_to, _G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Li_pos_secondary, node_Li_pos_to,  -_G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Lr_pos_secondary, node_Lr_pos_secondary, _G, Ylin_val,
# 				Ylin_row, Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Lr_pos_to, node_Lr_pos_secondary, -_G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Li_pos_secondary, node_Li_pos_secondary, _G, Ylin_val,
# 				Ylin_row, Ylin_col, idx_Y)
# 			#
#             Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
# 				node_Li_pos_to, node_Li_pos_secondary, -_G, Ylin_val, Ylin_row,
# 				Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_ground_neutral(self,
                             neg_real_node,
                             neg_imag_node,
                             gnd_real_V,
                             gnd_imag_V,
                             Ylin_val,
                             Ylin_row,
                             Ylin_col,
                             idx_Y):
        # Stamp Vng = 0 + j0
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(gnd_real_V, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col,
                                                            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(gnd_imag_V, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col,
                                                            idx_Y)
        # #
        # # Stamp Current Sources
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(neg_real_node, gnd_real_V, 1, Ylin_val, Ylin_row, Ylin_col,
                                                            idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(neg_imag_node, gnd_imag_V, 1, Ylin_val, Ylin_row, Ylin_col,
                                                            idx_Y)
        

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
	
    def stamp_ground_neutral_dual(self, neg_real_node, neg_imag_node, gnd_real_Lr,
							 gnd_imag_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y):
		# Stamp Vng = 0 + j0
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			neg_real_node, gnd_real_Lr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			neg_imag_node, gnd_imag_Li, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
		# Stamp Current Sources
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			gnd_real_Lr, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(
			gnd_imag_Li, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
		#

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_linear(self,
                     Ylin_val,
                     Ylin_row,
                     Ylin_col,
                     Jlin_row,
                     Jlin_val,
                     idx_Y,
                     idx_J,
                     stamped_ground):

        (nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to,
         nodes_Vr_secondary, nodes_Vi_secondary, 
         nodes_Vr_primary, nodes_Vi_primary, L_dict) = self.append_nodes()

        if self.stamp_dual:
            nodes_Lr_from = L_dict['nodes_Lr_from']
            nodes_Li_from = L_dict['nodes_Li_from']
            nodes_Lr_to = L_dict['nodes_Lr_to']
            nodes_Li_to = L_dict['nodes_Li_to']
            nodes_Lr_primary = L_dict['nodes_Lr_primary']
            nodes_Li_primary = L_dict['nodes_Li_primary']
            
        (_A, _B, _C) = (0, 1, 2)
        _N = 3

        try:

            # Option 1 - Stamp gwye-gwye
            if self.connect_type == self._GWYE_GWYE:
                _phases_list = []
                if self.phases & 0x01 == int(0x01):
                    _phases_list.append((_A, _N))

                if self.phases & 0x02 == int(0x02):
                    _phases_list.append((_B, _N))

                if self.phases & 0x04 == int(0x04):
                    _phases_list.append((_C, _N))

                for phaseSet in _phases_list:
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]

                    ##########Stamp Primary Circuit##################
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
                        Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        nodes_Vr_from[_pos], nodes_Vi_from[_pos],
                        nodes_Vr_from[_neg], nodes_Vi_from[_neg],
                        nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                        nodes_Vr_secondary[_winding],
                        nodes_Vi_secondary[_winding],
                        nodes_Vr_primary[_winding],
                        nodes_Vi_primary[_winding], phaseSet[0])
					
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary_dual(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodes_Lr_from[_pos], nodes_Li_from[_pos],
								nodes_Lr_from[_neg], nodes_Li_from[_neg],
								nodes_Lr_to[_neg], nodes_Li_to[_neg],
								nodes_Lr_to[_winding],
								nodes_Li_to[_winding],
								nodes_Lr_primary[_winding],
								nodes_Li_primary[_winding], phaseSet[0])
                    
                    ##########Stamp Secondary Circuit##################
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
                        Ylin_val, Ylin_row, Ylin_col, idx_Y,
                        nodes_Vr_to[_pos], nodes_Vi_to[_pos],
                        nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                        nodes_Vr_secondary[_winding], nodes_Vi_secondary[_winding],
                        nodes_Vr_primary[_winding],
                        nodes_Vi_primary[_winding], phaseSet[0])
						
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary_dual(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodes_Lr_to[_pos], nodes_Li_to[_pos],
								nodes_Lr_to[_neg], nodes_Li_to[_neg],
								nodes_Lr_primary[_winding],
								nodes_Li_primary[_winding], phaseSet[0])
                # #########Stamp Neutral to Ground##################
                _neg = _N
                # One for set of three phase
                # Stamp primary
                if self.from_node not in stamped_ground:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                        nodes_Vr_from[_neg], nodes_Vi_from[_neg],
                        self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary,
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    stamped_ground.add(self.from_node)
					
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral_dual(
							nodes_Lr_from[_neg], nodes_Li_from[_neg], self.nodeGnd_Lr_primary,
							self.nodeGnd_Li_primary, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        stamped_ground.add(self.from_node)
                
                # Stamp secondary
                if self.to_node not in stamped_ground:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(
                        nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                        self.nodeGnd_Vr_secondary, self.nodeGnd_Vi_secondary,
                        Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    stamped_ground.add(self.to_node)

                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral_dual(
							nodes_Lr_to[_neg], nodes_Li_to[_neg], self.nodeGnd_Lr_secondary,
							self.nodeGnd_Li_secondary, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        stamped_ground.add(self.to_node)
                        
            elif self.connect_type == self._DELTA_DELTA:
                _phases = [(_A, _B), (_B, _C), (_C, _A)]
                for phaseSet in _phases:
                    # In delta - delta get a tuple of AB, BC, CA
                    # #########Stamp Primary Circuit##################
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(Ylin_val, Ylin_row, Ylin_col,
                                                                             idx_Y, nodes_Vr_from[_pos],
                                                                             nodes_Vi_from[_pos],
                                                                             nodes_Vr_from[_neg],
                                                                             nodes_Vi_from[_neg],
                                                                             nodes_Vr_to[_neg],
                                                                             nodes_Vi_to[_neg],
                                                                             nodes_Vr_secondary[_winding],
                                                                             nodes_Vi_secondary[_winding],
                                                                             nodes_Vr_primary[_winding],
                                                                             nodes_Vi_primary[_winding],
                                                                             phaseSet[0])
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary_dual(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodes_Lr_from[_pos], nodes_Li_from[_pos],
								nodes_Lr_from[_neg], nodes_Li_from[_neg],
								nodes_Lr_to[_neg], nodes_Li_to[_neg],
								nodes_Lr_to[_winding],
								nodes_Li_to[_winding],
								nodes_Lr_primary[_winding],
								nodes_Li_primary[_winding], phaseSet[0])
                        
                    # #########Stamp Secondary Circuit##################
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(Ylin_val, Ylin_row, Ylin_col,
                                                                               idx_Y,
                                                                               nodes_Vr_to[_pos],
                                                                               nodes_Vi_to[_pos], nodes_Vr_to[_neg],
                                                                               nodes_Vi_to[_neg],
                                                                               nodes_Vr_secondary[_winding], nodes_Vi_secondary[_winding],
                                                                               nodes_Vr_primary[_winding],
                                                                               nodes_Vi_primary[_winding],
                                                                               phaseSet[0])
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary_dual(
								Ylin_val, Ylin_row, Ylin_col, idx_Y,
								nodes_Lr_to[_pos], nodes_Li_to[_pos],
								nodes_Lr_to[_neg], nodes_Li_to[_neg],
								nodes_Lr_primary[_winding],
								nodes_Li_primary[_winding], phaseSet[0])

        except ValueError:
            print('Illegal connect_type for the regulator')

        return Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J, stamped_ground


class VariableRegulator(NonlinearElement):
    """A power control element that regulates the voltage at a particular bus.

    Attributes:
        name (str): Name of the regulator.
        ID  (str): Regulator ID.
        from_node (str): From node that the element is connected to.
        to_node (str): To node that the element is connected to.
        connected_transformer (str): Name of the transformer that the regulator is connected to.
        connected_winding (int): Number of the winding that the regulator is monitoring.
        reg_type (str): Indicates if the regulator is Type A or B.
        connect_type (int): Type 1 or 2. Indicates Wye connection or Delta Connection.
        tapA (float): Phase A tap value of the regulator.
        tapB (float): Phase B tap value of the regulator.
        tapC (float): Phase C tap value of the regulator.
        phases (int): The phases present in the regulator.
        tap_positions (list): A collection of the phase taps.
        node_Vr_meas (dict): Collects the real ammeter measurement nodes.
        node_Vi_meas (dict): Collects the imaginary ammeter measurement nodes.
        aR (list): Calculates and collects the effective regulator ratios.
        tapA_change_count: The number of phase A tap changes made.
        tapB_change_count: The number of phase B tap changes made.
        tapC_change_count: The number of phase C tap changes made.
        G (float): Regulator loss conductance.
        B (float): Regulator loss susceptance.
    """
    # Create variable names
    _GWYE_GWYE = 1
    _DELTA_DELTA = 2
    _ids = count(0)
    _regulator_locations = dict()

    def __init__(self,
                 name,
                 ID,
                 phases,
                 from_node,
                 to_node,
                 connect_type,
                 reg_type,
                 band_center,
                 bandwidth,
                 rated_power,
                 lowstep,
                 highstep,
                 num_tap_steps,
                 pt_ratio,
                 pt_phase,
                 ct_phase,
                 connected_transformer,
                 connected_winding):
        """

        Args:
            name (str): Name of the regulator.
            ID  (str): Regulator ID.
            phases (int): The phases present in the regulator.
            from_node (str): From node that the element is connected to.
            to_node (str): To node that the element is connected to.
            connect_type (int): Type 1 or 2. Indicates Wye connection or Delta Connection.
            reg_type (str): Indicates if the regulator is Type A or B.
            band_center (float): the band center voltage of the regulator
            bandwidth (float): specifies how far the voltage can deviate from the band center
            rated_power (float): rated power of the regulator
            lowstep (int): minimum number of tap steps to lower by from neutral position
            highstep (int): maximum number of tap steps to raise by from neutral position
            num_tap_steps (int): total number of possible tap steps
            pt_ratio (float): The turns ratio used for the power transducer with a line-drop compensator.
            pt_phase (str): The phase of the power transducer used to monitor the voltage
            ct_phase (str): The phase of the current transformer
            connected_transformer (object): The transformer the voltage regulator is connected to
            connected_winding (int): the winding of the connected transformer being monitored by the regulator
        """

        super(VariableRegulator, self).__init__()

        # Basic Properties
        self.name = name
        self.ID = ID
        self.from_node = from_node
        self.to_node = to_node
        self.connected_transformer = connected_transformer
        self.connected_winding = connected_winding
        self.reg_type = reg_type
        self.connect_type = connect_type
        self.phases = phases if phases else (pt_phase | ct_phase)
        # grounded or not
        if self.connect_type == 1:
            self.isGnd = True
        else:
            self.isGnd = False

        # Regulator Settings
        self.pt_ratio = pt_ratio

        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        aR_max = 1.1
        aR_min = 0.9
        self.aR = [1.0, 1.0, 1.0]
        self.aR_step = (aR_max - aR_min) / num_tap_steps
        self.band_center = band_center / self.pt_ratio
        self.bandwidth = bandwidth / self.pt_ratio
        self.Vreg_step = 120.0 * self.aR_step
        self.Vmax = (self.band_center + (self.bandwidth / 2))
        self.Vmin = (self.band_center - (self.bandwidth / 2))
        Va_reg = self.band_center
        Vb_reg = self.band_center
        Vc_reg = self.band_center
        if self.phases & 0x01 == int(0x01):
            Va_reg = self.band_center
        if self.phases & 0x02 == int(0x02):
            Vb_reg = self.band_center
        if self.phases & 0x04 == int(0x04):
            Vc_reg = self.band_center
        Vmag_reg = [Va_reg, Vb_reg, Vc_reg]
        self.raise_taps, self.lower_taps = self.check_regulator(Vmag_reg)
        if not self.lower_taps[0] and not self.raise_taps[0]:
            flag_tapA = 'hold'
        elif self.lower_taps[0]:
            flag_tapA = 'lower'
        else:
            flag_tapA = 'raise'

        if not self.lower_taps[1] and not self.raise_taps[1]:
            flag_tapB = 'hold'
        elif self.lower_taps[1]:
            flag_tapB = 'lower'
        else:
            flag_tapB = 'raise'

        if not self.lower_taps[2] and not self.raise_taps[2]:
            flag_tapC = 'hold'
        elif self.lower_taps[0]:
            flag_tapC = 'lower'
        else:
            flag_tapC = 'raise'

        self.rated_power = rated_power
        self.highstep = highstep
        self.lowstep = lowstep

        # Compensator Settings
        self.tapA = self.calculate_taps(Va_reg, flag_tapA)
        self.tapB = self.calculate_taps(Vb_reg, flag_tapB)
        self.tapC = self.calculate_taps(Vc_reg, flag_tapC)

        self.aR_control = RegulatorRatioControl(self.name, self.band_center, self.Vmax, self.Vmin, num_tap_steps,
                                                self.Vreg_step)

        self.tap_positions = [self.tapA, self.tapB, self.tapC]
        self.initial_taps = [self.tapA, self.tapB, self.tapC]
        self.compute_effective_regulator_ratios()

        # Assign Extra Nodes
        self.nodeA_Vr_from = None
        self.nodeA_Vr_to = None
        self.nodeA_Vi_from = None
        self.nodeA_Vi_to = None

        self.nodeB_Vr_from = None
        self.nodeB_Vr_to = None
        self.nodeB_Vi_from = None
        self.nodeB_Vi_to = None

        self.nodeC_Vr_from = None
        self.nodeC_Vr_to = None
        self.nodeC_Vi_from = None
        self.nodeC_Vi_to = None

        self.nodeN_Vr_from = None
        self.nodeN_Vi_from = None
        self.nodeN_Vr_to = None
        self.nodeN_Vi_to = None

        self.nodeA_Vr_secondary = None
        self.nodeA_Vi_secondary = None
        self.nodeA_Vr_primary = None
        self.nodeA_Vi_primary = None

        self.nodeB_Vr_secondary = None
        self.nodeB_Vi_secondary = None
        self.nodeB_Vr_primary = None
        self.nodeB_Vi_primary = None

        self.nodeC_Vr_secondary = None
        self.nodeC_Vi_secondary = None
        self.nodeC_Vr_primary = None
        self.nodeC_Vi_primary = None

        self.nodeGnd_Vr_primary = None
        self.nodeGnd_Vi_primary = None
        self.nodeGnd_Vr_secondary = None
        self.nodeGnd_Vi_secondary = None

        self.nodeA_Vr_loadCenter = None
        self.nodeA_Vi_loadCenter = None
        self.nodeA_Vmag2_loadCenter = None
        self.nodeA_V_diff = None

        self.nodeB_Vr_loadCenter = None
        self.nodeB_Vi_loadCenter = None
        self.nodeB_Vmag2_loadCenter = None
        self.nodeB_V_diff = None

        self.nodeC_Vr_loadCenter = None
        self.nodeC_Vi_loadCenter = None
        self.nodeC_Vmag2_loadCenter = None
        self.nodeC_V_diff = None

        self.nodeA_aR = None
        self.nodeB_aR = None
        self.nodeC_aR = None

        self.node_Ir_self = None
        self.node_Ii_self = None
        self.node_Ir_mutual1 = None
        self.node_Ii_mutual1 = None
        self.node_Ir_mutual2 = None
        self.node_Ii_mutual2 = None

        self.node_set = []
        self.node_Vr_meas = {'A': 0, 'B': 0, 'C': 0}
        self.node_Vi_meas = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_diff = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_aR = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_loadCenter = {'A': [], 'B': [], 'C': []}
        self.nodes_V_to = {'A': [], 'B': [], 'C': []}
        self.nodes_V_source = {'A': [], 'B': [], 'C': []}

        # Default loss terms
        self.G = 1e9
        self.B = 0

        # Count the nos. of changes
        self.taps_changed = False

        # Initialize Voltages
        self.Va_mag = 0
        self.Vb_mag = 0
        self.Vc_mag = 0
        self.Va_norm = 0
        self.Vb_norm = 0
        self.Vc_norm = 0
        if self.connect_type == self._DELTA_DELTA:
            self.Vab_mag = 0
            self.Vbc_mag = 0
            self.Vca_mag = 0

        # Initialize Currents
        self.Ia_mag = 0.0
        self.Ia_ang = 0.0
        self.Ib_mag = 0.0
        self.Ib_ang = 0.0
        self.Ic_mag = 0.0
        self.Ic_ang = 0.0

        # Initialiaze Effective Regulator Ratios
        self.aR_A = self.aR[0]
        self.aR_B = self.aR[1]
        self.aR_C = self.aR[2]

        self.isTriplex = False
        self.all_phases = {0: "A", 1: "B", 2: "C", 3: "N"}

    def assign_nodes(self,
                     node_key,
                     node_index_,
                     node):
        # There are likely to be 12 additional nodes for the XFMRs
        # 6 additional rows for the voltage source eqns. (No Phase Shifters
        # 12 additional rows for the voltage sources eqn (w Phase Shifters
        if self.phases & 0x1 == 1:  # Check for phase A
            self.nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr
            self.nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi
            self.nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
            self.nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi

            self.nodeA_aR = node_index_.__next__()
            self.nodes_aR["A"] = self.nodeA_aR
            self.nodeA_Vr_loadCenter = node_index_.__next__()
            self.nodeA_Vi_loadCenter = node_index_.__next__()
            self.nodeA_Vmag2_loadCenter = node_index_.__next__()
            self.nodeA_V_diff = node_index_.__next__()
            self.nodes_V_loadCenter["A"] = [self.nodeA_Vr_loadCenter, self.nodeA_Vi_loadCenter,
                                            self.nodeA_Vmag2_loadCenter]
            self.nodes_V_source["A"] = [self.nodeA_Vr_from, self.nodeA_Vi_from]
            self.nodes_V_to["A"] = [self.nodeA_Vr_to, self.nodeA_Vi_to]
            self.nodes_V_diff["A"] = self.nodeA_V_diff
            Nodes.Tap_index.extend([self.nodeA_aR, self.nodeA_Vr_loadCenter,
                                    self.nodeA_Vi_loadCenter, self.nodeA_V_diff])
            self.node_set.extend([self.nodeA_aR])

        if self.phases & 0x2 == 2:  # Check for phase B
            self.nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr
            self.nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi
            self.nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
            self.nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi

            self.nodeB_aR = node_index_.__next__()
            self.nodes_aR["B"] = self.nodeB_aR
            self.nodeB_Vr_loadCenter = node_index_.__next__()
            self.nodeB_Vi_loadCenter = node_index_.__next__()
            self.nodeB_Vmag2_loadCenter = node_index_.__next__()
            self.nodeB_V_diff = node_index_.__next__()
            self.nodes_V_loadCenter["B"] = [self.nodeB_Vr_loadCenter, self.nodeB_Vi_loadCenter,
                                            self.nodeB_Vmag2_loadCenter]
            self.nodes_V_source["B"] = [self.nodeB_Vr_from, self.nodeB_Vi_from]
            self.nodes_V_to["B"] = [self.nodeB_Vr_to, self.nodeB_Vi_to]
            self.nodes_V_diff["B"] = self.nodeB_V_diff
            Nodes.Tap_index.extend([self.nodeB_aR, self.nodeB_Vr_loadCenter,
                                    self.nodeB_Vi_loadCenter, self.nodeB_V_diff])
            self.node_set.extend([self.nodeB_aR])

        if self.phases & 0x4 == 4:  # Check for phase C
            self.nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
            self.nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi
            self.nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
            self.nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi

            self.nodeC_aR = node_index_.__next__()
            self.nodes_aR["C"] = self.nodeC_aR
            self.nodeC_Vr_loadCenter = node_index_.__next__()
            self.nodeC_Vi_loadCenter = node_index_.__next__()
            self.nodeC_Vmag2_loadCenter = node_index_.__next__()
            self.nodeC_V_diff = node_index_.__next__()
            self.nodes_V_loadCenter["C"] = [self.nodeC_Vr_loadCenter, self.nodeC_Vi_loadCenter,
                                            self.nodeC_Vmag2_loadCenter]
            self.nodes_V_source["C"] = [self.nodeC_Vr_from, self.nodeC_Vi_from]
            self.nodes_V_to["C"] = [self.nodeC_Vr_to, self.nodeC_Vi_to]
            self.nodes_V_diff["C"] = self.nodeC_V_diff
            Nodes.Tap_index.extend([self.nodeC_aR, self.nodeC_Vr_loadCenter,
                                    self.nodeC_Vi_loadCenter, self.nodeC_V_diff])
            self.node_set.extend([self.nodeC_aR])

        self.nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
        self.nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi
        self.nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
        self.nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi
        #############################################################
        # Add ground voltage source if needed
        if self.connect_type == self._GWYE_GWYE or self.connect_type == self._DELTA_DELTA:  # Wye-Wye
            # Add for the three phase extra node in the secondary
            # Add for the three phases extra node for the voltage sources (primarySource)
            # circuit (for both real and imaginary)
            if self.phases & 0x1 == 1:  # Check for phase A
                self.nodeA_Vr_secondary = node_index_.__next__()
                self.nodeA_Vi_secondary = node_index_.__next__()
                self.nodeA_Vr_primary = node_index_.__next__()
                self.nodeA_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["A"] = self.nodeA_Vr_primary
                self.node_Vi_meas["A"] = self.nodeA_Vi_primary
                Nodes.I_index.extend([self.nodeA_Vr_primary, self.nodeA_Vi_primary])
                self.node_set.extend([self.nodeA_Vr_secondary, self.nodeA_Vi_secondary, self.nodeA_Vr_primary,
                                      self.nodeA_Vi_primary])
            if self.phases & 0x2 == 2:  # Check for phase B
                self.nodeB_Vr_secondary = node_index_.__next__()
                self.nodeB_Vi_secondary = node_index_.__next__()
                self.nodeB_Vr_primary = node_index_.__next__()
                self.nodeB_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["B"] = self.nodeB_Vr_primary
                self.node_Vi_meas["B"] = self.nodeB_Vi_primary
                Nodes.I_index.extend([self.nodeB_Vr_primary, self.nodeB_Vi_primary])
                self.node_set.extend([self.nodeB_Vr_secondary, self.nodeB_Vi_secondary, self.nodeB_Vr_primary,
                                      self.nodeB_Vi_primary])
            if self.phases & 0x4 == 4:  # Check for phase C
                self.nodeC_Vr_secondary = node_index_.__next__()
                self.nodeC_Vi_secondary = node_index_.__next__()
                self.nodeC_Vr_primary = node_index_.__next__()
                self.nodeC_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["C"] = self.nodeC_Vr_primary
                self.node_Vi_meas["C"] = self.nodeC_Vi_primary
                Nodes.I_index.extend([self.nodeC_Vr_primary, self.nodeC_Vi_primary])
                self.node_set.extend([self.nodeC_Vr_secondary, self.nodeC_Vi_secondary, self.nodeC_Vr_primary,
                                      self.nodeC_Vi_primary])
            # Check for ground and add further nodes
            if self.isGnd:
                if (self.from_node,
                    self.to_node) not in self._regulator_locations:
                    self.nodeGnd_Vr_primary = node_index_.__next__()
                    self.nodeGnd_Vi_primary = node_index_.__next__()
                    self.nodeGnd_Vr_secondary = node_index_.__next__()
                    self.nodeGnd_Vi_secondary = node_index_.__next__()
                else:
                    self.nodeGnd_Vr_primary = self._regulator_locations[(
                        self.from_node, self.to_node)][0]
                    self.nodeGnd_Vi_primary = self._regulator_locations[(
                        self.from_node, self.to_node)][1]
                    self.nodeGnd_Vr_secondary = self._regulator_locations[(
                        self.from_node, self.to_node)][2]
                    self.nodeGnd_Vi_secondary = self._regulator_locations[(
                        self.from_node, self.to_node)][3]
        else:
            print('Incorrect regulator type for assigning nodes')

        return node_index_

    def check_connections(self):
        # check if regulator connections already exist
        if (self.from_node, self.to_node) not in self._regulator_locations:
            if self.connect_type == 1:
                self._regulator_locations[(self.from_node, self.to_node)] = [
                    self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary, self.nodeGnd_Vr_secondary,
                    self.nodeGnd_Vi_secondary
                ]
            else:
                pass

    def initialize_regulator(self,
                             Vinit):

        if self.phases & 0x1 == 1:  # Check for phase
            Vinit[self.nodeA_aR] = self.aR[0]
            Vinit[self.nodeA_Vr_loadCenter] = self.band_center
            Vinit[self.nodeA_Vi_loadCenter] = 0
            Vinit[self.nodeA_Vmag2_loadCenter] = self.band_center

        if self.phases & 0x2 == 2:  # Check for phase B
            Vbr = self.band_center * np.cos(np.deg2rad(-120))
            Vbi = self.band_center * np.sin(np.deg2rad(-120))
            Vinit[self.nodeB_aR] = self.aR[1]
            Vinit[self.nodeB_Vr_loadCenter] = Vbr
            Vinit[self.nodeB_Vi_loadCenter] = Vbi
            Vinit[self.nodeB_Vmag2_loadCenter] = np.sqrt(Vbr ** 2 + Vbi ** 2)

        if self.phases & 0x4 == 4:  # Check for phase C
            Vcr = self.band_center * np.cos(np.deg2rad(120))
            Vci = self.band_center * np.sin(np.deg2rad(120))
            Vinit[self.nodeC_aR] = self.aR[2]
            Vinit[self.nodeC_Vr_loadCenter] = Vcr
            Vinit[self.nodeC_Vi_loadCenter] = Vci
            Vinit[self.nodeC_Vmag2_loadCenter] = np.sqrt(Vcr ** 2 + Vci ** 2)

        return Vinit

    def compute_effective_regulator_ratios(self):
        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        if self.reg_type == 'A':  # Type A
            if self.raise_taps[0]:
                self.aR[0] = (1 + self.aR_step * self.tap_positions[0]) ** -1
            else:
                self.aR[0] = (1 - self.aR_step * self.tap_positions[0]) ** -1

            if self.raise_taps[1]:
                self.aR[1] = (1 + self.aR_step * self.tap_positions[1]) ** -1
            else:
                self.aR[1] = (1 - self.aR_step * self.tap_positions[1]) ** -1

            if self.raise_taps[2]:
                self.aR[2] = (1 + self.aR_step * self.tap_positions[2]) ** -1
            else:
                self.aR[2] = (1 - self.aR_step * self.tap_positions[2]) ** -1

        else:  # Type B
            if self.raise_taps[0]:
                self.aR[0] = 1 - self.aR_step * self.tap_positions[0]
            else:
                self.aR[0] = 1 + self.aR_step * self.tap_positions[0]

            if self.raise_taps[1]:
                self.aR[1] = 1 - self.aR_step * self.tap_positions[1]
            else:
                self.aR[1] = 1 + self.aR_step * self.tap_positions[1]

            if self.raise_taps[2]:
                self.aR[2] = 1 - self.aR_step * self.tap_positions[2]
            else:
                self.aR[2] = 1 + self.aR_step * self.tap_positions[2]

    def calculate_taps(self,
                       Vload,
                       change_type):
        if change_type == 'raise':
            Vlim = self.Vmin
            tap = np.around(np.abs(Vlim - Vload) / self.Vreg_step, 0)
        elif change_type == 'lower':
            Vlim = self.Vmax
            tap = np.around(np.abs(Vlim - Vload) / self.Vreg_step, 0)
        else:
            tap = np.around(np.abs(self.band_center - Vload) / self.Vreg_step, 0)

        return tap

    def check_regulator(self,
                        Vabc):

        self.Va_mag = np.abs(Vabc[0])
        self.Vb_mag = np.abs(Vabc[1])
        self.Vc_mag = np.abs(Vabc[2])

        # Set raise or lower tap flags depending on connection type
        raise_tapA = False
        lower_tapA = False
        raise_tapB = False
        lower_tapB = False
        raise_tapC = False
        lower_tapC = False
        if self.phases & 0x01 == int(0x01):
            raise_tapA = True if self.Va_mag < self.Vmin else False
            lower_tapA = True if self.Va_mag > self.Vmax else False
        if self.phases & 0x02 == int(0x02):
            # self.Vb_norm = self.Vb_mag / self.pt_ratio  # normalize voltage for 120 V base
            raise_tapB = True if self.Vb_mag < self.Vmin else False
            lower_tapB = True if self.Vb_mag > self.Vmax else False
        if self.phases & 0x04 == int(0x04):
            # self.Vc_norm = self.Vc_mag / self.pt_ratio  # normalize voltage for 120 V base
            raise_tapC = True if self.Vc_mag < self.Vmin else False
            lower_tapC = True if self.Vc_mag > self.Vmax else False
        raise_taps = [raise_tapA, raise_tapB, raise_tapC]
        lower_taps = [lower_tapA, lower_tapB, lower_tapC]

        return raise_taps, lower_taps

    @staticmethod
    def calc_G_B(r,
                 x):
        if r == 0:
            G = None
        else:
            G = r / (r ** 2 + x ** 2)
        if x == 0:
            B = None
        else:
            B = -x / (r ** 2 + x ** 2)
        return G, B

    def append_nodes(self):
        # Put together the set of nodes
        nodes_Vr_from = []
        nodes_Vi_from = []
        nodes_Vr_to = []
        nodes_Vi_to = []
        nodes_Vr_secondary = []
        nodes_Vi_secondary = []
        nodes_taps = []
        nodes_V_loadCenters = []
        nodes_Vr_relays = []
        nodes_Vi_relays = []
        nodes_Ir_comp = []
        nodes_Ii_comp = []
        nodes_Vr_load = []
        nodes_Vi_load = []
        Ir_sec = []
        Ii_sec = []

        # Find the from (Primary) bus nodes to stamp
        if self.phases & 0x01 == int(0x01):
            nodes_Vr_from.append(self.nodeA_Vr_from)
            nodes_Vi_from.append(self.nodeA_Vi_from)
            nodes_Vr_to.append(self.nodeA_Vr_to)
            nodes_Vi_to.append(self.nodeA_Vi_to)
            nodes_Vr_secondary.append(self.nodeA_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeA_Vi_secondary)

            Ir_sec.append(self.nodeA_Vr_primary)
            Ii_sec.append(self.nodeA_Vi_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            nodes_V_loadCenters.append(-1)
            nodes_Vr_relays.append(-1)
            nodes_Vi_relays.append(-1)
            nodes_Vr_load.append(-1)
            nodes_Vi_load.append(-1)
            nodes_taps.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)

        if self.phases & 0x02 == int(0x02):
            nodes_Vr_from.append(self.nodeB_Vr_from)
            nodes_Vi_from.append(self.nodeB_Vi_from)
            nodes_Vr_to.append(self.nodeB_Vr_to)
            nodes_Vi_to.append(self.nodeB_Vi_to)
            nodes_Vr_secondary.append(self.nodeB_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeB_Vi_secondary)

            Ir_sec.append(self.nodeB_Vr_primary)
            Ii_sec.append(self.nodeB_Vi_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            nodes_V_loadCenters.append(-1)
            nodes_Vr_relays.append(-1)
            nodes_Vi_relays.append(-1)
            nodes_Vr_load.append(-1)
            nodes_Vi_load.append(-1)
            nodes_taps.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)

        if self.phases & 0x04 == int(0x04):
            nodes_Vr_from.append(self.nodeC_Vr_from)
            nodes_Vi_from.append(self.nodeC_Vi_from)
            nodes_Vr_to.append(self.nodeC_Vr_to)
            nodes_Vi_to.append(self.nodeC_Vi_to)
            nodes_Vr_secondary.append(self.nodeC_Vr_secondary)
            nodes_Vi_secondary.append(self.nodeC_Vi_secondary)

            Ir_sec.append(self.nodeC_Vr_primary)
            Ii_sec.append(self.nodeC_Vi_primary)
        else:
            nodes_Vr_from.append(-1)
            nodes_Vi_from.append(-1)
            nodes_Vr_to.append(-1)
            nodes_Vi_to.append(-1)
            nodes_Vr_secondary.append(-1)
            nodes_Vi_secondary.append(-1)
            nodes_V_loadCenters.append(-1)
            nodes_Vr_relays.append(-1)
            nodes_Vi_relays.append(-1)
            nodes_Vr_load.append(-1)
            nodes_Vi_load.append(-1)
            nodes_taps.append(-1)
            Ir_sec.append(-1)
            Ii_sec.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)

        if self.connect_type == self._GWYE_GWYE:
            nodes_Vr_from.append(self.nodeN_Vr_from)
            nodes_Vi_from.append(self.nodeN_Vi_from)
            nodes_Vr_to.append(self.nodeN_Vr_to)
            nodes_Vi_to.append(self.nodeN_Vi_to)

        return nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
               nodes_Vr_secondary, nodes_Vi_secondary, nodes_taps, nodes_Vr_relays, \
               nodes_Vi_relays, nodes_V_loadCenters, nodes_Vr_load, nodes_Vi_load, nodes_Ir_comp, \
               nodes_Ii_comp, Ir_sec, Ii_sec

    def get_nodes(self,
                  node_key,
                  node,
                  V):
        """
        Find the indices of the solution vector and the values at those indices.
        :param node: The vector of all system nodes.
        :param V: The solution vector.
        :return: Three dictionaries which hold the phase and sequence nodes.
        """
        # Find the from bus nodes
        nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr
        nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi
        nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr
        nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi
        nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
        nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi

        # Find the to bus nodes
        nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
        nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi
        nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
        nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi
        nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
        nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi

        # Find the from and to node voltages
        (Var_from, Vai_from, Vbr_from, Vbi_from, Vcr_from, Vci_from) = \
            (V[nodeA_Vr_from], V[nodeA_Vi_from], V[nodeB_Vr_from], V[nodeB_Vi_from], V[nodeC_Vr_from],
             V[nodeC_Vi_from])

        (Var_to, Vai_to, Vbr_to, Vbi_to, Vcr_to, Vci_to) = \
            (V[nodeA_Vr_to], V[nodeA_Vi_to], V[nodeB_Vr_to], V[nodeB_Vi_to], V[nodeC_Vr_to],
             V[nodeC_Vi_to])

        # Find the relay voltages and Taps
        Va_source = 0.0
        Vb_source = 0.0
        Vc_source = 0.0
        Va_loadCenter = 0.0
        Vb_loadCenter = 0.0
        Vc_loadCenter = 0.0
        Va_mag2_loadCenter = 0.0
        Vb_mag2_loadCenter = 0.0
        Vc_mag2_loadCenter = 0.0
        Ia_prim = 0.0
        Ib_prim = 0.0
        Ic_prim = 0.0
        Ia_mag_prim = 0.0
        Ib_mag_prim = 0.0
        Ic_mag_prim = 0.0
        Va_neg_from = 0.0
        Vb_neg_from = 0.0
        Vc_neg_from = 0.0
        Va_neg_to = 0.0
        Vb_neg_to = 0.0
        Vc_neg_to = 0.0
        Va_diff = 0.0
        Vb_diff = 0.0
        Vc_diff = 0.0
        aR_A = self.aR_A
        aR_B = self.aR_B
        aR_C = self.aR_C

        Va_pos_sec = 0.0
        Vb_pos_sec = 0.0
        Vc_pos_sec = 0.0

        Va_neg_sec = 0.0
        Vb_neg_sec = 0.0
        Vc_neg_sec = 0.0
        Vn_neg_sec = 0.0

        if self.phases & 0x1 == 1:
            Va_source = complex(V[self.nodeA_Vr_from], V[self.nodeA_Vi_from])
            Va_loadCenter = complex(V[self.nodeA_Vr_loadCenter], V[self.nodeA_Vi_loadCenter])
            Va_mag2_loadCenter = float(V[self.nodeA_Vmag2_loadCenter])
            Ia_prim = complex(V[self.nodeA_Vr_primary], V[self.nodeA_Vi_primary])
            Va_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
            Va_diff = float(V[self.nodeA_V_diff])
            Va_neg_from = complex(V[self.nodeN_Vr_from], V[self.nodeN_Vi_from])
            aR_A = V[self.nodeA_aR][0]
            Va_neg_sec = complex(V[nodeA_Vr_to], V[nodeA_Vi_to])
            Va_pos_sec = complex(V[self.nodeA_Vr_secondary], V[self.nodeA_Vi_secondary])

        if self.phases & 0x2 == 2:
            Vb_source = complex(V[self.nodeB_Vr_from], V[self.nodeB_Vi_from])
            Vb_loadCenter = complex(V[self.nodeB_Vr_loadCenter], V[self.nodeB_Vi_loadCenter])
            Vb_mag2_loadCenter = float(V[self.nodeB_Vmag2_loadCenter])
            Ib_prim = complex(V[self.nodeB_Vr_primary], V[self.nodeB_Vi_primary])
            Vb_neg_from = complex(V[self.nodeN_Vr_from], V[self.nodeN_Vi_from])
            Vb_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
            aR_B = V[self.nodeB_aR][0]
            Vb_diff = float(V[self.nodeB_V_diff])
            Vb_neg_sec = complex(V[nodeB_Vr_to], V[nodeB_Vi_to])
            Vb_pos_sec = complex(V[self.nodeB_Vr_secondary], V[self.nodeB_Vi_secondary])

        if self.phases & 0x4 == 4:
            Vc_source = complex(V[self.nodeC_Vr_from], V[self.nodeC_Vi_from])
            Vc_loadCenter = complex(V[self.nodeC_Vr_loadCenter], V[self.nodeC_Vi_loadCenter])
            Vc_mag2_loadCenter = float(V[self.nodeC_Vmag2_loadCenter])
            Ic_prim = complex(V[self.nodeC_Vr_primary], V[self.nodeC_Vi_primary])
            Vc_neg_from = complex(V[self.nodeN_Vr_from], V[self.nodeN_Vi_from])
            Vc_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
            aR_C = V[self.nodeC_aR][0]
            Vc_diff = float(V[self.nodeC_V_diff])
            Vc_neg_sec = complex(V[nodeC_Vr_to], V[nodeC_Vi_to])
            Vc_pos_sec = complex(V[self.nodeC_Vr_secondary], V[self.nodeC_Vi_secondary])

        if self.connect_type == self._GWYE_GWYE:
            Vn_neg_sec = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])

        self.aR = [aR_A, aR_B, aR_C]

        Vabc = {
            'Vr_from': [Var_from, Vbr_from, Vcr_from],
            'Vi_from': [Vai_from, Vbi_from, Vci_from],
            'node_Vr_from': [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from],
            'node_Vi_from': [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from],
            'To': [complex(Var_to, Vai_to), complex(Vbr_to, Vbi_to), complex(Vcr_to, Vci_to)],
            'node_Vr_to': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to],
            'node_Vi_to': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to],
            'Vsource': [Va_source, Vb_source, Vc_source],
            'Vlc': [Va_loadCenter, Vb_loadCenter, Vc_loadCenter],
            'Vlc2': [Va_mag2_loadCenter, Vb_mag2_loadCenter, Vc_mag2_loadCenter],
            'Iprim': [Ia_prim, Ib_prim, Ic_prim],
            'Vneg_from': [Va_neg_from, Vb_neg_from, Vc_neg_from],
            'Vneg_to': [Va_neg_to, Vb_neg_to, Vc_neg_to],
            'aR': [aR_A, aR_B, aR_C],
            'Imag_prim': [Ia_mag_prim, Ib_mag_prim, Ic_mag_prim],
            'Vneg_sec': [Va_neg_sec, Vb_neg_sec, Vc_neg_sec, Vn_neg_sec],
            'Vpos_sec': [Va_pos_sec, Vb_pos_sec, Vc_pos_sec],
            'Vdiff': [Va_diff, Vb_diff, Vc_diff]

        }

        Vabc = SimpleNamespace(**Vabc)

        return Vabc

    def update_values(self,
                      V):
        # # Update Effective Regulator Ratios # #
        self.aR_A = np.around(V[self.nodeA_aR][0], 3) if hasattr(self, 'nodeA_aR') else 0.0
        self.aR_B = np.around(V[self.nodeB_aR][0], 3) if hasattr(self, 'nodeB_aR') else 0.0
        self.aR_C = np.around(V[self.nodeC_aR][0], 3) if hasattr(self, 'nodeC_aR') else 0.0

        m = (1.1 - 0.9) / (-16 - 16)
        b = 0.9 - (m * 16)

        self.tapA = int((self.aR_A - b) / m)
        self.tapB = int((self.aR_B - b) / m)
        self.tapC = int((self.aR_C - b) / m)

        # # Update Regulator Currents # #
        deg = True
        Ia = complex(V[self.nodeA_Vr_primary], V[self.nodeA_Vi_primary]) if hasattr(self, 'nodeA_Vr_primary') else 0
        self.Ia_mag = np.around(np.abs(Ia), 2)
        self.Ia_ang = np.around(np.angle(Ia, deg=deg), 2)
        Ib = complex(V[self.nodeB_Vr_primary], V[self.nodeB_Vi_primary]) if hasattr(self, 'nodeB_Vr_primary') else 0
        self.Ib_mag = np.around(np.abs(Ib), 2)
        self.Ib_ang = np.around(np.angle(Ib, deg=deg), 2)
        Ic = complex(V[self.nodeC_Vr_primary], V[self.nodeC_Vi_primary]) if hasattr(self, 'nodeC_Vr_primary') else 0
        self.Ic_mag = np.around(np.abs(Ic), 2)
        self.Ic_ang = np.around(np.angle(Ic, deg=deg), 2)

    def stamp_primary(self,
                      Ylin_val,
                      Ylin_row,
                      Ylin_col,
                      idx_Y,
                      node_Vr_pos_from,
                      node_Vi_pos_from,
                      node_Vr_neg_from,
                      node_Vi_neg_from,
                      node_Vr_neg_secondary,
                      node_Vi_neg_secondary,
                      node_Vr_pos_secondary,
                      node_Vi_pos_secondary,
                      node_Vr_primary,
                      node_Vi_primary,
                      phase,
                      homotopy=False,
                      h_factor=None):

        # Get Effective Regulator Ratio of the Corresponding Phase
        if not homotopy:
            aR = self.aR[phase]
        else:
            aR = self.aR[phase] + h_factor * (1 - self.aR[phase])
        # aR = self.aR[phase] + h_factor * (1 - self.aR[phase])

        # # PRIMARY VOLTAGES # #
        # Real #
        # Vr_primary = Vr_pos_from - aR * Vr_secondary_pos  + aR * Vr_secondary_neg - Vr_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected

        idx_Y = self.stamp_Y(node_Vr_primary, node_Vr_pos_from, 1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vr_primary, node_Vr_pos_secondary, -aR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vr_primary, node_Vr_neg_secondary, aR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vr_primary, node_Vr_neg_from, -1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        # Imaginary #
        # Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        idx_Y = self.stamp_Y(node_Vi_primary, node_Vi_pos_from, 1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vi_primary, node_Vi_pos_secondary, -aR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vi_primary, node_Vi_neg_from, -1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vi_primary, node_Vi_neg_secondary, aR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)

        # PRIMARY CURRENTS
        # Real #
        idx_Y = self.stamp_Y(node_Vr_pos_from, node_Vr_primary, 1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vr_neg_from, node_Vr_primary, -1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        # Imaginary #
        idx_Y = self.stamp_Y(node_Vi_pos_from, node_Vi_primary, 1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vi_neg_from, node_Vi_primary, -1, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        return idx_Y

    def stamp_secondary(self,
                        Ylin_val,
                        Ylin_row,
                        Ylin_col,
                        idx_Y,
                        node_Vr_pos_to,
                        node_Vi_pos_to,
                        node_Vr_neg_secondary,
                        node_Vi_neg_secondary,
                        node_Vr_pos_secondary,
                        node_Vi_pos_secondary,
                        node_Vr_primary,
                        node_Vi_primary,
                        phase,
                        homotopy=False,
                        h_factor=None,
                        B_factor=None,
                        G_factor=None):

        # Calculate conductance and susceptance of regulator if homotopy
        if not homotopy:
            # # Get Effective Regulator Ratio of the Corresponding Phase # #
            aR = self.aR[phase]
            _B = self.B
            _G = self.G
        else:
            aR = self.aR[phase] + h_factor * (1 - self.aR[phase])
            # aR = self.aR[phase]
            # _B = B_factor * self.B
            # _G = G_factor * self.G
            _B = self.B
            _G = self.G

        if aR:
            # # SECONDARY CURRENT SOURCES # #
            # Real #
            # Ir_pos_secondary = - aR * Ir_primary
            # Ir_neg_secondary = aR * Ir_primary
            idx_Y = self.stamp_Y(node_Vr_pos_secondary, node_Vr_primary, -aR, Ylin_val,
                                 Ylin_row, Ylin_col, idx_Y)
            idx_Y = self.stamp_Y(node_Vr_neg_secondary, node_Vr_primary, aR, Ylin_val,
                                 Ylin_row, Ylin_col, idx_Y)

            # Imaginary #
            # Ii_pos_secondary = - aR * Ii_primary
            # Ii_neg_secondary = aR * Ii_primary
            idx_Y = self.stamp_Y(node_Vi_pos_secondary, node_Vi_primary, -aR, Ylin_val,
                                 Ylin_row, Ylin_col, idx_Y)
            idx_Y = self.stamp_Y(node_Vi_neg_secondary, node_Vi_primary, aR, Ylin_val,
                                 Ylin_row, Ylin_col, idx_Y)

        if _G:
            # Branch stamping for the conductance
            #
            idx_Y = self.stamp_Y(
                node_Vr_pos_to, node_Vr_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vr_pos_to, node_Vr_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vi_pos_to, node_Vi_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vi_pos_to, node_Vi_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vr_pos_secondary, node_Vr_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vr_pos_secondary, node_Vr_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vi_pos_secondary, node_Vi_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(
                node_Vi_pos_secondary, node_Vi_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        return idx_Y

    def stamp_aR(self,
                 phase,
                 Vlc,
                 Ynlin_val,
                 Ynlin_row,
                 Ynlin_col,
                 Jnlin_val,
                 Jnlin_row,
                 idx_Y,
                 idx_J,
                 homotopy=False,
                 h_factor=None):

        # Vr_lc, Vi_lc = Vlc.real, Vlc.imag
        Vlc_mag = Vlc

        if not homotopy:
            aR = self.aR[phase]
        else:
            aR = self.aR[phase] + h_factor * (1 - self.aR[phase])

        # Compute partials and historical aR value
        (dFaR, FaR_k, flag_aR_max,
         flag_aR_min, flag_aR_lim) = self.aR_control.compute_linearized_aR(Vlc_mag, aR, phase)

        # Collect nodes
        node_Vmag_loadCenter = self.nodes_V_loadCenter[self.all_phases[phase]][2]
        node_aR = self.nodes_aR[self.all_phases[phase]]

        if flag_aR_max and flag_aR_lim:
            aR_max = 1.1
            idx_J = self.stamp_J(node_aR, aR_max, Jnlin_val, Jnlin_row,
                                 idx_J)

        elif flag_aR_min and flag_aR_lim:
            aR_min = 0.9
            idx_J = self.stamp_J(node_aR, aR_min, Jnlin_val, Jnlin_row,
                                 idx_J)
        else:
            idx_Y = self.stamp_Y(node_aR, node_Vmag_loadCenter, dFaR.dV[phase], Ynlin_val, Ynlin_row, Ynlin_col,
                                 idx_Y)

            FaR = FaR_k - dFaR.daR * aR - dFaR.dV[phase] * Vlc

            idx_J = self.stamp_J(node_aR, -FaR, Jnlin_val, Jnlin_row, idx_J)

        idx_Y = self.stamp_Y(node_aR, node_aR, dFaR.daR, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

        return idx_Y, idx_J

    def stamp_aR_primary(self,
                         Ylin_val,
                         Ylin_row,
                         Ylin_col,
                         idx_Y,
                         node_Vr_primary,
                         node_Vi_primary,
                         V_secondary_pos,
                         V_secondary_neg,
                         phase):

        # Collect the aR node and the secondary voltages
        node_aR = self.nodes_aR[self.all_phases[phase]]

        Vr_secondary_pos, Vi_secondary_pos = V_secondary_pos.real, V_secondary_pos.imag
        Vr_secondary_neg, Vi_secondary_neg = V_secondary_neg.real, V_secondary_neg.imag

        # Real Stamps #
        # Vr_primary = Vr_pos_from - aR * Vr_secondary_pos  + aR * Vr_secondary_neg - Vr_neg_from
        dVr_primary_daR = -Vr_secondary_pos + Vr_secondary_neg

        idx_Y = self.stamp_Y(node_Vr_primary, node_aR, dVr_primary_daR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)

        # Imaginary Stamps#
        # Vi_primary = Vi_neg_from - aR * Vi_secondary_pos  + aR * Vi_secondary_neg - Vi_neg_from
        dVi_primary_daR = -Vi_secondary_pos + Vi_secondary_neg
        idx_Y = self.stamp_Y(node_Vi_primary, node_aR, dVi_primary_daR, Ylin_val,
                             Ylin_row, Ylin_col, idx_Y)

        return idx_Y

    def stamp_aR_secondary(self,
                           Ynlin_val,
                           Ynlin_row,
                           Ynlin_col,
                           idx_Y,
                           node_Vr_pos_secondary,
                           node_Vi_pos_secondary,
                           node_Vr_neg_secondary,
                           node_Vi_neg_secondary,
                           I_primary,
                           phase):

        # Collect the aR node and the primary currents
        node_aR = self.nodes_aR[self.all_phases[phase]]
        Ir_primary, Ii_primary = I_primary.real, I_primary.imag

        # Real #
        # Ir_pos_secondary = - aR * Ir_primary
        # Ir_neg_secondary = aR * Ir_primary
        dIr_pos_sec_daR = -Ir_primary
        dIr_neg_sec_daR = Ir_primary

        idx_Y = self.stamp_Y(node_Vr_pos_secondary, node_aR, dIr_pos_sec_daR,
                             Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vr_neg_secondary, node_aR, dIr_neg_sec_daR,
                             Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

        # Imaginary #
        # Ii_pos_secondary = - aR * Ii_primary
        # Ii_neg_secondary = aR * Ii_primary
        dIi_pos_sec_daR = - Ii_primary
        dIi_neg_sec_daR = Ii_primary

        idx_Y = self.stamp_Y(node_Vi_pos_secondary, node_aR, dIi_pos_sec_daR,
                             Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

        idx_Y = self.stamp_Y(node_Vi_neg_secondary, node_aR, dIi_neg_sec_daR,
                             Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

        return idx_Y

    def stamp_ground_neutral(self,
                             neg_real_node,
                             neg_imag_node,
                             gnd_real_V,
                             gnd_imag_V,
                             Ylin_val,
                             Ylin_row,
                             Ylin_col,
                             idx_Y):
        # Stamp Vng = 0 + j0
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(gnd_real_V, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col,
                                                           idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(gnd_imag_V, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col,
                                                           idx_Y)
        #
        # Stamp Current Sources
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(neg_real_node, gnd_real_V, 1, Ylin_val, Ylin_row, Ylin_col,
                                                           idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_Y(neg_imag_node, gnd_imag_V, 1, Ylin_val, Ylin_row, Ylin_col,
                                                           idx_Y)
        #

        return idx_Y

    def stamp_output_voltage_control(self,
                                     Ylin_val,
                                     Ylin_row,
                                     Ylin_col,
                                     Jlin_val,
                                     Jlin_row,
                                     idx_Y,
                                     idx_J,
                                     Vmag_lc,
                                     Vlc,
                                     phase):

        node_Vr_substation = self.nodes_V_source[self.all_phases[phase]][0]
        node_Vi_substation = self.nodes_V_source[self.all_phases[phase]][1]
        node_Vr_loadCenter = self.nodes_V_loadCenter[self.all_phases[phase]][0]
        node_Vi_loadCenter = self.nodes_V_loadCenter[self.all_phases[phase]][1]
        node_Vmag_loadcenter = self.nodes_V_loadCenter[self.all_phases[phase]][2]

        # Ir_self, Ii_self = I_self.real, I_self.imag
        # Ir_mutual1, Ii_mutual1 = I_mutual1.real, I_mutual1.imag
        # Vr_source, Vi_source = Vsource.real, Vsource.imag
        Vr_lc, Vi_lc = Vlc.real, Vlc.imag

        # - Load Center (LC) Voltage Partials - #
        dVrLC_dVrLC = -1
        dVrLC_dVrSub = 1 / self.pt_ratio

        # _Vr_lc = Vr_hist_lc - dVrLC_dVrLC * Vr_lc - dVrLC_dVrSub * Vr_source - dVrLC_dIrSelf * Ir_self - \
        #          dVrLC_dIiSelf * Ii_self - dVrLC_dIrMutual1 * Ir_mutual1 - dVrLC_dIiMutual1 * Ii_mutual1 - \
        #          dVrLC_dIrMutual2 * Ir_mutual2 - dVrLC_dIiMutual2 * Ii_mutual1

        idx_Y = self.stamp_Y(node_Vr_loadCenter, node_Vr_loadCenter, dVrLC_dVrLC, Ylin_val, Ylin_row, Ylin_col,
                             idx_Y)
        idx_Y = self.stamp_Y(node_Vr_loadCenter, node_Vr_substation, dVrLC_dVrSub, Ylin_val, Ylin_row,
                             Ylin_col,
                             idx_Y)

        dViLC_dViLC = -1
        dViLC_dViSub = 1 / self.pt_ratio

        # _Vi_lc = Vi_hist_lc - dViLC_dViLC * Vi_lc - dViLC_dViSub * Vi_source - dViLC_dIrSelf * Ir_self - \
        #          dViLC_dIrSelf * Ir_self - dViLC_dIrMutual1 * Ir_mutual1 - dViLC_dIiMutual1 * Ii_mutual1 - \
        #          dViLC_dIrMutual2 * Ir_mutual2 - dViLC_dIiMutual2 * Ii_mutual2

        idx_Y = self.stamp_Y(node_Vi_loadCenter, node_Vi_loadCenter, dViLC_dViLC, Ylin_val, Ylin_row, Ylin_col,
                             idx_Y)
        idx_Y = self.stamp_Y(node_Vi_loadCenter, node_Vi_substation, dViLC_dViSub, Ylin_val, Ylin_row,
                             Ylin_col,
                             idx_Y)

        Vmag_lc_hist = Vr_lc ** 2 + Vi_lc ** 2 - Vmag_lc ** 2
        dVmag_dVr = 2 * Vr_lc
        dVmag_dVi = 2 * Vi_lc
        dVmag_dVmag = -2 * Vmag_lc

        _Vmag_lc = Vmag_lc_hist - dVmag_dVr * Vr_lc - dVmag_dVi * Vi_lc - dVmag_dVmag * Vmag_lc

        idx_Y = self.stamp_Y(node_Vmag_loadcenter, node_Vr_loadCenter, dVmag_dVr, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vmag_loadcenter, node_Vi_loadCenter, dVmag_dVi, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)
        idx_Y = self.stamp_Y(node_Vmag_loadcenter, node_Vmag_loadcenter, dVmag_dVmag, Ylin_val, Ylin_row,
                             Ylin_col, idx_Y)

        idx_J = self.stamp_J(node_Vmag_loadcenter, -_Vmag_lc, Jlin_val, Jlin_row, idx_J)

        return idx_Y, idx_J

    def stamp_voltage_difference_from_limit(self,
                                            Ynlin_Val,
                                            Ynlin_row,
                                            Ynlin_col,
                                            Jnlin_val,
                                            Jnlin_row,
                                            idx_Y,
                                            idx_J,
                                            Vmag_lc,
                                            Vdiff,
                                            phase):

        node_Vmag_loadCenter = self.nodes_V_loadCenter[self.all_phases[phase]][2]
        node_V_diff = self.nodes_V_diff[self.all_phases[phase]]

        # Collect the kth iteration load center voltage
        # Vr_lc, Vi_lc = Vlc.real, Vlc.imag

        # Determine the value of the limit voltage
        # The object is to always maintain the load center or relay voltage between the compensator band limits

        if Vmag_lc < self.Vmin:
            V_limit = self.Vmin
        elif Vmag_lc > self.Vmax:
            V_limit = self.Vmax
        else:
            V_limit = self.band_center

        # # Voltage Difference From Limit Stamps # #

        V_diff_hist = V_limit - Vmag_lc - Vdiff
        dVdiff_dVdiff = -1
        dVdiff_dVmag_lc = -1

        _Vdiff = V_diff_hist - dVdiff_dVmag_lc * Vmag_lc - dVdiff_dVdiff * Vdiff

        idx_Y = self.stamp_Y(node_V_diff, node_V_diff, dVdiff_dVdiff, Ynlin_Val, Ynlin_row, Ynlin_col, idx_Y)
        idx_Y = self.stamp_Y(node_V_diff, node_Vmag_loadCenter, dVdiff_dVmag_lc, Ynlin_Val, Ynlin_row, Ynlin_col,
                             idx_Y)

        idx_J = self.stamp_J(node_V_diff, - _Vdiff, Jnlin_val, Jnlin_row, idx_J)

        return idx_Y, idx_J

    def stamp_linear_ground(self,
                            Ynlin_val,
                            Ynlin_row,
                            Ynlin_col,
                            idx_Y,
                            stamped_ground):

        _neg = 3

        nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
        nodes_Vr_secondary, nodes_Vi_secondary, nodes_taps, nodes_Vr_relay, \
        nodes_Vi_relay, nodes_V_loadCenter, nodes_Vr_load, nodes_Vi_load, nodes_Ir_comp, \
        nodes_Ii_comp, nodes_Vr_primary, nodes_Vi_primary, L_dict = self.append_nodes()

        if self.from_node not in stamped_ground:
            idx_Y = self.stamp_ground_neutral(
                nodes_Vr_from[_neg], nodes_Vi_from[_neg], self.nodeGnd_Vr_primary,
                self.nodeGnd_Vi_primary, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            stamped_ground.add(self.from_node)
        # # Stamp secondary
        #
        if self.to_node not in stamped_ground:
            idx_Y = self.stamp_ground_neutral(
                nodes_Vr_to[_neg], nodes_Vi_to[_neg], self.nodeGnd_Vr_secondary,
                self.nodeGnd_Vi_secondary, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            stamped_ground.add(self.to_node)

        return idx_Y, stamped_ground

    def stamp_nonlinear(self,
                        node_key,
                        node,
                        V,
                        Ynlin_val,
                        Ynlin_row,
                        Ynlin_col,
                        Jnlin_val,
                        Jnlin_row,
                        idx_Y,
                        idx_J,
                        homotopy_option=False,
                        h_factor=None,
                        G_homotopy=None,
                        B_homotopy=None):

        V_k = self.get_nodes(node_key, node, V)

        nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
        nodes_Vr_secondary, nodes_Vi_secondary, nodes_taps, nodes_Vr_relay, \
        nodes_Vi_relay, nodes_V_loadCenter, nodes_Vr_load, nodes_Vi_load, nodes_Ir_comp, \
        nodes_Ii_comp, nodes_Vr_primary, nodes_Vi_primary, L_dict = self.append_nodes()

        (_A, _B, _C) = (0, 1, 2)
        _N = 3

        try:

            # Option 1 - Stamp gwye-gwye
            if self.connect_type == self._GWYE_GWYE:
                _phases_list = []
                if self.phases & 0x01 == int(0x01):
                    _phases_list.append((_A, _N))

                if self.phases & 0x02 == int(0x02):
                    _phases_list.append((_B, _N))

                if self.phases & 0x04 == int(0x04):
                    _phases_list.append((_C, _N))

                # # Stamp Effective Regulator Ratio Control Circuit and the Load Center Voltage # #
                for phaseSet in _phases_list:
                    _phase = phaseSet[0]

                    # # Stamp Voltage Difference # #
                    idx_Y, idx_J = self.stamp_output_voltage_control(Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val,
                                                                     Jnlin_row, idx_Y, idx_J,
                                                                     V_k.Vlc2[_phase], V_k.Vlc[_phase], _phase)

                    idx_Y, idx_J = self.stamp_voltage_difference_from_limit(Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val,
                                                                            Jnlin_row, idx_Y, idx_J, V_k.Vlc2[_phase],
                                                                            V_k.Vdiff[_phase], _phase)
                    # # Stamp Effective Regulator Ratio # #
                    if not homotopy_option:
                        idx_Y, idx_J = self.stamp_aR(_phase, V_k.Vlc2[_phase], Ynlin_val, Ynlin_row, Ynlin_col,
                                                     Jnlin_val,
                                                     Jnlin_row, idx_Y, idx_J)
                    else:
                        idx_Y, idx_J = self.stamp_aR(_phase, V_k.Vlc2[_phase], Ynlin_val, Ynlin_row, Ynlin_col,
                                                     Jnlin_val,
                                                     Jnlin_row, idx_Y, idx_J, homotopy_option, h_factor)

                for phaseSet in _phases_list:
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]

                    ##########Stamp Primary Circuit##################
                    # from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]

                    #
                    # ##########Stamp Secondary Circuit##################
                    idx_Y = self.stamp_primary(Ynlin_val, Ynlin_row, Ynlin_col,
                                               idx_Y, nodes_Vr_from[_pos],
                                               nodes_Vi_from[_pos],
                                               nodes_Vr_from[_neg],
                                               nodes_Vi_from[_neg],
                                               nodes_Vr_to[_neg],
                                               nodes_Vi_to[_neg],
                                               nodes_Vr_secondary[_winding],
                                               nodes_Vi_secondary[_winding],
                                               nodes_Vr_primary[_winding],
                                               nodes_Vi_primary[_winding],
                                               phaseSet[0])

                    if not homotopy_option:
                        idx_Y = self.stamp_secondary(Ynlin_val, Ynlin_row, Ynlin_col,
                                                     idx_Y,
                                                     nodes_Vr_to[_pos],
                                                     nodes_Vi_to[_pos], nodes_Vr_to[_neg],
                                                     nodes_Vi_to[_neg],
                                                     nodes_Vr_secondary[_winding],
                                                     nodes_Vi_secondary[_winding],
                                                     nodes_Vr_primary[_winding],
                                                     nodes_Vi_primary[_winding],
                                                     phaseSet[0])
                    else:
                        G_factor = G_homotopy * h_factor
                        B_factor = B_homotopy * h_factor
                        idx_Y = self.stamp_secondary(Ynlin_val, Ynlin_row, Ynlin_col,
                                                     idx_Y,
                                                     nodes_Vr_to[_pos],
                                                     nodes_Vi_to[_pos], nodes_Vr_to[_neg],
                                                     nodes_Vi_to[_neg],
                                                     nodes_Vr_secondary[_winding],
                                                     nodes_Vi_secondary[_winding],
                                                     nodes_Vr_primary[_winding],
                                                     nodes_Vi_primary[_winding],
                                                     phaseSet[0], homotopy_option,
                                                     h_factor, B_factor, G_factor)

            # _, _, _, idx_Y = self.stamp_aR_primary(Ynlin_val, Ynlin_row, Ynlin_col,
            #                                        idx_Y, nodes_Vr_primary[_winding],
            #                                        nodes_Vi_primary[_winding],
            #                                        V_k.Vpos_sec[_winding],
            #                                        V_k.Vneg_sec[_neg], _winding)
            #
            # _, _, _, idx_Y = self.stamp_aR_secondary(Ynlin_val, Ynlin_row, Ynlin_col,
            #                                          idx_Y,
            #                                          nodes_Vr_secondary[_winding],
            #                                          nodes_Vi_secondary[_winding],
            #                                          nodes_Vr_to[_neg],
            #                                          nodes_Vi_to[_neg],
            #                                          V_k.Iprim[_winding],
            #                                          _winding)

            elif self.connect_type == self._DELTA_DELTA:
                _phases = [(_A, _B), (_B, _C), (_C, _A)]
                for phaseSet in _phases:
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]

                    # Update Tap Positions and Effective Regulator Ratios #
                    if self.phases & 0x01 == int(0x01):
                        self.tapA = V.Tap[0]
                    if self.phases & 0x02 == int(0x02):
                        self.tapB = V.Tap[1]
                    if self.phases & 0x04 == int(0x04):
                        self.tapC = V.Tap[2]
                    self.tap_positions = [self.tapA, self.tapB, self.tapC]
                    self.compute_effective_regulator_ratios()

        except ValueError:
            print('Illegal connect_type for the regulator')

        return idx_Y, idx_J
