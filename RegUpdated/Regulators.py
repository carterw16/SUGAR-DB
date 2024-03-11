"""Creates regulators for distribution system analysis.

  Author(s): Amrit Pandey, Naeem Turner-Bandele, Elizabeth Foster
  Created Date: 04-11-2017
  Updated Date: 2-26-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu, emfoster@andrew.cmu.edu
  Status: Implementing infeasibility

  self consist of an autotransformer and a load tap changing (LTC) mechanism. They are used to regulate the
  voltage on a particular bus. This class  specifies the core properties of regulator objects and the stamps for
  regulators. See W.H. Kersting Distribution System Analysis for details on how regulators work.
"""

from __future__ import division

from itertools import count
from types import SimpleNamespace

import numpy as np

from classes.Elements import LinearElement
from classes.Nodes import Nodes


class Regulators(LinearElement):
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
    _DELTA_DELTA = 2  # Validated
    _ids = count(0)
    _regulator_locations = dict()

    def __init__(self, node_index_, node, name, ID, phases, from_node, to_node,
                 connect_type, reg_type, nom_voltage, band_center, bandwidth, rated_power,
                 fixed_taps, tapA, tapB, tapC, lowstep, highstep, num_tap_steps,
                 aR_max, aR_min, pt_ratio, pt_phase, ct_prim, ct_phase, connected_transformer,
                 connected_winding, tap_ctrl_scheme, tap_ctrl_fcn, tap_ctrl_smoothing, 
                 compensator_resistance = 1, compensator_reactance = 1, stamp_dual = False):
                 

        super(Regulators, self).__init__()
        
        # Basic Properties
        self.name = name
        self.ID = ID
        self.from_node = from_node
        self.to_node = to_node
        self.connected_transformer = connected_transformer
        self.connected_winding = connected_winding
        self.reg_type = reg_type
        self.connect_type = connect_type
        self.Vnom = nom_voltage
        self.phases = phases if phases else (self.pt_phase | self.ct_phase)
        
        self.stamp_dual = stamp_dual
        
        # grounded or not
        if self.connect_type == 1:
            self.isGnd = True
        else:
            self.isGnd = False
        
        # Regulator Settings
        self.band_center = band_center
        self.bandwidth = bandwidth
        self.pt_ratio = pt_ratio
        self.ct_prim = ct_prim
        # Regulator Control Phase(s)
        self.pt_phase = pt_phase
        self.ct_phase = ct_phase
        # default value chosen based on W.H. Kersting, Distribution System Modeling and Analysis guidance
        self.ct_sec = 5  
        
        # Tap Values
        self.fixed_taps = fixed_taps
        if self.fixed_taps: 
            self.tapA = tapA
            self.tapB = tapB
            self.tapC = tapC
        
        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        self.aR_max = aR_max
        self.aR_min = aR_min
        self.aR = [1.0, 1.0, 1.0]
        self.aR_step = (self.aR_max - self.aR_min) / num_tap_steps
        Vnorm = (self.Vnom / np.sqrt(3)) / self.pt_ratio
        self.Vbase = Vnorm if np.isclose(Vnorm, 120.0,
                                         1e-1) else self.Vnom / self.pt_ratio
        self.Vreg_step = 120.0 * self.aR_step
        self.Vmax = (self.band_center + (self.bandwidth / 2))
        self.Vmin = (self.band_center - (self.bandwidth / 2))
        Va_reg = self.band_center
        Vb_reg = self.band_center
        Vc_reg = self.band_center
        if self.phases & 0x01 == int(0x01):
            Va_reg = self.Vbase
        if self.phases & 0x02 == int(0x02):
            Vb_reg = self.Vbase
        if self.phases & 0x04 == int(0x04):
            Vc_reg = self.Vbase
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

        if not self.fixed_taps:
            self.tapA = self.calculate_taps(Va_reg, flag_tapA)
            self.tapB = self.calculate_taps(Vb_reg, flag_tapB)
            self.tapC = self.calculate_taps(Vc_reg, flag_tapC)

        self.rated_power = rated_power
        self.highstep = highstep
        self.lowstep = lowstep

        self.tap_positions = [self.tapA, self.tapB, self.tapC]
        self.initial_taps = [self.tapA, self.tapB, self.tapC]
        self.compute_effective_regulator_ratios()
        
        # Compensator Settings
        if not self.fixed_taps:
            self.Y_series = np.zeros((3, 3), dtype=complex)
            self.compensator_R = compensator_resistance / self.ct_sec if compensator_resistance else 0
            self.compensator_X = compensator_reactance / self.ct_sec if compensator_reactance else 0
            self.compensator_Z = complex(self.compensator_R, self.compensator_X)
            self.G_comp, self.B_comp = self.calc_G_B(self.compensator_R,
                                                     self.compensator_X)
            self.ct = self.ct_prim / self.ct_sec
            self.tap_control_scheme = tap_ctrl_scheme
            self.tap_control_function = tap_ctrl_fcn
            self.smoothing_factor = tap_ctrl_smoothing

        # Assign Extra Nodes
        self.node_Ir_self = -1
        self.node_Ii_self = -1
        self.node_Ir_mutual1 = -1
        self.node_Ii_mutual1 = -1
        self.node_Ir_mutual2 = -1
        self.node_Ii_mutual2 = -1
        self.node_set = []
        self.node_Vr_meas = {'A': 0, 'B': 0, 'C': 0}
        self.node_Vi_meas = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_diff = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_aR = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_Tap = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_ImagPrim = {'A': 0, 'B': 0, 'C': 0}
        self.nodes_V_loadCenter = {'A': [], 'B': [], 'C': []}
        self.nodes_V_load = {'A': [], 'B': [], 'C': []}
        self.nodes_V_to = {'A': [], 'B': [], 'C': []}
        self.nodes_V_source = {'A': [], 'B': [], 'C': []}
        self.nodes_I_comp = {'A': [0, 0, 0], 'B': [0, 0, 0], 'C': [0, 0, 0]}
        if self.stamp_dual:
            self.node_Lr_meas = {'A': 0, 'B': 0, 'C': 0}
            self.node_Li_meas = {'A': 0, 'B': 0, 'C': 0}

        self.assign_nodes(node_index_, node)
        

        # Default loss terms
        self.G = 1e9
        self.B = 0

        # check if regulator connections already exist
        if (self.from_node, self.to_node) not in self._regulator_locations:
            if self.connect_type == 1:
                self._regulator_locations[(self.from_node, self.to_node)] = [
                    self.nodeGnd_Vr_primary, self.nodeGnd_Vi_primary,
                    self.nodeGnd_Vr_secondary, self.nodeGnd_Vi_secondary
                ]
            else:
                pass

        # Count the nos. of changes
        self.tapA_change_count = count(0)
        self.tapB_change_count = count(0)
        self.tapC_change_count = count(0)
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

        self.update_lin_stamp = False
        self.isTriplex = False
        self.all_phases = {0: "A", 1: "B", 2: "C", 3: "N"}
        ## Temporary Test Parameters ##
        self.Vrelay = [113.5, 115.0, 112.75]
        self.Vdrop = [
            complex(6.88, 6.529),
            complex(2.863, -7.611),
            complex(-10.525, 2.604)
        ]
        self.Vreg = [
            complex(120.08, 0),
            complex(-60, -104.025),
            complex(-60, 104.025)
        ]
        self.Icomp = [
            complex(3.752643550898216, -0.036653988983750094),
            complex(-2.4162429679078747, -0.13354845573934707),
            complex(-0.35342170731556494, 4.452853672813877)
        ]
        
    def initialize_regulator(self, Vinit):

        if not self.fixed_taps:
            if self.phases & 0x1 == 1:  # Check for phase
                Vinit[self.nodeA_aR] = self.aR[0]
                Vinit[self.nodeA_Tap] = self.tapA
                Vinit[self.nodeA_Vr_loadCenter] = self.band_center
                Vinit[self.nodeA_Vi_loadCenter] = 0

            if self.phases & 0x2 == 2:  # Check for phase B
                Vinit[self.nodeB_aR] = self.aR[1]
                Vinit[self.nodeB_Tap] = self.tapB
                Vinit[self.nodeB_Vr_loadCenter] = self.band_center * np.cos(
                    np.deg2rad(-120))
                Vinit[self.nodeB_Vi_loadCenter] = self.band_center * np.sin(
                    np.deg2rad(-120))
            if self.phases & 0x4 == 4:  # Check for phase C
                Vinit[self.nodeC_aR] = self.aR[2]
                Vinit[self.nodeC_Tap] = self.tapC
                Vinit[self.nodeC_Vr_loadCenter] = self.band_center * np.cos(
                    np.deg2rad(120))
                Vinit[self.nodeC_Vi_loadCenter] = self.band_center * np.sin(
                    np.deg2rad(120))
        return Vinit
    
    def compute_effective_regulator_ratios(self):
        # Calculate aR
        # For raise position use minus sign and for lower use + sign
        if self.reg_type == 'A':  # Type A
            if self.raise_taps[0]:
                self.aR[0] = (1 + self.aR_step * self.tap_positions[0])**-1
            else:
                self.aR[0] = (1 - self.aR_step * self.tap_positions[0])**-1

            if self.raise_taps[1]:
                self.aR[1] = (1 + self.aR_step * self.tap_positions[1])**-1
            else:
                self.aR[1] = (1 - self.aR_step * self.tap_positions[1])**-1

            if self.raise_taps[2]:
                self.aR[2] = (1 + self.aR_step * self.tap_positions[2])**-1
            else:
                self.aR[2] = (1 - self.aR_step * self.tap_positions[2])**-1

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
                
                
                
    def calculate_taps(self, Vload, change_type):
        if change_type == 'raise':
            Vlim = self.Vmin
            tap = np.around(np.abs(Vlim - Vload) / self.Vreg_step, 0)
        elif change_type == 'lower':
            Vlim = self.Vmax
            tap = np.around(np.abs(Vlim - Vload) / self.Vreg_step, 0)
        else:
            tap = np.around(
                np.abs(self.band_center - Vload) / self.Vreg_step, 0)

        return tap 
    
    
    def check_regulator(self, Vabc):

        self.Va_mag = np.abs(Vabc[0])
        self.Vb_mag = np.abs(Vabc[1])
        self.Vc_mag = np.abs(Vabc[2])

        # compute voltage limits
        self.band_center + self.bandwidth / 2
        self.band_center - self.bandwidth / 2

        # Initialize raise and lower taps
        raise_taps = []
        lower_taps = []

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
    
    def normalize_voltages(self, Vmag):
        return Vmag / self.pt_ratio

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
    
    def assign_nodes(self, node_index_, node):
        # There are likely to be 12 additional nodes for the XFMRs
        # 6 additional rows for the voltage source eqns. (No Phase Shifters
        # 12 additional rows for the voltage sources eqn (w Phase Shifters
        if self.phases & 0x1 == 1:  # Check for phase A
            self.nodeA_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vr
            self.nodeA_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vi
            self.nodeA_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vr
            self.nodeA_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vi
            
            if self.stamp_dual:
                self.nodeA_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Lr
                self.nodeA_Li_from = node[Nodes.nodeKey[self.from_node]].nodeA_Li
                self.nodeA_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Lr
                self.nodeA_Li_to = node[Nodes.nodeKey[self.to_node]].nodeA_Li
 
            if not self.fixed_taps:
                self.nodeA_Tap = node_index_.__next__()
                self.nodeA_aR = node_index_.__next__()
                self.nodes_aR["A"] = self.nodeA_aR
                self.nodes_Tap["A"] = self.nodeA_Tap
                self.nodeA_Vr_loadCenter = node_index_.__next__()
                self.nodeA_Vi_loadCenter = node_index_.__next__()
                self.nodeA_Vmag2_loadCenter = node_index_.__next__()
                self.nodeA_V_diff = node_index_.__next__()
                self.nodes_V_loadCenter["A"] = [
                    self.nodeA_Vr_loadCenter, self.nodeA_Vi_loadCenter,
                    self.nodeA_Vmag2_loadCenter
                ]
                self.nodes_V_source["A"] = [
                    self.nodeA_Vr_from, self.nodeA_Vi_from
                ]
                self.nodes_V_to["A"] = [self.nodeA_Vr_to, self.nodeA_Vi_to]
                self.nodes_V_diff["A"] = self.nodeA_V_diff
                Nodes.Tap_index.extend([
                    self.nodeA_Tap, self.nodeA_aR, self.nodeA_Vr_loadCenter,
                    self.nodeA_Vi_loadCenter, self.nodeA_V_diff
                ])
                self.node_set.extend([self.nodeA_Tap, self.nodeA_aR])
                
                # TO DO: @emfoster@andrew.cmu.edu
                # Add stamp dual features for changing taps?
        
        if self.phases & 0x2 == 2:  # Check for phase B
            self.nodeB_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vr
            self.nodeB_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vi
            self.nodeB_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vr
            self.nodeB_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vi
            
            if self.stamp_dual:
                self.nodeB_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Lr
                self.nodeB_Li_from = node[Nodes.nodeKey[self.from_node]].nodeB_Li
                self.nodeB_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Lr
                self.nodeB_Li_to = node[Nodes.nodeKey[self.to_node]].nodeB_Li            
            
            if not self.fixed_taps:
                self.nodeB_aR = node_index_.__next__()
                self.nodes_aR["B"] = self.nodeB_aR
                self.nodeB_Tap = node_index_.__next__()
                self.nodes_Tap["B"] = self.nodeB_Tap
                self.nodeB_Vr_loadCenter = node_index_.__next__()
                self.nodeB_Vi_loadCenter = node_index_.__next__()
                self.nodeB_Vmag2_loadCenter = node_index_.__next__()
                self.nodeB_V_diff = node_index_.__next__()
                self.nodes_V_loadCenter["B"] = [
                    self.nodeB_Vr_loadCenter, self.nodeB_Vi_loadCenter,
                    self.nodeB_Vmag2_loadCenter
                ]
                self.nodes_V_source["B"] = [
                    self.nodeB_Vr_from, self.nodeB_Vi_from
                ]
                self.nodes_V_to["B"] = [self.nodeB_Vr_to, self.nodeB_Vi_to]
                self.nodes_V_diff["B"] = self.nodeB_V_diff
                Nodes.Tap_index.extend([
                    self.nodeB_Tap, self.nodeB_aR, self.nodeB_Vr_loadCenter,
                    self.nodeB_Vi_loadCenter, self.nodeB_V_diff
                ])
                self.node_set.extend([self.nodeB_Tap, self.nodeB_aR])
                
                # TO DO: @emfoster@andrew.cmu.edu
                # Add stamp dual features for changing taps?
        
        if self.phases & 0x4 == 4:  # Check for phase C
            self.nodeC_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vr
            self.nodeC_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vi
            self.nodeC_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vr
            self.nodeC_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vi
            
            if self.stamp_dual:
                self.nodeC_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Lr
                self.nodeC_Li_from = node[Nodes.nodeKey[self.from_node]].nodeC_Li
                self.nodeC_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Lr
                self.nodeC_Li_to = node[Nodes.nodeKey[self.to_node]].nodeC_Li
                
            if not self.fixed_taps:
                self.nodeC_aR = node_index_.__next__()
                self.nodes_aR["C"] = self.nodeC_aR
                self.nodeC_Tap = node_index_.__next__()
                self.nodes_Tap["C"] = self.nodeC_Tap
                self.nodeC_Vr_loadCenter = node_index_.__next__()
                self.nodeC_Vi_loadCenter = node_index_.__next__()
                self.nodeC_Vmag2_loadCenter = node_index_.__next__()
                self.nodeC_V_diff = node_index_.__next__()
                self.nodes_V_loadCenter["C"] = [
                    self.nodeC_Vr_loadCenter, self.nodeC_Vi_loadCenter,
                    self.nodeC_Vmag2_loadCenter
                ]
                self.nodes_V_source["C"] = [
                    self.nodeC_Vr_from, self.nodeC_Vi_from
                ]
                self.nodes_V_to["C"] = [self.nodeC_Vr_to, self.nodeC_Vi_to]
                self.nodes_V_diff["C"] = self.nodeC_V_diff
                Nodes.Tap_index.extend([
                    self.nodeC_Tap, self.nodeC_aR, self.nodeC_Vr_loadCenter,
                    self.nodeC_Vi_loadCenter, self.nodeC_V_diff
                ])
                self.node_set.extend([self.nodeC_Tap, self.nodeC_aR])
                
                 # TO DO: @emfoster@andrew.cmu.edu
                # Add stamp dual features for changing taps?

        self.nodeN_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vr
        self.nodeN_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeN_Vi
        self.nodeN_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vr
        self.nodeN_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeN_Vi
        
        if self.stamp_dual:    
            self.nodeN_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeN_Lr
            self.nodeN_Li_from = node[Nodes.nodeKey[self.from_node]].nodeN_Li
            self.nodeN_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeN_Lr
            self.nodeN_Li_to = node[Nodes.nodeKey[self.to_node]].nodeN_Li           
        #############################################################
        # Add ground voltage source if needed
        if self.connect_type == self._GWYE_GWYE or self.connect_type == self._DELTA_DELTA:  # Wye-Wye
            # Add for the three phase extra node in the secondary
            # Add for the three phases extra node for the voltage sources
            # circuit (for both real and imaginary)
            if self.phases & 0x1 == 1:  # Check for phase A
                self.nodeA_Vr_secondary = node_index_.__next__()
                self.nodeA_Vi_secondary = node_index_.__next__()
                self.nodeA_Vr_primary = node_index_.__next__()
                self.nodeA_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["A"] = self.nodeA_Vr_primary
                self.node_Vi_meas["A"] = self.nodeA_Vi_primary

                
                if self.stamp_dual:
                    self.nodeA_Lr_secondary = node_index_.__next__()
                    self.nodeA_Li_secondary = node_index_.__next__()
                    self.nodeA_Lr_primary = node_index_.__next__()
                    self.nodeA_Li_primary = node_index_.__next__()
                    self.node_Lr_meas["A"] = self.nodeA_Lr_primary
                    self.node_Li_meas["A"] = self.nodeA_Li_primary
                    
                    Nodes.I_index.extend(
                        [self.nodeA_Vr_primary, self.nodeA_Vi_primary, 
                         self.nodeA_Lr_primary, self.nodeA_Li_primary])
                    self.node_set.extend([
                        self.nodeA_Vr_secondary, self.nodeA_Vi_secondary,
                        self.nodeA_Vr_primary, self.nodeA_Vi_primary,
                        self.nodeA_Lr_secondary, self.nodeA_Li_secondary,
                        self.nodeA_Lr_primary, self.nodeA_Li_primary
                        ])
                    
                else:
                    Nodes.I_index.extend(
                        [self.nodeA_Vr_primary, self.nodeA_Vi_primary])
                    self.node_set.extend([
                        self.nodeA_Vr_secondary, self.nodeA_Vi_secondary,
                        self.nodeA_Vr_primary, self.nodeA_Vi_primary
                    ])
                    
            if self.phases & 0x2 == 2:  # Check for phase B
                self.nodeB_Vr_secondary = node_index_.__next__()
                self.nodeB_Vi_secondary = node_index_.__next__()
                self.nodeB_Vr_primary = node_index_.__next__()
                self.nodeB_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["B"] = self.nodeB_Vr_primary
                self.node_Vi_meas["B"] = self.nodeB_Vi_primary
                
                if self.stamp_dual:
                    self.nodeB_Lr_secondary = node_index_.__next__()
                    self.nodeB_Li_secondary = node_index_.__next__()
                    self.nodeB_Lr_primary = node_index_.__next__()
                    self.nodeB_Li_primary = node_index_.__next__()
                    self.node_Lr_meas["B"] = self.nodeB_Lr_primary
                    self.node_Li_meas["B"] = self.nodeB_Li_primary
                    
                    Nodes.I_index.extend(
                        [self.nodeB_Vr_primary, self.nodeB_Vi_primary, 
                         self.nodeB_Lr_primary, self.nodeB_Li_primary])
                    self.node_set.extend([
                        self.nodeB_Vr_secondary, self.nodeB_Vi_secondary,
                        self.nodeB_Vr_primary, self.nodeB_Vi_primary,
                        self.nodeB_Lr_secondary, self.nodeB_Li_secondary,
                        self.nodeB_Lr_primary, self.nodeB_Li_primary
                        ])
                    
                else:
                    Nodes.I_index.extend(
                        [self.nodeB_Vr_primary, self.nodeB_Vi_primary])
                    self.node_set.extend([
                        self.nodeB_Vr_secondary, self.nodeB_Vi_secondary,
                        self.nodeB_Vr_primary, self.nodeB_Vi_primary
                    ])
                    
            if self.phases & 0x4 == 4:  # Check for phase C
                self.nodeC_Vr_secondary = node_index_.__next__()
                self.nodeC_Vi_secondary = node_index_.__next__()
                self.nodeC_Vr_primary = node_index_.__next__()
                self.nodeC_Vi_primary = node_index_.__next__()
                self.node_Vr_meas["C"] = self.nodeC_Vr_primary
                self.node_Vi_meas["C"] = self.nodeC_Vi_primary
                
                if self.stamp_dual:
                    self.nodeC_Lr_secondary = node_index_.__next__()
                    self.nodeC_Li_secondary = node_index_.__next__()
                    self.nodeC_Lr_primary = node_index_.__next__()
                    self.nodeC_Li_primary = node_index_.__next__()
                    self.node_Lr_meas["C"] = self.nodeC_Lr_primary
                    self.node_Li_meas["C"] = self.nodeC_Li_primary
                    
                    Nodes.I_index.extend(
                        [self.nodeC_Vr_primary, self.nodeC_Vi_primary, 
                         self.nodeC_Lr_primary, self.nodeC_Li_primary])
                    self.node_set.extend([
                        self.nodeC_Vr_secondary, self.nodeC_Vi_secondary,
                        self.nodeC_Vr_primary, self.nodeC_Vi_primary,
                        self.nodeC_Lr_secondary, self.nodeC_Li_secondary,
                        self.nodeC_Lr_primary, self.nodeC_Li_primary
                        ])
                    
                else:
                    Nodes.I_index.extend(
                        [self.nodeC_Vr_primary, self.nodeC_Vi_primary])
                    self.node_set.extend([
                        self.nodeC_Vr_secondary, self.nodeC_Vi_secondary,
                        self.nodeC_Vr_primary, self.nodeC_Vi_primary
                    ])
            
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
                    self.nodeVGnd_rp = self._regulator_locations[(
                        self.from_node, self.to_node)][0]
                    self.nodeVGnd_ip = self._regulator_locations[(
                        self.from_node, self.to_node)][1]
                    self.nodeVGnd_rs = self._regulator_locations[(
                        self.from_node, self.to_node)][2]
                    self.nodeVGnd_is = self._regulator_locations[(
                        self.from_node, self.to_node)][3]
                    
                    if self.stamp_dual:
                        self.nodeGnd_Lr_p = self._regulator_locations[(
                            self.from_node, self.to_node)][0]
                        self.nodeGnd_Li_p = self._regulator_locations[(
                            self.from_node, self.to_node)][1]
                        self.nodeGnd_Lr_s = self._regulator_locations[(
                            self.from_node, self.to_node)][2]
                        self.nodeGnd_Li_s = self._regulator_locations[(
                            self.from_node, self.to_node)][3]
                
        else:
            print('Incorrect regulator type for assigning nodes')
        return None
    
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
        nodes_I_comp = []
        nodes_Vr_load = []
        nodes_Vi_load = []
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
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeA_Lr_from)
                nodes_Li_from.append(self.nodeA_Li_from)
                nodes_Lr_to.append(self.nodeA_Lr_to)
                nodes_Li_to.append(self.nodeA_Li_to)
                nodes_Lr_secondary.append(self.nodeA_Lr_secondary)
                nodes_Li_secondary.append(self.nodeA_Li_secondary)
                nodes_Lr_primary.append(self.nodeA_Lr_primary)
                nodes_Li_primary.append(self.nodeA_Li_primary)
            
            if not self.fixed_taps:
                pass
                # nodes_taps.append(self.nodeA_Tap)
                # nodes_V_loadCenters.append(self.nodeA_Vi_compRelay)
                # nodes_Vr_relays.append(self.nodeA_Vr_compRelay)
                # nodes_Vi_relays.append(self.nodeA_Vi_compRelay)
                # nodes_Vr_load.append(self.nodeA_Vr_compReg)
                # nodes_Vi_load.append(self.nodeA_Vi_compReg)
                # nodes_Vr_compDrop.append(self.nodeA_Vr_compDrop)
                # nodes_Vi_compDrop.append(self.nodeA_Vi_compDrop)
                # nodes_I_comp.append(self.nodeA_I_comp)
                # nodes_Ir_comp.append(self.nodeA_Ir_comp)
                # nodes_Ii_comp.append(self.nodeA_Ii_comp)
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
            nodes_I_comp.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)
            
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
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeB_Lr_from)
                nodes_Li_from.append(self.nodeB_Li_from)
                nodes_Lr_to.append(self.nodeB_Lr_to)
                nodes_Li_to.append(self.nodeB_Li_to)
                nodes_Lr_secondary.append(self.nodeB_Lr_secondary)
                nodes_Li_secondary.append(self.nodeB_Li_secondary)
                nodes_Lr_primary.append(self.nodeB_Lr_primary)
                nodes_Li_primary.append(self.nodeB_Li_primary)
                
            if not self.fixed_taps:
                pass
                # nodes_taps.append(self.nodeB_Tap)
                # nodes_V_loadCenters.append(self.nodeB_V_compRelay)
                # nodes_Vr_relays.append(self.nodeB_Vr_compRelay)
                # nodes_Vi_relays.append(self.nodeB_Vi_compRelay)
                # nodes_Vr_load.append(self.nodeB_Vr_compReg)
                # nodes_Vi_load.append(self.nodeB_Vi_compReg)
                # nodes_Vr_compDrop.append(self.nodeB_Vr_compDrop)
                # nodes_Vi_compDrop.append(self.nodeB_Vi_compDrop)
                # nodes_I_comp.append(self.nodeB_I_comp)
                # nodes_Ir_comp.append(self.nodeB_Ir_comp)
                # nodes_Ii_comp.append(self.nodeB_Ii_comp)
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
            nodes_I_comp.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)
            
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
            
            if self.stamp_dual:
                nodes_Lr_from.append(self.nodeC_Lr_from)
                nodes_Li_from.append(self.nodeC_Li_from)
                nodes_Lr_to.append(self.nodeC_Lr_to)
                nodes_Li_to.append(self.nodeC_Li_to)
                nodes_Lr_secondary.append(self.nodeC_Lr_secondary)
                nodes_Li_secondary.append(self.nodeC_Li_secondary)
                nodes_Lr_primary.append(self.nodeC_Lr_primary)
                nodes_Li_primary.append(self.nodeC_Li_primary)
                
            if not self.fixed_taps:
                pass
                # nodes_taps.append(self.nodeC_Tap)
                # # nodes_V_loadCenters.append(self.nodeC_V_compRelay)
                # nodes_Vr_relays.append(self.nodeC_Vr_compRelay)
                # nodes_Vi_relays.append(self.nodeC_Vi_compRelay)
                # nodes_Vr_load.append(self.nodeC_Vr_compReg)
                # nodes_Vi_load.append(self.nodeC_Vi_compReg)
                # nodes_Vr_compDrop.append(self.nodeC_Vr_compDrop)
                # nodes_Vi_compDrop.append(self.nodeC_Vi_compDrop)
                # nodes_I_comp.append(self.nodeC_I_comp)
                # nodes_Ir_comp.append(self.nodeC_Ir_comp)
                # nodes_Ii_comp.append(self.nodeC_Ii_comp)
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
            nodes_I_comp.append(-1)
            nodes_Ir_comp.append(-1)
            nodes_Ii_comp.append(-1)
            
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
                      'nodes_Lr_secondary': nodes_Lr_secondary,
                      'nodes_Li_secondary': nodes_Li_secondary,
                      'nodes_Lr_primary': nodes_Lr_primary,
                      'nodes_Li_primary': nodes_Li_primary}
        else:
            L_dict = None
        return nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
               nodes_Vr_secondary, nodes_Vi_secondary, nodes_taps, nodes_Vr_relays, \
               nodes_Vi_relays, nodes_V_loadCenters, nodes_Vr_load, nodes_Vi_load, nodes_Ir_comp, \
               nodes_Ii_comp, nodes_I_comp, Ir_sec, Ii_sec, L_dict
               
    def get_nodes(self, node, V):
        """
		Find the indices of the solution vector and the values at those indices.
		:param node: The vector of all system nodes.
		:param V: The solution vector.
		:return: Three dictionaries which hold the phase and sequence nodes.
		"""
        # Find the from bus nodes
        nodeA_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vr
        nodeA_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeA_Vi
        nodeB_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vr
        nodeB_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeB_Vi
        nodeC_Vr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vr
        nodeC_Vi_from = node[Nodes.nodeKey[self.from_node]].nodeC_Vi

        # Find the to bus nodes
        nodeA_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vr
        nodeA_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeA_Vi
        nodeB_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vr
        nodeB_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeB_Vi
        nodeC_Vr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vr
        nodeC_Vi_to = node[Nodes.nodeKey[self.to_node]].nodeC_Vi

        # Find the from and to node voltages
        (Var_from, Vai_from, Vbr_from, Vbi_from, Vcr_from, Vci_from) = \
            (V[nodeA_Vr_from], V[nodeA_Vi_from], V[nodeB_Vr_from], V[nodeB_Vi_from], V[nodeC_Vr_from],
             V[nodeC_Vi_from])

        (Var_to, Vai_to, Vbr_to, Vbi_to, Vcr_to, Vci_to) = \
            (V[nodeA_Vr_to], V[nodeA_Vi_to], V[nodeB_Vr_to], V[nodeB_Vi_to], V[nodeC_Vr_to],
             V[nodeC_Vi_to])
        
        if self.stamp_dual:
            nodeA_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeA_Lr
            nodeA_Li_from = node[Nodes.nodeKey[self.from_node]].nodeA_Li
            nodeB_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeB_Lr
            nodeB_Li_from = node[Nodes.nodeKey[self.from_node]].nodeB_Li
            nodeC_Lr_from = node[Nodes.nodeKey[self.from_node]].nodeC_Lr
            nodeC_Li_from = node[Nodes.nodeKey[self.from_node]].nodeC_Li
    
            # Find the to bus nodes
            nodeA_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeA_Lr
            nodeA_Li_to = node[Nodes.nodeKey[self.to_node]].nodeA_Li
            nodeB_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeB_Lr
            nodeB_Li_to = node[Nodes.nodeKey[self.to_node]].nodeB_Li
            nodeC_Lr_to = node[Nodes.nodeKey[self.to_node]].nodeC_Lr
            nodeC_Li_to = node[Nodes.nodeKey[self.to_node]].nodeC_Li

        # Find the relay voltages and Taps
        TapA = self.tapA
        TapB = self.tapB
        TapC = self.tapC
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
        Iar_comp = 0.0
        Iai_comp = 0.0
        Ibr_comp = 0.0
        Ibi_comp = 0.0
        Icr_comp = 0.0
        Ici_comp = 0.0
        Va_pos = 0.0
        Vb_pos = 0.0
        Vc_pos = 0.0
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
        Va_pos_to = 0.0
        Vb_pos_to = 0.0
        Vc_pos_to = 0.0

        if not self.fixed_taps:
            if self.phases & 0x1 == 1:
                TapA = np.around(V[self.nodeA_Tap][0])
                Va_source = complex(V[self.nodeA_Vr_from],
                                    V[self.nodeA_Vi_from])
                Va_loadCenter = complex(V[self.nodeA_Vr_loadCenter],
                                        V[self.nodeA_Vi_loadCenter])
                Va_mag2_loadCenter = V[self.nodeA_Vmag2_loadCenter]
                # Va_load = complex(V[self.nodeA_Vr_load], V[self.nodeA_Vi_load])
                Ia_prim = complex(V[self.nodeA_Vr_primary],
                                  V[self.nodeA_Vi_primary])
                Va_pos = complex(V[self.nodeA_Vr_secondary],
                                 V[self.nodeA_Vi_secondary])
                Va_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
                Va_diff = V[self.nodeA_V_diff]
                Va_neg_from = complex(V[self.nodeN_Vr_from],
                                      V[self.nodeN_Vi_from])
                aR_A = V[self.nodeA_aR][0]
                Va_pos_to = complex(V[self.nodeA_Vr_secondary],
                                    V[self.nodeA_Vi_secondary])
                # Iar_comp = V[self.nodeA_Ir_comp][0]
                # Iai_comp = V[self.nodeA_Ii_comp][0]
                np.abs(complex(Iar_comp, Iai_comp))
            if self.phases & 0x2 == 2:
                TapB = np.around(V[self.nodeB_Tap][0])
                Vb_source = complex(V[self.nodeB_Vr_from],
                                    V[self.nodeB_Vi_from])
                Vb_loadCenter = complex(V[self.nodeB_Vr_loadCenter],
                                        V[self.nodeB_Vi_loadCenter])
                Vb_mag2_loadCenter = V[self.nodeB_Vmag2_loadCenter]
                # Vb_load = complex(V[self.nodeB_Vr_load], V[self.nodeB_Vi_load])
                Ib_prim = complex(V[self.nodeB_Vr_primary],
                                  V[self.nodeB_Vi_primary])
                Vb_pos = complex(V[self.nodeB_Vr_secondary],
                                 V[self.nodeB_Vi_secondary])
                Vb_neg_from = complex(V[self.nodeN_Vr_from],
                                      V[self.nodeN_Vi_from])
                Vb_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
                aR_B = V[self.nodeB_aR][0]
                Vb_diff = V[self.nodeB_V_diff]
                Vb_pos_to = complex(V[self.nodeB_Vr_secondary],
                                    V[self.nodeB_Vi_secondary])
                # Ibr_comp = V[self.nodeB_Ir_comp][0]
                # Ibi_comp = V[self.nodeB_Ii_comp][0]
                np.abs(complex(Ibr_comp, Ibi_comp))

            if self.phases & 0x4 == 4:
                TapC = np.around(V[self.nodeC_Tap][0])
                Vc_source = complex(V[self.nodeC_Vr_from],
                                    V[self.nodeC_Vi_from])
                Vc_loadCenter = complex(V[self.nodeC_Vr_loadCenter],
                                        V[self.nodeC_Vi_loadCenter])
                Vc_mag2_loadCenter = V[self.nodeC_Vmag2_loadCenter]
                # Vc_load = complex(V[self.nodeC_Vr_load], V[self.nodeC_Vi_load])
                Ic_prim = complex(V[self.nodeC_Vr_primary],
                                  V[self.nodeC_Vi_primary])
                Vc_pos = complex(V[self.nodeC_Vr_secondary],
                                 V[self.nodeC_Vi_secondary])
                Vc_neg_from = complex(V[self.nodeN_Vr_from],
                                      V[self.nodeN_Vi_from])
                Vc_neg_to = complex(V[self.nodeN_Vr_to], V[self.nodeN_Vi_to])
                aR_C = V[self.nodeC_aR][0]
                Vc_diff = V[self.nodeC_V_diff]
                Vc_pos_to = complex(V[self.nodeC_Vr_secondary],
                                    V[self.nodeC_Vi_secondary])
                # Icr_comp = V[self.nodeC_Ir_comp][0]
                # Ici_comp = V[self.nodeC_Ii_comp][0]
                np.abs(complex(Icr_comp, Ici_comp))

        self.aR = [aR_A, aR_B, aR_C]
        self.tap_positions = [TapA, TapB, TapC]

        if self.stamp_dual:
            Vabc = {
                'Vr_from': [Var_from, Vbr_from, Vcr_from],
                'Vi_from': [Vai_from, Vbi_from, Vci_from],
                'node_Vr_from': [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from],
                'node_Vi_from': [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from],
                'To': [
                    complex(Var_to, Vai_to),
                    complex(Vbr_to, Vbi_to),
                    complex(Vcr_to, Vci_to)
                ],
                'node_Vr_to': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to],
                'node_Vi_to': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to],
                'Tap': [TapA, TapB, TapC],
                'Vsource': [Va_source, Vb_source, Vc_source],
                'Vlc': [Va_loadCenter, Vb_loadCenter, Vc_loadCenter],
                'Vlc2': [
                    Va_mag2_loadCenter, Vb_mag2_loadCenter, Vc_mag2_loadCenter
                ],
                'Iprim': [Ia_prim, Ib_prim, Ic_prim],
                'Vpos': [Va_pos, Vb_pos, Vc_pos],
                'Vneg_from': [Va_neg_from, Vb_neg_from, Vc_neg_from],
                'Vneg_to': [Va_neg_to, Vb_neg_to, Vc_neg_to],
                'aR': [aR_A, aR_B, aR_C],
                'Imag_prim': [Ia_mag_prim, Ib_mag_prim, Ic_mag_prim],
                'Vpos_to': [Va_pos_to, Vb_pos_to, Vc_pos_to],
                'Vdiff': [Va_diff, Vb_diff, Vc_diff],
                'Lr_from': [V[nodeA_Lr_from], V[nodeB_Lr_from], V[nodeC_Lr_from]],
                'Li_from': [V[nodeA_Li_from], V[nodeB_Li_from], V[nodeC_Li_from]],
                'node_Lr_from': [nodeA_Lr_from, nodeB_Lr_from, nodeC_Lr_from],
                'node_Li_from': [nodeA_Li_from, nodeB_Li_from, nodeC_Li_from],
                'Lr_to': [V[nodeA_Lr_to], V[nodeB_Lr_to], V[nodeC_Lr_to]],
                'Li_to': [V[nodeA_Li_to], V[nodeB_Li_to], V[nodeC_Li_to]],
                'node_Lr_to': [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to],
                'node_Li_to': [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to]
            }
        else:
            Vabc = {
                'Vr_from': [Var_from, Vbr_from, Vcr_from],
                'Vi_from': [Vai_from, Vbi_from, Vci_from],
                'node_Vr_from': [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from],
                'node_Vi_from': [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from],
                'To': [
                    complex(Var_to, Vai_to),
                    complex(Vbr_to, Vbi_to),
                    complex(Vcr_to, Vci_to)
                ],
                'node_Vr_to': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to],
                'node_Vi_to': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to],
                'Tap': [TapA, TapB, TapC],
                'Vsource': [Va_source, Vb_source, Vc_source],
                'Vlc': [Va_loadCenter, Vb_loadCenter, Vc_loadCenter],
                'Vlc2': [
                    Va_mag2_loadCenter, Vb_mag2_loadCenter, Vc_mag2_loadCenter
                ],
                'Iprim': [Ia_prim, Ib_prim, Ic_prim],
                'Vpos': [Va_pos, Vb_pos, Vc_pos],
                'Vneg_from': [Va_neg_from, Vb_neg_from, Vc_neg_from],
                'Vneg_to': [Va_neg_to, Vb_neg_to, Vc_neg_to],
                'aR': [aR_A, aR_B, aR_C],
                'Imag_prim': [Ia_mag_prim, Ib_mag_prim, Ic_mag_prim],
                'Vpos_to': [Va_pos_to, Vb_pos_to, Vc_pos_to],
                'Vdiff': [Va_diff, Vb_diff, Vc_diff]
            }

        Vabc = SimpleNamespace(**Vabc)

        return Vabc 
    
    def update_values(self, V):

        # # Update Taps # #
        self.tapA = np.around(V[self.nodeA_Tap][0], 0) if hasattr(
            self, 'nodeA_Tap') else self.tapA
        self.tapB = np.around(V[self.nodeB_Tap][0], 0) if hasattr(
            self, 'nodeB_Tap') else self.tapB
        self.tapC = np.around(V[self.nodeC_Tap][0], 0) if hasattr(
            self, 'nodeC_Tap') else self.tapC

        # # Update Effective Regulator Ratios # #
        self.aR_A = np.around(V[self.nodeA_aR][0], 5) if hasattr(
            self, 'nodeA_aR') else np.around(self.aR_A, 5)
        self.aR_B = np.around(V[self.nodeB_aR][0], 5) if hasattr(
            self, 'nodeB_aR') else np.around(self.aR_B, 5)
        self.aR_C = np.around(V[self.nodeC_aR][0], 5) if hasattr(
            self, 'nodeC_aR') else np.around(self.aR_C, 5)

        # # Update Regulator Currents # #
        deg = True
        Ia = complex(V[self.nodeA_Vr_primary],
                     V[self.nodeA_Vi_primary]) if hasattr(
                         self, 'nodeA_Vr_primary') else 0
        self.Ia_mag = np.around(np.abs(Ia), 2)
        self.Ia_ang = np.around(np.angle(Ia, deg=deg), 2)
        Ib = complex(V[self.nodeB_Vr_primary],
                     V[self.nodeB_Vi_primary]) if hasattr(
                         self, 'nodeB_Vr_primary') else 0
        self.Ib_mag = np.around(np.abs(Ib), 2)
        self.Ib_ang = np.around(np.angle(Ib, deg=deg), 2)
        Ic = complex(V[self.nodeC_Vr_primary],
                     V[self.nodeC_Vi_primary]) if hasattr(
                         self, 'nodeC_Vr_primary') else 0
        self.Ic_mag = np.around(np.abs(Ic), 2)
        self.Ic_ang = np.around(np.angle(Ic, deg=deg), 2)

        print('Regulator %s has taps: %d, %d, %d' %
              (self.name, self.tapA, self.tapB, self.tapC))
        
    def stamp_primary(self, Ylin_val, Ylin_row, Ylin_col, idx_Y, 
                     node_Vr_pos_from, node_Vi_pos_from, node_Vr_neg_from, 
                     node_Vi_neg_from,  node_Vr_neg_secondary, 
                     node_Vi_neg_secondary, node_Vr_pos_secondary,
                     node_Vi_pos_secondary, node_Vr_primary, node_Vi_primary,
                     phase):

        # # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]      
        # PRIMARY VOLTAGES
        # Real #
        # Vr_primary = Vr_pos_from - tr * Vr_secondary_pos  + tr * Vr_secondary_neg - Vr_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_primary, node_Vr_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_primary, node_Vr_pos_secondary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_primary, node_Vr_neg_secondary, aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_primary, node_Vr_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        # Imaginary #
        # Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_primary, node_Vi_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_primary, node_Vi_pos_secondary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_primary, node_Vi_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_primary, node_Vi_neg_secondary, aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)

        # PRIMARY CURRENTS
        # Real #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_pos_from, node_Vr_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vr_neg_from, node_Vr_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        # Imaginary #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_pos_from, node_Vi_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Vi_neg_from, node_Vi_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        if not self.fixed_taps:
            pass

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    def stamp_primary_dual(self, Ylin_val, Ylin_row, Ylin_col, idx_Y, 
                     node_Lr_pos_from, node_Li_pos_from, node_Lr_neg_from, 
                     node_Li_neg_from,  node_Lr_neg_secondary, 
                     node_Li_neg_secondary, node_Lr_pos_secondary,
                     node_Li_pos_secondary, node_Lr_primary, node_Li_primary,
                     phase):

        # # Get Effective Regulator Ratio of the Corresponding Phase # #
        aR = self.aR[phase]      
        # DUAL OF PRIMARY VOLTAGES - TRANSPOSE OF LINEAR STAMPS
        # Real #
        # Vr_primary = Vr_pos_from - tr * Vr_secondary_pos  + tr * Vr_secondary_neg - Vr_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_pos_from, node_Lr_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_pos_secondary, node_Lr_primary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_neg_secondary, node_Lr_primary, aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_neg_from, node_Lr_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        # Imaginary #
        # Vi_primary = Vi_neg_from - tr * Vi_secondary_pos  + tr * Vi_secondary_neg - Vi_neg_from
        # neg (negative) could be the neutral node or it could be another phase node if delta connected
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_pos_from, node_Li_primary, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_pos_secondary, node_Li_primary, -aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_neg_from, node_Li_primary, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_neg_secondary, node_Li_primary,  aR, Ylin_val, Ylin_row,
            Ylin_col, idx_Y)

        # PRIMARY CURRENTS
        # Real #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_primary, node_Lr_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Lr_primary, node_Lr_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        # Imaginary #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_primary, node_Li_pos_from, 1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            node_Li_primary, node_Li_neg_from, -1, Ylin_val, Ylin_row, Ylin_col,
            idx_Y)

        if not self.fixed_taps:
            pass

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_secondary(self, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                       node_Vr_pos_to, node_Vi_pos_to,
                       node_Vr_neg_secondary, node_Vi_neg_secondary,
                       node_Vr_pos_secondary, node_Vi_pos_secondary,
                       node_Vr_primary, node_Vi_primary, phase, 
                       homotopy = None, B_factor = None, G_factor = None): 
         # Calculate conductance and susceptance of regulator if homotopy
        if not homotopy:
            # # Get Effective Regulator Ratio of the Corresponding Phase # #
            aR = self.aR[phase]
            _B = self.B
            _G = self.G
        else:
            aR = 1.0
            _B = B_factor * self.B
            _G = G_factor * self.G

        if aR:
            # # SECONDARY CURRENT SOURCES # #
            # Real #
            # Ir_pos_secondary = - aR * Ir_primary
            # Ir_neg_secondary = aR * Ir_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_secondary, node_Vr_primary, -aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_neg_secondary, node_Vr_primary, aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

            # Imaginary #
            # Ii_pos_secondary = - aR * Ii_primary
            # Ii_neg_secondary = aR * Ii_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_secondary, node_Vi_primary, -aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_neg_secondary, node_Vi_primary, aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        # Branch stamping for the susceptances
        if _B:
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_to, node_Vi_pos_to, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_to, node_Vi_pos_secondary, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_to, node_Vr_pos_to, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_to, node_Vr_pos_secondary, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_secondary, node_Vi_pos_secondary, -_B, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_secondary, node_Vi_pos_to, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_secondary, node_Vr_pos_secondary, _B, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_secondary, node_Vr_pos_to, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        # Adding resistance if existent
        if _G:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_to, node_Vr_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_to, node_Vr_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_to, node_Vi_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_to, node_Vi_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_secondary, node_Vr_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vr_pos_secondary, node_Vr_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_secondary, node_Vi_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Vi_pos_secondary, node_Vi_pos_to, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_secondary_dual(self, Ylin_val, Ylin_row, Ylin_col, idx_Y,
                       node_Lr_pos_to, node_Li_pos_to,
                       node_Lr_neg_secondary, node_Li_neg_secondary,
                       node_Lr_pos_secondary, node_Li_pos_secondary,
                       node_Lr_primary, node_Li_primary, phase, 
                       homotopy = None, B_factor = None, G_factor = None): 
         # Calculate conductance and susceptance of regulator if homotopy
        if not homotopy:
            # # Get Effective Regulator Ratio of the Corresponding Phase # #
            aR = self.aR[phase]
            _B = self.B
            _G = self.G
        else:
            aR = 1.0
            _B = B_factor * self.B
            _G = G_factor * self.G

        if aR:
            # # SECONDARY CURRENT SOURCES # #
            # Real #
            # Ir_pos_secondary = - aR * Ir_primary
            # Ir_neg_secondary = aR * Ir_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_primary, node_Lr_pos_secondary, -aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
               node_Lr_primary,  node_Lr_neg_secondary, aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

            # Imaginary #
            # Ii_pos_secondary = - aR * Ii_primary
            # Ii_neg_secondary = aR * Ii_primary
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_primary, node_Li_pos_secondary, -aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_primary, node_Li_neg_secondary, aR, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        # Branch stamping for the susceptances
        if _B:
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_to, node_Lr_pos_to, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_secondary, node_Lr_pos_to, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_to, node_Li_pos_to, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_secondary, node_Li_pos_to, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_secondary, node_Lr_pos_secondary, -_B, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_to,node_Lr_pos_secondary, _B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_secondary, node_Li_pos_secondary, _B, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_to, node_Li_pos_secondary, -_B, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        # Adding resistance if existent
        if _G:
            # Branch stamping for the conductance
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_to, node_Lr_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_secondary,node_Lr_pos_to,  -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_to, node_Li_pos_to, _G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_secondary, node_Li_pos_to,  -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_secondary, node_Lr_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Lr_pos_to, node_Lr_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_secondary, node_Li_pos_secondary, _G, Ylin_val,
                Ylin_row, Ylin_col, idx_Y)
            #
            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
                node_Li_pos_to, node_Li_pos_secondary, -_G, Ylin_val, Ylin_row,
                Ylin_col, idx_Y)

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
    
    @staticmethod
    def stamp_ground_neutral(self, neg_real_node, neg_imag_node, gnd_real_V,
                             gnd_imag_V, Ylin_val, Ylin_row, Ylin_col, idx_Y):
        # Stamp Vng = 0 + j0
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            gnd_real_V, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            gnd_imag_V, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        # Stamp Current Sources
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            neg_real_node, gnd_real_V, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            neg_imag_node, gnd_imag_V, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #

        return Ylin_val, Ylin_row, Ylin_col, idx_Y
   
    @staticmethod
    def stamp_ground_neutral_dual(self, neg_real_node, neg_imag_node, gnd_real_Lr,
                             gnd_imag_Li, Ylin_val, Ylin_row, Ylin_col, idx_Y):
        # Stamp Vng = 0 + j0
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            neg_real_node, gnd_real_Lr, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            neg_imag_node, gnd_imag_Li, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        # Stamp Current Sources
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            gnd_real_Lr, neg_real_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #
        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stampY(
            gnd_imag_Li, neg_imag_node, 1, Ylin_val, Ylin_row, Ylin_col, idx_Y)
        #

        return Ylin_val, Ylin_row, Ylin_col, idx_Y

    def stamp_linear(self, node, V, Ylin_val, Ylin_row, Ylin_col, idx_Y, 
                     stamped_ground, manual_init, homotopy_option = False,
                     h_factor = None, G_homotopy = None, B_homotopy = None):

        self.get_nodes(node, V)
        
        nodes_Vr_from, nodes_Vi_from, nodes_Vr_to, nodes_Vi_to, \
        nodes_Vr_secondary, nodes_Vi_secondary, nodes_taps, nodes_Vr_relay, \
        nodes_Vi_relay, nodes_V_loadCenter, nodes_Vr_load, nodes_Vi_load, nodes_Ir_comp, \
        nodes_Ii_comp, nodes_I_comp, nodes_Vr_primary, nodes_Vi_primary, L_dict = self.append_nodes()

        if self.stamp_dual:
            nodes_Lr_from = L_dict['nodes_Lr_from']
            nodes_Li_from = L_dict['nodes_Li_from']
            nodes_Lr_to = L_dict['nodes_Lr_to']
            nodes_Li_to = L_dict['nodes_Li_to']
            nodes_Lr_secondary = L_dict['nodes_Lr_secondary']
            nodes_Li_secondary = L_dict['nodes_Li_secondary']
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
                    ##########Stamp Primary Circuit##################
                    #from(+)|from(-)|to(+)|to(-)|extraNode|voltageSource
                    _pos = phaseSet[0]
                    _neg = phaseSet[1]
                    _winding = phaseSet[0]
                    
                    if not self.fixed_taps and manual_init:     
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
                                nodes_Lr_secondary[_winding], 
                                nodes_Li_secondary[_winding],
                                nodes_Lr_primary[_winding], 
                                nodes_Li_primary[_winding], phaseSet[0])
                        
                        if not homotopy_option:
                            ##########Stamp Secondary Circuit##################
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodes_Vr_to[_pos], nodes_Vi_to[_pos],
                                nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                                nodes_Vr_secondary[_winding],
                                nodes_Vi_secondary[_winding],
                                nodes_Vr_primary[_winding],
                                nodes_Vi_primary[_winding], phaseSet[0])
        
                            if self.stamp_dual: 
                                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary_dual(
                                    Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                    nodes_Lr_to[_pos], nodes_Li_to[_pos],
                                    nodes_Lr_to[_neg], nodes_Li_to[_neg],
                                    nodes_Lr_secondary[_winding],
                                    nodes_Li_secondary[_winding],
                                    nodes_Lr_primary[_winding],
                                    nodes_Li_primary[_winding], phaseSet[0])
                        else:
                            G_factor = G_homotopy * h_factor
                            B_factor = B_homotopy * h_factor
                            Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodes_Vr_to[_pos], nodes_Vi_to[_pos],
                                nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                                nodes_Vr_secondary[_winding],
                                nodes_Vi_secondary[_winding],
                                nodes_Vr_primary[_winding],
                                nodes_Vi_primary[_winding], phaseSet[0],
                                homotopy_option, B_factor, G_factor)
                            
                            if self.stamp_dual:
                                Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary_dual(
                                    Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                    nodes_Lr_to[_pos], nodes_Li_to[_pos],
                                    nodes_Lr_to[_neg], nodes_Li_to[_neg],
                                    nodes_Lr_secondary[_winding],
                                    nodes_Li_secondary[_winding],
                                    nodes_Lr_primary[_winding],
                                    nodes_Li_primary[_winding], phaseSet[0],
                                    homotopy_option, B_factor, G_factor)
                            
                # #########Stamp Neutral to Ground##################
                _neg = _N
                # One for set of three phase
                # Stamp primary
                if self.from_node not in stamped_ground:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(self,
                        nodes_Vr_from[_neg], nodes_Vi_from[_neg], self.nodeGnd_Vr_primary,
                        self.nodeGnd_Vi_primary, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    stamped_ground.add(self.from_node)
                    
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral_dual(self,
                            nodes_Lr_from[_neg], nodes_Li_from[_neg], self.nodeGnd_Lr_primary,
                            self.nodeGnd_Li_primary, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                        stamped_ground.add(self.from_node)
                # Stamp secondary
                if self.to_node not in stamped_ground:
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral(self,
                        nodes_Vr_to[_neg], nodes_Vi_to[_neg], self.nodeGnd_Vr_secondary,
                        self.nodeGnd_Vi_secondary, Ylin_val, Ylin_row, Ylin_col, idx_Y)
                    stamped_ground.add(self.to_node)
                    
                    if self.stamp_dual:
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_ground_neutral_dual(self,
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
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_primary(
                            Ylin_val, Ylin_row, Ylin_col, idx_Y,
                            nodes_Vr_from[_pos], nodes_Vi_from[_pos], 
                            nodes_Vr_from[_neg], nodes_Vi_from[_neg], 
                            nodes_Vr_to[_neg], nodes_Vi_to[_neg],
                            nodes_Vr_secondary[_winding], 
                            nodes_Vi_secondary[_winding],
                            nodes_Vr_primary[_winding], 
                            nodes_Vi_primary[_winding], phaseSet[0])

                    # #########Stamp Secondary Circuit##################
                    Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary(
                                Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                nodes_Vr_to[_pos], nodes_Vi_to[_pos],
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
                                nodes_Lr_secondary[_winding], 
                                nodes_Li_secondary[_winding],
                                nodes_Lr_primary[_winding], 
                                nodes_Li_primary[_winding], phaseSet[0])
    
                        # #########Stamp Secondary Circuit##################
                        Ylin_val, Ylin_row, Ylin_col, idx_Y = self.stamp_secondary_dual(
                                    Ylin_val, Ylin_row, Ylin_col, idx_Y,
                                    nodes_Lr_to[_pos], nodes_Li_to[_pos],
                                    nodes_Lr_to[_neg], nodes_Li_to[_neg],
                                    nodes_Lr_secondary[_winding],
                                    nodes_Li_secondary[_winding],
                                    nodes_Lr_primary[_winding],
                                    nodes_Li_primary[_winding], phaseSet[0],
                                    homotopy_option, B_factor, G_factor)
                    
        except ValueError:
            print('Illegal connect_type for the regulator')

        return Ylin_val, Ylin_row, Ylin_col, idx_Y, stamped_ground

    def adjust_effective_regulator_ratio(self, tap, change_type):
        if change_type == 'raise':
            if self.reg_type == 'A':
                aR = (1 + 0.00625 * tap)**-1
            else:
                aR = 1 - (0.00625 * tap)
        else:
            if self.reg_type == 'A':
                aR = (1 - 0.00625 * tap)**-1
            else:
                aR = 1 + (0.00625 * tap)
        return aR

