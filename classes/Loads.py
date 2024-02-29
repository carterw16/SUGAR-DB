"""Provides Load, Constant, and ZIP classes for distribution system analysis.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-11-2017
  Updated Date: 10-14-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

  Load defines the core properties and operations that apply to all loads. Constant is for loads that are constant
   power. ZIP defines loads that are dependent on fractional quantities of impedance (Z), current (I), and power (P).
   Both Constant and ZIP have their own stamp_nonlinear functions which add the linearized elements of their nonlin
   components.
"""

# import built-in modules
from __future__ import division

from itertools import count
from types import SimpleNamespace

import numpy as np
# import third-party modules
from termcolor import colored

# import custom modules
from classes.Elements import NonlinearElement
from classes.GlobalVars import connected_nodes


class Load(NonlinearElement):
    _ids = count(0)

    def __init__(self,
                 ID,
                 name,
                 parent,
                 phases,
                 is_delta=False,
                 load_class=0,
                 stamp_dual = False):
        """
        Initializes the load class.
        :param ID: The load ID.
        :param name: The name of the load.
        :param parent: The bus (node) that the load is located on.
        :param phases: Indicates the number of phases present.
        :param load_class: Specifies whether load is of type 0 - Unknown, 1-Residential, 2-Commercial,
        3-Industrial, or 4-Agricultural.
        :param voltageA: The voltage on phase A of the load.
        :param voltageB: The voltage on phase B of the load.
        :param voltageC: The voltage on phase C of the load.
        :param voltageAB: The voltage on the delta-phase AB of the load.
        :param voltageBC: The voltage on the delta-phase BC of the load.
        :param voltageCA: The voltage on the delta-phase CA of the load.
        """

        super(Load, self).__init__()

        if parent:
            self.ID = parent  # Load bus number
        else:
            self.ID = ID  # Load bus number

        self.name = name  # Name of the load
        self.Type = "Constant"
        self.phases = int(phases) if phases else None
        self.isDelta = is_delta
        if self.isDelta:
            self.isGnd = False
        else:
            self.isGnd = True

        self.load_phases = [0x01, 0x02, 0x04]
        self.stamp_dual = stamp_dual

        # Assign the load numbers
        self.ids = Load._ids.__next__()
        self.load_class = int(
            load_class
        ) if load_class else None  # [0 - Unknown, 1-R, 2-C, 3-I, 4-A]

        # Load Currents
        self.Ia_mag = 0.0
        self.Ia_ang = 0.0
        self.Ib_mag = 0.0
        self.Ib_ang = 0.0
        self.Ic_mag = 0.0
        self.Ic_ang = 0.0

        # Final Load Powers
        self.Pa = 0.0
        self.Pb = 0.0
        self.Pc = 0.0
        self.Qa = 0.0
        self.Qb = 0.0
        self.Qc = 0.0

        # Add the load id to connected nodes
        connected_nodes.add(self.ID)

    def get_nodes(self, node_key, node, V):
        """
        Find the indices of the solution vector and the values at those indices.
        :param node: The vector of all system nodes.
        :param V: The solution vector.
        :return: Two dictionaries which hold the from and to node indices and their respective values.
        """
        nodeA_Vr_from = node[node_key[self.ID]].nodeA_Vr
        nodeA_Vi_from = node[node_key[self.ID]].nodeA_Vi
        nodeB_Vr_from = node[node_key[self.ID]].nodeB_Vr
        nodeB_Vi_from = node[node_key[self.ID]].nodeB_Vi
        nodeC_Vr_from = node[node_key[self.ID]].nodeC_Vr
        nodeC_Vi_from = node[node_key[self.ID]].nodeC_Vi

        if not self.isDelta:
            # Assuming neutral is grounded
            # Find the nodes to stamp

            nodeA_Vr_to = node[node_key[self.ID]].nodeN_Vr
            nodeA_Vi_to = node[node_key[self.ID]].nodeN_Vi
            nodeB_Vr_to = node[node_key[self.ID]].nodeN_Vr
            nodeB_Vi_to = node[node_key[self.ID]].nodeN_Vi
            nodeC_Vr_to = node[node_key[self.ID]].nodeN_Vr
            nodeC_Vi_to = node[node_key[self.ID]].nodeN_Vi
        else:
            # Assuming neutral is not grounded
            # Find the nodes to stamp

            nodeA_Vr_to = node[node_key[self.ID]].nodeA_Vr
            nodeA_Vi_to = node[node_key[self.ID]].nodeA_Vi
            nodeB_Vr_to = node[node_key[self.ID]].nodeB_Vr
            nodeB_Vi_to = node[node_key[self.ID]].nodeB_Vi
            nodeC_Vr_to = node[node_key[self.ID]].nodeC_Vr
            nodeC_Vi_to = node[node_key[self.ID]].nodeC_Vi

        if self.stamp_dual:
            nodeA_Lr_from = node[node_key[self.ID]].nodeA_dual_eq_var_r
            nodeB_Lr_from = node[node_key[self.ID]].nodeB_dual_eq_var_r
            nodeC_Lr_from = node[node_key[self.ID]].nodeC_dual_eq_var_r
            
            nodeA_Li_from = node[node_key[self.ID]].nodeA_dual_eq_var_i
            nodeB_Li_from = node[node_key[self.ID]].nodeB_dual_eq_var_i
            nodeC_Li_from = node[node_key[self.ID]].nodeC_dual_eq_var_i
          
            if not self.isDelta:
                nodeA_Lr_to = node[node_key[self.ID]].nodeN_dual_eq_var_r
                nodeB_Lr_to = node[node_key[self.ID]].nodeN_dual_eq_var_r
                nodeC_Lr_to = node[node_key[self.ID]].nodeN_dual_eq_var_r
                
                nodeA_Li_to = node[node_key[self.ID]].nodeN_dual_eq_var_i
                nodeB_Li_to = node[node_key[self.ID]].nodeN_dual_eq_var_i
                nodeC_Li_to = node[node_key[self.ID]].nodeN_dual_eq_var_i
                
            else:
                nodeA_Lr_to = node[node_key[self.ID]].nodeA_dual_eq_var_r
                nodeB_Lr_to = node[node_key[self.ID]].nodeB_dual_eq_var_r
                nodeC_Lr_to = node[node_key[self.ID]].nodeC_dual_eq_var_r
                
                nodeA_Li_to = node[node_key[self.ID]].nodeA_dual_eq_var_i
                nodeB_Li_to = node[node_key[self.ID]].nodeB_dual_eq_var_i
                nodeC_Li_to = node[node_key[self.ID]].nodeC_dual_eq_var_i
                
            L_from = {'LR':[V[nodeA_Lr_from], V[nodeB_Lr_from], V[nodeC_Lr_from]],
                'nodeLR': [nodeA_Lr_from, nodeB_Lr_from, nodeC_Lr_from],
                'LI':[V[nodeA_Li_from], V[nodeB_Li_from], V[nodeC_Li_from]],
                'nodeLI': [nodeA_Li_from, nodeB_Li_from, nodeC_Li_from]}
            L_to = {'LR':[V[nodeA_Lr_to], V[nodeB_Lr_to], V[nodeC_Lr_to]],
                    'nodeLR': [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to],
                    'LI':[V[nodeA_Li_to], V[nodeB_Li_to], V[nodeC_Li_to]],
                    'nodeLI': [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to]}
            L_from = SimpleNamespace(**L_from)
            L_to = SimpleNamespace(**L_to)
            
        else:
            L_from = None
            L_to = None
        # Find the node voltages

        (Vr_nodeA_from, Vi_nodeA_from, Vr_nodeB_from, Vi_nodeB_from, Vr_nodeC_from, Vi_nodeC_from) = \
         (V[nodeA_Vr_from], V[nodeA_Vi_from], V[nodeB_Vr_from], V[nodeB_Vi_from], V[nodeC_Vr_from],
          V[nodeC_Vi_from])

        (Vr_nodeA_to, Vi_nodeA_to, Vr_nodeB_to, Vi_nodeB_to, Vr_nodeC_to, Vi_nodeC_to) = \
         (V[nodeA_Vr_to], V[nodeA_Vi_to], V[nodeB_Vr_to], V[nodeB_Vi_to], V[nodeC_Vr_to],
          V[nodeC_Vi_to])

        # Add node voltages and nodes to dictionaries

        V_from = {
            'VR': [Vr_nodeA_from, Vr_nodeB_from, Vr_nodeC_from],
            'VI': [Vi_nodeA_from, Vi_nodeB_from, Vi_nodeC_from],
            'nodeVR': [nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from],
            'nodeVI': [nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from]
        }

        V_to = {
            'VR': [Vr_nodeA_to, Vr_nodeB_to, Vr_nodeC_to],
            'VI': [Vi_nodeA_to, Vi_nodeB_to, Vi_nodeC_to],
            'nodeVR': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to],
            'nodeVI': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]
        }
        
        V_from = SimpleNamespace(**V_from)
        V_to = SimpleNamespace(**V_to)

        return V_from, V_to, L_from, L_to

    def _load_calculations(self, _cP, _cQ, _cG, _cB, _cIr, _cIi, 
                            node_Vr_from, node_Vi_from, node_Vr_to, node_Vi_to, 
                            Vr_node_from, Vi_node_from, Vr_node_to, Vi_node_to, 
                            Ynlin_val, Ynlin_row, Ynlin_col, 
                            Jnlin_val, Jnlin_row, idx_Y, idx_J, 
                            node_Lr_from, node_Li_from, node_Lr_to, node_Li_to,
                            Lr_from, Li_from, Lr_to, Li_to):

        """
				Stamps the partial derivatives and their phase changes into the admittance (Y) matrix.

				:param _cP: Constant real power.
				:param _cQ: Constant reactive power.
				:param _cG: Constant conductance.
				:param _cB: Constant susceptance.
				:param _cIr: Constant real current.
				:param _cIi: Constant imaginary current.
				:param node_Vr_from: From node of the real circuit.
				:param node_Vi_from: From node of the imaginary circuit.
				:param node_Vr_to: To node of the real circuit.
				:param node_Vi_to: To node of the imaginary circuit.
				:param Vr_node_from: Voltage of the from node in the real circuit.
				:param Vi_node_from: Voltage of the from node in the imaginary circuit.
				:param Vr_node_to: Voltage of the to node in the real circuit.
				:param Vi_node_to: Voltage of the to node in the imaginary circuit.
				:param Ynlin_val: The value of the partial derivative to be stamped into the admittance matrix.
				:param Ynlin_row: The admittance matrix row.
				:param Ynlin_col: The admittance matrix column.
				:param Jnlin_val: The value of the partial derivative to be stamped into the excitation vector.
				:param Jnlin_row: The row of the excitation vector.
				:param _isGnd: Indicated whether load is grounded.
				:return: None
			"""

        if not self.isGnd:
            # Calculate partials for PQ loads
            Vr_across = Vr_node_from - Vr_node_to
            Vi_across = Vi_node_from - Vi_node_to
            
            if self.stamp_dual:
                Lr_across = Lr_from - Lr_to
                Li_across = Li_from - Lr_to
            else:
                Lr_across = 0
                Li_across = 0

            # Calculate Partials from the voltage difference between from and to node
            (_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVil, _dIildVrl, Lr_hist, Li_hist, \
                dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi, 
                dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi) = \
                Constant.calc_partials(self, _cP, _cQ, Vr_across, Vi_across, Lr_across, Li_across)

            # 4 real and 4 imaginary stamps for the from circuit
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from, _dIrldVrl + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)  # G
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_to, -_dIrldVrl - _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from, _dIrldVil - _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_to, -_dIrldVil + _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from, _dIildVil + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_to, -_dIildVil - _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from, _dIildVrl + _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_to, -_dIildVrl - _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            # 4 real and imaginary circuits for the to circuit
            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_from, -_dIrldVrl - _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_to, _dIrldVrl + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_from, -_dIrldVil + _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_to, _dIrldVil - _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_from, -_dIildVil - _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_to, _dIildVil + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_from, -_dIildVrl - _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_to, _dIildVrl + _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

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
                # LR_from
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
 
        else:
            # Calculate partials for PQ loads
            Vr_across = Vr_node_from
            Vi_across = Vi_node_from
            
            if self.stamp_dual:
                Lr_across = Lr_from 
                Li_across = Li_from 
            else:
                Lr_across = 0
                Li_across = 0

            # Calculate Partials from the voltage difference between from and to node
            (_Irl_hist, _Iil_hist, _dIrldVrl, _dIrldVil, _dIildVil, _dIildVrl, Lr_hist, Li_hist,
                dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi, 
                dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi) = \
                Constant.calc_partials(self, _cP, _cQ, Vr_across, Vi_across, Lr_across, Li_across)

            # 4 real and 4 imaginary stamps for the from circuit
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from, _dIrldVrl + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)  # G
            #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from, _dIrldVil - _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from, _dIildVil + _cG,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)
            #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from, _dIildVrl + _cB,
                                Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

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


class Constant(Load):

    def __init__(self,
                 ID,
                 name,
                 parent,
                 phases, 
                 is_delta=False,
                 load_class=0,
                 constant_powerA=0,
                 constant_powerB=0,
                 constant_powerC=0,
                 constant_currentA=0,
                 constant_currentB=0,
                 constant_currentC=0,
                 constant_impedanceA=None,
                 constant_impedanceB=None,
                 constant_impedanceC=None,
                 constant_powerAB=0,
                 constant_powerBC=0,
                 constant_powerCA=0,
                 nom_voltage = 208,
                 stamp_dual = False):
        """Initialize the Constant Load class.
		:param ID: The load ID.
		:param name: The name of the load.
		:param parent: The bus (node) that the load is located on.
		:param phases: Indicates the number of phases present.
		:param load_class: Specifies whether load is of type 0 - Unknown, 1-Residential, 2-Commercial,
		3-Industrial, or 4-Agricultural.
		:param constant_powerA: The constant power value of the load at phase A in a wye-connection.
		:param constant_powerB: The constant power value of the load at phase B in a wye-connection.
		:param constant_powerC: The constant power value of the load at phase C in a wye-connection.
		:param constant_currentA: The constant current value of the load at phase A in a wye-connection.
		:param constant_currentB: The constant current value of the load at phase B in a wye-connection.
		:param constant_currentC: The constant current value of the load at phase C in a wye-connection.
		:param constant_impedanceA: The constant impedance value of the load at phase A in a wye-connection.
		:param constant_impedanceB: The constant impedance value of the load at phase B in a wye-connection.
		:param constant_impedanceC: The constant impedance value of the load at phase C in a wye-connection.
		:param constant_powerAB: The constant power value of the load at phase AB in a delta-connection.
		:param constant_powerBC: The constant power value of the load at phase BC in a delta-connection.
		:param constant_powerCA: The constant power value of the load at phase CA in a delta-connection.
		"""
        # Make the Const_Load class inherit all the methods and properties of the Load class.
        super(Constant, self).__init__(ID, name, parent, phases,
                                       is_delta, load_class, stamp_dual = stamp_dual)
        self.nom_voltage = nom_voltage
        # Const_Load properties
        self.cP_A = constant_powerA  # [VA]
        self.cP_B = constant_powerB  # [VA]
        self.cP_C = constant_powerC  # [VA]
        self.cP_AB = constant_powerAB  # [VA]
        self.cP_BC = constant_powerBC  # [VA]
        self.cP_CA = constant_powerCA  # [VA]
        self.cI_A = constant_currentA  # [Amp]
        self.cI_B = constant_currentB  # [Amp]
        self.cI_C = constant_currentC  # [Amp]
        self.cZ_A = constant_impedanceA  # [ohms]
        self.cZ_B = constant_impedanceB  # [ohms]
        self.cZ_C = constant_impedanceC  # [ohms]

        (self.cG_A, self.cB_A) = self.calc_G_B(float(self.cZ_A.real),
                                               float(self.cZ_A.imag))
        (self.cG_B, self.cB_B) = self.calc_G_B(float(self.cZ_B.real),
                                               float(self.cZ_B.imag))
        (self.cG_C, self.cB_C) = self.calc_G_B(float(self.cZ_C.real),
                                               float(self.cZ_C.imag))

        # Check for consistency between phase and load values
        if self.cP_A or self.cZ_A or self.cI_A or self.cP_AB or self.cP_CA:
            if self.phases & 0x01 != int(0x01):
                print(
                    colored(
                        'Warning: Load %s value for phase A given when phase not connected',
                        'yellow') % self.name)
        if self.cP_B or self.cZ_B or self.cI_B or self.cP_AB or self.cP_BC:
            if self.phases & 0x02 != int(0x02):
                print(
                    colored(
                        'Warning: Load %s value for phase B given when phase not connected',
                        'yellow') % self.name)
        if self.cP_C or self.cZ_C or self.cI_C or self.cP_BC or self.cP_CA:
            if self.phases & 0x04 != int(0x04):
                print(
                    colored(
                        'Warning: Load %s value for phase C given when phase not connected',
                        'yellow') % self.name)

        # Collect Power Parameters into a list
        self._cP = [self.cP_A.real, self.cP_B.real, self.cP_C.real]
        self._cQ = [self.cP_A.imag, self.cP_B.imag, self.cP_C.imag]
        self._cG = [self.cG_A, self.cG_B, self.cG_C]
        self._cB = [self.cB_A, self.cB_B, self.cB_C]
        self._cIr = [self.cI_A.real, self.cI_B.real, self.cI_C.real]
        self._cIi = [self.cI_A.imag, self.cI_B.imag, self.cI_C.imag]

        self._cP_LL = [self.cP_AB.real, self.cP_BC.real, self.cP_CA.real]
        self._cQ_LL = [self.cP_AB.imag, self.cP_BC.imag, self.cP_CA.imag]

    @staticmethod
    def calc_G_B(r, x):
        """
		Calculate the conductance (G) and susceptance (B) params that will be used for phase shifting.
		:param r: resistance
		:param x: admittance
		:return: G and B if at least r or x exists.
		"""
        if r or x:
            G = r / (r**2 + x**2)
            B = -x / (r**2 + x**2)
            return G, B
        else:
            return 0, 0

    def calc_partials(self, P, Q, Vr, Vi, Lr, Li):
        """
		Calculate partial derivatives for the constant load.
		:param P: Real power of the load.
		:param Q: Reactive power of the load.
		:param Vr: Real voltage of the load.
		:param Vi: Imaginary voltage of the load.
		:return:
		"""
        QVrVi = Q * Vr * Vi
        PVrVi = P * Vr * Vi
        Vr2Vi2 = Vr * Vr - Vi * Vi
        Vmag2 = Vr * Vr + Vi * Vi
        Vmag4 = Vmag2 * Vmag2
        Irl_hist = (P * Vr + Q * Vi) / Vmag2 if Vmag2 != 0 else 0
        Iil_hist = (P * Vi - Q * Vr) / Vmag2 if Vmag2 != 0 else 0
        dIrldVrl = -(P * Vr2Vi2 + 2.0 * QVrVi) / Vmag4  if Vmag4 != 0 else 0# dIr/dVr
        dIrldVil = (Q * Vr2Vi2 - 2.0 * PVrVi) / Vmag4  if Vmag4 != 0 else 0# dIr/dVi
        dIildVil = (P * Vr2Vi2 + 2.0 * QVrVi) / Vmag4  if Vmag4 != 0 else 0# dIi/dVi
        dIildVrl = dIrldVil  # dIi/dVi
        
        if self.stamp_dual:
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
            
        return (Irl_hist, Iil_hist, dIrldVrl, dIrldVil, dIildVil, dIildVrl, Lr_hist, Li_hist,
                dLr_dVr, dLr_dVi, dLr_dLr, dLr_dLi, 
                dLi_dVr, dLi_dVi, dLi_dLr, dLi_dLi)

    def stamp_constant_phase_load(self, node_key, node, V, lf, Ynlin_val,
                                  Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
                                  idx_Y, idx_J):

        V_from, V_to, L_from, L_to = self.get_nodes(node_key, node, V)

        nodeA_Vr_to = V_to.nodeVR[0]
        nodeA_Vi_to = V_to.nodeVI[0]
        nodeB_Vr_to = V_to.nodeVR[1]
        nodeB_Vi_to = V_to.nodeVI[1]
        nodeC_Vr_to = V_to.nodeVR[2]
        nodeC_Vi_to = V_to.nodeVI[2]

        node_Vr_from = [V_from.nodeVR[0], V_from.nodeVR[1], V_from.nodeVR[2]]
        node_Vi_from = [V_from.nodeVI[0], V_from.nodeVI[1], V_from.nodeVI[2]]

        Vr_from = [V_from.VR[0], V_from.VR[1], V_from.VR[2]]
        Vi_from = [V_from.VI[0], V_from.VI[1], V_from.VI[2]]
        
        if self.stamp_dual:
            nodeA_Lr_to = L_to.nodeLR[0]
            nodeA_Li_to = L_to.nodeLI[0]
            nodeB_Lr_to = L_to.nodeLR[1]
            nodeB_Li_to = L_to.nodeLI[1]
            nodeC_Lr_to = L_to.nodeLR[2]
            nodeC_Li_to = L_to.nodeLI[2]
            
            node_Lr_from = [L_from.nodeLR[0], L_from.nodeLR[1], L_from.nodeLR[2]]
            node_Li_from = [L_from.nodeLI[0], L_from.nodeLI[1], L_from.nodeLI[2]]
            
            Lr_from = [L_from.LR[0], L_from.LR[1], L_from.LR[2]]
            Li_from = [L_from.LI[0], L_from.LI[1], L_from.LI[2]]
        else:
            node_Lr_from = [0, 0, 0]
            node_Li_from = [0, 0, 0]
            
            Lr_from = [0, 0, 0]
            Li_from = [0, 0, 0]
            
            node_Lr_to = [0, 0, 0]
            node_Li_to = [0, 0, 0]
            
            Lr_to = [0, 0, 0]
            Li_to = [0, 0, 0]
            
        if self.isDelta:
            node_Vr_to = [nodeB_Vr_to, nodeC_Vr_to, nodeA_Vr_to]
            node_Vi_to = [nodeB_Vi_to, nodeC_Vi_to, nodeA_Vi_to]

            Vr_to = [V_to.VR[1], V_to.VR[2], V_to.VR[0]]
            Vi_to = [V_to.VI[1], V_to.VI[2], V_to.VI[0]]
            
            if self.stamp_dual:
                node_Lr_to = [nodeB_Lr_to, nodeC_Lr_to, nodeA_Lr_to]
                node_Li_to = [nodeB_Li_to, nodeC_Li_to, nodeA_Li_to]
                
                Lr_to = [L_to.LR[1], L_to.LR[2], L_to.LR[0]]
                Li_to = [L_to.LI[1], L_to.LI[2], L_to.LI[0]]
               
        else:
            node_Vr_to = [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to]
            node_Vi_to = [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]

            Vr_to = [V_to.VR[0], V_to.VR[1], V_to.VI[2]]
            Vi_to = [V_to.VI[0], V_to.VI[1], V_to.VI[2]]
            
            if self.stamp_dual:
                node_Lr_to = [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to]
                node_Li_to = [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to]
                
                Lr_to = [L_to.LR[0], L_to.LR[1], L_to.LR[2]]
                Li_to = [L_to.LI[0], L_to.LI[1], L_to.LI[2]]                

        for _phase in range(3):
            # scale initial loading by load factor or power stepping factor if necessary
            _cP = self._cP[_phase] * lf
            _cQ = self._cQ[_phase] * lf
            _cG = self._cG[_phase] * lf
            _cB = self._cB[_phase] * lf
            _cIr = self._cIr[_phase] * lf
            _cIi = self._cIi[_phase] * lf

            # stamp phase load if constant power params present
            if _cP or _cQ or _cG or _cB or _cIr or _cIi:
                if self.load_phases[_phase] & self.phases == int(
                        self.load_phases[_phase]):
                    idx_Y, idx_J = self._load_calculations(
                        _cP, _cQ, _cG, _cB, _cIr, _cIi, node_Vr_from[_phase],
                        node_Vi_from[_phase], node_Vr_to[_phase],
                        node_Vi_to[_phase], Vr_from[_phase], Vi_from[_phase],
                        Vr_to[_phase], Vi_to[_phase], Ynlin_val, Ynlin_row,
                        Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J,
                        node_Lr_from[_phase], node_Li_from[_phase], node_Lr_to[_phase], node_Li_to[_phase],
                        Lr_from[_phase], Li_from[_phase], Lr_to[_phase], Li_to[_phase])

        return idx_Y, idx_J

    def stamp_constant_delta_load(self, node_key, node, V, lf, Ynlin_val,
                                  Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
                                  idx_Y, idx_J):

        V_from, V_to, L_from, L_to = self.get_nodes(node_key, node, V)

        nodeA_Vr_to = V_to.nodeVR[0]
        nodeA_Vi_to = V_to.nodeVI[0]
        nodeB_Vr_to = V_to.nodeVR[1]
        nodeB_Vi_to = V_to.nodeVI[1]
        nodeC_Vr_to = V_to.nodeVR[2]
        nodeC_Vi_to = V_to.nodeVI[2]

        node_Vr_from = [V_from.nodeVR[0], V_from.nodeVR[1], V_from.nodeVR[2]]
        node_Vi_from = [V_from.nodeVI[0], V_from.nodeVI[1], V_from.nodeVI[2]]
        node_Vr_to = [nodeB_Vr_to, nodeC_Vr_to, nodeA_Vr_to]
        node_Vi_to = [nodeB_Vi_to, nodeC_Vi_to, nodeA_Vi_to]

        Vr_from = [V_from.VR[0], V_from.VR[1], V_from.VR[2]]
        Vi_from = [V_from.VI[0], V_from.VI[1], V_from.VI[2]]
        Vr_to = [V_to.VR[1], V_to.VR[2], V_to.VR[0]]
        Vi_to = [V_to.VI[1], V_to.VI[2], V_to.VI[0]]
        
        if self.stamp_dual:
            nodeA_Lr_to = L_to.nodeLR[0]
            nodeA_Li_to = L_to.nodeLI[0]
            nodeB_Lr_to = L_to.nodeLR[1]
            nodeB_Li_to = L_to.nodeLI[1]
            nodeC_Lr_to = L_to.nodeLR[2]
            nodeC_Li_to = L_to.nodeLI[2]
            
            node_Lr_from = [L_from.nodeLR[0], L_from.nodeLR[1], L_from.nodeLR[2]]
            node_Li_from = [L_from.nodeLI[0], L_from.nodeLI[1], L_from.nodeLI[2]]

            Lr_from = [L_from.LR[0], L_from.LR[1], L_from.LR[2]]
            Li_from = [L_from.LI[0], L_from.LI[1], L_from.LI[2]]
            
            node_Lr_to = [nodeB_Lr_to, nodeC_Lr_to, nodeA_Lr_to]
            node_Li_to = [nodeB_Li_to, nodeC_Li_to, nodeA_Li_to]
            
            Lr_to = [L_to.LR[1], L_to.LR[2], L_to.LR[0]]
            Li_to = [L_to.LI[1], L_to.LI[2], L_to.LI[0]]         
            

        else:
            node_Lr_from = None
            node_Li_from = None
            
            Lr_from = None
            Li_from = None
            
            node_Lr_to = None
            node_Li_to = None
            
            Lr_to = None
            Li_to = None       
            

        for _idx in range(3):
            _cP_LL_x = self._cP_LL[_idx] * lf 
            _cQ_LL_x = self._cQ_LL[_idx] * lf 

            if _cP_LL_x or _cQ_LL_x:
                idx_Y, idx_J = self._load_calculations(
                    _cP_LL_x, _cQ_LL_x, 0, 0, 0, 0, node_Vr_from[_idx],
                    node_Vi_from[_idx], node_Vr_to[_idx], node_Vi_to[_idx],
                    Vr_from[_idx], Vi_from[_idx], Vr_to[_idx], Vi_to[_idx],
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
                    idx_Y, idx_J,
                    node_Lr_from[_idx], node_Li_from[_idx], node_Lr_to[_idx], node_Li_to[_idx],
                        Lr_from[_idx], Li_from[_idx], Lr_to[_idx], Li_to[_idx])
        
        return idx_Y, idx_J

    def stamp_nonlinear(self,
                        node_key,
                        node,
                        V,
                        lf,
                        Ynlin_val,
                        Ynlin_row,
                        Ynlin_col,
                        Jnlin_val,
                        Jnlin_row,
                        idx_Y,
                        idx_J,
                        enable_CM=False):
        """
        Stamp the nonlinear linearized partials for the Constant Load.

        Args:
            Ynlin_val: Stores the admittance matrix value that is stamped. Type float.
            Ynlin_row: Stores the admittance matrix row number where the stamp is made. Type int.
            Ynlin_col: Stores the admittance matrix column number where the stamp is made. Type int.
            Jnlin_val: Stores the excitation vector value that is stamped. Type float.
            Jnlin_row: Stores the excitation vector row number where the stamp is made. Type int.
            node: The vector containing a list of all nodes.
            idx_Y: Admittance matrix stamp index. Type int.
            idx_J: Excitation vector stamp index. Type int.
            V: The solution vector.
            lf: Load factor. Type float.
            enable_CM: Turn current measures on/off. Disabled for now.

        Returns:
            Excitation vector stamp index and admittance matrix stamp index. 
        """

        idx_Y, idx_J = self.stamp_constant_phase_load(node_key, node, V, lf,
                                                      Ynlin_val, Ynlin_row,
                                                      Ynlin_col, Jnlin_val,
                                                      Jnlin_row, idx_Y, idx_J)

        if any(self._cP_LL) > 0 or any(self._cQ_LL) > 0:
            idx_Y, idx_J = self.stamp_constant_delta_load(
                node_key, node, V, lf, Ynlin_val, Ynlin_row, Ynlin_col,
                Jnlin_val, Jnlin_row, idx_Y, idx_J)


        return idx_Y, idx_J
    
    def calc_residual(self, Vsol, node_key, node, res_eqn, lf):
        V_from, V_to, L_from, L_to = self.get_nodes(node_key, node, Vsol)
        
        nodeA_Vr_to = V_to.nodeVR[0]
        nodeA_Vi_to = V_to.nodeVI[0]
        nodeB_Vr_to = V_to.nodeVR[1]
        nodeB_Vi_to = V_to.nodeVI[1]
        nodeC_Vr_to = V_to.nodeVR[2]
        nodeC_Vi_to = V_to.nodeVI[2]

        node_Vr_from = [V_from.nodeVR[0], V_from.nodeVR[1], V_from.nodeVR[2]]
        node_Vi_from = [V_from.nodeVI[0], V_from.nodeVI[1], V_from.nodeVI[2]]
        
        Vr_from = [V_from.VR[0], V_from.VR[1], V_from.VR[2]]
        Vi_from = [V_from.VI[0], V_from.VI[1], V_from.VI[2]]
        
        if self.isDelta:
            node_Vr_to = [nodeB_Vr_to, nodeC_Vr_to, nodeA_Vr_to]
            node_Vi_to = [nodeB_Vi_to, nodeC_Vi_to, nodeA_Vi_to]
            
            Vr_to = [V_to.VR[1], V_to.VR[2], V_to.VR[0]]
            Vi_to = [V_to.VI[1], V_to.VI[2], V_to.VI[0]]
        
        else:
            node_Vr_to = [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to]
            node_Vi_to = [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to]
            
            Vr_to = [V_to.VR[0], V_to.VR[1], V_to.VR[2]]
            Vi_to = [V_to.VI[0], V_to.VI[1], V_to.VI[2]]
        
        if self.stamp_dual:
            nodeA_Lr_to = L_to.nodeLR[0]
            nodeA_Li_to = L_to.nodeLI[0]
            nodeB_Lr_to = L_to.nodeLR[1]
            nodeB_Li_to = L_to.nodeLI[1]
            nodeC_Lr_to = L_to.nodeLR[2]
            nodeC_Li_to = L_to.nodeLI[2]
            
            node_Lr_from = [L_from.nodeLR[0], L_from.nodeLR[1], L_from.nodeLR[2]]
            node_Li_from = [L_from.nodeLI[0], L_from.nodeLI[1], L_from.nodeLI[2]]

            Lr_from = [L_from.LR[0], L_from.LR[1], L_from.LR[2]]
            Li_from = [L_from.LI[0], L_from.LI[1], L_from.LI[2]]
            
            if self.isDelta:
                node_Lr_to = [nodeB_Lr_to, nodeC_Lr_to, nodeA_Lr_to]
                node_Li_to = [nodeB_Li_to, nodeC_Li_to, nodeA_Li_to]
                
                Lr_to = [L_to.LR[1], L_to.LR[2], L_to.LR[0]]
                Li_to = [L_to.LI[1], L_to.LI[2], L_to.LI[0]]  
            else:
                node_Lr_to = [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to]
                node_Li_to = [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to]
                
                Lr_to = [L_to.LR[0], L_to.LR[1], L_to.LR[2]]
                Li_to = [L_to.LI[0], L_to.LI[1], L_to.LI[2]]  

        else:
            node_Lr_from = [0, 0, 0]
            node_Li_from = [0, 0, 0]
            
            Lr_from = [0, 0, 0]
            Li_from = [0, 0, 0]
            
            node_Lr_to = [0, 0, 0]
            node_Li_to = [0, 0, 0]
            
            Lr_to = [0, 0, 0]
            Li_to = [0, 0, 0] 
        
        
        for _phase in range(3):                        
            if not self.isDelta: 
                if self.isGnd:
                    _cP = self._cP[_phase] * lf
                    _cQ = self._cQ[_phase] * lf 
                    _cG = self._cG[_phase] * lf
                    _cB = self._cB[_phase] * lf
                    _cIr = self._cIr[_phase] * lf
                    _cIi = self._cIi[_phase] * lf
                       
                    # Equality
                    ipq_from_denom = (Vr_from[_phase])**2 + (Vi_from[_phase])**2
                    ipq_from_to_r = _cP*(Vr_from[_phase]) + _cQ*(Vi_from[_phase])
                    ipq_from_to_i = _cP*(Vi_from[_phase]) - _cQ*(Vr_from[_phase]) 
                    
                    if self.stamp_dual:
                        # Stationarity
                        dipq_dvr_from_r = -(_cP * (Vr_from[_phase] * Vr_from[_phase] - Vi_from[_phase] * Vi_from[_phase]) + 2.0 * _cQ * Vr_from[_phase] * Vi_from[_phase])  # dIr/dVr
                        dipq_dvi_from_r = (_cQ * (Vr_from[_phase] * Vr_from[_phase] - Vi_from[_phase] * Vi_from[_phase]) - 2.0 * _cP * Vr_from[_phase] * Vi_from[_phase])  # dIr/dVi
                        dipq_dvi_from_i = (_cP * (Vr_from[_phase] * Vr_from[_phase] - Vi_from[_phase] * Vi_from[_phase]) + 2.0 * _cQ * Vr_from[_phase] * Vi_from[_phase])  # dIi/dVi
                        dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                        dipq_from_denom = ipq_from_denom**2       
                    
                    if abs(ipq_from_denom) > 1e-6:
                        res_eqn[node_Vr_from[_phase]] += _cIr + ipq_from_to_r/ipq_from_denom + _cG * Vr_from[_phase] - _cB * Vi_from[_phase]
                        res_eqn[node_Vi_from[_phase]] += _cIi + ipq_from_to_i/ipq_from_denom + _cB * Vr_from[_phase] + _cG * Vi_from[_phase]
                        
                        if self.stamp_dual:
                            res_eqn[node_Lr_from[_phase]] += dipq_dvr_from_r/dipq_from_denom * Lr_from[_phase] + dipq_dvi_from_r/dipq_from_denom * Li_from[_phase]
                            res_eqn[node_Li_from[_phase]] += dipq_dvr_from_i/dipq_from_denom * Lr_from[_phase] + dipq_dvi_from_i/dipq_from_denom * Li_from[_phase]

                else:                    
                    _cP = self._cP[_phase] * lf
                    _cQ = self._cQ[_phase] * lf 
                    _cG = self._cG[_phase] * lf
                    _cB = self._cB[_phase] * lf
                    _cIr = self._cIr[_phase] * lf
                    _cIi = self._cIi[_phase] * lf
                    
                    ipq_from_denom = (Vr_from[_phase] - Vr_to[_phase])**2 + (Vi_from[_phase] - Vi_to[_phase])**2
                    ipq_from_to_r = _cP*(Vr_from[_phase]-Vr_to[_phase]) + _cQ*(Vi_from[_phase] - Vi_to[_phase])
                    ipq_from_to_i = _cP*(Vi_from[_phase]-Vi_to[_phase]) - _cQ*(Vr_from[_phase] - Vr_to[_phase])
                    
                    ipq_to_denom = (Vr_to[_phase] - Vr_from[_phase])**2 + (Vi_to[_phase] - Vi_from[_phase])**2
                    ipq_to_from_r = _cP*(Vr_to[_phase] - Vr_from[_phase]) + _cQ*(Vi_to[_phase] - Vi_from[_phase])
                    ipq_to_from_i = _cP*(Vi_to[_phase] - Vi_from[_phase]) - _cQ*(Vr_to[_phase] - Vr_from[_phase])
                                    
                    if self.stamp_dual:
                        # Stationarity
                        dipq_dvr_from_r = -(_cP * ((Vr_from[_phase] - Vr_to[_phase])**2 - (Vi_from[_phase] - Vi_to[_phase])**2) + 2.0 * _cQ * Vr_from[_phase] * Vi_from[_phase])  # dIr/dVr
                        dipq_dvi_from_r = (_cQ * ((Vr_from[_phase] - Vr_to[_phase])**2 - (Vi_from[_phase] - Vi_to[_phase])**2) - 2.0 * _cP * Vr_from[_phase] * Vi_from[_phase])  # dIr/dVi
                        dipq_dvi_from_i = (_cP * ((Vr_from[_phase] - Vr_to[_phase])**2 - (Vi_from[_phase] - Vi_to[_phase])**2) + 2.0 * _cQ * Vr_from[_phase] * Vi_from[_phase])  # dIi/dVi
                        dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                        dipq_from_denom = ipq_from_denom**2     

                        # Stationarity
                        # dipq_dvr_to_r = -(_cP * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) + 2.0 * _cQ * Vr_to[_phase] * Vi_to[_phase])  # dIr/dVr
                        # dipq_dvi_to_r = (_cQ * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) - 2.0 * _cP * Vr_to[_phase] * Vi_to[_phase])  # dIr/dVi
                        # dipq_dvi_to_i = (_cP * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) + 2.0 * _cQ * Vr_to[_phase] * Vi_to[_phase])  # dIi/dVi
                        # dipq_dvr_to_i = dipq_dvi_to_r  # dIi/dVi
                        # dipq_to_denom = ipq_to_denom**2                       
                    
                    if abs(ipq_from_denom) > 1e-6:
                        res_eqn[node_Vr_from[_phase]] += _cIr + ipq_from_to_r/ipq_from_denom + _cG * (Vr_from[_phase] - Vr_to[_phase]) - _cB * (Vi_from[_phase] - Vi_to[_phase])
                        res_eqn[node_Vi_from[_phase]] += _cIi + ipq_from_to_i/ipq_from_denom + _cB * (Vr_from[_phase] - Vr_to[_phase]) + _cG * (Vi_from[_phase] - Vi_to[_phase])
                        res_eqn[node_Vr_to[_phase]] += -_cIr + ipq_to_from_r/ipq_to_denom - _cG * (Vr_from[_phase] - Vr_to[_phase]) + _cB * (Vi_from[_phase] - Vi_to[_phase])
                        res_eqn[node_Vi_to[_phase]] += -_cIi + ipq_to_from_i/ipq_to_denom - _cB * (Vr_from[_phase] - Vr_to[_phase]) - _cG * (Vi_from[_phase] - Vi_to[_phase])
                        
                        if self.stamp_dual:
                            res_eqn[node_Lr_from[_phase]] += dipq_dvr_from_r/dipq_from_denom * (Lr_from[_phase] - Lr_to[_phase]) + dipq_dvi_from_r/dipq_from_denom * (Li_from[_phase] - Li_to[_phase])
                            res_eqn[node_Li_from[_phase]] += dipq_dvr_from_i/dipq_from_denom * (Lr_from[_phase] - Lr_to[_phase]) + dipq_dvi_from_i/dipq_from_denom * (Li_from[_phase] - Li_to[_phase])
                    
                            res_eqn[node_Lr_to[_phase]] += -(dipq_dvr_from_r/dipq_from_denom * (Lr_from[_phase] - Lr_to[_phase]) + dipq_dvi_from_r/dipq_from_denom * (Li_from[_phase] - Li_to[_phase]))#dipq_dvr_to_r/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_r/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase])
                            res_eqn[node_Li_to[_phase]] += -(dipq_dvr_from_i/dipq_from_denom * (Lr_from[_phase] - Lr_to[_phase]) + dipq_dvi_from_i/dipq_from_denom * (Li_from[_phase] - Li_to[_phase]))#dipq_dvr_to_i/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_i/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase]) 
                             
            elif self.isDelta:                
                _cP = self._cP[_phase] * lf
                _cQ = self._cQ[_phase] * lf 
                _cG = self._cG[_phase] * lf
                _cB = self._cB[_phase] * lf 
                _cIr = self._cIr[_phase] * lf
                _cIi = self._cIi[_phase] * lf
                
                ipq_from_denom = (Vr_from[_phase] - Vr_to[_phase])**2 + (Vi_from[_phase] - Vi_to[_phase])**2
                ipq_from_to_r = _cP*(Vr_from[_phase]-Vr_to[_phase]) + _cQ*(Vi_from[_phase] - Vi_to[_phase])
                ipq_from_to_i = _cP*(Vi_from[_phase]-Vi_to[_phase]) - _cQ*(Vr_from[_phase] - Vr_to[_phase])
                
                ipq_to_denom =  ipq_from_denom #(Vr_to[_phase] - Vr_from[_phase])**2 + (Vi_to[_phase] - Vi_from[_phase])**2
                ipq_to_from_r = -ipq_from_to_r# _cP*(Vr_to[_phase] - Vr_from[_phase]) + _cQ*(Vi_to[_phase] - Vi_from[_phase])
                ipq_to_from_i = -ipq_from_to_i# _cP*(Vi_to[_phase] - Vi_from[_phase]) - _cQ*(Vr_to[_phase] - Vr_from[_phase])
                
                if self.stamp_dual:
                        Vr_across = (Vr_from[_phase] - Vr_to[_phase])
                        Vi_across = (Vi_from[_phase] - Vi_to[_phase])
                        # Stationarity
                        dipq_dvr_from_r = -(_cP * (Vr_across**2 - Vi_across**2) + 2.0 * _cQ * Vr_across * Vi_across)  # dIr/dVr
                        dipq_dvi_from_r = (_cQ * (Vr_across**2 - Vi_across**2) - 2.0 * _cP * Vr_across * Vi_across)  # dIr/dVi
                        dipq_dvi_from_i = (_cP * (Vr_across**2 - Vi_across**2) + 2.0 * _cQ * Vr_across * Vi_across)  # dIi/dVi
                        dipq_dvr_from_i = dipq_dvi_from_r  # dIi/dVi
                        dipq_from_denom = ipq_from_denom**2     

                        # # Stationarity
                        # dipq_dvr_to_r = -(_cP * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) + 2.0 * _cQ * Vr_to[_phase] * Vi_to[_phase])  # dIr/dVr
                        # dipq_dvi_to_r = (_cQ * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) - 2.0 * _cP * Vr_to[_phase] * Vi_to[_phase])  # dIr/dVi
                        # dipq_dvi_to_i = (_cP * ((Vr_to[_phase] - Vr_from[_phase])**2 - (Vi_to[_phase] - Vi_from[_phase])**2) + 2.0 * _cQ * Vr_to[_phase] * Vi_to[_phase])  # dIi/dVi
                        # dipq_dvr_to_i = dipq_dvi_to_r  # dIi/dVi
                        # dipq_to_denom = ipq_to_denom**2  
                        
                if abs(ipq_from_denom) > 1e-6:
                    res_eqn[node_Vr_from[_phase]] +=_cIr + ipq_from_to_r/ipq_from_denom + _cG * (Vr_from[_phase] - Vr_to[_phase]) - _cB * (Vi_from[_phase] - Vi_to[_phase])
                    res_eqn[node_Vi_from[_phase]] +=_cIi + ipq_from_to_i/ipq_from_denom + _cB * (Vr_from[_phase] - Vr_to[_phase]) + _cG * (Vi_from[_phase] - Vi_to[_phase])
                    
                    res_eqn[node_Vr_to[_phase]] += -_cIr + ipq_to_from_r/ipq_to_denom - _cG * (Vr_from[_phase] - Vr_to[_phase]) + _cB * (Vi_from[_phase] - Vi_to[_phase])
                    res_eqn[node_Vi_to[_phase]] += -_cIi + ipq_to_from_i/ipq_to_denom - _cB * (Vr_from[_phase] - Vr_to[_phase]) - _cG * (Vi_from[_phase] - Vi_to[_phase])
                    
                    if self.stamp_dual:
                        Lr_across = Lr_from[_phase] - Lr_to[_phase]
                        Li_across = Li_from[_phase] - Li_to[_phase]
                        
                        res_eqn[node_Lr_from[_phase]] += dipq_dvr_from_r/dipq_from_denom * Lr_across + dipq_dvi_from_r/dipq_from_denom * Li_across
                        res_eqn[node_Li_from[_phase]] += dipq_dvr_from_i/dipq_from_denom * Lr_across + dipq_dvi_from_i/dipq_from_denom * Li_across
                
                        res_eqn[node_Lr_to[_phase]] += -(dipq_dvr_from_r/dipq_from_denom * Lr_across + dipq_dvi_from_r/dipq_from_denom * Li_across)#dipq_dvr_to_r/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_r/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase])
                        res_eqn[node_Li_to[_phase]] += -(dipq_dvr_from_i/dipq_from_denom * Lr_across + dipq_dvi_from_i/dipq_from_denom * Li_across)#dipq_dvr_to_i/dipq_to_denom * (- Lr_from[_phase] + Lr_to[_phase]) + dipq_dvi_to_i/dipq_to_denom * (- Li_from[_phase] + Li_to[_phase]) 
                
        return res_eqn
    
class Exponential(Load):

    def __init__(self,
                 ID,
                 name,
                 parent,
                 phases,
                 is_delta,
                 nominal_voltage,
                 power_A,
                 power_B,
                 power_C,
                 power_AB,
                 power_BC,
                 power_CA,
                 CVRwatts,
                 CVRvars,
                 load_class=0):
        """The exponential load model is a general load model used to quantify CVR impact.

		The exponential load model highlights the voltage dependency of P and Q using exponential parameters
		for the Conservation Voltage Reduction (CVR) factor.
		The governing equations for exponential loads are:
		P = P0 * (V/V0) ** CVRwatts
		Q = Q0 * (V/V0) ** CVRvars
		where V0, P0, and Q0 are the initial load voltage magnitude, real power, and active power. V is the load
		voltage magnitude of the current iteration.

		Attributes:
			nominal_voltage: The rated voltage of the load. Type=float.
			power_A: Initial phase A power of the load. Type=Complex
			power_B: Initial phase B power of the load. Type=Complex
			power_C: Initial phase C power of the load. Type=Complex
			power_AB: Initial phase AB power of the load. Type=Complex
			power_BC: Initial phase BC power of the load. Type=Complex
			power_CA: Initial phase CA power of the load. Type=Complex
			CVRwatts: defines the relationship between voltage and active power. Type=list.
			CVRvars: defines the relationship between voltage and reactive power. Type=list.
		"""

        super(Exponential, self).__init__(ID, name, parent, phases,
                                          is_delta, load_class)

        if self.isDelta:
            self.P0 = [power_AB.real, power_BC.real, power_CA.real]
            self.Q0 = [power_AB.imag, power_BC.imag, power_CA.imag]
        else:
            self.P0 = [power_A.real, power_B.real, power_C.real]
            self.Q0 = [power_A.imag, power_B.imag, power_C.imag]
        self.CVRwatts = CVRwatts
        self.CVRvars = CVRvars
        self.V0 = nominal_voltage

    def calc_partials(self, Vr, Vi, P0, Q0, CVRwatts, CVRvars, lamR, lamI):
        """ Calculate partial derivatives for the exponential load.

		Args:
		  Vr: The real component of the exponential load voltage. Type float.
		  Vi: The imaginary component of the exponential load voltage. Type float.
		  P0: Initial real power of the load. Type float.
		  Q0: Initial reactive power of the load. Type float.
		  CVRwatts: defines the relationship between voltage and active power. Type=float.
		  CVRvars: defines the relationship between voltage and reactive power. Type=float.
		Returns:
			Historical currents and partial derivatives for the exponential load with respect to the load's real and
			imaginary components.
		"""

        V0 = self.V0

        QVrVi = Q0 * Vr * Vi * V0**CVRvars
        PVrVi = P0 * Vr * Vi * V0**CVRwatts
        Vi2Vr2 = Vi**2 - Vr**2
        Vmag2 = Vr**2 + Vi**2
        Vmag4 = Vmag2 * Vmag2
        V0cvr = V0**(-CVRvars - CVRwatts)
        Ir_hist = (((P0 * Vr * V0**CVRwatts) +
                    (Q0 * Vi * V0**CVRvars)) * V0cvr) / Vmag2
        Ii_hist = (((P0 * Vi * V0**CVRwatts) -
                    (Q0 * Vr * V0**CVRvars)) * V0cvr) / Vmag2

        dIrdVr = (
            (P0 *
             (V0**CVRwatts) * Vi2Vr2 - 2 * QVrVi) * V0cvr) / Vmag4  # dIr/dVr
        dIrdVi = (
            (-Q0 *
             (V0**CVRvars) * Vi2Vr2 - 2 * PVrVi) * V0cvr) / Vmag4  # dIr/dVi
        dIidVi = (
            (-P0 *
             (V0**CVRwatts) * Vi2Vr2 + 2.0 * QVrVi) * V0cvr) / Vmag4  # dIi/dVi
        dIidVr = dIrdVi  # dIr/dVi
    
            
        return (Ir_hist, Ii_hist, dIrdVr, dIrdVi, dIidVi, dIidVr)

    def stamp_phase_load(self, node_Vr_from, node_Vi_from, node_Vr_to,
                         node_Vi_to, Vr_from, Vi_from, Vr_to, Vi_to, P0, Q0,
                         CVRwatts, CVRvars, Ynlin_val, Ynlin_row, Ynlin_col,
                         Jnlin_val, Jnlin_row, idx_Y, idx_J, lamR, lamI, stamp_dual):
        """"Stamp the phase load circuits.

			Stamps the from and to phase load circuits. The to circuit is stamped only if the load is delta connected.

			Args:
				node_Vr_from: The real from voltage node. Type int.
				node_Vi_from: The imaginary from voltage node. Type int.
				node_Vr_to: The real to voltage node. Type int.
				node_Vi_to: The imaginary to voltage node. Type int.
				Vr_from: The real from voltage of a particular phase. Type float.
				Vi_from: The imaginary from voltage of a particular phase. Type float.
				Vr_to: The real to voltage of a particular phase. Type float.
				Vi_to: The imaginary to voltage of a particular phase. Type float.
				P0: Per-phase initial real power.
				Q0: Per-phase initial imaginary power.
				CVRwatts: defines the relationship between voltage and active power. Type=float.
				CVRvars: defines the relationship between voltage and reactive power. Type=float.
				Ynlin_val: Stores the admittance matrix value that is stamped. Type float.
				Ynlin_row: Stores the admittance matrix row number where the stamp is made. Type int.
				Ynlin_col: Stores the admittance matrix column number where the stamp is made. Type int.
				Jnlin_val: Stores the excitation vector value that is stamped. Type float.
				Jnlin_row: Stores the excitation vector row number where the stamp is made. Type int.
				idx_Y: Admittance matrix stamp index. Type int.
				idx_J: Excitation vector stamp index. Type int.

			Returns:
				Excitation vector stamp index and admittance matrix stamp index.

		"""

        if not self.isGnd:
            # calculate line to line voltage
            Vr = Vr_from - Vr_to
            Vi = Vi_from - Vi_to

            # compute partials and historical currents
            (Ir_hist, Ii_hist, dIr_dVr, dIr_dVi, dIi_dVr, dIi_dVi, lamR_hist, lamI_hist,
                dlamR_dVr, dlamR_dVi, dlamR_dlamR, dlamR_dlamI, dlamI_dVr, dlamI_dVi,
                dlamI_dlamR, dlamI_dlamI)  = \
                    self.calc_partials(Vr, Vi, P0, Q0, CVRwatts, CVRvars, lamR, lamI, stamp_dual)

            # computer independent current source currnets
            Ir = Ir_hist - dIr_dVr * Vr - dIr_dVi * Vi
            Ii = Ii_hist - dIi_dVi * Vi - dIi_dVr * Vr

            # #  From Circuit Stamps # #
            # - Real Stamps - #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from, dIr_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_to, -dIr_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from, dIr_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_to, -dIr_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Imaginary Stamps - #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from, dIi_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_to, -dIi_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from, dIi_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_to, -dIi_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Independent Current Source Stamps - #
            idx_J = self.stamp_J(node_Vr_from, -Ir, Jnlin_val, Jnlin_row, idx_J)

            idx_J = self.stamp_J(node_Vi_from, -Ii, Jnlin_val, Jnlin_row, idx_J)

            # # To Circuit Stamps # #
            # - Real Stamps - #
            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_from, -dIr_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_to, node_Vr_to, dIr_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_from, -dIr_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_to, node_Vi_to, dIr_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Imaginary Stamps - #
            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_from, -dIi_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_to, node_Vi_to, dIi_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_from, dIi_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_to, node_Vr_to, dIi_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Independent Current Source Stamps - #
            idx_J = self.stamp_J(node_Vr_to, -Ir, Jnlin_val, Jnlin_row, idx_J)

            idx_J = self.stamp_J(node_Vi_to, -Ii, Jnlin_val, Jnlin_row, idx_J)

        else:
            # calculate phase voltages
            Vr = Vr_from
            Vi = Vi_from
            

            # compute partials and historical currents
            (Ir_hist, Ii_hist, dIr_dVr, dIr_dVi, dIi_dVr, dIi_dVi, lamR_hist, lamI_hist,
                dlamR_dVr, dlamR_dVi, dlamR_dlamR, dlamR_dlamI, dlamI_dVr, dlamI_dVi,
                dlamI_dlamR, dlamI_dlamI)  = \
                    self.calc_partials(Vr, Vi, P0, Q0, CVRwatts, CVRvars, lamR, lamI, stamp_dual) 

            # compute independent current source currnets
            Ir = Ir_hist - dIr_dVr * Vr - dIr_dVi * Vi
            Ii = Ii_hist - dIi_dVi * Vi - dIi_dVr * Vr

            # #  From Circuit Stamps # #
            # - Real Stamps - #
            idx_Y = self.stamp_Y(node_Vr_from, node_Vr_from, dIr_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vr_from, node_Vi_from, dIr_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Imaginary Stamps - #
            idx_Y = self.stamp_Y(node_Vi_from, node_Vi_from, dIi_dVi, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            idx_Y = self.stamp_Y(node_Vi_from, node_Vr_from, dIi_dVr, Ynlin_val, Ynlin_row, Ynlin_col, idx_Y)

            # - Independent Current Source Stamps - #
            idx_J = self.stamp_J(node_Vr_from, -Ir, Jnlin_val, Jnlin_row, idx_J)

            idx_J = self.stamp_J(node_Vi_from, -Ii, Jnlin_val, Jnlin_row, idx_J)    

        return idx_Y, idx_J

    def stamp_nonlinear(self,
                        node_key,
                        node,
                        V,
                        lf,
                        Ynlin_val,
                        Ynlin_row,
                        Ynlin_col,
                        Jnlin_val,
                        Jnlin_row,
                        idx_Y,
                        idx_J,
                        stamp_dual,
                        enable_CM=False):
        """" Stamp the exponential load partials.

        Args:
            Ynlin_val: Stores the admittance matrix value that is stamped. Type float.
            Ynlin_row: Stores the admittance matrix row number where the stamp is made. Type int.
            Ynlin_col: Stores the admittance matrix column number where the stamp is made. Type int.
            Jnlin_val: Stores the excitation vector value that is stamped. Type float.
            Jnlin_row: Stores the excitation vector row number where the stamp is made. Type int.
            node: The vector containing a list of all nodes.
            idx_Y: Admittance matrix stamp index. Type int.
            idx_J: Excitation vector stamp index. Type int.
            V: The solution vector.
            lf: Load factor. Type float.
            enable_CM: Turn current measures on/off. Disabled for now. Type bool.

        Returns:
            Excitation vector stamp index and admittance matrix stamp index.

        """
        # Get nodes and voltages to use
        V_from, V_to = self.get_nodes(node_key, node, V)

        for _phase in range(3):

            # collect from and to nodes of each phase
            node_Vr_from = V_from.nodeVR[_phase]
            node_Vi_from = V_from.nodeVI[_phase]
            node_Vr_to = V_to.nodeVR[_phase]
            node_Vi_to = V_to.nodeVI[_phase]

            # collect from and to voltage values of each phase
            Vr_from = V_from.VR[_phase]
            Vi_from = V_from.VI[_phase]
            Vr_to = V_to.VR[_phase]
            Vi_to = V_to.VI[_phase]

            # stamp the phase load if that phase is present
            if self.load_phases[_phase] & self.phases == int(
                    self.load_phases[_phase]):
                # scale initial loading by load factor or power stepping factor if necessary
                P0 = self.P0[_phase] * lf
                Q0 = self.Q0[_phase] * lf
                # retrieve phase CVR value
                CVRwatts = self.CVRwatts[_phase]
                CVRvars = self.CVRvars[_phase]

                self.stamp_phase_load(node_Vr_from, node_Vi_from, node_Vr_to,
                                      node_Vi_to, Vr_from, Vi_from, Vr_to,
                                      Vi_to, P0, Q0, CVRwatts, CVRvars,
                                      Ynlin_val, Ynlin_row, Ynlin_col,
                                      Jnlin_val, Jnlin_row, idx_Y, idx_J)

        return idx_Y, idx_J       
        

class ZIP(Load):

    def __init__(self,
                 ID,
                 name,
                 parent,
                 phases,
                 is_delta=False,
                 load_class=0,
                 base_powerA=0,
                 base_powerB=0,
                 base_powerC=0,
                 percent_Z=None,
                 percent_I=None,
                 percent_P=None):
        """
		Initialize the ZIP load class.
		:param ID: The load ID.
		:param name: The name of the load.
		:param parent: The bus (node) that the load is located on.
		:param phases: Indicates the number of phases present.
		:param load_class: Specifies whether load is of type 0 - Unknown, 1-Residential, 2-Commercial,
		3-Industrial, or 4-Agricultural.
		:param voltageA: The voltage on phase A of the load.
		:param voltageB: The voltage on phase B of the load.
		:param voltageC: The voltage on phase C of the load.
		:param voltageAB: The voltage on the delta-phase AB of the load.
		:param voltageBC: The voltage on the delta-phase BC of the load.
		:param voltageCA: The voltage on the delta-phase CA of the load.
		:param base_powerA: Nominal power on phase A before ZIP fractionals.
		:param base_powerB: Nominal power on phase B before ZIP fractionals.
		:param base_powerC: Nominal power on phase C before ZIP fractionals.
		:param power_pfA: Load power factor of the phase A constant power load portion.
		:param current_pfA: Load power factor of the phase A constant current load portion.
		:param impedance_pfA: Load power factor of the phase A constant impedance load portion.
		:param power_pfB:  Load power factor of the phase B constant power load portion.
		:param current_pfB: Load power factor of the phase B constant current load portion.
		:param impedance_pfB: Load power factor of the phase B constant impedance load portion.
		:param power_pfC: Load power factor of the phase C constant power load portion.
		:param current_pfC: Load power factor of the phase C constant current load portion.
		:param impedance_pfC: Load power factor of the phase C constant impedance load portion.
		:param power_fracA: Phase A base constant power fraction.
		:param current_fracA: Phase A base constant current fraction.
		:param impedance_fracA: Phase A base constant impedance fraction.
		:param power_fracB: Phase B base constant power fraction.
		:param current_fracB: Phase B base constant current fraction.
		:param impedance_fracB: Phase B base constant impedance fraction.
		:param power_fracC: Phase C base constant power fraction.
		:param current_fracC: Phase C base constant current fraction.
		:param impedance_fracC: Phase C base constant impedance fraction.
		"""

        super(ZIP, self).__init__(ID, name, parent, phases,
                                  is_delta, load_class)

        self.Type = "ZIP"

        self.base_powerA = float(base_powerA.real) if base_powerA else 0  # [VA]
        self.base_powerB = float(base_powerB.real) if base_powerB else 0  # [VA]
        self.base_powerC = float(base_powerC.real) if base_powerC else 0  # [VA]

        _len_Percent_Z = len(percent_Z)
        _len_Percent_I = len(percent_I)
        _len_Percent_P = len(percent_P)

        self.Z_pA = float(percent_Z[0][0]) if _len_Percent_Z >= 1 else 0
        self.Z_pB = float(percent_Z[1][0]) if _len_Percent_Z >= 2 else 0
        self.Z_pC = float(percent_Z[2][0]) if _len_Percent_Z == 3 else 0
        self.Z_qA = float(percent_Z[0][1]) if _len_Percent_Z >= 1 else 0
        self.Z_qB = float(percent_Z[1][1]) if _len_Percent_Z >= 2 else 0
        self.Z_qC = float(percent_Z[2][1]) if _len_Percent_Z == 3 else 0

        self.I_pA = float(percent_I[0][0]) if _len_Percent_I >= 1 else 0
        self.I_pB = float(percent_I[1][0]) if _len_Percent_I >= 2 else 0
        self.I_pC = float(percent_I[2][0]) if _len_Percent_I == 3 else 0
        self.I_qA = float(percent_I[0][1]) if _len_Percent_I >= 1 else 0
        self.I_qB = float(percent_I[1][1]) if _len_Percent_I >= 2 else 0
        self.I_qC = float(percent_I[2][1]) if _len_Percent_I == 3 else 0

        self.P_pA = float(percent_P[0][0]) if _len_Percent_P >= 1 else 0
        self.P_pB = float(percent_P[1][0]) if _len_Percent_P >= 2 else 0
        self.P_pC = float(percent_P[2][0]) if _len_Percent_P == 3 else 0
        self.P_qA = float(percent_P[0][1]) if _len_Percent_P >= 1 else 0
        self.P_qB = float(percent_P[1][1]) if _len_Percent_P >= 2 else 0
        self.P_qC = float(percent_P[2][1]) if _len_Percent_P == 3 else 0

        self._base_power = np.array(
            [self.base_powerA, self.base_powerB, self.base_powerC])
        self._Zr_perc = np.array([self.Z_pA, self.Z_pB, self.Z_pC])
        self._Ir_perc = np.array([self.I_pA, self.I_pB, self.I_pC])
        self._Pr_perc = np.array([self.P_pA, self.P_pB, self.P_pC])
        self._Zi_perc = np.array([self.Z_qA, self.Z_qB, self.Z_qC])
        self._Ii_perc = np.array([self.I_qA, self.I_qB, self.I_qC])
        self._Pi_perc = np.array([self.P_qA, self.P_qB, self.P_qC])

    def stamp_nonlinear(self,
                        node_key,
                        node,
                        V,
                        lf,
                        Ynlin_val,
                        Ynlin_row,
                        Ynlin_col,
                        Jnlin_val,
                        Jnlin_row,
                        idx_Y,
                        idx_J):
        """
		Stamp the nonlinear linearized ZIP load partials.

		:param Ynlin_val: Admittance matrix value to be stamped.
		:param Ynlin_row: Admittance matrix row.
		:param Ynlin_col: Admittance matrix column.
		:param Jnlin_val: Excitation vector partial derivative value to be stamp.
		:param Jnlin_row: Excitation vector row.
		:param node: The vector containing a list of all nodes.
		:param V: The solution vector.
		:param lf: Load factor.
		:return:
		"""
        V_from, V_to = super(ZIP, self).get_nodes(node_key, node, V)

        _base_power = self._base_power * lf
        _Zr_perc = self._Zr_perc
        _Ir_perc = self._Ir_perc
        _Pr_perc = self._Pr_perc
        _Zi_perc = self._Zi_perc
        _Ii_perc = self._Ii_perc
        _Pi_perc = self._Pi_perc

        for _phase_idx in range(3):

            if self.load_phases[_phase_idx] & self.phases == int(
                    self.load_phases[_phase_idx]):
                _cP = _base_power[_phase_idx] * _Pr_perc[_phase_idx]
                _cQ = _base_power[_phase_idx] * _Pi_perc[_phase_idx]
                _cG = _base_power[_phase_idx] * _Zr_perc[_phase_idx]
                _cB = _base_power[_phase_idx] * _Zi_perc[_phase_idx]
                _cIp = _base_power[_phase_idx] * _Ir_perc[_phase_idx]
                _cIq = _base_power[_phase_idx] * _Ii_perc[_phase_idx]

                node_Vr_from = V_from.nodeVR[_phase_idx]
                node_Vi_from = V_from.nodeVI[_phase_idx]
                node_Vr_to = V_to.nodeVR[_phase_idx]
                node_Vi_to = V_to.nodeVI[_phase_idx]

                Vr_from = V_from.VR[_phase_idx]
                Vi_from = V_from.VI[_phase_idx]
                Vr_to = V_to.VR[_phase_idx]
                Vi_to = V_to.VI[_phase_idx]

                if _cP or _cQ or _cG or _cB or _cIp or _cIq:
                    idx_Y, idx_J = self._load_calculations(
                        _cP, _cQ, _cG, _cB, _cIp, _cIq, node_Vr_from,
                        node_Vi_from, node_Vr_to, node_Vi_to, Vr_from, Vi_from,
                        Vr_to, Vi_to, Ynlin_val, Ynlin_row, Ynlin_col,
                        Jnlin_val, Jnlin_row, idx_Y, idx_J)

        return idx_Y, idx_J
