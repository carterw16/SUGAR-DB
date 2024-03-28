"""The class for lines for Three-phase power flow
File Created: 04/11/2017 by Amrit Pandey
File Amended: 11/28/2022 by Elizabeth Foster
Changes: Reintroduce older line calculations and implement file structure"""

from __future__ import division

import logging
from itertools import count
from types import MethodType, SimpleNamespace
import numpy as np

from classes.Elements import LinearElement
from classes.GlobalVars import _PHASE

from types import MethodType

from .LinesStamp import (stamp_linear, stamp_phase_ckt, stamp_tx_stepping,
                         stamp_neutral, stamp_linear_triplex)

# A = 1 | B = 2 | C = 4 | D = 7 | N - 8 | G - 16 | S - 32


class Lines(LinearElement):
    impedance_factor = 1
    _ids = count(0)
    
    @staticmethod
    def calc_Zii(r_i, GMRi):
        try:
            # return Zii in ohms/mile
            Zii = r_i + 0.09530 + 1j * 0.12134 * (np.math.log(GMRi ** -1) + 7.93402)
        except TypeError:
            Zii = -1000
        return Zii

    @staticmethod
    def calc_Zij(Dij):
        # return Zij in ohms/mile
        Zij = 0.09530 + 1j * 0.12134 * (np.math.log(Dij ** -1) + 7.93402)
        return Zij
    
    @staticmethod
    def calc_underground_Zii(r_i, GMRi):
        # returns Zii in ohms/mile
        Zii = r_i + 0.09327 + 1j * 0.12134 * (np.math.log(GMRi**-1) + 7.95153)
        return Zii

    @staticmethod
    def calc_underground_Zij(Dij):
        # returns Zij in ohms/mile
        Zij = 0.09327 + 1j * 0.12134 * (np.math.log(Dij**-1) + 7.95153)
        return Zij
    
    @staticmethod
    def kron_reduction(Zij, Zin, Znj, Znn):
        #  _r = Zij.shape[0]
        #  _c = Zij.shape[1]
        if hasattr(Znn, "__len__"):
            Znn_inv = np.linalg.inv(Znn)
        else:
            Znn_inv = Znn**-1
        _tempZ = np.dot(np.dot(Zin, Znn_inv), Znj)
        Zabc = Zij - _tempZ
        return Zabc

    @staticmethod
    def setNonetoZero(x):
        if isinstance(x, type(None)):
            x = 0
            return x
        else:
            return x

    @staticmethod
    def removeRowColZero(Zmatrix):
        # Convert Matrix to Ymatrix
        # Find the rows of zero
        _rowZeroIdx = np.all(Zmatrix == 0, axis=1)
        Zmatrix_red = Zmatrix[~_rowZeroIdx]
        # Find the cols of zero
        _colZeroIdx = np.all(Zmatrix_red == 0, axis=0)
        Zmatrix_red = Zmatrix_red[:, ~_colZeroIdx]
        return Zmatrix_red, _rowZeroIdx, _colZeroIdx

    def __init__(self,
                 ID,
                 name,
                 type,
                 phases,
                 nominal_voltage,
                 from_node,
                 to_node,
                 length,
                 impedance_matrix,
                 capacitance_matrix,
                 spacing=None,
                 unused_phases=None,
                 conductor_phases=None,
                 stamp_dual = False):

        super(Lines, self).__init__()
        self.stamp_linear = MethodType(stamp_linear, self)
        self.stamp_neutral = MethodType(stamp_neutral, self)
        self.stamp_linear_triplex = MethodType(stamp_linear_triplex, self)
        self.stamp_phase_ckt = MethodType(stamp_phase_ckt, self)
        self.stamp_tx_stepping = MethodType(stamp_tx_stepping, self)
        
        self.name = name
        self.phases = [
            key for key, val in _PHASE.items() if val == int(phases)
        ].pop()
        self._phases = _PHASE[self.phases]
        self.type = type
        self.from_node = from_node
        self.to_node = to_node
        self.id = self._ids.__next__()
        self.ID = ID
        self.freq = 60
        self.length = float(length)  # [miles]
        self.nominal_voltage = float(
            nominal_voltage) if nominal_voltage else None
        self._Zmatrix = impedance_matrix
        self._Cmatrix = capacitance_matrix
        self.conductor_phases = conductor_phases if conductor_phases else None
        self._conductor_phases = _PHASE[
            self.conductor_phases] if conductor_phases else self._phases
        self.isTriplex = True if self.type == "triplex" else False
        self.spacing = spacing
        self.unused_phases = unused_phases
        self.hasShunt = False
        self.Zmatrix = None
        self.Ymatrix = np.zeros((3, 3), dtype=complex)
        self.Yshunt = np.zeros((3, 3), dtype=complex)
        self.neutral_exists = False
        self.stamp_dual = stamp_dual

        # # Link Cython Methods to Class # #
        self.stamp_linear = MethodType(stamp_linear, self)
        self.stamp_linear_triplex = MethodType(stamp_linear_triplex, self)
        self.stamp_phase_ckt = MethodType(stamp_phase_ckt, self)
        self.stamp_tx_stepping = MethodType(stamp_tx_stepping, self)

        # Initialize line spacings
        self.D_AN = 0.0
        self.D_BN = 0.0
        self.D_CN = 0.0
        self.D_NE = 0.0
        self.D_AB = 0.0
        self.D_AC = 0.0
        self.D_BC = 0.0
        self.D_AE = 0.0
        self.D_BE = 0.0
        self.D_CE = 0.0
        self.init_spacings()

        # Initialize line current magnitudes and angles
        self.Ia_mag = 0
        self.Ia_ang = 0
        self.Ib_mag = 0
        self.Ib_ang = 0
        self.Ic_mag = 0
        self.Ic_ang = 0

    def augment_matrices(self):
        _complexZero = np.zeros((1, 1), dtype=complex)
        for ele in self.unused_phases:
            if ele == 'A':
                self._Zmatrix = np.insert(self._Zmatrix,
                                          0,
                                          _complexZero,
                                          axis=0)
                self._Zmatrix = np.insert(self._Zmatrix,
                                          0,
                                          _complexZero,
                                          axis=1)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          0,
                                          _complexZero,
                                          axis=0)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          0,
                                          _complexZero,
                                          axis=1)
            elif ele == 'B':
                self._Zmatrix = np.insert(self._Zmatrix,
                                          1,
                                          _complexZero,
                                          axis=0)
                self._Zmatrix = np.insert(self._Zmatrix,
                                          1,
                                          _complexZero,
                                          axis=1)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          1,
                                          _complexZero,
                                          axis=0)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          1,
                                          _complexZero,
                                          axis=1)
            elif ele == 'C':
                self._Zmatrix = np.insert(self._Zmatrix,
                                          2,
                                          _complexZero,
                                          axis=0)
                self._Zmatrix = np.insert(self._Zmatrix,
                                          2,
                                          _complexZero,
                                          axis=1)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          2,
                                          _complexZero,
                                          axis=0)
                self._Cmatrix = np.insert(self._Cmatrix,
                                          2,
                                          _complexZero,
                                          axis=1)

    def calcInverse(self):
        num_nnz = np.count_nonzero(self.Zmatrix)
        if num_nnz == 1:
            self.Zmatrix[
                self.Zmatrix == 0 +
                0j] = np.inf  # Division in the next step will result in 0 + 0j in Y
            Ymatrix = np.divide(np.ones(np.shape(self.Zmatrix), dtype=complex),
                                self.Zmatrix)
            return Ymatrix
        # Convert Matrix to Ymatrix
        # Find the rows of zero
        _rowZeroIdx = np.all(self.Zmatrix == 0, axis=1)
        Zmatrix_red = self.Zmatrix[~_rowZeroIdx]
        # Find the cols of zero
        _colZeroIdx = np.all(Zmatrix_red == 0, axis=0)
        Zmatrix_red = Zmatrix_red[:, ~_colZeroIdx]
        # Find the inverse of reduced matrix
        Ymatrix_red = np.linalg.inv(Zmatrix_red)
        # Remap to original format
        if len(Ymatrix_red) == 2:
            _rowIdx = list(_rowZeroIdx).index(True)
            _colIdx = list(_colZeroIdx).index(True)
            _complexZero = np.zeros((1, 1), dtype=complex)
            Ymatrix = np.insert(Ymatrix_red, _rowIdx, _complexZero, axis=0)
            Ymatrix = np.insert(Ymatrix, _colIdx, _complexZero, axis=1)
        else:
            Ymatrix = Ymatrix_red
        return Ymatrix
    
    def calcInverse_DSS(self):
        num_nnz = np.count_nonzero(self.Zmatrix)
        if num_nnz == 1:
            self.Zmatrix[
                self.Zmatrix == 0 +
                0j] = np.inf  # Division in the next step will result in 0 + 0j in Y
            Ymatrix = np.divide(np.ones(np.shape(self.Zmatrix), dtype=complex),
                                self.Zmatrix)
            return Ymatrix
        # Convert Matrix to Ymatrix
        # Find the rows of zero
        _rowZeroIdx = np.all(self.Zmatrix == 0, axis=1)
        Zmatrix_red = self.Zmatrix[~_rowZeroIdx]
        # Find the cols of zero
        _colZeroIdx = np.all(Zmatrix_red == 0, axis=0)
        Zmatrix_red = Zmatrix_red[:, ~_colZeroIdx]
        # Find the inverse of reduced matrix
        Ymatrix_red = np.linalg.inv(Zmatrix_red)
        # Remap to original format
        if len(Ymatrix_red) < len(self.phases):
            _rowIdx = list(_rowZeroIdx).index(True)
            _colIdx = list(_colZeroIdx).index(True)
            _complexZero = np.zeros((1, 1), dtype=complex)
            Ymatrix = np.insert(Ymatrix_red, _rowIdx, _complexZero, axis=0)
            Ymatrix = np.insert(Ymatrix, _colIdx, _complexZero, axis=1)
        else:
            Ymatrix = Ymatrix_red
        return Ymatrix
    
    def init_spacings(self):
        zero = float(0)
        neg_one = float(-1)
        default = float(4)
        if hasattr(self.spacing, '_distance_AN'):
            if float(getattr(self.spacing, '_distance_AN')) == zero:
                logging.debug(
                    'Spacing AN for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_AN = default
            else:
                self.D_AN = float(self.spacing._distance_AN) 
        else:
            self.D_AN = neg_one

        if hasattr(self.spacing, '_distance_BN'):
            if float(getattr(self.spacing, '_distance_BN')) == zero:
                logging.debug(
                    'Spacing BN for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_BN = default
            else:
                self.D_BN = float(self.spacing._distance_BN) 
        else:
            self.D_BN = neg_one

        if hasattr(self.spacing, '_distance_CN'):
            if float(getattr(self.spacing, '_distance_CN')) == zero:
                logging.debug(
                    'Spacing CN for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_CN = default
            else:
                self.D_CN = float(self.spacing._distance_CN) 
        else:
            self.D_CN = neg_one
        if hasattr(self.spacing, '_distance_AB'):
            if float(getattr(self.spacing, '_distance_AB')) == zero:
                logging.debug(
                    'Spacing AB for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_AB = default
            else:
                self.D_AB = float(self.spacing._distance_AB)
        else:
            self.D_AB = neg_one
        if hasattr(self.spacing, '_distance_AC'):
            if float(getattr(self.spacing, '_distance_AC')) == zero:
                logging.debug(
                    'Spacing AC for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_AC = default
            else:
                self.D_AC = float(self.spacing._distance_AC)
        else:
            self.D_AC = neg_one

        if hasattr(self.spacing, '_distance_BC'):
            if float(getattr(self.spacing, '_distance_BC')) == zero:
                logging.debug(
                    'Spacing AC for %d is given as 0. Please Correct. Using Default'
                    % self.id)
                self.D_BC = default
            else:
                self.D_BC = float(self.spacing._distance_BC) 
        else:
            self.D_BC = neg_one

        self.D_NE = float(self.spacing._distance_NE) if hasattr(
            self.spacing, '_distance_NE') else float(25)
        self.D_AE = float(self.spacing._distance_AE) if hasattr(
            self.spacing, '_distance_AE') else float(25)
        self.D_BE = float(self.spacing._distance_BE) if hasattr(
            self.spacing, '_distance_BE') else float(25)
        self.D_CE = float(self.spacing._distance_CE) if hasattr(
            self.spacing, '_distance_CE') else float(25)

    def get_nodes(self, node, node_key):
        """
		Find the indices of the solution vector and the values at those indices.
		:param node: The vector of all system nodes.
		:param V: The solution vector.
		:return: Two dictionaries which hold the from/to phase nodes for the transmission lines.
		"""

        if self.isTriplex:
            # Find the from bus nodes to stamp
            node1_Vr_from = node[node_key[self.from_node]].node1_Vr
            node1_Vi_from = node[node_key[self.from_node]].node1_Vi
            node2_Vr_from = node[node_key[self.from_node]].node2_Vr
            node2_Vi_from = node[node_key[self.from_node]].node2_Vi
            nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
            nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi

            # Find the to bus nodes to stamp
            node1_Vr_to = node[node_key[self.to_node]].node1_Vr
            node1_Vi_to = node[node_key[self.to_node]].node1_Vi
            node2_Vr_to = node[node_key[self.to_node]].node2_Vr
            node2_Vi_to = node[node_key[self.to_node]].node2_Vi
            nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
            nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi
            
            if self.stamp_dual:
                node1_Lr_from = node[node_key[self.from_node]].node1_dual_eq_var_r
                node1_Li_from = node[node_key[self.from_node]].node1_dual_eq_var_i
                node2_Lr_from = node[node_key[self.from_node]].node2_dual_eq_var_r
                node2_Li_from = node[node_key[self.from_node]].node2_dual_eq_var_i
                nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
                nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i
                
                node1_Lr_to = node[node_key[self.to_node]].node1_dual_eq_var_r
                node1_Li_to = node[node_key[self.to_node]].node1_dual_eq_var_i
                node2_Lr_to = node[node_key[self.to_node]].node2_dual_eq_var_r
                node2_Li_to = node[node_key[self.to_node]].node2_dual_eq_var_i
                nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
                nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i
                
                nodes_from = {
                    'VR': [node1_Vr_from, node2_Vr_from, nodeN_Vr_from],
                    'VI': [node1_Vi_from, node2_Vi_from, nodeN_Vi_from],
                    'LR': [node1_Lr_from, node2_Lr_from, nodeN_Lr_from],
                    'LI': [node1_Li_from, node2_Li_from, nodeN_Li_from]
                }
                nodes_to = {
                    'VR': [node1_Vr_to, node2_Vr_to, nodeN_Vr_to],
                    'VI': [node1_Vi_to, node2_Vi_to, nodeN_Vi_to],
                    'LR': [node1_Lr_to, node2_Lr_to, nodeN_Lr_to],
                    'LI': [node1_Li_to, node2_Li_to, nodeN_Li_to]
                }
                
            else:
                nodes_from = {
                    'VR': [node1_Vr_from, node2_Vr_from, nodeN_Vr_from],
                    'VI': [node1_Vi_from, node2_Vi_from, nodeN_Vi_from]
                }
                nodes_to = {
                    'VR': [node1_Vr_to, node2_Vr_to, nodeN_Vr_to],
                    'VI': [node1_Vi_to, node2_Vi_to, nodeN_Vi_to]
                }
            
            nodes_from = SimpleNamespace(**nodes_from)
            nodes_to = SimpleNamespace(**nodes_to)

        else:
            if self.line_parameters['file type'] == 'DSS':
                if self.stamp_dual:
                    nodes_from = {
                        'VR': [],
                        'VI': [],
                        'LR': [],
                        'LI': []
                    }
        
                    nodes_to = {
                        'VR': [],
                        'VI': [],
                        'LR': [],
                        'LI': []
                    }
                else:
                    nodes_from = {
                        'VR': [],
                        'VI': []
                    }
        
                    nodes_to = {
                        'VR': [],
                        'VI': []
                    } 

                if 'A' in self.phases:
                    nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr 
                    nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi 

                    nodes_from['VR'].append(nodeA_Vr_from)
                    nodes_from['VI'].append(nodeA_Vi_from)
                    
                    nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
                    nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi

                    nodes_to['VR'].append(nodeA_Vr_to)
                    nodes_to['VI'].append(nodeA_Vi_to)
                    if self.stamp_dual:
                        nodeA_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
                        nodeA_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i

                        nodes_from['LR'].append(nodeA_Lr_from)
                        nodes_from['LI'].append(nodeA_Li_from)

                        nodeA_Lr_to = node[node_key[self.to_node]].nodeA_dual_eq_var_r
                        nodeA_Li_to = node[node_key[self.to_node]].nodeA_dual_eq_var_i

                        nodes_to['LR'].append(nodeA_Lr_to)
                        nodes_to['LI'].append(nodeA_Li_to)
                
                if 'B' in self.phases:
                    nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr 
                    nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi 

                    nodes_from['VR'].append(nodeB_Vr_from)
                    nodes_from['VI'].append(nodeB_Vi_from)

                    nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
                    nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi

                    nodes_to['VR'].append(nodeB_Vr_to)
                    nodes_to['VI'].append(nodeB_Vi_to)
                    if self.stamp_dual:
                        nodeB_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
                        nodeB_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
                        
                        nodes_from['LR'].append(nodeB_Lr_from)
                        nodes_from['LI'].append(nodeB_Li_from)

                        nodeB_Lr_to = node[node_key[self.to_node]].nodeB_dual_eq_var_r
                        nodeB_Li_to = node[node_key[self.to_node]].nodeB_dual_eq_var_i

                        nodes_to['LR'].append(nodeB_Lr_to)
                        nodes_to['LI'].append(nodeB_Li_to)
                
                if 'C' in self.phases:
                    nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
                    nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi

                    nodes_from['VR'].append(nodeC_Vr_from)
                    nodes_from['VI'].append(nodeC_Vi_from)

                    nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
                    nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi

                    nodes_to['VR'].append(nodeC_Vr_to)
                    nodes_to['VI'].append(nodeC_Vi_to)
                    if self.stamp_dual:
                        nodeC_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
                        nodeC_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i

                        nodes_from['LR'].append(nodeC_Lr_from)
                        nodes_from['LI'].append(nodeC_Li_from)

                        nodeC_Lr_to = node[node_key[self.to_node]].nodeC_dual_eq_var_r
                        nodeC_Li_to = node[node_key[self.to_node]].nodeC_dual_eq_var_i

                        nodes_to['LR'].append(nodeC_Lr_to)
                        nodes_to['LI'].append(nodeC_Li_to)
                
                nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
                nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi 
                
                nodes_from['VR'].append(nodeN_Vr_from)
                nodes_from['VI'].append(nodeN_Vi_from)

                nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
                nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi

                nodes_to['VR'].append(nodeN_Vr_to)
                nodes_to['VI'].append(nodeN_Vi_to)
                if self.stamp_dual:
                    nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
                    nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i

                    nodes_from['LR'].append(nodeN_Lr_from)
                    nodes_from['LI'].append(nodeN_Li_from)
                    
                    nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
                    nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i

                    nodes_to['LR'].append(nodeN_Lr_to)
                    nodes_to['LI'].append(nodeN_Li_to)

            else:
                nodeA_Vr_from = node[node_key[self.from_node]].nodeA_Vr 
                nodeA_Vi_from = node[node_key[self.from_node]].nodeA_Vi 
                nodeB_Vr_from = node[node_key[self.from_node]].nodeB_Vr 
                nodeB_Vi_from = node[node_key[self.from_node]].nodeB_Vi 
                nodeC_Vr_from = node[node_key[self.from_node]].nodeC_Vr
                nodeC_Vi_from = node[node_key[self.from_node]].nodeC_Vi

                if self.stamp_dual:
                    nodeA_Lr_from = node[node_key[self.from_node]].nodeA_dual_eq_var_r
                    nodeA_Li_from = node[node_key[self.from_node]].nodeA_dual_eq_var_i
                    nodeB_Lr_from = node[node_key[self.from_node]].nodeB_dual_eq_var_r
                    nodeB_Li_from = node[node_key[self.from_node]].nodeB_dual_eq_var_i
                    nodeC_Lr_from = node[node_key[self.from_node]].nodeC_dual_eq_var_r
                    nodeC_Li_from = node[node_key[self.from_node]].nodeC_dual_eq_var_i
                    nodeN_Lr_from = node[node_key[self.from_node]].nodeN_dual_eq_var_r
                    nodeN_Li_from = node[node_key[self.from_node]].nodeN_dual_eq_var_i
                    
                    nodeA_Lr_to = node[node_key[self.to_node]].nodeA_dual_eq_var_r
                    nodeA_Li_to = node[node_key[self.to_node]].nodeA_dual_eq_var_i
                    nodeB_Lr_to = node[node_key[self.to_node]].nodeB_dual_eq_var_r
                    nodeB_Li_to = node[node_key[self.to_node]].nodeB_dual_eq_var_i
                    nodeC_Lr_to = node[node_key[self.to_node]].nodeC_dual_eq_var_r
                    nodeC_Li_to = node[node_key[self.to_node]].nodeC_dual_eq_var_i
                    nodeN_Lr_to = node[node_key[self.to_node]].nodeN_dual_eq_var_r
                    nodeN_Li_to = node[node_key[self.to_node]].nodeN_dual_eq_var_i
            
                
                nodeN_Vr_from = node[node_key[self.from_node]].nodeN_Vr
                nodeN_Vi_from = node[node_key[self.from_node]].nodeN_Vi

                # Find the to bus nodes to stamp
                nodeA_Vr_to = node[node_key[self.to_node]].nodeA_Vr
                nodeA_Vi_to = node[node_key[self.to_node]].nodeA_Vi
                nodeB_Vr_to = node[node_key[self.to_node]].nodeB_Vr
                nodeB_Vi_to = node[node_key[self.to_node]].nodeB_Vi
                nodeC_Vr_to = node[node_key[self.to_node]].nodeC_Vr
                nodeC_Vi_to = node[node_key[self.to_node]].nodeC_Vi
                nodeN_Vr_to = node[node_key[self.to_node]].nodeN_Vr
                nodeN_Vi_to = node[node_key[self.to_node]].nodeN_Vi

                if self.stamp_dual:
                    nodes_from = {
                        'VR': [
                            nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from, nodeN_Vr_from
                        ],
                        'VI': [
                            nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from, nodeN_Vi_from
                        ],
                        'LR': [
                            nodeA_Lr_from, nodeB_Lr_from, nodeC_Lr_from, nodeN_Lr_from
                        ],
                        'LI':[
                            nodeA_Li_from, nodeB_Li_from, nodeC_Li_from, nodeN_Li_from]
                    }
        
                    nodes_to = {
                        'VR': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to, nodeN_Vr_to],
                        'VI': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to, nodeN_Vi_to],
                        'LR': [nodeA_Lr_to, nodeB_Lr_to, nodeC_Lr_to, nodeN_Lr_to],
                        'LI': [nodeA_Li_to, nodeB_Li_to, nodeC_Li_to, nodeN_Li_to]
                    }
                    
                else:
                    nodes_from = {
                        'VR': [
                            nodeA_Vr_from, nodeB_Vr_from, nodeC_Vr_from, nodeN_Vr_from
                        ],
                        'VI': [
                            nodeA_Vi_from, nodeB_Vi_from, nodeC_Vi_from, nodeN_Vi_from
                        ]
                    }
        
                    nodes_to = {
                        'VR': [nodeA_Vr_to, nodeB_Vr_to, nodeC_Vr_to, nodeN_Vr_to],
                        'VI': [nodeA_Vi_to, nodeB_Vi_to, nodeC_Vi_to, nodeN_Vi_to]
                    }

                nodes_from = SimpleNamespace(**nodes_from)
                nodes_to = SimpleNamespace(**nodes_to)

        return nodes_from, nodes_to

