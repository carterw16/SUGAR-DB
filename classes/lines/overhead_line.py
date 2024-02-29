import numpy as np

from classes.GlobalVars import _PHASE
from classes.lines.base import Lines


class OverheadLines(Lines):

    def __init__(self, ID, name, type, phases, nominal_voltage, from_node,
                 to_node, line_parameters, length, spacing, unused_phases, diameter, conductor_phases, 
                 gmr, resistance, impedance_matrix, capacitance_matrix, stamp_dual = False):
        
        super(OverheadLines,
              self).__init__(ID, name, type, phases, nominal_voltage, from_node,
                             to_node, length, impedance_matrix,
                             capacitance_matrix, spacing, unused_phases, conductor_phases, stamp_dual)
        self.line_parameters = line_parameters

        if line_parameters['file type'] == 'DSS':
            self.length = line_parameters['length']
            if line_parameters['given Z'] == True:
                self._Zmatrix = line_parameters['impedance matrix']
                self._Cmatrix = line_parameters['capacitance matrix']
            
                self.neutral_exists = True if len(
                    self._Zmatrix) > 3 and self.phases.find('N') else False
                
                # Set single phase variable if line is single phase
                if self._conductor_phases in [_PHASE['AN'], _PHASE['A']]:
                    self.single_phase = 0
                elif self._conductor_phases in [_PHASE['BN'], _PHASE['B']]:
                    self.single_phase = 1
                elif self._conductor_phases in [_PHASE['CN'], _PHASE['C']]:
                    self.single_phase = 2
                else:
                    self.single_phase = -1
                
                # if neutral exists, remove the neutral row and column
                if self.neutral_exists:
                    # check size of the capacitance and impedance matrices to ensure they are 3x3
                    if len(self._Zmatrix) == 4:
                        # if neutral included in matrix delete neutral row and column
                        self._Zmatrix = np.delete(self._Zmatrix, 3, axis=0)
                        self._Zmatrix = np.delete(self._Zmatrix, 3, axis=1)
                    if len(self._Cmatrix) == 4:
                        self._Cmatrix = np.delete(self._Cmatrix, 3, axis=0)
                        self._Cmatrix = np.delete(self._Cmatrix, 3, axis=1)
            
                # self.calc_overhead_Yshunt()
                if self._Cmatrix.size != 0 or len(self._Cmatrix) != 0:
                    self.hasShunt = True
                    _Y = 376.9911 * 1e-6 * 1j * self._Cmatrix
                    self._Yshunt = _Y
                    self.Yshunt = self._Yshunt * self.length

                # Calculate final values for impedance and shunt matrices
                # self.calc_overhead_Z()
                self.Zmatrix = self._Zmatrix * self.length * Lines.impedance_factor
                self.Ymatrix = self.calcInverse_DSS()  # calculate Ymatrix

        else:
            # What information do we need to calculate line parameters?
            # Name, ratings, gmr, resistance, diameter
            self.length = length 
            self.dia = diameter
            self.dia_N = self.dia[-1] if self.dia else None
            self.rad_N = self.dia_N / 24 if self.dia else None

            self.GMR = list(gmr.values())
            self.resistance = list(resistance.values())
            if self._conductor_phases in [
                    _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']
            ]:
                if self._conductor_phases == _PHASE['ABN']:
                    self._phases_id = [0, 1, 3]
                elif self._conductor_phases == _PHASE['BCN']:
                    self._phases_id = [1, 2, 3]
                else:
                    self._phases_id = [0, 2, 3]
            
            self.distances = [[0, self.D_AB, self.D_AC, self.D_AN],
                            [self.D_AB, 0, self.D_BC, self.D_BN],
                            [self.D_AC, self.D_BC, 0, self.D_CN],
                            [self.D_AN, self.D_BN, self.D_CN, 0]]
            self.distances = np.array(self.distances, dtype=float)
        

            # Initiliaze Radii
            if self.dia:
                if self._conductor_phases in [
                        _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']
                ]:
                    self.rad = [ele / 24 for ele in self.dia]
                elif self._conductor_phases in [
                        _PHASE['AN'], _PHASE['BN'], _PHASE['CN']
                ]:
                    non_zero_idx = [
                        idx for idx, val in enumerate(self.dia) if val != 0
                    ]
                    self.rad = self.dia[non_zero_idx[0]] / 24
                else:
                    self.rad = [ele / 24 for ele in self.dia if ele]

            # Set single phase variable if line is single phase
            if self._conductor_phases in [_PHASE['AN'], _PHASE['A']]:
                self.single_phase = 0
            elif self._conductor_phases in [_PHASE['BN'], _PHASE['B']]:
                self.single_phase = 1
            elif self._conductor_phases in [_PHASE['CN'], _PHASE['C']]:
                self.single_phase = 2
            else:
                self.single_phase = -1

            if self.unused_phases:
                self.unused_phases.sort()
                self.augment_matrices()
                self.neutral_exists = True if len(
                    self._Zmatrix) > 3 and self.phases.find('N') else False
                
                # if neutral exists, remove the neutral row and column
                if self.neutral_exists:
                    # check size of the capacitance and impedance matrices to ensure they are 3x3
                    if len(self._Zmatrix) == 4:
                        # if neutral included in matrix delete neutral row and column
                        self._Zmatrix = np.delete(self._Zmatrix, 3, axis=0)
                        self._Zmatrix = np.delete(self._Zmatrix, 3, axis=1)
                    if len(self._Cmatrix) == 4:
                        self._Cmatrix = np.delete(self._Cmatrix, 3, axis=0)
                        self._Cmatrix = np.delete(self._Cmatrix, 3, axis=1)

            self.calc_overhead_Yshunt()
            self.hasShunt = True

            # Calculate final values for impedance and shunt matrices
            self.calc_overhead_Z()
            self.Zmatrix = self._Zmatrix * self.length * Lines.impedance_factor
            self.Ymatrix = self.calcInverse()  # calculate Ymatrix
            self.Yshunt = self._Yshunt * self.length
    
    def calc_overhead_Z(self):
        # Three-phase calculations
        if self._conductor_phases in [_PHASE['ABC'], _PHASE['ABCN'], _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']]:
            if self._conductor_phases in [_PHASE['ABCN'], _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']]:
                _numCond = 4
            else:
                _numCond = 3
            skip_phase = {_PHASE['ABN']: 2, _PHASE['ACN']: 1, _PHASE['BCN']:0, _PHASE['ABC']: -1, _PHASE['ABCN']: -1}
            # Store Z Matrix as 2d list
            _Z = np.zeros((_numCond, _numCond), dtype=complex)  # Size 4 matrix
            for i in range(_numCond):
                if i != skip_phase[self._conductor_phases]:
                    for j in range(_numCond):
                        if i == j:
                            _Z[i, i] = self.calc_Zii(self.resistance[i], self.GMR[i])
                        else:
                            if j != skip_phase[self._conductor_phases]:
                                _D_ij = self.distances[i, j]
                                _Z[i, j] = self.calc_Zij(_D_ij)
            
            if self._phases == _PHASE['ABCN']:
                # Apply Kron's reduction if wye
                _Zij = _Z[:3, :3]
                _Zin = _Z[:3, 3, np.newaxis]  # singleton converted to 3x1
                _Znj = _Z[3, :3, np.newaxis].T  # singleton converted to 1x3
                _Znn = _Z[3, 3]
                self._Zmatrix = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)
            elif self._phases in [_PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']]:
                # Remove rows and columns of zeros
                _Z = _Z[~np.all(_Z == 0, axis=1)]  # For row
                _Z = _Z.T
                _Z = _Z[~np.all(_Z == 0, axis=1)]  # For column
                _Z = _Z.T
                # Apply Kron's reduction if wye
                # Slicing is not exclusive in last term careful
                # Get in the form [Z_ij Z_in; Z_nj Z_nn]
                _Zij = _Z[:2, :2]
                _Zin = _Z[:2, 2, np.newaxis]  # singleton converted to 3x1
                _Znj = _Z[2, :2, np.newaxis].T  # singleton converted to 1x3
                _Znn = _Z[2, 2]
                Zmatrix_2 = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)
                if self._phases == _PHASE['ABN']:
                    self._Zmatrix = np.array([[Zmatrix_2[0][0], Zmatrix_2[0][1], 0],
                                            [Zmatrix_2[1][0], Zmatrix_2[1][1], 0],
                                            [0, 0, 0]], dtype=complex)
                elif self._phases == _PHASE['BCN']:
                    self._Zmatrix = np.array([[0, 0, 0],
                                            [0, Zmatrix_2[0][0], Zmatrix_2[0][1]],
                                            [0, Zmatrix_2[1][0], Zmatrix_2[1][1]]], dtype=complex)
                elif self._phases == _PHASE['ACN']:
                    self._Zmatrix = np.array([[Zmatrix_2[0][0], 0, Zmatrix_2[0][1]],
                                            [0, 0, 0],
                                            [Zmatrix_2[1][0], 0, Zmatrix_2[1][1]]], dtype=complex)
            else:
                self._Zmatrix = _Z[:3, :3]

        # Single Phase Calculations
        elif self._phases in [_PHASE['AN'], _PHASE['BN'], _PHASE['CN']]:
            single_phase_index = {_PHASE['AN']:0, _PHASE['BN']:1, _PHASE['CN']:2}
            single_phase_spacing = {_PHASE['AN']:self.D_AN, _PHASE['BN']: self.D_BN, _PHASE['CN']:self.D_CN}
            # Single Phase Overhead Conductors
            self._Zmatrix = np.zeros((3, 3), dtype=complex)
            _Zij = self.calc_Zii(self.resistance[single_phase_index[self._phases]], self.GMR[single_phase_index[self._phases]])
            _Zin = self.calc_underground_Zij(single_phase_spacing[self._phases]) # D_spN just the single phase line spacing An, Bn, Cn, etc
            _Znj = self.calc_underground_Zij(single_phase_spacing[self._phases])
            _Znn = self.calc_Zii(self.resistance[-1], self.GMR[-1])
            _z = _Zij - (_Zin * _Znj) / _Znn  # Due to the values being scalar
            self._Zmatrix[self.single_phase, self.single_phase] = _z

        # Two phase one neutral calculation
        # elif self._phases in [_PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']]:
        #     # Store Z Matrix as 2d list
        #     _Z = np.zeros((4, 4), dtype=complex)  # Size 4 matrix
        #     for i in range(4):
        #         for j in range(4):
        #             if i == j:
        #                 if self._2phase_idx[i] & self._2phase_idx[j] == 1:
        #                     _Z[i, i] = self.calc_Zii(self.resistance[i], self.GMR[i])
        #             else:
        #                 if self._2phase_idx[i] & self._2phase_idx[j] == 1:
        #                     _D_ij = self.distances[i, j]
        #                     _Z[i, j] = self.calc_Zij(_D_ij)
        #     # Remove rows and columns of zeros
        #     _Z = _Z[~np.all(_Z == 0, axis=1)]  # For row
        #     _Z = _Z.T
        #     _Z = _Z[~np.all(_Z == 0, axis=1)]  # For column
        #     _Z = _Z.T
        #     # Apply Kron's reduction if wye
        #     # Slicing is not exclusive in last term careful
        #     # Get in the form [Z_ij Z_in; Z_nj Z_nn]
        #     _Zij = _Z[:2, :2]
        #     _Zin = _Z[:2, 2, np.newaxis]  # singleton converted to 3x1
        #     _Znj = _Z[2, :2, np.newaxis].T  # singleton converted to 1x3
        #     _Znn = _Z[2, 2]
        #     self.Zmatrix_2 = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)
        #     if self._phases == _PHASE['ABN']:
        #         self._Zmatrix = np.array([[self.Zmatrix_2[0][0], self.Zmatrix_2[0][1], 0],
        #                                    [self.Zmatrix_2[1][0], self.Zmatrix_2[1][1], 0],
        #                                    [0, 0, 0]], dtype=complex)
        #     elif self._phases == _PHASE['BCN']:
        #         self._Zmatrix = np.array([[0, 0, 0],
        #                                    [0, self.Zmatrix_2[0][0], self.Zmatrix_2[0][1]],
        #                                    [0, self.Zmatrix_2[1][0], self.Zmatrix_2[1][1]]], dtype=complex)
        #     elif self._phases == _PHASE['ACN']:
        #         self._Zmatrix = np.array([[self.Zmatrix_2[0][0], 0, self.Zmatrix_2[0][1]],
        #                                    [0, 0, 0],
        #                                    [self.Zmatrix_2[1][0], 0, self.Zmatrix_2[1][1]]], dtype=complex)


    def calc_overhead_Yshunt(self):
        if self._conductor_phases in [
                _PHASE['ABC'], _PHASE['ABCN'], _PHASE['ABCD']
        ]:
            _S = np.zeros((4, 4), dtype=complex)
            # Convert NoneType to 0
            [self.D_AE, self.D_BE, self.D_CE,
             self.D_NE] = map(self.setNonetoZero,
                              [self.D_AE, self.D_BE, self.D_CE, self.D_NE])
            _S[0, 0] = 2 * self.D_AE
            _S[1, 1] = 2 * self.D_BE
            _S[2, 2] = 2 * self.D_CE
            _S[3, 3] = 2 * self.D_NE
            for i in range(4):
                for j in range(4):
                    if (i == 0 and j == 1) or (j == 0 and i == 1):
                        _horizDist = self.D_AB
                        _vertDist = 2 * self.D_BE
                    elif (i == 0 and j == 2) or (j == 0 and i == 2):
                        _horizDist = self.D_AC
                        _vertDist = 2 * self.D_AE
                    elif (i == 0 and j == 3) or (j == 0 and i == 3):
                        _horizDist = (self.D_AN**2 -
                                      (self.D_AE - self.D_NE)**2)**0.5
                        _vertDist = self.D_AE + self.D_NE
                    elif (i == 1 and j == 2) or (j == 1 and i == 2):
                        _horizDist = self.D_BC
                        _vertDist = 2 * self.D_CE
                    elif (i == 1 and j == 3) or (j == 1 and i == 3):
                        _horizDist = (self.D_BN**2 -
                                      (self.D_BE - self.D_NE)**2)**0.5
                        _vertDist = self.D_BE + self.D_NE
                    elif (i == 2 and j == 3) or (j == 2 and i == 3):
                        _horizDist = (self.D_CN**2 -
                                      (self.D_CE - self.D_NE)**2)**0.5
                        _vertDist = self.D_CE + self.D_NE
                    if i != j:  # No diagonal terms exist here
                        _S[i, j] = (_horizDist**2 + _vertDist**2)**0.5
            # Multiply the diagonal elements by half as we are adding transpose to it
            if self._conductor_phases == _PHASE['ABC']:
                self.rad.append(-1)
            _dist_shunt_mat = [[
                self.rad[0] * 0.5, self.D_AB, self.D_AC, self.D_AN
            ], [0, self.rad[1] * 0.5, self.D_BC, self.D_BN],
                               [0, 0, self.rad[2] * 0.5, self.D_CN],
                               [0, 0, 0, self.rad[3] * 0.5]]
            _dist_shunt_mat = np.array(_dist_shunt_mat, dtype=float)
            _dist_shunt_mat += np.transpose(_dist_shunt_mat)
            _temp_mat = np.divide(_S, _dist_shunt_mat)
            Pprim = 11.17689 * np.log(_temp_mat)
            # Apply Kron's reduction if wye
            # Slicing is not exclusive in last term careful
            # There could be a problem if diagonal terms are zero
            if self._conductor_phases in [_PHASE['ABCN'], _PHASE['ABCD']]:
                _Pij = Pprim[:3, :3]
                _Pin = Pprim[:3, 3, np.newaxis]  # singleton to 3x1
                _Pnj = Pprim[3, :3, np.newaxis].T  # singleton to 1x3
                _Pnn = Pprim[3, 3]  # scalar
                Pred = self.kron_reduction(_Pij, _Pin, _Pnj, _Pnn)
            else:
                Pred = Pprim[0:3, 0:3]
            _C = np.linalg.inv(Pred)
            _Y = 376.9911 * 1e-6 * 1j * _C
            self._Yshunt = _Y
        elif self._conductor_phases in [
                _PHASE['AN'], _PHASE['BN'], _PHASE['CN'], _PHASE['A'],
                _PHASE['B'], _PHASE['C']
        ]:
            if self._conductor_phases in [
                    _PHASE['A'], _PHASE['B'], _PHASE['C']
            ]:
                self.rad = self.rad[self.single_phase]
            _S = np.zeros((2, 2), dtype=float)
            # Distance Calculations
            _D_E = [self.D_AE, self.D_BE, self.D_CE]
            _D_N = [self.D_AN, self.D_BN, self.D_CN]
            _horizDist = (_D_N[self.single_phase]**2 -
                          (_D_E[self.single_phase] - self.D_NE)**2)**0.5
            _vertDist = _D_E[self.single_phase] + self.D_NE
            # S matrix calculation
            _S[0, 1] = (_horizDist**2 + _vertDist**2)**0.5
            _S[1, 0] = _S[0, 1]
            _S[0, 0] = _D_E[self.single_phase]
            _S[1, 1] = self.D_NE
            # distance shunt matrix
            _dist_shunt_mat = np.array([[self.rad, _D_N[self.single_phase]],
                                        [_D_N[self.single_phase], self.rad_N]])

            _temp_mat = np.divide(_S, _dist_shunt_mat)
            Pprim = 11.17689 * np.log(_temp_mat)

            _Pij = Pprim[0, 0]
            _Pin = Pprim[0, 1]
            _Pnj = Pprim[1, 0]
            _Pnn = Pprim[1, 1]
            Pred = self.kron_reduction(_Pij, _Pin, _Pnj, _Pnn)
            _C = Pred**-1
            _Y = 376.9911 * 1e-6 * 1j * _C
            self._Yshunt = np.zeros((3, 3), dtype=complex)
            self._Yshunt[self.single_phase, self.single_phase] = _Y
        else:
            _S = np.zeros((4, 4), dtype=complex)
            # Convert NoneType to 0
            [self.D_AE, self.D_BE, self.D_CE,
             self.D_NE] = map(self.setNonetoZero,
                              [self.D_AE, self.D_BE, self.D_CE, self.D_NE])
            if self._conductor_phases & 0x1 == 1:
                _S[0, 0] = 2 * self.D_AE
            else:
                self.D_AB = -1
                self.D_AC = -1
                self.D_AN = -1
            if self._conductor_phases & 0x2 == 2:
                _S[1, 1] = 2 * self.D_BE
            else:
                self.D_AB = -1
                self.D_BC = -1
                self.D_BN = -1
            if self._conductor_phases & 0x4 == 4:
                _S[2, 2] = 2 * self.D_CE
            else:
                self.D_BC = -1
                self.D_AC = -1
                self.D_CN = -1
            if self._conductor_phases & 0x8 == 8:
                _S[3, 3] = 2 * self.D_NE
            else:
                self.D_AN = -1
                self.D_BN = -1
                self.D_CN = -1
            for i in range(4):
                for j in range(4):
                    flag = 0
                    if (i == 0 and j == 1) or (j == 0 and i == 1):
                        if self.D_AB != -1:
                            _horizDist = self.D_AB
                            _vertDist = 2 * self.D_BE
                            flag = 1
                    elif (i == 0 and j == 2) or (j == 0 and i == 2):
                        if self.D_AC != -1:
                            _horizDist = self.D_AC
                            _vertDist = 2 * self.D_AE
                            flag = 1
                    elif (i == 0 and j == 3) or (j == 0 and i == 3):
                        if self.D_AN != -1:
                            _horizDist = (self.D_AN**2 -
                                          (self.D_AE - self.D_NE)**2)**0.5
                            _vertDist = self.D_AE + self.D_NE
                            flag = 1
                    elif (i == 1 and j == 2) or (j == 1 and i == 2):
                        if self.D_BC != -1:
                            _horizDist = self.D_BC
                            _vertDist = 2 * self.D_CE
                            flag = 1
                    elif (i == 1 and j == 3) or (j == 1 and i == 3):
                        if self.D_BN != -1:
                            _horizDist = (self.D_BN**2 -
                                          (self.D_BE - self.D_NE)**2)**0.5
                            _vertDist = self.D_BE + self.D_NE
                            flag = 1
                    elif (i == 2 and j == 3) or (j == 2 and i == 3):
                        if self.D_CN != -1:
                            _horizDist = (self.D_CN**2 -
                                          (self.D_CE - self.D_NE)**2)**0.5
                            _vertDist = self.D_CE + self.D_NE
                            flag = 1
                    if i != j:  # No diagonal terms exist here
                        if flag == 1:
                            _S[i, j] = (_horizDist**2 + _vertDist**2)**0.5
            # Multiple the diag by half as we are adding with transpose
            _dist_shunt_mat = [
                [self.rad[0] * 0.5, self.D_AB, self.D_AC, self.D_AN],
                [0, self.rad[1] * 0.5, self.D_BC, self.D_BN],
                [0, 0, self.rad[2] * 0.5, self.D_CN],
                [0, 0, 0, self.rad[3] * 0.5]
            ]
            _dist_shunt_mat = np.array(_dist_shunt_mat, dtype=float)
            _dist_shunt_mat += np.transpose(_dist_shunt_mat)
            _dist_shunt_mat[_dist_shunt_mat == -1] = 0
            # Reduce the dimensions of _S and _dist_shunt_mat
            (_S_red, _rowIdx, _colIdx) = self.removeRowColZero(_S)
            (_dist_shunt_mat_red, _rowIdx,
             _colIdx) = self.removeRowColZero(_dist_shunt_mat)
            _temp_mat = np.divide(_S_red, _dist_shunt_mat_red)
            Pprim = 11.17689 * np.log(_temp_mat)
            # Apply Kron's reduction if wye
            # Slicing is not exclusive in last term careful
            # There could be a problem if diagonal terms are zero
            if self._conductor_phases & 0x8 == 8:
                _Pij = Pprim[:2, :2]
                _Pin = Pprim[:2, 2, np.newaxis]  # singleton to 3x1
                _Pnj = Pprim[2, :2, np.newaxis].T  # singleton to 1x3
                _Pnn = Pprim[2, 2]  # scalar
                Pred = self.kron_reduction(_Pij, _Pin, _Pnj, _Pnn)
            else:
                Pred = Pprim[0:2, 0:2]
            _C = np.linalg.inv(Pred)
            _Y = 376.9911 * 1e-6 * 1j * _C
            # Map into 3x3 structure
            _rowIdx = list(_rowIdx).index(True)
            _colIdx = list(_colIdx).index(True)
            _Y = np.insert(_Y, _rowIdx, 0 + 1j * 0, axis=0)
            _Y = np.insert(_Y, _colIdx, 0 + 1j * 0, axis=1)
            self._Yshunt = _Y
        return None
