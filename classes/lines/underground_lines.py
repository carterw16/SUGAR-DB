from classes.lines.base import Lines
from classes.GlobalVars import _PHASE

import numpy as np
import logging

class UndergroundLines(Lines):
    def __init__(self,
                 ID,
                 name,
                 type,
                 phases,
                 nominal_voltage,
                 from_node,
                 to_node,
                 line_parameters,
                 length,
                 impedance_matrix,
                 capacitance_matrix,
                 spacing,
                 unused_phases,
                 conductor_phases,
                 outer_diameter=None,
                 neutral_diameter=None,
                 neutral_strands=None,
                 conductor_diameter=None,
                 conductor_resistance=None,
                 conductor_gmr=None,
                 neutral_resistance=None,
                 neutral_gmr=None,
                 stamp_dual = False):

        super(UndergroundLines,
              self).__init__(ID, name, type, phases, nominal_voltage, from_node,
                             to_node, length, impedance_matrix,
                             capacitance_matrix, spacing, unused_phases,
                             conductor_phases, stamp_dual)
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
                self.Zmatrix = self._Zmatrix * self.length * Lines.impedance_factor

                if self._Cmatrix.size != 0 or len(self._Cmatrix) != 0:
                    self._Yshunt = self._Cmatrix * 2 * np.pi * self.freq * 1e-9 * 1j
                    self.hasShunt = True
                    self.Yshunt = self._Yshunt * self.length
                
                self.Ymatrix = self.calcInverse_DSS()  # calculate Ymatrix
        else:
            self.length = length
            # Initialize UG Line Parameters
            self.conductor_resistance = [ele for ele in conductor_resistance.values() if ele
                                ] if conductor_resistance else None
            self.conductor_gmr = [ele for ele in conductor_gmr.values() if ele
                                ] if conductor_gmr else None
            self.neutral_gmr = [ele for ele in neutral_gmr if ele
                            ] if neutral_gmr else None
            self.neutral_resistance = [ele for ele in neutral_resistance if ele
                                    ] if neutral_resistance else None
            self.outer_diameter = [ele for ele in outer_diameter if ele
                                ] if outer_diameter else None
            self.neutral_diameter = [ele for ele in neutral_diameter if ele
                                    ] if neutral_diameter else None
            self.neutral_strands = [ele for ele in neutral_strands if ele
                                ] if neutral_strands else None
            self.conductor_diameter = [ele for ele in conductor_diameter if ele
                                    ] if conductor_diameter else None
            self.D_ppN = 0

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

            # Initialize Underground Distances and Parameters
            self.init_underground_params()
            self.calc_underground_Z()
            self._Zmatrix = self.Zmatrix_U
            self.Zmatrix = self._Zmatrix * self.length * Lines.impedance_factor

            if self._Cmatrix.size != 0 or len(self._Cmatrix) != 0:
                self._Yshunt = self._Cmatrix * 2 * np.pi * self.freq * 1e-9 * 1j
                self.hasShunt = True
            else:
                self.calc_underground_Yshunt()
                self._Yshunt = self.Yshunt_U
                self.hasShunt = True

            # Calculate final values for impedance and shunt matrices
            self.Ymatrix = self.calcInverse()  # calculate Ymatrix
            self.Yshunt = self._Yshunt * self.length

    @staticmethod
    def calcEquiNeutDist(D_phaseN, D_phase2N, neutral_strands, ID):
        if D_phaseN == -1:
            logging.debug('PHASE TO NEUTRAL DISTANCE NOT SPECIFIED FOR %s' % ID)
            return -1
        D = (D_phaseN ** neutral_strands - D_phase2N ** neutral_strands) \
            ** (1 / neutral_strands)
        return D

    def init_underground_params(self):
        if self._conductor_phases in [
                _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']
        ]:
            if self._conductor_phases == _PHASE['ABN']:
                self._phases_id = [0, 1, 3]
            elif self._conductor_phases == _PHASE['BCN']:
                self._phases_id = [1, 2, 3]
            else:
                self._phases_id = [0, 2, 3]

            # 2 CONCENTRIC NEUTRAL PARAMETERS
            self.concN_RES = [
                self.neutral_resistance[ele] / self.neutral_strands[ele]
                for ele in range(2)
            ]
            # radius of circle passing through center of concentric metal strands
            self.nR = [
                (self.outer_diameter[ele] - self.neutral_diameter[ele]) / 24
                for ele in range(2)
            ]
            _tempVar = [
                self.conductor_gmr[ele] * self.neutral_strands[ele]
                for ele in range(2)
            ]
            self.concN_GMR = []
            for ele in range(2):
                self.concN_GMR.append(
                    (_tempVar[ele] *
                     self.nR[ele]**(self.neutral_strands[ele] - 1))**(
                         1 / self.neutral_strands[ele]))

            # The spacing between conductors and their concentric neutral is equal to resistance
            self.DppN = 3 * [0]
            _j = 0  # Mapping three to 2
            for _i in range(3):
                if _i in self._phases_id:
                    self.DppN[_i] = self.nR[_j]
                    _j += 1
                else:
                    self.DppN[_i] = -1
            [self.D_AAN, self.D_BBN,
             self.D_CCN] = [self.DppN[0], self.DppN[1], self.DppN[2]]
            _idx_phase = 0
            if 0 in self._phases_id:
                self.D_NAN = self.calcEquiNeutDist(
                    self.D_AN, self.D_AAN, self.neutral_strands[_idx_phase],
                    self.ID)
                self.D_BAN = self.calcEquiNeutDist(
                    self.D_AB, self.D_AAN, self.neutral_strands[_idx_phase],
                    self.ID) if 1 in self._phases_id else -1
                self.D_CAN = self.calcEquiNeutDist(
                    self.D_AC, self.D_AAN, self.neutral_strands[_idx_phase],
                    self.ID) if 2 in self._phases_id else -1
                _idx_phase += 1
            else:
                [self.D_NAN, self.D_BAN, self.D_CAN] = [-1, -1, -1]
            if 1 in self._phases_id:
                self.D_ABN = self.calcEquiNeutDist(
                    self.D_AB, self.D_BBN, self.neutral_strands[_idx_phase],
                    self.ID) if 0 in self._phases_id else -1
                self.D_NBN = self.calcEquiNeutDist(
                    self.D_BN, self.D_BBN, self.neutral_strands[_idx_phase],
                    self.ID)
                self.D_CBN = self.calcEquiNeutDist(
                    self.D_BC, self.D_BBN, self.neutral_strands[_idx_phase],
                    self.ID) if 2 in self._phases_id else -1
                _idx_phase += 1
            else:
                [self.D_NBN, self.D_ABN, self.D_CBN] = [-1, -1, -1]
            if 2 in self._phases_id:
                self.D_ACN = self.calcEquiNeutDist(
                    self.D_AC, self.D_CCN, self.neutral_strands[_idx_phase],
                    self.ID) if 0 in self._phases_id else -1
                self.D_BCN = self.calcEquiNeutDist(
                    self.D_BC, self.D_CCN, self.neutral_strands[_idx_phase],
                    self.ID) if 1 in self._phases_id else -1
                self.D_NCN = self.calcEquiNeutDist(
                    self.D_CN, self.D_CCN, self.neutral_strands[_idx_phase],
                    self.ID)
                _idx_phase += 1
            else:
                [self.D_NCN, self.D_ACN, self.D_BCN] = [-1, -1, -1]
            # Create 7x7 distance matrix
            self.distances_underground = self.create_distances_underground()

        elif self._conductor_phases in [
                _PHASE['AN'], _PHASE['BN'], _PHASE['CN'], _PHASE['A'],
                _PHASE['B'], _PHASE['C']
        ]:
            if self._conductor_phases in [_PHASE['AN'], _PHASE['A']]:
                self.single_phase = 0
                self.D_spN = self.D_AN
            elif self._conductor_phases in [_PHASE['BN'], _PHASE['B']]:
                self.single_phase = 1
                self.D_spN = self.D_BN
            else:
                self.single_phase = 2
                self.D_spN = self.D_CN
            self.concN_RES = self.neutral_resistance[0] / self.neutral_strands[0]
            self.nR = (self.outer_diameter[0] - self.neutral_diameter[0]) / 24
            _tempVar = self.neutral_gmr[0] * self.neutral_strands[
                0] * self.nR**(self.neutral_strands[0] - 1)
            self.concN_GMR = _tempVar**(1 / self.neutral_strands[0])

            # Line Spacing - It's defined at start
            self.D_ppN = self.nR
            if self._conductor_phases in [
                    _PHASE['AN'], _PHASE['BN'], _PHASE['CN']
            ]:
                self.distances_underground = np.array(
                    [[0, self.D_ppN, self.D_spN], [self.D_ppN, 0, self.D_spN],
                     [self.D_spN, self.D_spN, 0]])
            else:
                self.distances_underground = np.array([[0, self.D_ppN],
                                                       [self.D_ppN, 0]])

        else:
            self.nR = [
                (self.outer_diameter[ele] - self.neutral_diameter[ele]) / 24
                for ele in range(3)
            ]
            self.concN_RES = [
                self.neutral_resistance[ele] / self.neutral_strands[ele]
                for ele in range(3)
            ]
            _tempVar = [
                self.neutral_gmr[ele] * self.neutral_strands[ele]
                for ele in range(3)
            ]
            self.concN_GMR = []
            for ele in range(3):
                self.concN_GMR.append(
                    (_tempVar[ele] *
                     self.nR[ele]**(self.neutral_strands[ele] - 1))**(
                         1 / self.neutral_strands[ele]))

            # Line Spacing Initialization
            # 7-wire system with neutral - Distance matrix is 7x7
            if self._conductor_phases == _PHASE['ABCN']:
                [self.D_AAN, self.D_BBN,
                 self.D_CCN] = [self.nR[0], self.nR[1], self.nR[2]]
                if self.neutral_strands[0]:
                    (self.D_NAN, self.D_BAN, self.D_CAN) = map(
                        self.calcEquiNeutDist,
                        [self.D_AN, self.D_AB, self.D_AC], 3 * [self.D_AAN],
                        3 * [self.neutral_strands[0]], 3 * [self.ID])
                if self.neutral_strands[1]:
                    (self.D_ABN, self.D_NBN, self.D_CBN) = map(
                        self.calcEquiNeutDist,
                        [self.D_AB, self.D_BN, self.D_BC], 3 * [self.D_BBN],
                        3 * [self.neutral_strands[1]], 3 * [self.ID])
                if self.neutral_strands[2]:
                    (self.D_ACN, self.D_BCN, self.D_NCN) = map(
                        self.calcEquiNeutDist,
                        [self.D_AC, self.D_BC, self.D_CN], 3 * [self.D_CCN],
                        3 * [self.neutral_strands[2]], 3 * [self.ID])
                # Create distance matrix
                self.distances_underground = self.create_distances_underground()
            # 6 wire system without neutral - Distance Matrix is 6x6
            elif self._conductor_phases == _PHASE['ABC']:
                # The spacing between conductors and their concentric neutral is equal to resistance
                [self.D_AAN, self.D_BBN,
                 self.D_CCN] = [self.nR[0], self.nR[1], self.nR[2]]
                if self.neutral_strands[0]:
                    (self.D_BAN,
                     self.D_CAN) = map(self.calcEquiNeutDist,
                                       [self.D_AB, self.D_AC], 2 * [self.D_AAN],
                                       2 * [self.neutral_strands[0]],
                                       2 * [self.ID])
                if self.neutral_strands[1]:
                    (self.D_ABN,
                     self.D_CBN) = map(self.calcEquiNeutDist,
                                       [self.D_AB, self.D_BC], 2 * [self.D_BBN],
                                       2 * [self.neutral_strands[1]],
                                       2 * [self.ID])
                if self.neutral_strands[2]:
                    (self.D_ACN,
                     self.D_BCN) = map(self.calcEquiNeutDist,
                                       [self.D_AC, self.D_BC], 2 * [self.D_CCN],
                                       2 * [self.neutral_strands[2]],
                                       2 * [self.ID])
                # Create distance matrix
                self.distances_underground = self.create_distances_underground_6wire(
                )

    def create_distances_underground(self):
        # D12 = D45 and D23 = D56 and so on
        distances_underground = np.array \
         ([[0, self.D_AB, self.D_AC, self.D_AAN, self.D_ABN, self.D_ACN, self.D_AN],
           [0, 0, self.D_BC, self.D_BAN, self.D_BBN, self.D_BCN, self.D_BN],
           [0, 0, 0, self.D_CAN, self.D_CBN, self.D_CCN, self.D_CN],
           [0, 0, 0, 0, self.D_AB, self.D_AC, self.D_NAN],
           [0, 0, 0, 0, 0, self.D_BC, self.D_NBN],
           [0, 0, 0, 0, 0, 0, self.D_NCN],
           [0, 0, 0, 0, 0, 0, 0]])
        distances_underground += np.transpose(distances_underground)
        return distances_underground

    def create_distances_underground_6wire(self):
        # D12 = D45 and D23 = D56 and so on
        distances_underground = np.array \
         ([[0, self.D_AB, self.D_AC, self.D_AAN, self.D_ABN, self.D_ACN],
           [0, 0, self.D_BC, self.D_BAN, self.D_BBN, self.D_BCN],
           [0, 0, 0, self.D_CAN, self.D_CBN, self.D_CCN],
           [0, 0, 0, 0, self.D_AB, self.D_AC],
           [0, 0, 0, 0, 0, self.D_BC],
           [0, 0, 0, 0, 0, 0]])
        distances_underground += np.transpose(distances_underground)
        return distances_underground

    def calc_underground_Z(self):
        # When calculating for three phase conductor + 3 neutral + 1 additional neutral
        if self._conductor_phases in [_PHASE['ABCN']]:
            _Z = np.zeros(
                (7, 7), dtype=complex
            )  # Represent the neutral for each individual conductor + 1
            _R = [
                self.conductor_resistance[0:-1] + self.concN_RES +
                [self.neutral_resistance[-1]]
            ][0]
            _GMR = [
                self.conductor_gmr[0:-1] + self.concN_GMR +
                [self.neutral_gmr[-1]]
            ][0]
            for i in range(7):
                for j in range(7):
                    if i == j:
                        _Z[i, j] = self.calc_underground_Zii(_R[i], _GMR[i])
                    else:
                        _D_ij = self.distances_underground[i][j]
                        _Z[i, j] = self.calc_underground_Zij(_D_ij)
            # Evaluate Zij, Zin, Znj, Znn
            _Zij = _Z[:3, :3]
            _Zin = _Z[:3, 3:]
            _Znj = _Z[3:, :3]
            _Znn = _Z[3:, 3:]
            self.Zmatrix_U = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)

        # When calculating for three phase conductor + 3 neutral
        elif self._conductor_phases in [_PHASE['ABC']]:
            _Z = np.zeros(
                (6, 6), dtype=complex
            )  # Represent the neutral for each individual conductor
            _R = [self.conductor_resistance + self.concN_RES][0]
            _GMR = [self.conductor_gmr + self.concN_GMR][0]
            for i in range(6):
                for j in range(6):
                    if i == j:
                        _Z[i, j] = self.calc_underground_Zii(_R[i], _GMR[i])
                    else:
                        _D_ij = self.distances_underground[i][j]
                        _Z[i, j] = self.calc_underground_Zij(_D_ij)
            # Evaluate Zij, Zin, Znj, Znn
            _Zij = _Z[:3, :3]
            _Zin = _Z[:3, 3:]
            _Znj = _Z[3:, :3]
            _Znn = _Z[3:, 3:]
            self.Zmatrix_U = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)

        # Two conductors and neutral
        elif self._conductor_phases in [
                _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']
        ]:
            _Z = np.zeros(
                (7, 7), dtype=complex
            )  # Represent the neutral for each individual conductor + 1
            _R = [
                self.conductor_resistance + self.concN_RES +
                [self.neutral_resistance]
            ][0]  # These are 5 in this
            # case
            _GMR = [self.conductor_gmr + self.concN_GMR + [self.neutral_gmr]
                   ][0]  # These are 5 in this case
            self._phases_id_ex = [
                self._phases_id[0], self._phases_id[0] + 3, self._phases_id[1],
                self._phases_id[1] + 3, 6
            ]
            _temp_i = 0
            for i in range(7):
                for j in range(7):
                    if i == j:
                        if i in self._phases_id_ex:
                            _Z[i, j] = self.calc_underground_Zii(
                                _R[_temp_i], _GMR[_temp_i])
                            _temp_i += 1
                    else:
                        _D_ij = self.distances_underground[i][j]
                        if _D_ij != -1:
                            _Z[i, j] = self.calc_underground_Zij(_D_ij)
            # Remove all the rows with 0
            _Z = _Z[~np.all(_Z == 0, axis=1)]
            # Remove all the cols with 0
            _Z = _Z[:, ~np.all(_Z == 0, axis=0)]
            _Zij = _Z[:2, :2]
            _Zin = _Z[:2, 2:]
            _Znj = _Z[2:, :2]
            _Znn = _Z[2:, 2:]
            cond_Zmatrix_U = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)
            self.Zmatrix_U = np.zeros((3, 3), dtype=complex)
            self.Zmatrix_U[self._phases_id[0],
                           self._phases_id[0]] = cond_Zmatrix_U[0, 0]
            self.Zmatrix_U[self._phases_id[1],
                           self._phases_id[1]] = cond_Zmatrix_U[1, 1]
            self.Zmatrix_U[self._phases_id[0],
                           self._phases_id[1]] = cond_Zmatrix_U[0, 1]
            self.Zmatrix_U[self._phases_id[1],
                           self._phases_id[0]] = cond_Zmatrix_U[1, 0]

        # Single conductor and neutral
        elif self._conductor_phases in [
                _PHASE['AN'],
                _PHASE['BN'],
                _PHASE['CN'],
                _PHASE['A'],
                _PHASE['B'],
                _PHASE['C'],
        ]:
            # Additional neutral conductor present
            self.Zmatrix_U = np.zeros((3, 3), dtype=complex)
            if self._conductor_phases & 0x08 == int(8):
                _Z = np.zeros((3, 3), dtype=complex)
                _R = [
                    self.conductor_resistance[0], self.concN_RES,
                    self.neutral_resistance[-1]
                ]
                _GMR = [
                    self.conductor_gmr[0], self.concN_GMR, self.neutral_gmr[-1]
                ]
                _numCond = 3
            # No additional neutral
            else:
                _Z = np.zeros((2, 2), dtype=complex)
                _R = [self.conductor_resistance[0], self.concN_RES]
                _GMR = [self.conductor_gmr[0], self.concN_GMR]
                _numCond = 2

            for i in range(_numCond):
                for j in range(_numCond):
                    if i == j:
                        _Z[i, j] = self.calc_underground_Zii(_R[i], _GMR[i])
                    else:
                        _D_ij = self.distances_underground[i][j]
                        _Z[i, j] = self.calc_underground_Zij(_D_ij)
            # Evaluate Zij, Zin, Znj, Znn
            _Zij = _Z[:1, :1]
            _Zin = _Z[:1, 1:]
            _Znj = _Z[1:, :1]
            _Znn = _Z[1:, 1:]
            _z = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)
            self.Zmatrix_U[self.single_phase, self.single_phase] = _z[0][0]
        # Send a warning message
        else:
            logging.debug('INCORRECT PHASES IN THE LINE-CONFIGURATION')

    def calc_underground_Yshunt(self):
        # Convert nR in feet to inches (nR - ft / dia - inches / neutral_dia - inches)
        _Rb = self.nR[0] * 12 if hasattr(self.nR, '__len__') else (self.nR * 12)
        _conductor_rad = self.conductor_diameter[0] / 2 if hasattr(self.conductor_diameter, '__len__') else \
         self.conductor_diameter / 2
        _conductor_neutral_rad = self.neutral_diameter[0] / 2 \
         if hasattr(self.neutral_diameter, '__len__') else self.neutral_diameter / 2
        _k = self.neutral_strands[0] if hasattr(
            self.neutral_strands, '__len__') else self.neutral_strands
        _Vp1 = np.log(_Rb / _conductor_rad) - (1 / _k) * np.log(
            _k * _conductor_neutral_rad / _Rb)
        _y = 77.3619 * 1j * 1e-6 / _Vp1

        # Create the Yshunt matrix (uS/miles)
        self.Yshunt_U = np.zeros((3, 3), dtype=complex)
        if self._conductor_phases in [_PHASE['ABC'], _PHASE['ABCN']]:
            np.fill_diagonal(self.Yshunt_U, _y)
        elif self._conductor_phases in [
                _PHASE['AN'], _PHASE['BN'], _PHASE['CN'], _PHASE['A'],
                _PHASE['B'], _PHASE['C']
        ]:
            self.Yshunt_U[self.single_phase, self.single_phase] = _y
        elif self._conductor_phases in [
                _PHASE['ABN'], _PHASE['BCN'], _PHASE['ACN']
        ]:
            self.Yshunt_U[self._phases_id[0], self._phases_id[0]] = _y
            self.Yshunt_U[self._phases_id[1], self._phases_id[1]] = _y
        else:
            logging.debug("Incorrect phases - Underground Shunt")