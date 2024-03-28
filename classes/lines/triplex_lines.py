from classes.lines.base import Lines

import numpy as np

class TriplexLines(Lines):
    def __init__(self, ID, name, type, phases, nominal_voltage, from_node,
                 to_node, line_parameters, length, impedance_matrix, capacitance_matrix, 
                 diameter, resistance, gmr, insulation_thickness,
                 stamp_dual = False):
        
        super(TriplexLines, self).__init__(ID, name, type, phases, nominal_voltage, from_node,
                             to_node, length, impedance_matrix,
                             capacitance_matrix, None, None, None, stamp_dual = stamp_dual)

        self.line_parameters = line_parameters
        self.length = length # Note: we are assuming length in feet from test file
        self.diameter = diameter
        self.resistance = resistance
        self.gmr = gmr
        self.insulation_thickness = insulation_thickness

        self.D_12 = (self.diameter[0] + 2 * self.insulation_thickness[0]) * (1 / 12)
        self.D_13 = (self.diameter[1] + self.insulation_thickness[1]) * (1 / 12)
        self.D_23 = self.D_13
        
        self.triplex_distance_matrix = np.array([[0, self.D_12, self.D_13],
                                                 [self.D_12, 0, self.D_23],
                                                 [self.D_13, self.D_23, 0]])

        self.calc_Triplex_Zmatrix()

        # Calculate final values for impedance and shunt matrices
        self.Zmatrix = self._Zmatrix * self.length * self.impedance_factor
        self.Ymatrix = np.linalg.inv(self.Zmatrix)
    
    def calc_Triplex_Zmatrix(self):
        # Phase parameters
        self.conductor_gmr = [self.gmr['A'],
                              self.gmr['B'],
                              self.gmr['N']]
        self.conductor_resistance = [self.resistance['A'],
                              self.resistance['B'],
                              self.resistance['N']]
        
        # Calculate 3x3 impedance matrix
        # Equations come from Kersting
        _Z = np.zeros((3, 3), dtype=complex)
        for i in range(3):
            for j in range(3):
                if i == j:
                    _Z[i, j] = self.conductor_resistance[i] + 0.09530 + 1j * 0.12134 * (
                            np.log(1 / self.conductor_gmr[i]) + 7.93402)
                else:
                    _Z[i, j] = 0.09530 + 1j * 0.12134 * (np.log(1 / self.triplex_distance_matrix[i, j]) + 7.93402)
        # Reduce the impedance matrix to 2x2
        _Zij = _Z[:2, :2]
        _Zin = _Z[:2, 2, np.newaxis]  # singleton to 2x1
        _Znj = _Z[2, :2, np.newaxis].T  # singleton to 1x2
        _Znn = _Z[2, 2]
        self._Zmatrix = self.kron_reduction(_Zij, _Zin, _Znj, _Znn)

