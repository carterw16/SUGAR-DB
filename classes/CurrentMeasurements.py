"""Measure the current of a powerflow element.

  Author(s): Naeem Turner-Bandele
  Created Date: 06-12-2020
  Updated Date: 10-14-2021
  Email: nturnerb@cmu.edu
  Status: Development

  Creates a current measurement device for a specified powerflow deliver element. After running the power flow,
  current measurements are given at the output. Currently implemented for lines only.

  Typical usage example:

  ammeter = CurrentMeasurements()
"""

from itertools import count

import numpy as np

from classes.Elements import NonlinearElement
from classes.GlobalVars import _PHASE
from classes.Nodes import Nodes
from classes.PTR import type_condition
from classes.CurrentMeasurementStamps import stamp_line_currents
from types import MethodType


class CurrentMeasurements(NonlinearElement):
    _ids = count(0)
    curr_meas_key = {}
    triplex_lines = False

    def __init__(self,
                 from_node,
                 to_node,
                 Imax,
                 name,
                 id,
                 obj_id,
                 obj_type,
                 phases,
                 Y_series,
                 Y_shunt=None):
        super(CurrentMeasurements, self).__init__()
        self.from_node = from_node
        self.to_node = to_node
        self.name = name if len(name) != 0 else obj_id
        self.Imax = Imax
        self.enable = True
        self.id = self._ids.__next__()
        self.type = obj_type
        self.isTriplex = True if self.type == 'Triplex Line' or self.type == 'XfmrCT' else False
        self.phases = phases if type(phases) is not int else [
            key for key, val in _PHASE.items() if val == int(phases)
        ].pop()
        self.all_phases = {
            0: "A",
            1: "B",
            2: "C",
            3: "N"
        } if not self.isTriplex else {
            0: "1",
            1: "2"
        }

        # Initialize Series and Shunt Admittances
        self.Y_series = Y_series
        self.Y_shunt = Y_shunt
        self.hasShunt = True if self.Y_shunt.any() else False

        # Initialize and Assign Nodes for Current Measurement
        self.node_Imag2 = -1
        self.node_Ir = -1
        self.node_Ii = -1
        self.node_self_Vr_from = -1
        self.node_self_Vr_to = -1
        self.node_self_Vi_from = -1
        self.node_self_Vi_to = -1
        self.node_mutual1_Vr_from = -1
        self.node_mutual1_Vr_to = -1
        self.node_mutual1_Vi_from = -1
        self.node_mutual1_Vi_to = -1
        self.node_mutual2_Vr_from = -1
        self.node_mutual2_Vr_to = -1
        self.node_mutual2_Vi_from = -1
        self.node_mutual2_Vi_to = -1

        if self.isTriplex:
            CurrentMeasurements.triplex_lines = self.isTriplex
            self.I1_mag = 0.0
            self.I1_ang = 0.0

            self.I2_mag = 0.0
            self.I2_ang = 0.0

            self.node1_Imag2 = -1
            self.node2_Imag2 = -1
            self.node1_Ir = -1
            self.node1_Ii = -1
            self.node2_Ir = -1
            self.node2_Ii = -1
        else:
            self.Ia_mag = 0.0
            self.Ia_ang = 0.0

            self.Ib_mag = 0.0
            self.Ib_ang = 0.0

            self.Ic_mag = 0.0
            self.Ic_ang = 0.0
            self.nodeA_Imag2 = -1
            self.nodeB_Imag2 = -1
            self.nodeC_Imag2 = -1
            self.nodeA_Ir = -1
            self.nodeA_Ii = -1
            self.nodeB_Ir = -1
            self.nodeB_Ii = -1
            self.nodeC_Ir = -1
            self.nodeC_Ii = -1

        self.nodes_voltages = {
            "A": [],
            "B": [],
            "C": []
        } if not self.isTriplex else {
            "1": [],
            "2": []
        }
        self.nodes_Imag2 = {
            "A": -1,
            "B": -1,
            "C": -1
        } if not self.isTriplex else {
            "1": -1,
            "2": -1
        }
        self.nodes_I = {
            "A": [],
            "B": [],
            "C": []
        } if not self.isTriplex else {
            "1": [],
            "2": []
        }

        CurrentMeasurements.curr_meas_key[(self.from_node, self.to_node,
                                           id)] = self.id

        # Link Cython Methods to Class #
        self.stamp_line_currents = MethodType(stamp_line_currents, self)

    def assign_nodes(self, node_index_, node_key, node):
        node_from = node[node_key[self.from_node]]
        node_to = node[node_key[self.to_node]]

        if not self.isTriplex:
            self.nodes_voltages["A"] = [
                node_from.nodeA_Vr, node_from.nodeA_Vi, node_to.nodeA_Vr,
                node_to.nodeA_Vi
            ]
            self.nodeA_Imag2 = node_index_.__next__()
            self.nodes_Imag2["A"] = self.nodeA_Imag2
            self.nodeA_Ir = node_index_.__next__()
            self.nodeA_Ii = node_index_.__next__()
            self.nodes_I["A"] = [self.nodeA_Ir, self.nodeA_Ii]

            self.nodes_voltages["B"] = [
                node_from.nodeB_Vr, node_from.nodeB_Vi, node_to.nodeB_Vr,
                node_to.nodeB_Vi
            ]
            self.nodeB_Imag2 = node_index_.__next__()
            self.nodes_Imag2["B"] = self.nodeB_Imag2
            self.nodeB_Ir = node_index_.__next__()
            self.nodeB_Ii = node_index_.__next__()
            self.nodes_I["B"] = [self.nodeB_Ir, self.nodeB_Ii]

            self.nodes_voltages["C"] = [
                node_from.nodeC_Vr, node_from.nodeC_Vi, node_to.nodeC_Vr,
                node_to.nodeC_Vi
            ]
            self.nodeC_Imag2 = node_index_.__next__()
            self.nodes_Imag2["C"] = self.nodeC_Imag2
            self.nodeC_Ir = node_index_.__next__()
            self.nodeC_Ii = node_index_.__next__()
            self.nodes_I["C"] = [self.nodeC_Ir, self.nodeC_Ii]
        else:
            self.nodes_voltages["1"] = [
                node_from.node1_Vr, node_from.node1_Vi, node_to.node1_Vr,
                node_to.node1_Vi
            ]
            self.node1_Imag2 = node_index_.__next__()
            self.nodes_Imag2["1"] = self.node1_Imag2
            self.node1_Ir = node_index_.__next__()
            self.node1_Ii = node_index_.__next__()
            self.nodes_I["1"] = [self.node1_Ir, self.node1_Ii]

            self.nodes_voltages["2"] = [
                node_from.node2_Vr, node_from.node2_Vi, node_to.node2_Vr,
                node_to.node2_Vi
            ]
            self.node2_Imag2 = node_index_.__next__()
            self.nodes_Imag2["2"] = self.node2_Imag2
            self.node2_Ir = node_index_.__next__()
            self.node2_Ii = node_index_.__next__()
            self.nodes_I["2"] = [self.node2_Ir, self.node2_Ii]

        if self.isTriplex:
            Nodes.I_index.extend(
                [self.node1_Ir, self.node1_Ii, self.node2_Ir, self.node2_Ii])
        else:
            Nodes.I_index.extend([
                self.nodeA_Ir, self.nodeA_Ii, self.nodeB_Ir, self.nodeB_Ii,
                self.nodeC_Ir, self.nodeC_Ii
            ])
        return node_index_

    def calc_currents(self, V):
        if self.isTriplex:
            self.I1_mag = np.around(
                np.abs(V[self.node1_Ir] + 1j * V[self.node1_Ii]).item(), 2)
            self.I1_ang = np.around(
                np.angle(V[self.node1_Ir] + 1j * V[self.node1_Ii],
                         deg=True).item(), 2)

            self.I2_mag = np.around(
                np.abs(V[self.node2_Ir] + 1j * V[self.node2_Ii]).item(), 2)
            self.I2_ang = np.around(
                np.angle(V[self.node2_Ir] + 1j * V[self.node2_Ii],
                         deg=True).item(), 2)

        else:
            self.Ia_mag = np.around(
                np.abs(V[self.nodeA_Ir] + 1j * V[self.nodeA_Ii]).item(), 2)
            self.Ia_ang = np.around(
                np.angle(V[self.nodeA_Ir] + 1j * V[self.nodeA_Ii],
                         deg=True).item(), 2)

            self.Ib_mag = np.around(
                np.abs(V[self.nodeB_Ir] + 1j * V[self.nodeB_Ii]).item(), 2)
            self.Ib_ang = np.around(
                np.angle(V[self.nodeB_Ir] + 1j * V[self.nodeB_Ii],
                         deg=True).item(), 2)

            self.Ic_mag = np.around(
                np.abs(V[self.nodeC_Ir] + 1j * V[self.nodeC_Ii]).item(), 2)
            self.Ic_ang = np.around(
                np.angle(V[self.nodeC_Ir] + 1j * V[self.nodeC_Ii],
                         deg=True).item(), 2)

    def stamp_nonlinear(self, V, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val,
                        Jnlin_row, idx_Y, idx_J):

        _phases = range(3) if not self.isTriplex else range(2)

        for _self in _phases:
            phases = {0, 1, 2} if not self.isTriplex else {0, 1}
            phases.remove(_self)
            _mutual1 = phases.pop()

            # Get currents and assign current nodes
            self.node_Ir = self.nodes_I[self.all_phases[_self]][0]
            self.node_Ii = self.nodes_I[self.all_phases[_self]][1]
            self.node_Imag2 = self.nodes_Imag2[self.all_phases[_self]]

            Ir = V[self.node_Ir]
            Ii = V[self.node_Ii]
            Imag2 = V[self.node_Imag2]

            # Collect self nodes and voltages
            self.node_self_Vr_from = self.nodes_voltages[
                self.all_phases[_self]][0]
            self.node_self_Vi_from = self.nodes_voltages[
                self.all_phases[_self]][1]
            self.node_self_Vr_to = self.nodes_voltages[
                self.all_phases[_self]][2]
            self.node_self_Vi_to = self.nodes_voltages[
                self.all_phases[_self]][3]

            V_self_from = complex(V[self.node_self_Vr_from],
                                  V[self.node_self_Vi_from])
            V_self_to = complex(V[self.node_self_Vr_to],
                                V[self.node_self_Vi_to])

            # Collect mutual 1 nodes and voltages
            self.node_mutual1_Vr_from = self.nodes_voltages[
                self.all_phases[_mutual1]][0]
            self.node_mutual1_Vi_from = self.nodes_voltages[
                self.all_phases[_mutual1]][1]
            self.node_mutual1_Vr_to = self.nodes_voltages[
                self.all_phases[_mutual1]][2]
            self.node_mutual1_Vi_to = self.nodes_voltages[
                self.all_phases[_mutual1]][3]

            V_mutual1_from = complex(V[self.node_mutual1_Vr_from],
                                     V[self.node_mutual1_Vi_from])
            V_mutual1_to = complex(V[self.node_mutual1_Vr_to],
                                   V[self.node_mutual1_Vi_to])

            # Collect mutual 2 nodes and voltages if they exist and then stamp
            if self.isTriplex:
                idx_Y, idx_J = self.stamp_line_currents(
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
                    idx_Y, idx_J, _self, V_self_from, V_self_to, _mutual1,
                    V_mutual1_from, V_mutual1_to, Imag2, Ir, Ii)
            else:
                _mutual2 = phases.pop()
                self.node_mutual2_Vr_from = self.nodes_voltages[
                    self.all_phases[_mutual2]][0]
                self.node_mutual2_Vi_from = self.nodes_voltages[
                    self.all_phases[_mutual2]][1]
                self.node_mutual2_Vr_to = self.nodes_voltages[
                    self.all_phases[_mutual2]][2]
                self.node_mutual2_Vi_to = self.nodes_voltages[
                    self.all_phases[_mutual2]][3]

                V_mutual2_from = complex(V[self.node_mutual2_Vr_from],
                                         V[self.node_mutual2_Vi_from])
                V_mutual2_to = complex(V[self.node_mutual2_Vr_to],
                                       V[self.node_mutual2_Vi_to])

                idx_Y, idx_J = self.stamp_line_currents(
                    Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row,
                    idx_Y, idx_J, _self, V_self_from, V_self_to, _mutual1,
                    V_mutual1_from, V_mutual1_to, Imag2, Ir, Ii, V_mutual2_from,
                    V_mutual2_to, _mutual2)

        return idx_Y, idx_J


def create_current_meas(enable_CM, overheadline, undergroundline, triplexline):
    curr_meas = []

    if enable_CM:
        curr_meas = initialize_ammeter(overheadline, curr_meas)
        curr_meas = initialize_ammeter(undergroundline, curr_meas)
        curr_meas = initialize_ammeter(triplexline, curr_meas)
    return curr_meas


def initialize_ammeter(obj, curr_meas):
    if obj:
        for ele in obj:
            obj_type = type_condition(ele)
            obj_id = ele.ID
            ammeter = CurrentMeasurements(ele.from_node, ele.to_node, 1,
                                          ele.name, ele.id, obj_id, obj_type,
                                          ele.phases,
                                          ele.Ymatrix, ele.Yshunt)
            curr_meas.append(ammeter)
    return curr_meas


def find_ammeter(ele, curr_meas, enable_CM):
    ammeter = []
    if enable_CM:
        from_node = ele.from_node
        to_node = ele.to_node
        obj_id = ele.id
        ammeter = curr_meas[CurrentMeasurements.curr_meas_key[(from_node,
                                                               to_node,
                                                               obj_id)]]

    return ammeter
