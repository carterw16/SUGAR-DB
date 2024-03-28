"""
Description : Class for the nodes in three-phase power flow
Date: 04/10/2017"""

from __future__ import division

import operator
from itertools import count

import numpy as np

bus_id_ = count(0)


def calcVrange(nodes, Vsol):
    Va_min = []
    Vb_min = []
    Vc_min = []
    # calculate the min and max angles and voltages on all nodes
    for nodeEle in nodes:
        nodeEle.calcMagAng(Vsol, True)

    for nodeEle in nodes:
        Vmin_idx, Va = min(enumerate([nodeEle.V_mag[0]]))
        Vmin_idx, Vb = min(enumerate([nodeEle.V_mag[1]]))
        try:
            Vmin_idx, Vc = min(enumerate([nodeEle.V_mag[2]]))
        except IndexError:
            Vc = 0
        Va_min.append(Va)
        Vb_min.append(Vb)
        Vc_min.append(Vc)

    Va_min = np.asarray(Va_min)

    Vmax_node, Vmax = max(enumerate([nodeEle.V_mag for nodeEle in nodes]),
                          key=operator.itemgetter(1))
    Vmin = min(Va_min[np.nonzero(Va_min)])
    Vmax_ang_node, Vmax_ang = max(enumerate(
        [nodeEle.V_ang for nodeEle in nodes]),
                                  key=operator.itemgetter(1))
    Vmin_ang_node, Vmin_ang = min(enumerate(
        [nodeEle.V_ang for nodeEle in nodes]),
                                  key=operator.itemgetter(1))
    return Vmax_node, max(Vmax), Vmin_ang_node, Vmin, Vmax_ang_node, max(
        Vmax_ang), Vmin_ang_node, min(Vmin_ang)


class Nodes:
    nodeKey = {}
    ExistTriplex = False
    voltage_index = []
    Vr_index = []
    Vi_index = []
    Vr_mag = []
    Vi_mag = []
    V_mag = []
    Q_index = []
    I_index = []
    Tap_index = []
    Q_mag = []

    @staticmethod
    def almostEqual(d1, d2, epsilon=10**-7):
        return abs(d2 - d1) < epsilon

    @staticmethod
    def _totuple(a):
        try:
            return tuple(i for i in a)
        except TypeError:
            return a

    @staticmethod
    def _convertline2phase(Vab, Vbc, Vca):
        _Veq = np.array([[1, -1, 0], [0, 1, -1], [-1, 0, 1]])
        _VLL = np.array([Vab, Vbc, Vca])
        _Vabc = np.linalg.solve(_Veq, _VLL)
        return Nodes._totuple(_Vabc)

    @staticmethod
    def sind(x):
        return np.sin(float(x) * np.pi / 180)

    @staticmethod
    def cosd(x):
        return np.cos(float(x) * np.pi / 180)

    def update_Vnom(self, Vsol):

        try:
            Va = abs(Vsol[self.nodeA_Vr] + 1j * Vsol[self.nodeA_Vi])
            err_nom = abs(Va - self.Vnom / np.sqrt(3))
            Vb = abs(Vsol[self.nodeB_Vr] + 1j * Vsol[self.nodeB_Vi])
            err_nom = abs(Vb - self.Vnom / np.sqrt(3))
        except AttributeError:
            Vc = abs(Vsol[self.nodeC_Vr] + 1j * Vsol[self.nodeC_Vi])
            err_nom = abs(Vc - self.Vnom / np.sqrt(3))

        try:
            self.Vnom = self.Vnom / np.sqrt(3)
        except err_nom > 0.10 * self.Vnom:
            self.Vnom = self.Vnom

    def __init__(self,
                 node_index_,
                 ID,
                 name,
                 phases,
                 nominal_voltage,
                 voltage_setpoint,
                 node_connect_type,
                 isTriplex=False,
                 bustype=0,
                 maxvoltage_error=1e-2,
                 busflags=0,
                 referencebus=None,
                 meanrepairtime=None,
                 degree=None,
                 stamp_dual=0,
                 voltage_type='LL'):
        self.ID = ID
        self.name = name
        self.phases = int(phases)
        self.voltage_type = voltage_type
        self.Vbase_LN = nominal_voltage * 1e-3 if voltage_type != 'LL' else (
            nominal_voltage * 1e-3) / np.sqrt(3)
        if voltage_type == 'LL':
            self.Vbase_LL = nominal_voltage * 1e-3 
        else:
            if isTriplex:
                self.Vbase_LL = nominal_voltage * 1e-3
            else:
                self.Vbase_LL = nominal_voltage * np.sqrt(3) * 1e-3
        self.Vnom = float(nominal_voltage)
        self.Vset = float(voltage_setpoint)
        # Unknown - 0 | PQ - 1 | PV -2 | Swing - 3
        self.bustype = int(bustype)
        self.stamp_dual = stamp_dual
        # Seperate triplex nodes from nodes (Assign ABCN to nodes)
        self.isTriplex = isTriplex
        Nodes.ExistTriplex = self.isTriplex
        self.print_slack = False
        self.V_unb = 0  # Voltage Unbalance
        # Check the node connect type
        # 0 - Delta 1 - Wye 2  - Triplex -1 - Unknown
        self.node_connect_type = node_connect_type
        # Check the degree of the nodes
        # self.degree = int(degree)
        if not self.isTriplex:
            self.Va = (self.Vnom * self.Vset) + 0 * 1j
            self.Vb = (self.Vnom * self.Vset) * (self.cosd(-120) +
                                                 self.sind(-120) * 1j)
            self.Vc = (self.Vnom * self.Vset) * (self.cosd(120) +
                                                 self.sind(120) * 1j)
            self.Vn = complex(0, 0)

            if self.node_connect_type == 0:
                self.Vab = self.Va - self.Vb
                self.Vbc = self.Vb - self.Vc  # [volts]
                self.Vca = self.Vc - self.Va
                self.Vn = complex(0, 0)

            self.assign_nodes(node_index_)

            # Collect the indices of the voltage vector
            Nodes.voltage_index.extend([
                self.nodeA_Vr, self.nodeA_Vi, self.nodeB_Vr, self.nodeB_Vi,
                self.nodeC_Vr, self.nodeC_Vi, self.nodeN_Vr, self.nodeN_Vi
            ])

            # Make sure the order is the same for the indices and normalizing magnitude
            Nodes.Vr_index.extend(
                [self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr])
            Nodes.Vi_index.extend(
                [self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi])
            Nodes.Vr_mag.extend(4 * [np.absolute(self.Va)])
            Nodes.Vi_mag.extend(4 * [np.absolute(self.Va)])

            Nodes.V_mag.extend(8 * [np.absolute(self.Va)])

        else:
            # Define V1, V2 and Vn [In triplex node Va,
            # Vb, Vc corresponds to V1, V2 and VN]
            if self.phases & 0x01 == 1:
                self.V1 = self.Vnom + 0 * 1j
            elif self.phases & 0x02 == 2:
                self.V1 = self.Vnom * self.cosd(120) + self.Vnom * self.sind(
                    120) * 1j
            elif self.phases & 0x04 == 4:
                self.V1 = self.Vnom * self.cosd(240) + self.Vnom * self.sind(
                    240) * 1j
            # Initialize hot wire 2
            self.V2 = -self.V1  # Is the negative of V1
            # Initialize neutral wire
            self.Vn = complex(0, 0)
            # Assign indices initially for 1, 2, N
            self.assign_triplex_nodes(node_index_)
            # Collect the indices of the voltage vector
            Nodes.voltage_index.extend([
                self.node1_Vr, self.node1_Vi, self.node2_Vr, self.node2_Vi,
                self.nodeN_Vr, self.nodeN_Vi
            ])
            self.node_set = [
                self.node1_Vr, self.node2_Vr, self.nodeN_Vr, self.node1_Vi,
                self.node2_Vi, self.nodeN_Vi
            ]
            Nodes.Vr_index.extend([self.node1_Vr, self.node2_Vr, self.nodeN_Vr])
            Nodes.Vi_index.extend([self.node1_Vi, self.node2_Vi, self.nodeN_Vi])
            Nodes.Vr_mag.extend(3 * [np.absolute(self.V1)])
            Nodes.Vi_mag.extend(3 * [np.absolute(self.V1)])
            Nodes.V_mag.extend(6 * [np.absolute(self.V1)])

        # Additional Parameters
        self.Ia = 0
        self.Ib = 0
        self.Ic = 0
        self.maxvoltage_error = maxvoltage_error
        self.busflags = int(
            busflags) if busflags else None  # Represents a PV bus in future
        self.referencebus = referencebus
        self.meanrepairtime = meanrepairtime
        # Bus id
        Nodes.nodeKey[self.ID] = bus_id_.__next__()

    def __str__(self):
        return 'Node name is :' + str(self.name)

    def assign_nodes(self, node_index_):
        # Assign indices initially for A, B, C, N
        self.nodeA_Vr = node_index_.__next__()
        self.nodeA_Vi = node_index_.__next__()
        self.nodeB_Vr = node_index_.__next__()
        self.nodeB_Vi = node_index_.__next__()
        self.nodeC_Vr = node_index_.__next__()
        self.nodeC_Vi = node_index_.__next__()
        self.nodeN_Vr = node_index_.__next__()
        self.nodeN_Vi = node_index_.__next__()

        # create dual voltage nodes if performing an optimization
        if self.stamp_dual:
            self.nodeA_Lr = node_index_.__next__()
            self.nodeA_Li = node_index_.__next__()
            self.nodeB_Lr = node_index_.__next__()
            self.nodeB_Li = node_index_.__next__()
            self.nodeC_Lr = node_index_.__next__()
            self.nodeC_Li = node_index_.__next__()
            self.nodeN_Lr = node_index_.__next__()
            self.nodeN_Li = node_index_.__next__()

        self.node_set = [
            self.nodeA_Vr, self.nodeB_Vr, self.nodeC_Vr, self.nodeN_Vr,
            self.nodeA_Vi, self.nodeB_Vi, self.nodeC_Vi, self.nodeN_Vi
        ]

        # Check for PV nodes
        if self.bustype == 2:
            self.nodeA_Veq = node_index_.__next__()
            self.nodeB_Veq = node_index_.__next__()
            self.nodeC_Veq = node_index_.__next__()
        else:
            self.nodeA_Veq = None
            self.nodeB_Veq = None
            self.nodeC_Veq = None
        # Add extra equations for slack bus
        if self.bustype == 3:
            self.slack_nodeA_Vr = node_index_.__next__()
            self.slack_nodeA_Vi = node_index_.__next__()
            self.slack_nodeB_Vr = node_index_.__next__()
            self.slack_nodeB_Vi = node_index_.__next__()
            self.slack_nodeC_Vr = node_index_.__next__()
            self.slack_nodeC_Vi = node_index_.__next__()

            # Create dual slack nodes if performing optimization
            if self.stamp_dual:
                self.slack_nodeA_Lr = node_index_.__next__()
                self.slack_nodeA_Li = node_index_.__next__()
                self.slack_nodeB_Lr = node_index_.__next__()
                self.slack_nodeB_Li = node_index_.__next__()
                self.slack_nodeC_Lr = node_index_.__next__()
                self.slack_nodeC_Li = node_index_.__next__()
            # Nodes.slackBus = self.ID

    def assign_triplex_nodes(self, node_index_):
        self.node1_Vr = node_index_.__next__()
        self.node1_Vi = node_index_.__next__()
        self.node2_Vr = node_index_.__next__()
        self.node2_Vi = node_index_.__next__()
        self.nodeN_Vr = node_index_.__next__()
        self.nodeN_Vi = node_index_.__next__()


        if self.stamp_dual:
            self.node1_Lr = node_index_.__next__()
            self.node1_Li = node_index_.__next__()
            self.node2_Lr = node_index_.__next__()
            self.node2_Li = node_index_.__next__()
            self.nodeN_Lr = node_index_.__next__()
            self.nodeN_Li = node_index_.__next__()   
            
            
            
    def calcMagAng(self, V, flag_Vrange=False):
        if not self.isTriplex:
            # Line to ground voltage
            _Vag = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi])
            _Vbg = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi])
            _Vcg = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi])
            # Line to neutral voltage
            _Van = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi]) - (
                V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
            _Vbn = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi]) - (
                V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
            _Vcn = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi]) - (
                V[self.nodeN_Vr] + 1j * V[self.nodeN_Vi])
            # Line to Line voltage
            _Vab = (V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi]) - (
                V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi])
            _Vbc = (V[self.nodeB_Vr] + 1j * V[self.nodeB_Vi]) - (
                V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi])
            _Vca = (V[self.nodeC_Vr] + 1j * V[self.nodeC_Vi]) - (
                V[self.nodeA_Vr] + 1j * V[self.nodeA_Vi])

            # Mag and Ang Line to Gnd
            Vag_mag = float(np.absolute(_Vag))
            Vbg_mag = float(np.absolute(_Vbg))
            Vcg_mag = float(np.absolute(_Vcg))
            Vag_ang = float(np.angle(_Vag, deg=True)) if Vag_mag != 0 else 0
            Vbg_ang = float(np.angle(_Vbg, deg=True)) if Vbg_mag != 0 else 0
            Vcg_ang = float(np.angle(_Vcg, deg=True)) if Vcg_mag != 0 else 0
            Vag_ang_rad = float(np.angle(_Vag,
                                         deg=False)) if Vag_mag != 0 else 0
            Vbg_ang_rad = float(np.angle(_Vbg,
                                         deg=False)) if Vbg_mag != 0 else 0
            Vcg_ang_rad = float(np.angle(_Vcg,
                                         deg=False)) if Vcg_mag != 0 else 0
            # Mag and Ang Line to Line
            Vab_mag = float(np.absolute(_Vab))
            Vbc_mag = float(np.absolute(_Vbc))
            Vca_mag = float(np.absolute(_Vca))
            Vab_ang = float(np.angle(_Vab, deg=True))
            Vbc_ang = float(np.angle(_Vbc, deg=True))
            Vca_ang = float(np.angle(_Vca, deg=True))
            #
            if self.phases & 0x08 == int(0x08):
                Vab_mag_pu = Vab_mag / (self.Vbase_LL * 1e3)
                Vbc_mag_pu = Vbc_mag / (self.Vbase_LL * 1e3)
                Vca_mag_pu = Vca_mag / (self.Vbase_LL * 1e3)
            else:
                Vab_mag_pu = Vab_mag / (self.Vbase_LL * 1e3)
                Vbc_mag_pu = Vbc_mag / (self.Vbase_LL * 1e3)
                Vca_mag_pu = Vca_mag / (self.Vbase_LL * 1e3)
            # Mag and Ang Line to Neutral
            Van_mag = float(np.absolute(_Van))
            Vbn_mag = float(np.absolute(_Vbn))
            Vcn_mag = float(np.absolute(_Vcn))
            Van_ang = float(np.angle(_Van, deg=True))
            Vbn_ang = float(np.angle(_Vbn, deg=True))
            Vcn_ang = float(np.angle(_Vcn, deg=True))
            Van_mag_pu = Van_mag / (self.Vbase_LN * 1e3)
            Vbn_mag_pu = Vbn_mag / (self.Vbase_LN * 1e3)
            Vcn_mag_pu = Vcn_mag / (self.Vbase_LN * 1e3)
            Vag_mag_pu = Vag_mag / self.Vnom
            Vbg_mag_pu = Vbg_mag / self.Vnom
            Vcg_mag_pu = Vcg_mag / self.Vnom

            # Store results
            # Line to neutral
            self.Vln = (np.around(Van_mag,
                                  2), np.around(Vbn_mag,
                                                2), np.around(Vcn_mag, 2))
            self.Vln_pu = (np.around(Van_mag_pu,
                                     5), np.around(Vbn_mag_pu,
                                                   5), np.around(Vcn_mag_pu, 5))
            self.Vll = (np.around(Vab_mag * 1e-3,
                                  5), np.around(Vbc_mag * 1e-3, 5),
                        np.around(Vca_mag * 1e-3, 5))
            self.Vll_pu = (np.around(Vab_mag_pu,
                                     5), np.around(Vbc_mag_pu,
                                                   5), np.around(Vca_mag_pu, 5))
            self.Vlg = (np.around(Vag_mag * 1e-3,
                                  5), np.around(Vbg_mag * 1e-3, 5),
                        np.around(Vcg_mag * 1e-3, 5))
            self.Vlg_pu = (np.around(Vag_mag_pu,
                                     5), np.around(Vbg_mag_pu,
                                                   5), np.around(Vcg_mag_pu, 5))
            self.Vbase_LL = np.around(self.Vbase_LL, 5)
            self.Vnom = np.around(self.Vnom, 5)

            if not flag_Vrange:
                self.V_mag = self.Vlg
                self.Vll_mag = self.Vll
            else:
                self.V_mag = self.Vlg_pu
                self.Vll_mag = self.Vll_pu

            self.V_ang = (np.around(Vag_ang,
                                    1), np.around(Vbg_ang,
                                                  1), np.around(Vcg_ang, 1))
            self.Vll_ang = (np.around(Vab_ang,
                                      5), np.around(Vbc_ang,
                                                    5), np.around(Vca_ang, 5))
            Vmag = np.asarray(self.V_mag)
            Vmax = np.max(Vmag) if Vmag.any() != 0 else 0
            Vavg = Vmag[Vmag != 0].mean() if Vmag.any() != 0 else 0

            self.V_unb = np.abs(Vmax - Vavg) / Vavg if Vmag.any() != 0 else 0.0
            self.V_unb = np.around(self.V_unb * 100, 3)
            # Calculate swing current
            if self.bustype == 3:
                self.Type = "Slack"
                currentA = -V[self.slack_nodeA_Vr] - 1j * V[self.slack_nodeA_Vi]
                currentB = -V[self.slack_nodeB_Vr] - 1j * V[self.slack_nodeB_Vi]
                currentC = -V[self.slack_nodeC_Vr] - 1j * V[self.slack_nodeC_Vi]

                self.Ia = currentA
                self.Ib = currentB
                self.Ic = currentC

                self.Ia_mag = float(np.absolute(currentA))
                self.Ib_mag = float(np.absolute(currentB))
                self.Ic_mag = float(np.absolute(currentC))
                self.Ia_ang = float(np.angle(currentA, deg=True))
                self.Ib_ang = float(np.angle(currentB, deg=True))
                self.Ic_ang = float(np.angle(currentC, deg=True))
                # Calculate Slack Bus Power
                self.Sa = currentA * self.Va
                self.Sb = currentB * self.Vb
                self.Sc = currentC * self.Vc
                self.S = np.absolute(self.Sa.real) + np.absolute(
                    self.Sb.real) + np.absolute(self.Sc.real) + 1j * (
                        np.absolute(self.Sa.imag) + np.absolute(self.Sb.imag) +
                        np.absolute(self.Sc.imag))
                if self.bustype == 3 and self.print_slack:
                    print('Slack bus power consumption is %f + %fi' %
                          (self.S.real, self.S.imag))
                    print(
                        'Slack Bus Voltages are : (%f, %f, %f, %f, %f, %f)' %
                        (Van_mag, Vbn_mag, Vcn_mag, Van_ang, Vbn_ang, Vcn_ang))
                    print('Slack Bus Currents are : (%f, %f, %f, %f, %f, %f)' %
                          (self.Ia_mag, self.Ib_mag, self.Ic_mag, self.Ia_ang,
                           self.Ib_ang, self.Ic_ang))
        else:
            # Line to ground voltage
            _V1g = (V[self.node1_Vr] + 1j * V[self.node1_Vi])
            _V2g = -(V[self.node2_Vr] + 1j * V[self.node2_Vi])
            _V12 = (V[self.node1_Vr] + 1j * V[self.node1_Vi]) - (
                V[self.node2_Vr] + 1j * V[self.node2_Vi])
            # Mag and Ang Line to Gnd
            V1g_mag = float(np.absolute(_V1g))
            V2g_mag = float(np.absolute(_V2g))
            V12_mag = float(np.absolute(_V12))
            V1g_pu = V1g_mag / self.Vnom
            V2g_pu = V2g_mag / self.Vnom
            V12_pu = V12_mag / self.Vbase_LL
            V1g_ang = float(np.angle(_V1g, deg=True)) if V1g_mag else 0
            V1g_ang_rad = float(np.angle(_V1g, deg=False)) if V1g_mag else 0
            V2g_ang = float(np.angle(_V2g, deg=True)) if V2g_mag else 0
            V2g_ang_rad = float(np.angle(_V2g, deg=False)) if V2g_mag else 0
            V12_ang = float(np.angle(_V12, deg=True)) if V12_mag else float(0)
            self.Vln = (np.around(V1g_mag * 1e-3,
                                  5), np.around(V2g_mag * 1e-3, 5))
            self.Vll = np.around(V12_mag * 1e-3, 5)
            self.Vln_pu = (V1g_pu, V2g_pu)
            self.Vll_pu = np.around(V12_pu * 1e-3, 5)
            self.V_mag = (V1g_mag, V2g_mag)
            self.V_ang = (np.around(V1g_ang, 5), np.around(V2g_ang, 5))
            self.Vll_ang = np.around(V12_ang, 5)
            self.V_unb = 0


class TriplexNodes:

    def __init__(self, name, ptr, bustype, current1, current12, current2,
                 currentN, impedance1, impedance12, impedance2,
                 maximum_voltage_error, mean_repair_time, nominal_voltage,
                 parent, phases, power1, power12, power2, reference_bus, shunt1,
                 shunt12, shunt2, voltage1, voltage12, voltage1N, voltage2,
                 voltage2N, voltageN):
        self.name = name
        self.ptr = int(ptr)
        self.bustype = int(bustype)
        self.parent = int(parent) if parent else None
        self.phases = int(phases)
        self.nominal_voltage = int(nominal_voltage)
        self.voltage1 = (voltage1) if voltage1 else None
        self.voltage2 = (voltage2) if voltage2 else None
        self.voltage12 = (voltage12) if voltage12 else None
        self.voltage1N = (voltage1N) if voltage1N else None
        self.voltage2N = (voltage2N) if voltage2N else None
        self.voltageN = (voltageN) if voltage12 else None
        self.reference_bus = reference_bus
        self.maximum_voltage_error = maximum_voltage_error
        self.mean_repair_time = mean_repair_time
        # Triplex Node Load Parameters
        self.current1 = (current1) if current1 else None
        self.current2 = (current2) if current2 else None
        self.current12 = (current12) if current12 else None
        self.currentN = (currentN) if currentN else None
        self.impedance1 = (impedance1) if impedance1 else None
        self.impedance2 = (impedance2) if impedance2 else None
        self.impedance12 = (impedance12) if impedance12 else None
        self.power1 = (power1) if power1 else None
        self.power2 = (power2) if power2 else None
        self.power12 = (power12) if power12 else None
        self.shunt1 = (shunt1) if shunt1 else None
        self.shunt2 = (shunt2) if shunt2 else None
        self.shunt12 = (shunt12) if shunt12 else None
