# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 17:53:27 2016

@author: Amrit
"""
from enum import Enum
from itertools import count

MVAbase = 100
f = 60
# Node Index

# Object Parser
object_parsed = set()
object_all = set()

# connected_nodes
connected_nodes = set()
_nodes_key = dict()


class PHASES(Enum):
    A = 0x1
    B = 0x2
    C = 0x4
    D = 0x7
    S = 0x20
    N = 0x8
    G = 0x10


class CONNECT_TYPE(Enum):
    WYE = 0
    DELTA = 1


# **** DO - NOT CHANGE THIS (LINE CONFIG USES IT)
_PHASE = dict()
_PHASE['ABC'] = 7
_PHASE['ACB'] = 7
_PHASE['A'] = 1
_PHASE['B'] = 2
_PHASE['C'] = 4
_PHASE['N'] = 8
_PHASE['AN'] = 9
_PHASE['BN'] = 10
_PHASE['CN'] = 12
_PHASE['ABCD'] = 71
_PHASE['ABCN'] = 15
_PHASE['ACBN'] = 15
_PHASE['NABC'] = 15
_PHASE['ACN'] = 13
_PHASE['ABN'] = 11
_PHASE['NAB'] = 11
_PHASE['NBC'] = 14
_PHASE['CBN'] = 14
_PHASE['BCN'] = 14
_PHASE['CB'] = 6
_PHASE['AB'] = 3
_PHASE['AC'] = 5
_PHASE['CA'] = 5
_PHASE['BC'] = 6
_PHASE['D'] = 0x40
_PHASE['AS'] = 33
_PHASE['BS'] = 34
_PHASE['CS'] = 36

_BUSTYPE = dict()
_BUSTYPE['PQ'] = 1
_BUSTYPE['PV'] = 2
_BUSTYPE['SWING'] = 3

_XFMR = dict()
_XFMR['UNKNOWN'] = 0
_XFMR['WYE_WYE'] = 1
_XFMR['GWYE_GWYE'] = 1
_XFMR['DELTA_DELTA'] = 2
_XFMR['DELTA_GWYE'] = 3
_XFMR['DELTA_WYE'] = 3
_XFMR['WYE_DELTA'] = 8
_XFMR['SINGLE_PHASE'] = 4
_XFMR['SINGLE_PHASE_CENTER_TAPPED'] = 5

_INSTALL_TYPE = dict()
_INSTALL_TYPE['UNKNOWN'] = 0
_INSTALL_TYPE['POLETOP'] = 1
_INSTALL_TYPE['PADMOUNT'] = 2
_INSTALL_TYPE['VAULT'] = 3

UN_PHASE = dict()
UN_PHASE[7] = (1, 2, 4)
UN_PHASE[9] = (1, 8)
UN_PHASE[10] = (2, 8)
UN_PHASE[12] = (4, 8)
UN_PHASE[11] = (1, 2, 8)
UN_PHASE[13] = (1, 4, 8)
UN_PHASE[14] = (2, 4, 8)
UN_PHASE[15] = (1, 2, 4, 8)
UN_PHASE[33] = (1, 0x20)
UN_PHASE[34] = (2, 0x20)
UN_PHASE[36] = (4, 0x20)
UN_PHASE[65] = (1, 0x40)
