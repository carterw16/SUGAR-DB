# -*- coding: utf-8 -*-
"""
Created on Fri Dec 08 10:36:44 2017

@author: Amrit
"""


class Links(object):

    def __init__(self, ID, name, phases, nominal_voltage, from_bus, to_bus,
                 parent, status):
        self.ID = ID
        self.name = name
        self.phases = int(phases) if phases else None
        self.Vnom = float(nominal_voltage) if nominal_voltage else None
        self.from_bus = from_bus
        self.to_bus = to_bus
        self.parent = int(parent) if parent else None
        self.status = status if status else None
