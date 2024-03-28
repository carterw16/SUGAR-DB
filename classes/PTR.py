from enum import Enum

from classes.Fuses import Fuses
from classes.GlobalVars import object_all, object_parsed
from classes.lines.base import Lines
from classes.Regulators import FixedTapRegulator, VariableRegulator
from classes.Switches import Switches
from classes.Transformers import Transformer


def type_condition(obj):
    obj_type = ""
    if isinstance(obj, Lines):
        if obj.type == 'overhead':
            obj_type = "Overhead Line"
        elif obj.type == 'underground':
            obj_type = 'Underground Line'
        else:
            obj_type = 'Triplex Line'
    elif isinstance(obj, Fuses):
        obj_type = "Fuse"
    elif isinstance(obj, Switches):
        obj_type = "Switch"
    elif isinstance(obj, FixedTapRegulator):
        obj_type = "Regulator"
    elif isinstance(obj, VariableRegulator):
        obj_type = "Regulator"
    elif isinstance(obj, Transformer):
        if obj.xfmr_type == 'Center Tap':
            obj_type = "XfmrCT"
        else:
            obj_type = "Xfmr"
    return obj_type


class OBJECT_TYPE(Enum):
    NODE = 0
    LINK = 1
    LINE = 2
    LINE_CONFIGURATION = 3
    LINE_SPACING = 4
    OVERHEAD_LINE = 5
    OVERHEAD_LINE_CONDUCTOR = 6
    UNDERGROUND_LINE = 7
    UNDERGROUND_LINE_CONDUCTOR = 8
    TRIPLEX_LINE = 9
    TRIPLEX_LINE_CONFIGURATION = 10
    TRIPLEX_LINE_CONDUCTOR = 11
    TRANSFORMER = 12
    TRANSFORMER_CONFIGURATION = 13
    LOAD = 14
    METER = 15
    TRIPLEX_NODE = 16
    TRIPLEX_METER = 17
    TRIPLEX_LOAD = 18
    REGULATOR = 19
    REGULATOR_CONFIGURATION = 20
    CAPACITOR = 21
    FUSE = 22
    SWITCH = 23
    RECLOSER = 24
    RELAY = 25
    SUBSTATION = 26
    UNKNOWN = -1
    UNKNOWN_BUG = 40


class PTR(object):
    # Object Dictionary Type
    object_name = dict()
    object_type = dict()

    @staticmethod
    def findindex(tuplelist, match):
        # type: (object, object) -> object
        for (idx, element) in enumerate(tuplelist):
            if match == element:
                return idx

    @staticmethod
    def check_object():
        PTR.object_not_parsed = object_all - object_parsed
        for ele in PTR.object_not_parsed:
            print('The following object %s with object type %s was not parsed' %
                  (PTR.object_name[ele], PTR.object_type[ele]))

    @staticmethod
    def parse_ptrs(casedata):
        # Parse all the ptr's
        ptr_data = casedata['Object']
        ptr_idx = dict()
        ptr_attr_list = ptr_data.dtype.names
        for ptr_attr in ptr_attr_list:
            ptr_idx[ptr_attr] = PTR.findindex(ptr_data.dtype.names, ptr_attr)
        ptr = []
        for ele in range(len(ptr_data)):
            ptr.append(
                PTR(ptr_data[ele][0][ptr_idx['name']],
                    ptr_data[ele][0][ptr_idx['ptr']],
                    ptr_data[ele][0][ptr_idx['parent']],
                    ptr_data[ele][0][ptr_idx['type']]))
        return ptr

    def __init__(self, name, ID, parent, type):
        self.name = name
        self.ID = int(ID)
        self.parent = int(parent) if parent else None
        self.type = OBJECT_TYPE(int(type))
        # Add all objects to a set
        if not self.parent:
            object_all.add(self.ID)
        else:
            object_all.add(self.parent)
        # Create object dictionaries
        PTR.object_name[self.ID] = self.name
        PTR.object_type[self.ID] = self.type
