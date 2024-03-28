""" Overcurrent protection device.

  Author(s): Amrit Pandey, Naeem Turner-Bandele
  Created Date: 04-10-2017
  Updated Date: 10-14-2021
  Email: amritanp@andrew.cmu.edu, nturnerb@cmu.edu
  Status: Development

"""

from itertools import count
from types import MethodType

import numpy as np

from .FuseStamp import stamp_linear


class Fuses:
    _GOOD = 1
    _BLOWN = 0
    _ids = count(0)

    def __init__(self, name, ID, from_node, to_node, phases,
                 phase_A_state, phase_B_state, phase_C_state, stamp_dual = False):
        super(Fuses, self).__init__()
        self.name = name
        self.id = self._ids.__next__()
        self.ID = ID
        self.from_node = from_node
        self.to_node = to_node
        self.phases = phases
        self.stamp_dual = stamp_dual
        # Closed : 1 Open : 0
        self.phase_A_state = phase_A_state
        self.phase_B_state = phase_B_state
        self.phase_C_state = phase_C_state

        # # Link Cython Methods to Class # #
        self.stamp_linear = MethodType(stamp_linear, self)

        # Initialize fuse currents #
        self.Ia_mag = 0
        self.Ia_ang = 0
        self.Ib_mag = 0
        self.Ib_ang = 0
        self.Ic_mag = 0
        self.Ic_ang = 0

    def assign_nodes(self, node_index_):
        if self.phases & 0x1 == 1:  # Check for phase A
            self.nodeA_Ir = node_index_.__next__()
            self.nodeA_Ii = node_index_.__next__()

        if self.phases & 0x2 == 2:  # Check for phase B
            self.nodeB_Ir = node_index_.__next__()
            self.nodeB_Ii = node_index_.__next__()
            
        if self.phases & 0x4 == 4:  # Check for phase C
            self.nodeC_Ir = node_index_.__next__()
            self.nodeC_Ii = node_index_.__next__()
        
        if self.stamp_dual:
            if self.phases & 0x1 == 1:  # Check for phase A
                self.nodeA_Lr_fuse = node_index_.__next__()
                self.nodeA_Li_fuse = node_index_.__next__()

            if self.phases & 0x2 == 2:  # Check for phase B
                self.nodeB_Lr_fuse = node_index_.__next__()
                self.nodeB_Li_fuse = node_index_.__next__()

            if self.phases & 0x4 == 4:  # Check for phase C
                self.nodeC_Lr_fuse = node_index_.__next__()
                self.nodeC_Li_fuse = node_index_.__next__()
        
        return node_index_

    def calc_currents(self, V):
        Ia = complex(V[self.nodeA_Ir], V[self.nodeA_Ii]) if hasattr(
            self, 'nodeA_Ir') else complex(0, 0)
        Ib = complex(V[self.nodeB_Ir], V[self.nodeB_Ii]) if hasattr(
            self, 'nodeB_Ir') else complex(0, 0)
        Ic = complex(V[self.nodeC_Ir], V[self.nodeC_Ii]) if hasattr(
            self, 'nodeC_Ir') else complex(0, 0)

        self.Ia_mag = np.around(np.abs(Ia), 2)
        self.Ib_mag = np.around(np.abs(Ib), 2)
        self.Ic_mag = np.around(np.abs(Ic), 2)
        self.Ia_ang = np.around(np.angle(Ia, deg=True), 2)
        self.Ib_ang = np.around(np.angle(Ib, deg=True), 2)
        self.Ic_ang = np.around(np.angle(Ic, deg=True), 2)
