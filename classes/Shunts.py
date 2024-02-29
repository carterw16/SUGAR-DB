# Check 1
from __future__ import division

from itertools import count

from classes.GlobalVars import MVAbase


class Shunts:
    _ids = count(0)
    Vnom = 1.0

    @staticmethod
    def MW_MVAR2pu(G_MW, B_MVAR):
        G_pu = (G_MW / MVAbase) / Shunts.Vnom**2
        B_pu = (B_MVAR / MVAbase) / Shunts.Vnom**2
        return (G_pu, B_pu)

    def __init__(self, Bus, G_MW, B_MVAR):
        self.Bus = Bus
        self.G_MW = G_MW
        self.B_MVAR = B_MVAR
        self.id = self._ids.__next__()
        (self.G_pu, self.B_pu) = Shunts.MW_MVAR2pu(self.G_MW, self.B_MVAR)

    def stamp(self, Vr_node, Vi_node, Ylin):
        Ylin[Vr_node, Vr_node] += self.G_pu
        Ylin[Vi_node, Vi_node] += self.G_pu
        Ylin[Vr_node, Vi_node] -= self.B_pu
        Ylin[Vi_node, Vr_node] += self.B_pu
        return Ylin
