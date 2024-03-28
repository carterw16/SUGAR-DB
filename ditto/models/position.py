from __future__ import absolute_import, division, print_function

from .base import DiTToHasTraits, Float


class Position(DiTToHasTraits):
    long = Float(help="""Decimal Longitude""")
    lat = Float(help="""Decimal Latitude""")
    elevation = Float(help="""Decimal elevation (meters)""")

    def build(self, model):
        self._model = model
