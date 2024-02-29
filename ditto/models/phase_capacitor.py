from __future__ import absolute_import, division, print_function

from .base import DiTToHasTraits, Float, Int, Unicode


class PhaseCapacitor(DiTToHasTraits):

    phase = Unicode(
        help="""The phase (A, B, C, N, s1, s2) of the capacitor section""",
        default_value=None,
    )
    var = Float(help="""The var rating of the capacitor phase""",
                default_value=None)
    switch = Int(
        help=
        """A boolean value to denote whether or not the capacitor is switched in. 1 means it's switched in, 0 means that it's not""",
        default_value=None,
    )
    sections = Int(
        help="""The maximum number of sections connected to this phase""",
        default_value=None,
    )
    normalsections = Int(
        help="""The normal number of sections connected to this phase""",
        default_value=None,
    )

    def build(self, model):
        self._model = model
