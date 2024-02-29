# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from .base import (Any, Bool, Complex, DiTToHasTraits, Float, Instance, Int,
                   List, Unicode, observe)


class PhaseLoad(DiTToHasTraits):
    phase = Unicode(help="""The phase (A, B, C, N, s1, s2) of the load""",
                    default_value=None)
    p = Float(
        help=
        """The active power of the load which is fixed. Positive values represent flow out of the node.""",
        default_value=None,
    )
    q = Float(
        help=
        """The reactive power of the load which is fixed. Positive values represent flow out of the node.""",
        default_value=None,
    )
    i = Complex(
        help="""The complex current of the load which is fixed.""",
        default_value=complex(0, 0),
    )
    z = Complex(
        help="""The complex impedance of the load which is fixed.""",
        default_value=complex(0, 0),
    )

    # Modification: Nicolas (August 2017)
    # opendss has 8 different load models. Without this field there is no way to capture
    # this information in DiTTo ==> Only constant P&Q and Zipload would be considered
    # Note: use_zip can probably be removed since it is equivalent to model=8
    model = Int(help="""opendss Load Model number.""", default_value=1)

    # Modification: Nicolas Gensollen (December 2017)
    # Drop flag is used if we created objects in the reader that we do not want to output.
    # This is much faster than looping over objects to remove them in a pre/post-processing step
    drop = Bool(
        help=
        """Set to 1 if the object should be dropped in the writing process. Otherwise leave 0.""",
        default_value=False,
    )

    cvrwatts = Float(
        help="This is the active power conservation voltage reduction ratio.",
        default_value=0.0)

    cvrvars = Float(
        help="This is the reactive power conservation voltage reduction ratio.",
        default_value=0.0)

    ppercentcurrent = Float(
        help=
        """This is the portion of active power load modeled as constant current.  Active portions of current, power and impedance should all add to 1. Used for ZIP models.""",
        default_value=None,
    )
    qpercentcurrent = Float(
        help=
        """ This is the portion of active power load modeled as constant impedance. Reactive portions of current, power and impedance should all add to 1. Used for ZIP models.""",
        default_value=None,
    )
    ppercentpower = Float(
        help=
        """This is the portion of active power load modeled as constant power. Active portions of current, power and impedance should all add to 1.  Used for ZIP models."""
    )
    qpercentpower = Float(
        help=
        """This is the portion of reactive power load modeled as constant current. Reactive portions of current, power and impedance should all add to 1. Used for ZIP models."""
    )
    ppercentimpedance = Float(
        help=
        """This is the portion of reactive power load modeled as  Active portions of current, power and impedance should all add to 1. constant impedance. Used for ZIP models.""",
        default_value=None,
    )
    qpercentimpedance = Float(
        help=
        """This is the portion of reactive power load modeled as constant impedance. Reactive portions of current, power and impedance should all add to 1. Used for ZIP models.""",
        default_value=None,
    )

    def build(self, model):
        self._model = model
