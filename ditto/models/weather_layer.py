from __future__ import absolute_import, division, print_function

from .base import Any, DiTToHasTraits, List, Unicode
from .position import Position


class WeatherLayer(DiTToHasTraits):

    name = Unicode(help="""Name of the weather object""")
    interval = Integer(
        help="""The time resolution (in seconds) for the measured data""")
    ghi = Any(help="""The input data for global horizontal irradiance""")
    temperature = Any(help="""The input data for temperature""")
    relative_humitidy = Any(help="""The input data for relative humidity""")
    surface_windspeed = Any(help="""The input data for surface windspeed""")
    positions = List(
        Instance(Position),
        help=
        """This parameter is a list of positional points describing the weather data. The positions are objects containing elements of long, lat and elevation (See Position object documentation).""",
    )

    def build(self, model):
        self._model = model
