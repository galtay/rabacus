""" A module for loading monochromatic spectra.  This module is only here to 
conform to the style of the other spectral shapes.  Normalization is handled 
through classes derived from :class:`~rabacus.rad_src.source.Source`. 
""" 

import numpy as np

__all__ = ['MonochromaticSpectrum']


class MonochromaticSpectrum:

    """ Provides a monochromatic shape. 

    Args: 
      `source`: a class derived from
      (:class:`~rabacus.rad_src.source.Source`)  

    .. seealso::

      :class:`~rabacus.rad_src.powerlaw.PowerlawSpectrum`, 
      :class:`~rabacus.rad_src.thermal.ThermalSpectrum`, 
      :class:`~rabacus.rad_src.hm12.HM12Spectrum`

    """

    def __init__( self, source ): 

        """ 
        """

        self.yvals = 1.0











