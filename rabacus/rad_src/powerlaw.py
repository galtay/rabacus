""" A module for loading powerlaw spectra.  This module is only concerned with 
the shape of the spectrum.  Normalization is handled through classes derived
from :class:`~rabacus.rad_src.source.Source`. """ 

import numpy as np

__all__ = ['PowerlawSpectrum']


class PowerlawSpectrum:

    """ Provides a powerlaw shape. 

    Args: 
      `source`: a class derived from
      (:class:`~rabacus.rad_src.source.Source`)  

      `alpha` (float): powerlaw index 

    .. seealso::

      :class:`~rabacus.rad_src.hm12.HM12Spectrum`,
      :class:`~rabacus.rad_src.thermal.ThermalSpectrum`, 
      :class:`~rabacus.rad_src.monochromatic.MonochromaticSpectrum`

    """

    def __init__( self, source, alpha ): 

        # attach input specific to thermal 
        #-----------------------------------------------------------------
        self.alpha = alpha

        # set shape using alpha
        #-----------------------------------------------------------------
        self.yvals = ( source.nu / source.nu[0] )**alpha
        self.yvals = self.yvals.magnitude

    








