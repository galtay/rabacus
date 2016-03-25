""" A module for loading blackbody spectra.  This module is only concerned with 
the shape of the spectrum.  Normalization is handled through classes derived
from :class:`~rabacus.rad_src.source.Source`. """ 

import numpy as np
from rabacus.utils import utils

__all__ = ['ThermalSpectrum']


class ThermalSpectrum:

    """ Provides a blackbody shape. 

    Args: 
      `source`: a class derived from
      (:class:`~rabacus.rad_src.source.Source`)  

      `T_eff` (float): blackbody temperature

    .. seealso::

      :class:`~rabacus.rad_src.hm12.HM12Spectrum`,
      :class:`~rabacus.rad_src.powerlaw.PowerlawSpectrum`, 
      :class:`~rabacus.rad_src.monochromatic.MonochromaticSpectrum`

    """

    def __init__( self, source, T_eff ): 

        # check input 
        #--------------------------------------------------------
        if hasattr( T_eff, 'units' ):
            T_eff.units = 'K'
        else:
            raise utils.NeedUnitsError('\nT_eff must have units of K\n')


        # attach input specific to thermal 
        #-----------------------------------------------------------------
        self.T_eff = T_eff


        # set shape using T_eff
        #-----------------------------------------------------------------
        t1 = 2 * source.PC.h * source.nu**3 / source.PC.c**2
        e1 = ( source.PC.h * source.nu / (source.PC.kb * T_eff) ).simplified
        t2 = np.exp( e1 )

        # the units of this quantity are the same as Inu
        # erg / ( s cm^2 Hz sr )

        self.yvals = t1 * 1.0 / (t2-1)
        self.yvals = self.yvals.magnitude










