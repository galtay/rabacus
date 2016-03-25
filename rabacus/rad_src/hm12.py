""" A module for loading the Haardt and Madau 2012 spectral model 
(see :class:`~rabacus.uv_bgnd.hm12.HM12_UVB_Table`).  
The returned spectrum will be normalized as in the HM12 model.  Further 
adjustments to normalization are handled through classes derived from 
:class:`~rabacus.rad_src.source.Source`. """ 

import numpy as np
from rabacus.uv_bgnd.hm12 import HM12_UVB_Table

__all__ = ['HM12Spectrum']


class HM12Spectrum:

    """ Provides the HM12 spectral shape. 

    Args: 
      `source`: a class derived from
      (:class:`~rabacus.rad_src.source.Source`)  

      `z` (float): redshift

    .. seealso:: 

      :class:`~rabacus.rad_src.powerlaw.PowerlawSpectrum`, 
      :class:`~rabacus.rad_src.thermal.ThermalSpectrum`, 
      :class:`~rabacus.rad_src.monochromatic.MonochromaticSpectrum`


    """

    def __init__( self, source, z ): 

        
        # attach input specific to HM12  
        #-----------------------------------------------------------------
        self.z = z
        self.tab = HM12_UVB_Table()

        Inu = self.tab.return_spectrum_lam( z, source.lam )

        # return the spectral shape
        #-----------------------------------------------------------------    
        self.yvals = Inu.magnitude




