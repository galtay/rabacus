""" Constants and functions relevant for helium atoms. """ 

import numpy as np

from rabacus.constants import physical
from rabacus.constants import units
import photo_xsection

__all__ = ['Helium']


class Helium:
    r""" Helium atom class. 

    Attributes:

      `U` (:class:`~rabacus.constants.units.Units`) 

      `PC` (:class:`~rabacus.constants.physical.PhysicalConstants`)
 
      `PX` (:class:`~rabacus.atomic.photo_xsection.PhotoXsections`) 


    """

    def __init__(self, px_fit_type='verner'):

        # attach input
        #--------------------------------------------------------
        self.px_fit_type = px_fit_type

        # create physical constants and units 
        #--------------------------------------------------------
        self.PC = physical.PhysicalConstants()
        self.U = units.Units()
        self.PX = photo_xsection.PhotoXsections( fit_type = px_fit_type )


        # Express in temperature units
        #--------------------------------------------------------
        self.T_He1 = self.PX.Eth_He1 / self.PC.kb
        self.T_He1.units = 'K'

        self.T_He2 = self.PX.Eth_He2 / self.PC.kb
        self.T_He2.units = 'K'





