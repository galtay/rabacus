""" Constants and functions relevant for hydrogen atoms. """ 

import numpy as np

from rabacus.constants import physical
from rabacus.constants import units
from rabacus.utils import utils
import photo_xsection

__all__ = ['Hydrogen']


class Hydrogen:
    r""" Hydrogen atom class. 

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
        self.T_H1 = self.PX.Eth_H1 / self.PC.kb
        self.T_H1.units = 'K'

        


    #===================================================================

    def _px_h1_analytic( self, ryd_in ):
        """ Returns the analytic hydrogenic photoionzation cross section
        with Z=1.  If the input is not a scalar float assume it is a 
        numpy array. 
        see http://adsabs.harvard.edu/abs/2009RvMP...81.1405M """

        if not hasattr( self, 'A0' ):
            fac1 = ( 2**9 * np.pi ) / ( 3 * np.exp(1.)**4 )
            fac2 = self.PC.alpha * np.pi * self.PC.Rbohr**2
            self.A0 = fac1 * fac2

        MAGIC = 10.0

        # if we have a single float
        if isinstance(ryd_in, float):
            ryd = ryd_in
            if ryd < 1.0:
                sigma = 0.0
            elif ryd == 1.0:
                sigma = self.A0
            else:
                eps = np.sqrt(ryd - 1)
                ieps = 1./eps
                efac = 4. - ( 4. * np.arctan(eps) ) * ieps
                num = np.exp(efac)
                den = 1. - np.exp(-2. * np.pi * ieps)        
                sigma = self.A0 * (ryd)**(-4) * num / den
            return sigma

        # if we have a numpy array
        else:
            ryd = ryd_in.copy()

            ilo = np.where( ryd < 1.0 )[0]
            if ilo.size > 0: ryd[ilo] = MAGIC   # anything > 0
        
            ione = np.where( ryd == 1.0 )[0]
            if ione.size > 0: ryd[ione] = 1.0 + 1.0e-10  

            eps = np.sqrt(ryd - 1)
            ieps = 1./eps
            efac = 4. - ( 4. * np.arctan(eps) ) * ieps
            num = np.exp(efac)
            den = 1. - np.exp(-2. * np.pi * ieps)        
            sigma = self.A0 * (ryd)**(-4) * num / den

            if ilo.size > 0: sigma[ilo] = 0.0
            if ione.size > 1: sigma[ione] = self.A0

            return sigma



    def _px_h1_powerlaw( self, ryd_in ):
        """ Returns the hydrogen I  photoionzation cross section.  
        If the input is not a scalar float assume it is a numpy array. 
        simple powerlaw fit """

        MAGIC = 10.0

        # if we have a single float
        if isinstance(ryd_in, float):
            ryd = ryd_in
            if ryd < 1.0:
                sigma = 0.0
            else:
                sigma = self.A0 * ryd**-3
            return sigma 

        # if we have a numpy array
        else:
            ryd = ryd_in.copy()
            ilo = np.where( ryd < 1.0 )[0]
            if ilo.size > 0: ryd[ilo] = MAGIC   # anything > 0
            sigma = self.A0 * ryd**-3
            if ilo.size > 0: sigma[ilo] = 0.0
            return sigma 



    #===================================================================

    def analytic_soltn_xH1(self, nH, y, k):
        """ Analytic equilibrium solution for the neutral fraction xH1 = nH1/nH 
        
        Args:

          `nH` (array): H number density

          `y` (array): electrons from elements heavier than hydrogen, 
          ne = nH * (xH2 + y)

          `k` (rates object): must contain the attributes ciH1, reH2, and H1i, 
          see :class:`~rabacus.atomic.chemistry.ChemistryRates`

        Returns: 

          `xH1` (array): neutral hydrogen fraction
           
        """
 
        RR = (k.ciH1 + k.reH2) * nH
        QQ = -( k.H1i + k.reH2 * nH + RR * (1+y) )
        PP = k.reH2 * nH * (1.0 + y) 
        
        dd = QQ*QQ - 4.0*RR*PP
 
        if utils.isiterable( dd ):
            if np.any( dd < 0.0 ):
                indx = np.where( dd < 0.0 )
                dd[indx] = dd[indx] * 0.0
        else:
            if dd < 0.0:
                dd = dd * 0.0 

        # QQ is always negative in this case
        #q = -0.5 * (QQ + np.sign(QQ) * np.sqrt(QQ*QQ - 4.0*RR*PP))
        q = -0.5 * ( QQ - np.sqrt(dd) )
        xH1 = PP / q
        
        return xH1


    #===================================================================

    def analytic_slab_soltn_tau(self, nH, k, xH1):
        """ Analytic inverse solution for monochromatic radiation incident onto 
        a constant density and temperature slab.  Given a neutral fraction this
        function returns the optical depth into the slab. 

        Args: 

          `nH` (array): H number density

          `k` (rates object): must contain the attributes ciH1, reH2, and H1i, 
          see :class:`~rabacus.atomic.chemistry.ChemistryRates`

          `xH1` (array): neutral hydrogen fraction

        Returns:

          `tauH` (array): depth into the slab for given neutral fraction, tau = 
          NH * sigma

        """

        y = 0.0

        xH1_0 = self.analytic_soltn_xH1( nH, y, k )
        xH1_c = k.reH2 / ( k.reH2 + k.ciH1 )

        t1 = ( 1.0 / xH1_0 - 1.0 / xH1 )

        num = xH1   * (1.0 - xH1_0)
        den = xH1_0 * (1.0 - xH1  ) 
        t2 = np.log( num / den )

        num = xH1   * (xH1_c - xH1_0)
        den = xH1_0 * (xH1_c - xH1  ) 
        t3 = (1.0 / xH1_c) * np.log( num / den )
                
        tauH = t1 + t2 + t3

        return tauH


    #===================================================================

    def analytic_soltn_H1i(self, nH, y, k, xH1):
        """ Analytic equilibrium solution for the photoionization rate given 
        the neutral fraction xH1 = nH1/nH 
        
        Args:

          `nH` (array): H number density

          `y` (array): electrons from elements heavier than hydrogen, 
          ne = nH * (xH2 + y)

          `k` (rates object): must contain the attributes ciH1 and reH2, 
          see :class:`~rabacus.atomic.chemistry.ChemistryRates`

          `xH1` (array): neutral hydrogen fraction


        Returns:

          `H1i` (array): hydrogen photoionization rate
           
        """
 
        t1 = k.reH2 * nH * ( 1.0-xH1+y ) * ( 1.0-xH1 ) / xH1
        t2 = k.ciH1 * nH * ( 1.0-xH1+y ) 

        H1i = t1 - t2

        return H1i


    #===================================================================

