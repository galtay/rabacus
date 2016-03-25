""" Handles photo-ionization cross-sections as described in Verner 96
http://adsabs.harvard.edu/abs/1996ApJ...465..487V """

import numpy as np
import os.path

from rabacus.constants import physical
from rabacus.constants import units
from rabacus.utils import utils


__all__ = ['PhotoXsections_Verner96']


class PhotoXsections_Verner96:

    """ Reads V96 photoionization cross-section table 

    Attributes:

      `U` (:class:`~rabacus.constants.units.Units`)
      
      `PC` (:class:`~rabacus.constants.physical.PhysicalConstants`)

      """ 

    def __init__( self ):

        # create physical constants and units 
        #--------------------------------------------------------
        self.U = units.Units()
        self.PC = physical.PhysicalConstants()

        # read data
        #-----------------------------------------------------
        local = os.path.dirname(os.path.realpath(__file__))
        fname = local + '/photo.dat'        
        dat = np.loadtxt( fname, unpack=True )
        self._dat = dat
        
        # organize data
        #-----------------------------------------------------
        self.Z = dat[0]        
        self.N = dat[1]
        self.Eth = dat[2] * self.U.eV
        self.Emax = dat[3] * self.U.eV
        self.E0 = dat[4] * self.U.eV
        self.sigma0 = dat[5] * 1.0e-18 * self.U.cm**2
        self.ya = dat[6]
        self.P = dat[7]
        self.yw = dat[8]
        self.y0 = dat[9]
        self.y1 = dat[10]



    def return_Eth( self, Z, N ):
        """ Returns threshold ionization energy for ions defined by `Z` and 
        `N`.

        Args:

          `Z` (int): atomic number (number of protons)        
          
          `N` (int): electron number (number of electrons)

        Returns:

          `Eth` (float): ionization energy

        """

        # Make sure input is OK
        #--------------------------------------------------------
        if utils.isiterable( Z ) or not isinstance(Z,int):
            raise utils.InputError, '\n Z must be an integer scalar'

        if utils.isiterable( N ) or not isinstance(N,int):
            raise utils.InputError, '\n N must be an integer scalar'

        if Z < 1 or Z > 26:
            raise utils.InputError, '\n We must have 1 <= Z <= 26'

        if N < 1 or N > Z:
            raise utils.InputError, '\n We must have 1 <= N <= Z'


        # calculate 
        #--------------------------------------------------------

        c1 = self.Z == Z
        c2 = self.N == N

        indx = np.where( c1 & c2 )

        indx = indx[0][0]

        Eth = self.Eth[indx]

        Eth.units = 'eV' 
        return Eth



    def return_fit( self, Z, N, E ):
        """ Returns a photo-ionization cross-section for an ion defined by 
        `Z` and `N` at energies `E`.

        Args:

          `Z` (int): atomic number (number of protons)        
          
          `N` (int): electron number (number of electrons)

          `E` (array): calculate cross-section at these energies

        Returns:

          `sigma` (array): photoionization cross-sections


        """


        # Make sure input is OK
        #--------------------------------------------------------

        if hasattr(E,'units'): 
            E.units = 'eV'
        else:
            raise utils.NeedUnitsError, '\n Input variable E must have units \n'

        if utils.isiterable( Z ) or not isinstance(Z,int):
            raise utils.InputError, '\n Z must be an integer scalar'

        if utils.isiterable( N ) or not isinstance(N,int):
            raise utils.InputError, '\n N must be an integer scalar'

        if Z < 1 or Z > 26:
            raise utils.InputError, '\n We must have 1 <= Z <= 26'

        if N < 1 or N > Z:
            raise utils.InputError, '\n We must have 1 <= N <= Z'


        # calculate fit
        #--------------------------------------------------------

        c1 = self.Z == Z
        c2 = self.N == N

        indx = np.where( c1 & c2 )

        indx = indx[0][0]

        Z = self.Z[indx]
        N = self.N[indx]
        Eth = self.Eth[indx]
        Emax = self.Emax[indx]
        E0 = self.E0[indx]
        sigma0 = self.sigma0[indx]
        ya = self.ya[indx]
        P = self.P[indx]
        yw = self.yw[indx]
        y0 = self.y0[indx]
        y1 = self.y1[indx]

        x = E / E0 - y0
        y = np.sqrt( x*x + y1*y1 )

        t1 = (x-1)*(x-1) + yw*yw
        t2 = y**(0.5*P - 5.5)
        t3 = (1+np.sqrt(y/ya))**(-P)
        
        F = t1 * t2 * t3

        sigma = sigma0 * F 
      
        # zero cross-section below threshold
        #--------------------------------------------------------
        if utils.isiterable( E ):
            indx = np.where( E < Eth )
            if indx[0].size > 0:
                sigma[indx] = 0.0 * self.U.cm**2
        else:
            if E < Eth:
                sigma = 0.0 * self.U.cm**2

        return sigma
