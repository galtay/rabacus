""" A general Jeans scale module. """


import numpy as np

from rabacus.constants import physical
from rabacus.constants import units


__all__ = ['Jeans']


class Jeans:
    r""" A Jeans scale class. 

    Provides access to Jeans scales functions.  Default values for Yp and fg 
    are taken from Planck Cosmological Parameters

    Yp = 0.248
    fg = Omega_b / Omega_m = 0.154

    Args: 

    Kwargs:
      
      `Yp` (float): helium mass fraction
      
      `fg` (float): gas fraction

      `gamma` (float): ratio of specific heats


    """ 

    def __init__( self, Yp=0.248, fg=0.154, gamma=5.0/3.0 ):

        self.Yp = Yp
        self.fg = fg
        self.gamma = gamma

        self.U = units.Units()
        self.PC = physical.PhysicalConstants()


    def mu( self, xH2, xHe2, xHe3 ):
        """ Mean molecular weight, i.e. the mean mass of an ion in atomic
        mass units, as a function of ionization state. 
        
        Args:

          `xH2` (float): nHII / nH
          
          `xHe2` (float): nHeII / nHe
          
          `xHe3` (float): nHeIII / nHe

        """         
        num = 4.0
        t1 = 3.0 * self.Yp
        t2 = 4.0 * xH2 * (1.0-self.Yp)
        t3 = (xHe2 + 2.0*xHe3) * self.Yp
        den = 4.0 - t1 + t2 + t3 
        mu = num / den
        return mu


    def t_dyn( self, nH ):
        """ Dynamical time 

        Args:

          `nH` (float): number density of hydrogen

        """ 
        num = ( 1.0-self.Yp ) * self.fg
        den = nH * self.PC.G * self.PC.m_p 
        t_dyn = np.sqrt( num / den )
        t_dyn.units = 's'
        return t_dyn


    def cs( self, T, mu ):
        """ Sound speed 

        Args:

          `T` (float): temperature 

          `mu` (float): mean molecular weight
          
        """ 
        cs = np.sqrt( self.gamma * self.PC.kb * T / ( mu * self.PC.m_p ) )
        cs.units = 'km/s'
        return cs


    def t_sc( self, L, T, mu ): 
        """ Sound crossing time. 

        Args:

          `L` (float): length scale
          
          `T` (float): temperature

          `mu` (float): mean molecular weight
          
        """ 
        cs = self.cs( T, mu )
        t_sc = L / cs
        t_sc.units = 's'
        return t_sc


    def L( self, nH, T, mu ):
        """ Jeans length 

        Args:

          `nH` (float): hydrogen number density
          
          `T` (float): temperature

          `mu` (float): mean molecular weight

        """ 
        cs = self.cs( T, mu )
        t_dyn = self.t_dyn( nH )
        L = cs * t_dyn
        L.units = 'kpc'
        return L


    def NH( self, nH, T, mu ):
        """ Jeans column density 

        Args:

          `nH` (float): hydrogen number density
          
          `T` (float): temperature

          `mu` (float): mean molecular weight

        """ 
        NH = nH * self.L( nH, T, mu )
        NH.units = 'cm^-2'
        return NH
