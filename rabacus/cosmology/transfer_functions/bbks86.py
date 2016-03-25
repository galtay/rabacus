""" The transfer function of BBKS86,
http://adsabs.harvard.edu/abs/1986ApJ...304...15B """

import numpy as np
#import matplotlib.pyplot as plt

#from ...constants import units as units
from rabacus.constants import units 

__all__ = ['TransferBBKS']





class TransferBBKS:

    """ Transfer function of BBKS86,
    http://adsabs.harvard.edu/abs/1986ApJ...304...15B

    This transfer function was designed for situations in which the baryon 
    density is much smaller than the cold dark matter density.  In symbols, 
    OmegaB << OmegaC.  Note that all quantities returned by this class
    are arbitrarily normalized.  

    Args: 

      `cpdict` (dict) A dictionary of cosmological parameters.  
      For example, see
      :class:`~rabacus.cosmology.parameters.planck.load.PlanckParameters`.


    The dictionray `cpdict` must include the following keys, 
       - ``omegac`` -> current CDM density in units of critical today
       - ``h``      -> Hubble parameter H0 = 100 h km/s/Mpc 
       - ``ns``     -> slope of primordial power spectrum 



    """ 

    def __init__( self, cpdict ): 

        # check input dictionary
        #--------------------------------------------------------
        need_keys = ['omegac', 'h', 'ns']
        for key in need_keys:
            if not cpdict.has_key( key ):
                raise InitError, ' key missing from cpdict \n' 

        # transfer variables
        #--------------------------------------------------------
        self.OmegaC = cpdict['omegac']
        self.h = cpdict['h']
        self.ns = cpdict['ns']
        self.Gamma = self.OmegaC * self.h
        self.cu = units.CosmoUnits(self.h, 1.0)


    def T_cdm( self, k ):
        """ CDM transfer function (Eq. G3). 

        Args: 

          `k` (real or array): wavenumber with units of inverse length

        Returns:
        
          `T_cdm` (real or array): cold dark matter transfer function at
          the scales `k`.  Normalization is arbitrary. 

        """

        # make sure we have k in 1/Mpch
        #----------------------------------------------
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'
        
        # construct the unitless quantity q
        #----------------------------------------------------------------
        q = ( k / (self.cu.h/self.cu.Mpc) ) / self.Gamma

        # calculate fit
        #----------------------------------------------------------------
        t1 = np.log( 1.0 + 2.34 * q )
        t2 = 2.34 * q
        t3 = 1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4    
        T = t1 / t2 * t3**(-1.0/4)

        return T


    def T_b( self, k ):
        """ Baryon transfer function (Eq. G4). 

        Args: 

          `k` (real or array): wavenumber with units of inverse length

        Returns:
        
          `T_b` (real or array): baryon transfer function at
          the scales `k`.  Normalization is arbitrary. 

        """ 

        # make sure we have k in 1/Mpch
        #----------------------------------------------
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        T_cdm = self.T_cdm( k )

        Rjr = 1.6 / np.sqrt( self.OmegaC ) * self.cu.kpc / self.cu.h
        kRjr2 = (k * Rjr)**2

        kRjr2 = kRjr2.simplified 

        T = T_cdm / (1.0 + kRjr2 / 2)
        return T


    def Pk_cdm( self, k ):
        """ CDM Power spectrum 

        Args: 

          `k` (real or array): wavenumber with units of inverse length

        Returns:
        
          `Pk_cdm` (real or array): CDM power spectrum at the scales `k`.  
          Normalization is arbitrary. 
        """ 

        # make sure we have k in 1/Mpch
        #----------------------------------------------
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        kk = np.array(k)
        T = self.T_cdm( k )
        Pk = kk**self.ns * T * T
        return Pk


    def Pk_b( self, k ):
        """ Baryon power spectrum 

        Args: 

          `k` (real or array): wavenumber with units of inverse length

        Returns:
        
          `Pk_b` (real or array): Baryon power spectrum at the scales `k`.  
          Normalization is arbitrary. 
        """ 

        # make sure we have k in 1/Mpch
        #----------------------------------------------
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        kk = np.array(k)
        T = self.T_b( k )
        Pk = kk**self.ns * T * T
        return Pk






#if __name__ == "__main__":
#
#    print 'Testing BBKS Transfer Function' 
#    print 'Note the units of k are h Mpc^-1'

#    cpdict = dict( omegam=0.2865, 
#                   omegab=0.04628,
#                   omegac=0.2865 - 0.04628,
#                   h=0.6932,
#                   ns=1.0 )


#    bbks = TransferBBKS( cpdict )

#    log_kmin = -4.0
#    log_kmax = 4.0
#    Nk = 1000

#    log_k = np.linspace( log_kmin, log_kmax, Nk )
#    k = 10**log_k / bbks.cu.Mpc / bbks.cu.h

#    T_cdm = bbks.T_cdm( k )
#    T_b   = bbks.T_b( k )

#    fig = plt.figure( figsize=(10,10) )
#    ax = fig.add_subplot( 111 )

#    ax.loglog( k, T_cdm  ) 
#    ax.loglog( k, T_b  ) 

#    ax.set_xlabel( 'k [h/Mpc]' )
#    ax.set_ylabel( 'T' )
