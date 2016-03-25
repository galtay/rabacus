r""" A general FLRW cosmology module.  Assumes a spatially flat universe. """ 

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d

from rabacus.constants import physical
from rabacus.constants import units
from rabacus.utils import utils



__all__ = ['Cosmology']



class Error(Exception):
    r"""Base class for exceptions in this module."""
    pass


class InitError(Error):
    r"""Exception raised for not having required keys in input. """
    def __init__(self, message ):
        message = message + '\n' + \
            'cpdict must contain the following keys ... \n' + \
            '"omegam" -> current matter density in units of critical \n' + \
            '"omegal" -> current lambda density in units of critical \n' + \
            '"omegab" -> current baryon density in units of critical \n' + \
            '"h"      -> Hubble parameter H0 = 100 h km/s/Mpc \n' + \
            '"sigma8  -> density variance in spheres w/ R = 8 Mpc/h \n' + \
            '"ns"     -> slope of primordial power spectrum \n' + \
            '"Yp"     -> primordial mass fraction of helium \n'
        Exception.__init__( self, message )



class Cosmology(object):
    r""" General Cosmology Class.  Assumes a spatialy flat universe.  

    Args: 

      `cpdict` (dict) A dictionary of cosmological parameters.  
      For example, see
      :class:`~rabacus.cosmology.parameters.planck.load.PlanckParameters`.

    The dictionray `cpdict` must include the following keys, 
       - ``omegam`` -> current matter density in units of critical today
       - ``omegal`` -> current lambda density in units of critical today
       - ``omegab`` -> current baryon density in units of critical today
       - ``h``      -> Hubble parameter H0 = 100 h km/s/Mpc 
       - ``sigma8`` -> amplitude of fluctuations in spheres w/ R = 8 Mpc/h
       - ``ns``     -> slope of primordial power spectrum 
       - ``Yp``     -> primordial mass fraction of helium


    Attributes:

      H0 (real): hubble parameter now, :math:`H_0`.

      OmegaB (real): baryon density / critical density now,
      :math:`\Omega_b` 
    
      OmegaC (real): cold dark matter density / critical density now,
      :math:`\Omega_c`

      OmegaL (real): dark energy density / critical density now,
      :math:`\Omega_{\Lambda}`

      OmegaM (real): matter density / critical density now,
      :math:`\Omega_m`

      Yp (real): primoridial helium mass fraction, 
      :math:`Y_p`

      cu (:class:`~rabacus.constants.units.CosmoUnits`): cosmological 
      units which are aware of the hubble parameter. 

      dH0 (real): hubble distance now, 
      :math:`d_{\rm H_0} = c / H_0`

      tH0 (real): hubble time now, 
      :math:`t_{\rm H_0} = 1 / H_0`

      rho_crit0 (real): critical density now, 
      :math:`\rho_{c,0} = (3 H_0^2) / (8 \pi G)`

      eps_crit0 (real): critical energy density now,
      :math:`\epsilon_{c,0} = \rho_{c,0} \, c^2`

      nH_crit0 (real): critical hydrogen number density now,
      :math:`n_{\rm H,c,0} = \rho_{c,0} \, \Omega_b \, (1-Y_p) \, /  \, 
      m_{\rm H}`

      

    """ 

    def __init__( self, cpdict, verbose=False, zlo=0.0, zhi=200.0, Nz=500 ):


        # check input dictionary
        #--------------------------------------------------------
        need_keys = ['omegam', 'omegal', 'omegab', 'h', 'sigma8', 'ns', 'Yp']
        for key in need_keys:
            if not cpdict.has_key( key ):
                raise InitError, ' key missing from cpdict \n' 
            
        # transfer input parametrs
        #--------------------------------------------------------
        self.OmegaM = cpdict['omegam']
        self.OmegaL = cpdict['omegal']
        self.OmegaB = cpdict['omegab']
        self.h = cpdict['h']
        self.sigma8 = cpdict['sigma8']
        self.ns = cpdict['ns']
        self.Yp = cpdict['Yp']

        self.OmegaC = self.OmegaM - self.OmegaB

        if cpdict.has_key('omegar'):
            self.OmegaR = cpdict['omegar']
        else:
            # if not provided we use OmegaM and zeq from Planck13 
            # column one of table 2 in http://arxiv.org/abs/1303.5076
            # OmegaR = OmegaM / (1+zeq), OmegaM=0.3175, zeq = 3402
            self.OmegaR = 0.3175 / (1.0+3402.)  

        self.cpdict = cpdict

        # make sure flat universe is specified
        #--------------------------------------------------------
        err = ( self.OmegaR + self.OmegaM + self.OmegaL ) - 1.0
        if err > 1.0e-3:
            raise InitError, ' OmegaR + OmegaM + OmegaL must be 1.0 \n' + \
                ' OmegaR: ' + str(self.OmegaR) + '\n' + \
                ' OmegaM: ' + str(self.OmegaM) + '\n' + \
                ' OmegaL: ' + str(self.OmegaL) + '\n' + \
                ' sum: ' + str(self.OmegaR + self.OmegaM + self.OmegaL) + '\n' 


        # create physical constants and cosmological units 
        #--------------------------------------------------------
        self.pc = physical.PhysicalConstants()
        self.cu = units.CosmoUnits( self.h, a=1.0 )

        # derive constants
        #--------------------------------------------------------

        # Hubble parameter now
        self.H0 = 1.0e2 * self.cu.h * self.cu.km / self.cu.s / self.cu.Mpc
        self.H0.units = self.cu.km / self.cu.s / (self.cu.Mpc / self.cu.h)

        # Hubble time now
        self.tH0 = 1.0 / self.H0
        self.tH0.units = self.cu.Myr/self.cu.h

        # Hubble distance now
        self.dH0 = self.pc.c / self.H0
        self.dH0.units = self.cu.Mpc/self.cu.h

        # Critical mass density now
        self.rho_crit0 = ( 3.0 * self.H0**2 ) / (8.0 * np.pi * self.pc.G)
        self.rho_crit0.units = (self.cu.Msun/self.cu.h) / \
            (self.cu.Mpc/self.cu.h)**3

        # Critical energy density now
        self.eps_crit0 = self.rho_crit0 * self.pc.c**2
        self.eps_crit0.units = self.cu.erg / self.cu.cm**3

        # Critical hydrogen number density now
        H_rho = self.rho_crit0 * self.OmegaB * (1.0-self.Yp)
        self.nH_crit0 = H_rho / ( self.pc.m_p + self.pc.m_e )
        self.nH_crit0.units = self.cu.cm**(-3)

        # normalize the growth functions such that ... 
        #    D1a(a=1) = 1.0
        #    D1z(z=0) = 1.0
        #---------------------------------------------------------
        self._D1a_Norm = 1.0 / quad( self._D1a_integrand, 0.0, 1.0 )[0]

        self._D1z_Norm = 1.0 / quad( self._D1z_integrand, 0.0, np.inf )[0]

        # make interpolating functions for ....
        # redshift and comoving distance,
        # redshift and lookback time
        # for redshift we use log(1+z)
        #--------------------------------------------------------
 
        # first tabulate values with respect to log(1+z)
        # then create interpolating functions 

        self._lopz = np.linspace( np.log10(1+zlo), 
                                  np.log10(1+zhi), 
                                  Nz )

        self._z = 10**self._lopz - 1
        self._a = 1.0 / (1.0 + self._z)

        self._Dc = self.Dcz( self._z ) 
        self._Dc.units = 'Mpc/hh'
        self._Dc2lopz = interp1d( self._Dc.magnitude, self._lopz )

        self._tL = self.tLz( self._z ) 
        self._tL.units = 'Myr/hh'
        self._tL2lopz = interp1d( self._tL.magnitude, self._lopz )


    def z_at_Dc(self, Dc):
        r""" Redshift at a given comoving distance, z(Dc).  Inverts the 
        function :func:`Dcz` by interpolating between tabulated values. 
        
        .. math:: 
          z(D_C) = {\rm Inverse}[ D_C(z) ]
        """ 
        if hasattr(Dc,'units'): 
            Dc.units = 'Mpc/hh'
        else:
            msg = '\n Input variable Dc must have units of length. \n'
            raise ValueError(msg)
        lopz = self._Dc2lopz( Dc )
        z = (10**lopz-1) 
        return z

    def z_at_tL(self, tL):
        r""" Redshift at a given lookback time, z(tL).  Inverts the 
        function :func:`tLz` by interpolating between tabulated values. 

        .. math:: 
          z(t_L) = {\rm Inverse}[ t_L(z) ]
        """ 
        if hasattr(tL,'units'): 
            tL.units = 'Myr/hh'
        else:
            msg = '\n Input variable tL must have units of time. \n'
            raise ValueError(msg)
        lopz = self._tL2lopz( tL )
        z = (10**lopz-1) 
        return z

    def tL_at_Dc(self, Dc):
        r""" Lookback time at a given comoving distance, tL(Dc).  First 
        finds a redshift `z` by inverting the function :func:`Dcz` using 
        tabulated values.  Second, calls :func:`tLz` to get a lookback 
        time from `z`.  
        
        .. math:: 
          z = {\rm Inverse}[ D_C(z) ] \\
          t_L(z) = \int_0^z \frac{dz'}{(1+z') H(z')} 
        """ 
        if hasattr(Dc,'units'): 
            Dc.units = 'Mpc/hh'
        else:
            msg = '\n Input variable Dc must have units of length. \n'
            raise ValueError(msg)
        lopz = self._Dc2lopz( Dc )
        z = (10**lopz-1) 
        tL = self.tLz(z)
        return tL

    def Dc_at_tL(self, tL):
        r""" Comoving distance at a given lookback time, Dc(tL). First 
        finds a redshift `z` by inverting the function :func:`tLz` using
        tabulated values.  Second, calls :func:`Dcz` to get a comoving
        distance from `z`. 

        .. math:: 
          z = {\rm Inverse}[ t_L(z) ] \\
          D_C(z) = d_{\rm H_0} \int_{0}^{z} \frac{dz'}{E(z')}
        """ 
        if hasattr(tL,'units'): 
            tL.units = 'Myr/hh'
        else:
            raise ValueError('\n Input variable tL must have units \n')
        lopz = self._tL2lopz( tL )
        z = (10**lopz-1) 
        Dc = self.Dcz(z)
        return Dc

    def Ea(self, a):
        r""" Helper function H(a) = H0 * E(a) 

        .. math::
          E(a) = ( \Omega_R a^{-4} + \Omega_M a^{-3} + \Omega_{\Lambda} )^{1/2} 
          
        """
        E = np.sqrt( self.OmegaR*a**(-4) + self.OmegaM*a**(-3) + self.OmegaL )
        return E

    def Ez(self, z):
        r""" Helper function H(z) = H0 * E(z) 

        .. math::
          E(z) = [ \Omega_R (1+z)^{4} + \Omega_M (1+z)^{3} + 
          \Omega_{\Lambda} ]^{1/2} 

        """
        E = np.sqrt( self.OmegaR*(1+z)**4 + self.OmegaM*(1+z)**3 + self.OmegaL )
        return E

    def iEa(self, a):
        r""" Reciprocal of :func:`Ea`.  

        .. math::
          E(a) = ( \Omega_R a^{-4} + \Omega_M a^{-3} + 
          \Omega_{\Lambda} )^{-1/2} 

        """
        iE = 1.0 / self.Ea(a)
        return iE

    def iEz(self, z):
        r""" Reciprocal of :func:`Ez`.  

        .. math::
          E(z) = [ \Omega_R (1+z)^{4} + \Omega_M (1+z)^{3} + 
          \Omega_{\Lambda} ]^{-1/2} 
        """
        iE = 1.0 / self.Ez(z)
        return iE

    def Ha(self, a):
        r""" Hubble parameter at scale factor `a`, H(a) = H0 * E(a) 
        
        .. math::
          H(a) = H_0 E(a) 
        """
        H = self.H0 * self.Ea(a)
        return H

    def Hz(self, z):
        r""" Hubble parameter at redshift `z`, H(z) = H0 * E(z) 

        .. math::
          H(z) = H_0 E(z) 
        """
        H = self.H0 * self.Ez(z)
        return H

    def tHa(self, a):
        r""" Hubble time at scale factor `a`, tH(a) 
        
        .. math:: 
          t_{\rm H(a)} = 1 / H(a)
        """
        t = 1.0 / self.Ha(a)
        t.units = 'Myr/hh'
        return t

    def tHz(self, z):
        r""" Hubble time at redshift `z`, tH(z) 

        .. math:: 
          t_{\rm H(z)} = 1 / H(z)
        """
        t = 1.0 / self.Hz(z)
        t.units = 'Myr/hh'
        return t

    def dHa(self, a):
        r""" Hubble distance at scale factor `a`, dH(a) 

        .. math::
          d_{\rm H(a)} = c / H(a)
        """
        d = self.pc.c / self.Ha(a)
        d.units = 'Mpc/hh'
        return d

    def dHz(self, z):
        r""" Hubble distance at redshift `z`, dH(z) 

        .. math::
          d_{\rm H(z)} = c / H(z)
        """
        d = self.pc.c / self.Hz(z)
        d.units = 'Mpc/hh'
        return d

    def rho_crita(self, a):
        r""" Critical mass density at scale factor `a`, rho_crit(a) 
        
        .. math:: 
          \rho_{c,a} = \rho_{c,0} \,  E(a)^2
        """ 
        rho = self.rho_crit0 * self.Ea(a)**2
        return rho
        
    def rho_critz(self, z):
        r""" Critical mass density at redshift `z`, rho_crit(z) 

        .. math:: 
          \rho_{c,z} = \rho_{c,0} \,  E(z)^2
        """ 
        rho = self.rho_crit0 * self.Ez(z)**2
        return rho

    def eps_crita(self, a):
        r""" Critical energy density at scale factor `a`, eps_crit(a) 

        .. math:: 
          \epsilon_{c,a} = \epsilon_{c,0} \,  E(a)^2
        """ 
        eps = self.eps_crit0 * self.Ea(a)**2
        return eps
        
    def eps_critz(self, z):
        r""" Critical energy density at redshift `z`, eps_crit(z) 

        .. math:: 
          \epsilon_{c,z} = \epsilon_{c,0} \,  E(z)^2
        """ 
        eps = self.eps_crit0 * self.Ez(z)**2
        return eps

    def nH_crita(self, a):
        r""" Critical hydrogen number density at scale factor `a`, nH_crit(a) 

        .. math:: 
          n_{\rm H,c,a} =  n_{\rm H,c,0}  \,  E(a)^2
        """ 
        nH = self.nH_crit0 * self.Ea(a)**2
        return nH
        
    def nH_critz(self, z):
        r""" Critical hydrogen number density at redshift `z`, nH_crit(z) 

        .. math:: 
          n_{\rm H,c,z} =  n_{\rm H,c,0}  \,  E(z)^2
        """ 
        nH = self.nH_crit0 * self.Ez(z)**2
        return nH

    def ta( self, a):
        r""" Time since a=0 or age of the Universe, t(a). 

        .. math::
          t(a) = \int_0^a \frac{da'}{a' H(a')} 

        """ 

        # quad works but takes a long time. 
        #-------------------------------------------------------------------
#        if utils.isiterable( a ):
#            t = np.array( [quad( self._dt_da, 0.0, aa )[0] for aa in a] )
#        else:
#            t = quad( self._dt_da, 0.0, a )[0]
#        t = t * self.cu.s
#        print 't = ', t.rescale("Myr/hh")

        # trap rule includes OmegaR and is fast and accurate to 6 decimals
        #-------------------------------------------------------------------
        N = 5000
        a_min = 1.0e-10
        if utils.isiterable( a ):
            t = np.zeros( a.size ) * self.cu.s
            for ii,aa in enumerate(a):
                xx = np.linspace( a_min, aa, N )
                yy = 1.0 / ( xx * self.Ha(xx) )
                t[ii] = utils.trap( xx, yy )
        else:
            xx = np.linspace( a_min, a, N )
            yy = 1.0 / ( xx * self.Ha(xx) )
            t = utils.trap( xx, yy )
        t.units = 'Myr/hh'

        # this is analytic but does not include OmegaR        
        #-------------------------------------------------------------------
#        pre = 2.0 / ( 3.0 * np.sqrt( self.OmegaL ) )
#        aeq = (self.OmegaM/self.OmegaL)**(1./3.)  
#        arg = (a/aeq)**(3./2.) + np.sqrt(1 + (a/aeq)**3)
#        t = pre * np.log(arg) * self.tH0

        return t 

    def _dt_da( self, a ):
        r""" Differential dt/da """ 
        Ha = self.Ha(a).rescale('1/s')
        Ha = Ha.magnitude
        dtda = 1.0 / (a * Ha)
        return dtda



    def tz( self, z):
        r""" Time since z=inf or age of the Universe t(z). 

        .. math::
          t(z) = \int_z^\infty \frac{dz'}{(1+z') H(z')} 

        """
        a = 1.0/(1.0+z)
        t = self.ta(a)
        return t 

    def tLa( self, a):
        r""" Lookback from a = 1 to a, tL(a). 

        .. math::
          t_{\rm L}(a) = \int_a^1 \frac{da'}{a' H(a')} 
        
        """ 
        t1 = self.ta(1.0)
        tL = t1 - self.ta(a)
        return tL 
    
    def tLz( self, z):
        r""" Lookback from z = 0 to z, tL(z). 

        .. math::
          t_{\rm L}(z) = \int_0^z \frac{dz'}{(1+z') H(z')} 

        """ 
        t1 = self.tz(0.0)
        tL = t1 - self.tz(z)
        return tL 

    def Dcz(self,z):
        r""" Comoving distance between z=0 and z, Dc(z) 
        
        .. math::
          D_C(z) = d_{\rm H_0} \int_{0}^{z} \frac{dz'}{E(z')} 
          
        """
        if utils.isiterable( z ):
            Dc = np.array( [quad( self.iEz, 0.0, zz )[0] for zz in z] )
            Dc = Dc * self.dH0
        else:
            Dc = quad( self.iEz, 0.0, z )[0] * self.dH0
        Dc.units = self.cu.Mpc/self.cu.h
        return Dc 

    def Dca(self,a):
        r""" Comoving distance between a=1 and a, Dc(a).  The scale factor
        `a` is converted to redshift `z` and then :func:`Dcz` is called.

        .. math::
          D_C(a) = d_{\rm H_0} \int_{a}^{1} \frac{da'}{a'^2 E(a')}         

        """
        z = 1.0/a - 1.0
        Dc = self.Dcz(z)
        return Dc 

    def dz2dDc(self,z1,z2):
        r""" Comoving distance between `z1` and `z2`.  Two calls to :func:`Dcz`
        are made and the results subtracted.  `z1` < `z2` produces positive
        comoving distance. 

        .. math::
          D_C(z1,z2) = d_{\rm H_0} \int_{z1}^{z2} \frac{dz'}{E(z')} 
        """
        if utils.isiterable(z1) and utils.isiterable(z2):
            if not len(z1) == len(z2):
                raise utils.InputError, "len(z1) /= len(z2)"
        Dc1 = self.Dcz(z1)
        Dc2 = self.Dcz(z2)
        Dc = Dc2 - Dc1
        return Dc

    def da2dDc(self,a1,a2):
        r""" Comoving distance between `a1` and `a2`. The scale factors
        `a1` and `a2` are converted to redshifts `z1` and `z2` and then 
        :func:`dz2dDc` is called.  `a1` < `a2` produces positive comoving
        distance. 

        .. math::
          D_C(a1,a2) = d_{\rm H_0} \int_{a1}^{a2} \frac{da'}{a'^2 E(a')} 
        """
        z2 = 1.0/a1 - 1.0
        z1 = 1.0/a2 - 1.0
        Dc = self.dz2dDc(z1,z2)
        return Dc 

    def da2dtL(self,a1,a2):
        r""" Lookback time between `a1` and `a2`.  Two calls to :func:`tLa`
        are made and the results subtracted.  `a1` < `a2` produces positive
        lookback time. 

        .. math::
          t_{\rm L}(a1,a2) = \int_{a1}^{a2} \frac{da'}{a' H(a')} 
        """  
        t1 = self.ta( a1 )
        t2 = self.ta( a2 )
        dt = t2 - t1
        return dt

    def dz2dtL(self,z1,z2):
        r""" Lookback time between `z1` and `z2`.  The redshifts `z1` and 
        `z2` are converted to scale factors and then :func:`da2dtL` is called.
        `z1` < `z2` produces positive lookback time. 

        .. math::
          t_{\rm L}(z1,z2) = \int_{z1}^{z2} \frac{dz'}{(1+z') H(z')} 
        """ 
        a2 = 1.0 / (1.0 + z1)
        a1 = 1.0 / (1.0 + z2)
        dt = self.da2dtL( a1, a2 )
        return dt

    def X(self,z1,z2):
        r""" Absorption distance between `z1` and `z2`,
        
        .. math::
          X(z1,z2) = \int_{z1}^{z2} dz' (1+z')^2 / E(z') 
        """ 
        X = quad( self._X_integrand, z1, z2 )
        return X[0]

    def _X_integrand( self, z ):
        r""" function to interate in order to find absorption distance. """ 
        dXdz = (1.0+z)*(1.0+z) * self.iEz(z)
        return dXdz

    # Cosmological growth functions
    #----------------------------------------------------------------
    def D1z( self, z ):
        r""" Linear growth function D1(z) 

        .. math::
          D_1(z) \propto H(z) \int_z^\infty \frac{1+z'}{H(z')^3} dz'

        The normalization is such that D1z(0) = 1

        """ 

        if utils.isiterable(z):
            int_u = [quad( self._D1z_integrand, zz, np.inf )[0] for zz in z]
        else:
            int_u = quad( self._D1z_integrand, z, np.inf )[0]

        D1 = self.Hz(z) / self.H0 * int_u * self._D1z_Norm
        return D1

    def _D1z_integrand( self, z ):
        r""" Integrand for D1z """ 
        return (1.0 + z) / self.Hz(z)**3
    


    def D1a(self, a):
        r""" Linear growth function D1(a).  

        .. math::
          D_1(a) \propto H(a) \int_0^a \frac{da'}{a'^3 H(a')^3}

        The normalization is such that D1a(1) = 1

          """ 

        if utils.isiterable(a):
            integral = [quad( self._D1a_integrand, 0.0, aa )[0] for aa in a]
        else:
            integral = quad( self._D1a_integrand, 0.0, a )[0]

        D1 = self.Ha(a) / self.H0 * integral * self._D1a_Norm
        return D1 

    def _D1a_integrand( self, a ):
        r""" Integrand for D1a """ 
        return 1.0 / ( a**3 * self.Ha(a)**3 )



    # Integrands
    #----------------------------------------------------------------





    # plotting
    #----------------------------------------------------------------
#    def _plot_D1a(self):
#        N = 100
#        aa = np.linspace( 1.0e-4, 1.0, N )
#        zz = 1.0/aa - 1.0
#        D1z = self.D1z(zz)
#        D1a = self.D1a(aa)
#
#        fig = plt.figure()
#        plt.plot( aa, D1a, color='green' )
#        plt.plot( aa, D1z, ls='--', color='red' )
#
#        plt.xlabel( 'a', fontsize=20 )
#        plt.ylabel( 'D1(a)', fontsize=20 )


#    def _plot_Dc_vs_z(self):
#        fig = plt.figure()
#        Dc = self._Dc[::8].copy()    
#        lopz = self._lopz[::8].copy()
#        plt.plot( lopz, Dc, 'ro' )    
#        z = self.zDc( Dc )
#        lopz = np.log10( 1 + z )
#        plt.plot( lopz, Dc, 'g-' )
#        plt.xlabel( r'$\log(1+z)$', fontsize=20 )
#        plt.ylabel( r'$D_C \, [\rm Mpc/h]$', fontsize=20 )


    def __str__(self):
        print 'Cosmology instance'
        print '  OmegaR:  %5.3f' % self.OmegaR
        print '  OmegaM:  %5.3f' % self.OmegaM
        print '  OmegaB:  %5.3f' % self.OmegaB
        print '  OmegaC:  %5.3f' % self.OmegaC
        print '  OmegaL:  %5.3f' % self.OmegaL
        print '  h:       %5.3f' % self.h
        print '  sigma_8: %5.3f' % self.sigma8
        print '  n_s:     %5.3f' % self.ns
        print '  Y_P:     %5.3f' % self.Yp
        print 
        print '  H0: ', self.H0
        print '  dH: ', self.dH0
        print '  tH: ', self.tH0
        print '  rho_crit_0: ', self.rho_crit0
        print '  eps_crit_0: ', self.eps_crit0
        return ''
