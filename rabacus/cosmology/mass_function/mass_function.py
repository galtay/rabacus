""" 
A halo mass function module.  See the following reference for a discussion, 
http://adsabs.harvard.edu/abs/2007ApJ...671.1160L
"""

import sys

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl



#import ..cosmology
#import bbks86
#from ....constants import physical
#import idle_hands.constants.physical_constants_units as pcu

import scipy.integrate
import scipy.interpolate

__all__ = ['MassFunction']




class MassFunction:

    """ A mass function class.  During initialization the normalization of 
    the power spectrum is set to match the sigma8 from cosmo. 

    Args: 

      `cosmo` (:class:`~rabacus.cosmology.general.Cosmology`): an instance 
      of the cosmology class. 

      `tf` (class): an instance of a transfer function class.  For example,
      :class:`~rabacus.cosmology.transfer_functions.bbks86.TransferBBKS`.

    """ 

    def __init__(self, cosmo, tf):

        """ cosmo is an instance of Cosmology and 
            tf is an instance of TransferFunction """

        self.cosmo = cosmo
        self.tf = tf
        self.cu = cosmo.cu

        self.power_norm = 1.0
        sig2 = self.sigma2_R( 8.0 * self.cu.Mpc / self.cu.h, z=0.0 )
        self.power_norm = self.cosmo.sigma8**2 / sig2         
        self.delta_c0 = 1.68647

        self.map_sig2_R()


    def calc_mf( self, z, fit='Warren06' ):
        """ Calculate mass function from high to low mass.  The only
        redshift dependence is in f(sigma) via a rescaling of sigma """ 

        N = self.Mw_arr.size
        self.dndM = np.zeros( N ) * (self.cu.h / self.cu.Mpc)**3
        self.dndlogM = np.zeros( N ) * (self.cu.h / self.cu.Mpc)**3

        rho_mean = self.cosmo.rho_crit0 * self.cosmo.OmegaM
        delta_c = self.delta_c0 / self.D1( z )

        for i in reversed(np.arange(N)):
            
            # go from high mass to low mass

            sigma = self.sig_arr[i]
            mass = self.Mw_arr[i]
            logM = np.log10(  np.array( mass ) )

            f = self.mult_func( sigma, z, fit=fit )

            dlnisigdM    = self.lnisig_of_Mw.derivatives(mass)[1]
            dlnisigdlogM = self.lnisig_of_logMw.derivatives(logM)[1]

            # http://adsabs.harvard.edu/abs/2008ApJ...688..709T Eq. 2
            dndM    = (rho_mean / mass) * f * dlnisigdM
            dndlogM = (rho_mean / mass) * f * dlnisigdlogM

            #print 'M,dndM,dndlogM: ', mass, dndM, dndlogM 
            
            self.dndM[i] = dndM
            self.dndlogM[i] = dndlogM

        print 'attribute dndM now contains dn/dM at z = ', str(z)
        print 'attribute dndlogM now contains dn/dlogM at z = ', str(z)



    def mult_func( self, sigma_in, z, fit='Warren06' ):

        r""" The multiplicity function :math:`f(\sigma)`.  This function 
        determines the shape of the mass function given the variation of 
        :math:`\sigma` with scale.  The variable `fit` determines the form
        of the multiplicity function.  A convenient variable is 
        :math:`\nu = \delta_{c,0} / \sigma`.  The following values for `fit`
        lead to the following multiplicity funcitons, 

          - ``PS``: `Press & Schechter 74 
            <http://adsabs.harvard.edu/abs/1974ApJ...187..425P>`_
        
            .. math::
              f = \sqrt{ \frac{2}{\pi} } \, \nu \exp( -\nu^2/2 )


          - ``ST``: `Sheth & Tormen 99 
            <http://adsabs.harvard.edu/abs/1999MNRAS.308..119S>`_

            .. math::
              f = A \sqrt{ \frac{2a}{\pi} } [1 + (a \nu^2)^{-p}]
              \nu \exp(-a \nu^2/2) 

              A = 0.3222, \, a = 0.75, \, p = 0.3


          - ``Jenkins01``: `Jenkins 01 
            <http://adsabs.harvard.edu/abs/2001MNRAS.321..372J>`_

            .. math::
              f = 0.315 \exp( -| \ln \sigma^{-1} + 0.61 |^{3.8} )


          - ``Warren06``: `Warren 06 
            <http://adsabs.harvard.edu/abs/2006ApJ...646..881W>`_

            .. math::
              f = A ( \sigma^{-a} + b ) \exp(-c/\sigma^2)

              A=0.7234, \, a=1.625, \, b=0.2538, \, c=1.1982


          - ``Tinker08``: `Tinker 08 
            <http://adsabs.harvard.edu/abs/2008ApJ...688..709T>`_

            .. math::
              f = A \left[ \left( \frac{\sigma}{b} \right)^{-a} 
              + 1 \right] \exp(-c/\sigma^2)

              \Delta=300

              A = 0.1 \, \log \Delta - 0.05
              
              a = 1.43 + ( \log \Delta - 2.3  )^{1.5}

              b = 1.00 + ( \log \Delta - 1.6  )^{-1.5}
        
              c = 1.20 + ( \log \Delta - 2.35 )^{1.6}


        """ 

        sigma = sigma_in * self.D1(z)
        nu = self.delta_c0 / sigma
        
        if fit == 'PS':
            # http://adsabs.harvard.edu/abs/2007ApJ...671.1160L Eq. 9
            f = np.sqrt(2/np.pi) * nu * np.exp(-0.5*nu*nu)

        elif fit == 'ST':
            # http://adsabs.harvard.edu/abs/2007ApJ...671.1160L Eq. 11
            A = 0.3222
            a = 0.75 # a = 0.707
            p = 0.3
            t1 = A * np.sqrt( 2*a/np.pi )
            t2 = 1 + (a*nu*nu)**(-p)
            t3 = nu * np.exp( -0.5*a*nu*nu )
            f = t1 * t2 * t3

        elif fit == "Jenkins01":
            # http://adsabs.harvard.edu/abs/2007ApJ...671.1160L Eq. 12
            t1 = np.abs( np.log(1.0/sigma) + 0.61 )
            f = 0.315 * np.exp( -t1**(3.8) )
            
        elif fit == "Warren06":
            # http://adsabs.harvard.edu/abs/2007ApJ...671.1160L Eq. 13
            A=0.7234
            a=1.625 
            b=0.2538 
            c=1.1982
            t1 = sigma**(-a) + b
            f = A * t1 * np.exp(1.0)**( -c / (sigma*sigma) )

        elif fit == "Tinker08":
            # http://adsabs.harvard.edu/abs/2008ApJ...688..709T Eq. 3
            Delta = 300.0
            logD = np.log10( Delta )

            if Delta < 1600.0:
                A = 0.1 * logD - 0.05
            else:
                A = 0.26

            a = 1.43 + ( logD - 2.3  )**(1.5)
            b = 1.00 + ( logD - 1.6  )**(-1.5)
            c = 1.20 + ( logD - 2.35 )**(1.6)
            t1 = (sigma/b)**(-a) + 1
            f = A * t1 * np.exp(1.0)**( -c / (sigma*sigma) )

        else:
            print 'fit not recognized'
            sys.exit(1)

        return f


    def map_sig2_R( self, nbins=200 ):
        """ Map out the relationship between sigma^2 and R.  In this 
        routine we map out the relationship at z=0 and assume that 
        different redshifts can be accomodated through a simple scaling
        with the growth function D1(z). """ 

        z = 0.0
        Rw_min = 8.0e-4 * self.cu.Mpc / self.cu.h
        Rw_max = 8.0e1  * self.cu.Mpc / self.cu.h
        logRw_min = np.log10( np.array(Rw_min) )
        logRw_max = np.log10( np.array(Rw_max) )

        rho_mean = self.cosmo.rho_crit0 * self.cosmo.OmegaM
        V_min = (4.0/3.0) * np.pi * Rw_min**3
        V_max = (4.0/3.0) * np.pi * Rw_max**3

        M_min = rho_mean * V_min
        M_max = rho_mean * V_max

        print 'min/max Rw correspond to ... '
        print 'Rw min/max: ', Rw_min, Rw_max
        print 'M  min/max: ', M_min, M_max
        print 'calculating relationship between R and sigma^2'
        print 'this may take a minute or so ... '

        N = nbins
        dlogRw = ( logRw_max - logRw_min ) / N

        self.logRw_arr = np.linspace( logRw_min, logRw_max, N )
        self.Rw_arr = 10**self.logRw_arr * self.cu.Mpc / self.cu.h
        self.Vw_arr = (4.0/3.0) * np.pi * self.Rw_arr**3
        self.Mw_arr = rho_mean * self.Vw_arr

        self.logVw_arr = np.log10( np.array( self.Vw_arr ) )
        self.logMw_arr = np.log10( np.array( self.Mw_arr ) )

        self.logsig_arr = np.zeros( N )
        self.sig_arr = np.zeros( N )
        self.isig_arr = np.zeros( N )
        self.lnisig_arr = np.zeros( N )

        for i in range(N):
            Rw = self.Rw_arr[i] 
            sig2 = self.sigma2_R( Rw, z )
            sig = np.sqrt(sig2)
            self.logsig_arr[i] = np.log10( sig )
            self.sig_arr[i] = sig

            #print 'i,Rw,sig2: ', i,Rw,sig
        
        self.isig_arr = 1.0 / self.sig_arr
        self.lnisig_arr = np.log( self.isig_arr )
        

        # create spline representations
        #----------------------------------------------------------
        fit = scipy.interpolate.UnivariateSpline( self.logsig_arr,
                                                  self.logRw_arr,
                                                  s=0)
        self.logRw_of_logsig = fit


        fit = scipy.interpolate.UnivariateSpline( self.logRw_arr,
                                                  self.logsig_arr,
                                                  s=0)
        self.logsig_of_logRw = fit


        fit = scipy.interpolate.UnivariateSpline( self.Mw_arr,
                                                  self.lnisig_arr,
                                                  s=0)
        self.lnisig_of_Mw = fit

        fit = scipy.interpolate.UnivariateSpline( self.logMw_arr,
                                                  self.lnisig_arr,
                                                  s=0)
        self.lnisig_of_logMw = fit

        fit = scipy.interpolate.UnivariateSpline( self.Mw_arr,
                                                  self.sig_arr,
                                                  s=0)
        self.sig_of_Mw = fit



    def delta_cz( self, z ):
        """ Redshift dependent critical collapse value 

        Args:
          `z` (real): redshift

        .. math:: 
          \delta_c(z) = \delta_{c,0} / D_1(z)

        """ 
        dz = self.delta_c0 / self.D1( z )
        return dz

    def Pk(self, k, z):
        """ Power spectrum from `tf` normalized to match sigma8 from 
        `cosmo`. """ 
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'
        Pk = self.tf.Pk_cdm(k)
        return Pk * self.power_norm * self.D1(z)**2


    def Delta2k(self, k, z):
        r""" Dimensionless power spectrum,

        Args: 

          `k` (real or array): wavenumber.

          `z` (real): redshift

        .. math::
          \Delta^2(k,z) = \frac{k^3}{2 \pi^2} P(k,z)

        """ 
        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'
        kk = np.array(k)
        Dk = kk * kk * kk / (2.0 * np.pi * np.pi) * self.Pk(k,z)
        return Dk 


    def D1(self,z):
        """ Linear growth function from `cosmo` """ 
        D = self.cosmo.D1z( z )
        return D

    def WkR(self, k, R):  
        r""" The fourier transform of a real space spherical top hat filter. 

        Args: 

          `k` (real or array): wavenumber.

          `R` (real): filter scale. 

        .. math::
          W =  \frac{ 3 j_1(x) }{ x }, \, x = k R \\
          j_1 = [\sin(x) - x \cos(x)] \, x^{-2}

        """

        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'

        y = k * R
        yi = 1.0 / y
        yi2 = yi * yi
        yi3 = yi2 * yi
        t1 = np.sin(y) * yi3
        t2 = np.cos(y) * yi2
        W = 3.0 * (t1-t2)
        return W

    def W2kR(self, k,R):  
        r""" The square of the fourier transform of a real space spherical top 
        hat filter.

        Args: 

          `k` (real or array): wavenumber.

          `R` (real): filter scale. 

        .. math::
          W^2 =  \left[ \frac{ 3 j_1(x) }{ x } \right]^2, \, x = k R \\
          j_1 = [\sin(x) - x \cos(x)] \, x^{-2}

        """

        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'

        W = self.WkR(k,R)
        W2 = W * W
        return W2


    def W2lnkR(self, lnk, R):
        r""" The square of the fourier transform of a real space spherical 
        top hat filter as a  function of the natural log of `k`. 

        Args: 

          `lnk` (real or array): natural log of wavenumber, ln( k [h/Mpc] ).

          `R` (real): filter scale. 
        
        """

        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'

        kk = np.exp(1.0)**lnk / self.cu.Mpc / self.cu.h
        W2 = self.W2kR( kk, R )
        return W2


    def dsig2_dk( self, k, R, z ):
        r""" Integrand for calculation of sigma^2.  Inputs k and R must
        have units. 

        Args: 

          `k` (real or array): wavenumber.

          `R` (real): filter scale. 

          `z` (real): redshift

        .. math::
          \frac{d\sigma^2}{dk} = \frac{ \Delta^2(k,z) W^2(k,R) }{ k }

        """ 

        if hasattr(k,'units'): 
            k.units = 'hh/Mpc'
        else:
            raise ValueError, '\n Input variable k must have units \n'

        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'

        kk = np.array(k)
        dsdk = self.Delta2k(k,z) * self.W2kR(k,R) / kk
        return dsdk


    def dsig2_dlnk( self, lnk, R, z ):
        r""" Input is ln(k) which is unitless but k must have units of h/Mpc 

        Args: 

          `lnk` (real or array): wavenumber, ln( k [h/Mpc] ).

          `R` (real): filter scale. 

          `z` (real): redshift

        .. math::
          \frac{d\sigma^2}{d \ln k} = \Delta^2(k,z) W^2(k,R)

        """
        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'
        
        kk = np.exp(1.0)**lnk / self.cu.Mpc / self.cu.h
        dsdlnk = self.Delta2k(kk,z) * self.W2kR(kk,R)
        return dsdlnk
    
    
    def dsig2_dlogk( self, logk, R, z ): 
        r""" Input is log10(k) which is unitless but k must have units of h/Mpc 

        Args: 

          `logk` (real or array): wavenumber, log10( k [h/Mpc] )

          `R` (real): filter scale. 

          `z` (real): redshift

        .. math::
          \frac{d\sigma^2}{d \log k} = \ln(10) \Delta^2(k,z) W^2(k,R)

        """
        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'
        
        kk = 10**logk / self.cu.Mpc / self.cu.h
        dsdlogk = np.log(10.) * self.Delta2k(kk,z) * self.W2kR(kk,R) 
        return dsdlogk
    





    def sigma2_R(self, R, z, method='romberg'):
        r""" Variance of density field smoothed on scale R 

        Args: 

          `R` (real): filter scale. 

          `z` (real): redshift

        .. math::
          \sigma^2(R) = \int \Delta^2(k,z) W^2(k,R) d \ln k

        """ 

        if hasattr(R,'units'): 
            R.units = 'Mpc/hh'
        else:
            raise ValueError, '\n Input variable R must have units \n'

        kw = np.array( 1.0/R )
        lnkw = np.log( kw )

        lnk_min = lnkw - 6.0
        lnk_max = lnkw + 6.0

        if method == 'romberg':
            sig2 = scipy.integrate.romberg( self.dsig2_dlnk, 
                                            lnk_min, 
                                            lnk_max,
                                            args=(R,z),
                                            tol=1.0e-12,
#                                            rtol=1.0e-12,
                                            divmax=20,
                                            vec_func=True )

        elif method == 'quad':
            sig2 = scipy.integrate.quad( self.dsig2_dlnk, 
                                         lnk_min, 
                                         lnk_max,
                                         args=(R,z),
                                         epsabs=1.0e-8,
                                         epsrel=1.0e-8 )[0]


        return sig2


    def sigma2_V(self, V, z):

        r""" Variance of density field smoothed on scale V 

        Args: 

          `V` (real): filter scale. 

          `z` (real): redshift

        .. math::

          R = \left( \frac{3 V}{4 \pi} \right)^{1/3},  \\

          \sigma^2(R) = \int \Delta^2(k,z) W^2(k,R) d \ln k

        """ 

        if hasattr(V,'units'): 
            V.units = 'Mpc^3/hh^3'
        else:
            raise ValueError, '\n Input variable V must have units \n'
        R = ( (3 * V) / (4 * np.pi) )**(1./3)
        total = self.sigma2_R(R,z)
        return total 


    def sigma2_M(self, M, z):

        r""" Variance of density field smoothed on scale V 

        Args: 

          `V` (real): filter scale. 

          `z` (real): redshift

        .. math::

          V = \frac{M}{\rho_{c,0} \, \Omega_M}, \, 
          R = \left( \frac{3 V}{4 \pi} \right)^{1/3},  \\

          \sigma^2(R) = \int \Delta^2(k,z) W^2(k,R) d \ln k

        """ 

        if hasattr(M,'units'): 
            M.units = 'Msun/hh'
        else:
            raise ValueError, '\n Input variable M must have units \n'
        rho_mean = self.cosmo.rho_crit0 * self.cosmo.OmegaM
        vol = M / rho_mean
        R = ( (3 * vol) / (4 * np.pi) )**(1./3)
        total = self.sigma2_R( R, z )
        return total 


    def n_eff_approx( self, vol, dlogM, z ):

        dlnM = dlogM * np.log(10.)

        mass = np.array(vol) * 10**(dlogM/2.) * self.cu.Msunh
        sig2 = self.sigma2_M( mass, z )
        sig = np.sqrt(sig2)
        t1 = 1.0 / sig

        mass = np.array(vol) * 10**(-dlogM/2.) * self.cu.Msunh
        sig2 = self.sigma2_M( mass, z )
        sig = np.sqrt(sig2)
        t2 = 1.0 / sig

        dlnisig = np.abs( np.log(t1) - np.log(t2) )

        neff = 6.0 * (dlnisig / dlnM) - 3.0

        return neff


    def frac_mass( self, sig, delta_c ):
        """ fraction of mass in halos per unit dln(1/sigma) """

        nu = delta_c / sig

        nu_prime = np.sqrt(0.707) * nu
        ln_isig = np.log( 1.0/sig )
        ln_gauss1 = np.exp( -(ln_isig-0.4)**2 / (2.0 * 0.6**2) )

        t1 = 0.3222 * np.sqrt(2.0/np.pi) * nu_prime
        t2 = np.exp( -1.08 * nu_prime**2 / 2.0 ) 
        t3 = (1.0 + 1.0/nu_prime**0.6 + 0.2*ln_gauss1)
        frac = t1 * t2 * t3
                
        return frac


#    def plot_W2kR( self ):
#        """ Plots the window function for R=8 Mpc/h """ 

#        Rw = 8.0 * self.CU.Mpc / self.CU.h
#        kw = np.array( 1.0/Rw )
#        lnkw = np.log( kw )
#        logkw = np.log10(kw)

#        N = 5000
#        lnk_min = lnkw - 2.0
#        lnk_max = lnkw + 4.0
#        kk_min = np.exp(1.0)**lnk_min
#        kk_max = np.exp(1.0)**lnk_max
#        logk_min = np.log10( kk_min )
#        logk_max = np.log10( kk_max )

#        lnk = np.linspace( lnk_min, lnk_max, N )
#        kk = np.exp(1.0)**lnk / self.CU.Mpc / self.CU.h
#        logk = np.log10( np.array(kk) )

#        fig = plt.figure( figsize=(10,8) )
#        fig.subplots_adjust( top=0.95, right = 0.95, 
#                             bottom=0.12, left=0.14 )
#        ax = fig.add_subplot(111)


#        WW = self.W2lnkR( lnk, Rw )
#        ax.plot( logk, np.log10(WW), lw=2.0, color='blue'  )
#        ax.plot( [logkw,logkw], [-10,10] )

#        ax.set_xlabel( r'$\log(k \, [\rm{h/Mpc}])$', fontsize=26 )

#        txt = r'$ \log \, [W_R(k)] $'
#        ax.set_ylabel( txt, fontsize=26 )

#        ax.set_xlim( (logk_min,logk_max) )
#        ax.set_ylim( (-6.0,0.5) )



#    def plot_dsig2_dk( self ):

#        """ Plot the sigma^2(R) integrand = Delta^2(k) [W(k,R)]^2 for several 
#        values of R """

#        mpl.rcParams['xtick.labelsize'] = 20
#        mpl.rcParams['ytick.labelsize'] = 20

#        z = 0.0
#        N = 100000

#        R_arr = [800.0, 80.0, 8.0, 0.8, 0.08] 
#        c_arr = ['red', 'orange', 'green', 'blue', 'purple']

#        fig = plt.figure( figsize=(10,8) )
#        fig.subplots_adjust( top=0.86, right = 0.95, 
#                             bottom=0.12, left=0.14 )
#        ax = fig.add_subplot(111)
        
        # first plot the unitless power spectrum
        #--------------------------------------------------
#        plt_lnk_min = -8.0
#        plt_lnk_max = 6.0
#        plt_kk_min = np.exp(1.0)**plt_lnk_min
#        plt_kk_max = np.exp(1.0)**plt_lnk_max
#        plt_logk_min = np.log10( plt_kk_min )
#        plt_logk_max = np.log10( plt_kk_max )
        
#        logk = np.linspace( plt_logk_min, plt_logk_max, N )
#        kk = 10**logk / self.CU.Mpc / self.CU.h
#        lnk = np.log( np.array(kk) )

#        D2k = self.Delta2k(kk, z)
#        ax.plot( logk, np.log10(D2k), lw=3.0, color='black'  )
        
#        for i in range(len(R_arr)):

#            color = c_arr[i]
#            Rw = R_arr[i] * self.CU.Mpc / self.CU.h
#            kw = np.array( 1.0/Rw )
#            logkw = np.log10( kw )
#            lnkw = np.log( kw )

#            lnk_min = lnkw - 4.0
#            lnk_max = lnkw + 8.0
#            kk_min = np.exp(1.0)**lnk_min
#            kk_max = np.exp(1.0)**lnk_max
#            logk_min = np.log10( kk_min )
#            logk_max = np.log10( kk_max )

#            logk = np.linspace( logk_min, logk_max, N )
#            kk = 10**logk
#            lnk = np.log( np.array(kk) )

#            dsig2_dlogk = self.dsig2_dlogk( logk, Rw, z )
#            dsig2_dlnk = self.dsig2_dlnk( lnk, Rw, z )

#            ax.plot( logk, np.log10(dsig2_dlnk), lw=1.5, color=color,
#                     alpha=0.5)

#            xx = (logkw - plt_logk_min) / (plt_logk_max - plt_logk_min) 

#            ax.arrow( xx, 0.98, 0.0, -0.05,
#                      transform=ax.transAxes, color=color )

#            txt = str(np.array(Rw)) + r' $ \rm {Mpc/h}$'
#            ax.text( 0.02, 0.82-0.08*i, txt,
#                     horizontalalignment='left',
#                     verticalalignment='center',
#                     transform=ax.transAxes, color=color,
#                     fontsize=20)

        
#        ax.grid()

#        ax2 = ax.twiny()
#        ax2.plot( [1,1], [1,1], alpha=0.0 )

#        xlo = np.log10( 1.0/np.exp(1.0)**plt_lnk_max )
#        xhi = np.log10( 1.0/np.exp(1.0)**plt_lnk_min )
#        ax2.set_xlim( xhi, xlo )

#        ax.set_xlim( (plt_logk_min, plt_logk_max) )
#        ax.set_ylim( (-7.4,1.7) )

#        ax.set_xlabel( r'$\log(k \, [\rm{h/Mpc}])$', fontsize=26 )

#        ax2.text( 0.5, 1.10, r'$\log(R \, [\rm{Mpc/h}])$',
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 transform=ax2.transAxes, color='black',
#                 fontsize=26)

#        txt = r'$ \log \, [\Delta^2(k) W_R(k)] = \log \,[$'
#        txt = txt + r'$ \frac{d \sigma^2_R}{ d \, \ln k }]  $'
#        ax.set_ylabel( txt, fontsize=26 )







