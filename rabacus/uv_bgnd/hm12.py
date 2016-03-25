""" Handles IO and standard spectra creation for files related to the Haardt 
and Madau 2012 UV Background http://adsabs.harvard.edu/abs/2012ApJ...746..125H.
Default units are ``eV`` and ``s`` for ease of comparison between this code 
and the results in the paper. """

import os.path
import numpy as np
from scipy.interpolate import interp1d

from rabacus.constants import physical
from rabacus.constants import units


__all__ = ['HM12_Photorates_Table', 'HM12_UVB_Table']


class HM12_Photorates_Table:

    """ Handles the tabulated HM12 photoionization and heating data from  
    Table 3 in http://adsabs.harvard.edu/abs/2012ApJ...746..125H.

    Attributes: 
      `u` (:class:`~rabacus.constants.units.Units`)
 
      `pc` (:class:`~rabacus.constants.physical.PhysicalConstants`) 

    """ 

    def __init__( self ):

        # create physical constants and units 
        #--------------------------------------------------------
        self.pc = physical.PhysicalConstants()
        self.u = units.Units()

        # read data
        #-----------------------------------------------------
        local = os.path.dirname(os.path.realpath(__file__))
        fname = local + '/hm12_dat/photorates.out'        
        dat = np.loadtxt( fname )
        self._dat = dat

        # organize data
        #-----------------------------------------------------
        self._z = dat[:,0]        
        self._H1i = dat[:,1]        # H1 photoionization
        self._H1h = dat[:,2]        # H1 photoheating
        self._He1i = dat[:,3]       # He1 photoionization
        self._He1h = dat[:,4]       # He1 photoheating
        self._He2i = dat[:,5]       # He2 photoionization
        self._He2h = dat[:,6]       # He2 photoheating
        self._comptonh = dat[:,7]   # compton heating

        # create interpolating functions
        #-----------------------------------------------------
        self.z = self._z.copy()
        self.l1pz = np.log10( 1.0 + self.z )
        self._H1i_fit = interp1d( self.l1pz, np.log10(self._H1i) )
        self._H1h_fit = interp1d( self.l1pz, np.log10(self._H1h) )
        self._He1i_fit = interp1d( self.l1pz, np.log10(self._He1i) )
        self._He1h_fit = interp1d( self.l1pz, np.log10(self._He1h) )
        self._He2i_fit = interp1d( self.l1pz, np.log10(self._He2i) )
        self._He2h_fit = interp1d( self.l1pz, np.log10(self._He2h) )
        self._comptonh_fit = interp1d( self.l1pz, np.log10(self._comptonh) )


    def H1i(self,z):
        """ H1 photoionization rate as a function of `z`. 
        
        Args: 
          `z` (float): redshift
           
        Returns:
          `H1i` (float): H1 photoionization rate

        """
        l1pz = np.log10( 1.0 + z )
        H1i = self._H1i_fit(l1pz) 
        H1i = 10**H1i / self.u.s
        return H1i

    def H1h(self,z):
        """ H1 photoheating rate as a function of `z`. 
        
        Args:
          `z` (float): redshift
           
        Returns:
          `H1h` (float): H1 photoheating rate

        """
        l1pz = np.log10( 1.0 + z )
        H1h = self._H1h_fit(l1pz) 
        H1h = 10**H1h * self.u.eV / self.u.s
        return H1h

    def He1i(self,z):
        """ He1 photoionization rate as a function of `z`. 
        
        Args: 
          `z` (float): redshift
           
        Returns:
          `He1i` (float): He1 photoionization rate

        """
        l1pz = np.log10( 1.0 + z )
        He1i = self._He1i_fit(l1pz) 
        He1i = 10**He1i / self.u.s
        return He1i

    def He1h(self,z):
        """ He1 photoheating rate as a function of `z`. 
        
        Args: 
          `z` (float): redshift
           
        Returns:
          `He1h` (float): He1 photoheating rate

        """
        l1pz = np.log10( 1.0 + z )
        He1h = self._He1h_fit(l1pz) 
        He1h = 10**He1h * self.u.eV / self.u.s
        return He1h

    def He2i(self,z):
        """ He2 photoionization rate as a function of `z`. 
        
        Args: 
          `z` (float): redshift
           
        Returns:
          `He2i` (float): He2 photoionization rate

        """
        l1pz = np.log10( 1.0 + z )
        He2i = self._He2i_fit(l1pz) 
        He2i = 10**He2i / self.u.s
        return He2i

    def He2h(self,z):
        """ He2 photoheating rate as a function of `z`. 
        
        Args:
          `z` (float): redshift
           
        Returns:
          `He2h` (float): He2 photoheating rate

        """
        l1pz = np.log10( 1.0 + z )
        He2h = self._He2h_fit(l1pz) 
        He2h = 10**He2h * self.u.eV / self.u.s
        return He2h

    def comptonh(self,z):
        """ Compton heating rate as a function of `z`. 
        
        Args:
          `z` (float): redshift
           
        Returns:
          `comptonh` (flaot): Compton heating rate

        """
        l1pz = np.log10( 1.0 + z )
        comptonh = self._comptonh_fit(l1pz) 
        comptonh = 10**comptonh * self.u.eV / self.u.s
        return comptonh


    


class HM12_UVB_Table:

    """ Handles the tabulated HM12 spectrum. 

    Attributes: 
      `U` (:class:`~rabacus.constants.units.Units`)
 
      `PC` (:class:`~rabacus.constants.physical.PhysicalConstants`) 

      `z` (array): 1-D array [Nz] redshift

      `Inu` (array): 2-D array [Nnu, Nz] specific intensity 
      [erg/(cm**2 s sr Hz)]

      `nu` (array): 1-D array [Nnu] photon frequency

      `lam` (array): 1-D array [Nnu] photon wavelength

      `E` (array): 1-D array [Nnu] photon energy


    """ 

    def __init__( self ):

        # create physical constants and units 
        #--------------------------------------------------------
        self.pc = physical.PhysicalConstants()
        self.u = units.Units()

        # read data
        #-----------------------------------------------------
        local = os.path.dirname(os.path.realpath(__file__))
        fname = local + '/hm12_dat/UVB.out'
        with open(fname, 'r') as f:
            dum = f.readlines()
        dum = dum[20].split()
        self.z = np.array( [float(i) for i in dum] )

        dat = np.loadtxt( fname, skiprows=21 )
        self._dat = dat

        Iunit = self.u.erg / (self.u.s * self.u.cm**2 * self.u.Hz * self.u.sr)
        self.Inu = dat[:,1:] * Iunit

        Inu_mag = self.Inu.magnitude
        isfinite = Inu_mag > 0
        Imin = Inu_mag[isfinite].min()
        
        self.logInu = Inu_mag.copy() 
        self.logInu[isfinite] = np.log10( Inu_mag[isfinite] )
        self.logInu[~isfinite] = np.log10( Imin )

        self.lam = dat[:,0] * self.u.angstrom
        self.nu = self.pc.c / self.lam
        self.E = self.nu * self.pc.h

        self.lam.units = 'cm'
        self.nu.units = 'Hz' 
        self.E.units = 'eV'



    def return_Inu_intrp( self, z ): 
        """ Generates an interpolating function giving log Inu for a given
        log lamda at the requested redshift. 

        Args:
          `z` (float): requested redshift 

        Returns:
          `intrpf` (function): interpolating function giving log Inu for a 
          given log lambda

        """ 

        assert isinstance( z, float )

        # find redshift indices
        #-----------------------------------------------------------------
        i_lo = np.argmin( np.abs( self.z - z ) )
        if self.z[i_lo] < z:
            i_hi = i_lo + 1
        else:
            i_hi = i_lo
            i_lo = i_lo - 1            
        z_frac = (z - self.z[i_lo]) / (self.z[i_hi] - self.z[i_lo])

        # create interpolating function
        #-----------------------------------------------------------------
        dlog_Inu = self.logInu[:,i_hi] - self.logInu[:,i_lo]  
        log_Inu = self.logInu[:,i_lo] + z_frac * dlog_Inu  

        self.lam.units = 'cm' 
        log_lam_cm = np.log10( self.lam.magnitude )
        intrpf = interp1d( log_lam_cm, log_Inu )

        return intrpf


    def return_spectrum_lam( self, z, lam ): 
        """ Interpolates the HM12 tabulated spectrum in redshift at the 
        requested wavelengths. 

        Args: 
          `z` (float): requested redshift 

          `lam` (array): requested wavelengths

        Returns:

          `Inu` (array): specific intensity at requested wavelengths

        """ 

        # get interpolating function
        #-----------------------------------------------------------------
        intrpf = self.return_Inu_intrp( z )

        # evaluate at requested wavelengths
        #-----------------------------------------------------------------
        lam.units = 'cm' 
        log_lam_cm = np.log10( lam.magnitude )
        log_Inu = intrpf( log_lam_cm )
        Inu = 10**log_Inu * self.Inu.units

        return Inu


    def return_spectrum_nu( self, z, nu ): 
        """ Interpolates the HM12 tabulated spectrum in redshift at the 
        requested frequencies. 

        Args: 
          `z` (float): requested redshift

          `nu` (array): requested frequencies

        Returns:
          `Inu` (array): specific intensity at requested frequencies

        """ 

        # get interpolating function
        #-----------------------------------------------------------------
        intrpf = self.return_Inu_intrp( z )

        # evaluate at requested wavelengths
        #-----------------------------------------------------------------
        nu.units = 'Hz'  # just to make sure we have an frequency
        lam = self.pc.c / nu
        lam.units = 'cm' 
        log_lam_cm = np.log10( lam.magnitude )
        log_Inu = intrpf( log_lam_cm )
        Inu = 10**log_Inu * self.Inu.units

        return Inu


    def return_spectrum_E( self, z, E ): 
        """ Interpolates the HM12 tabulated spectrum in redshift at the 
        requested energies. 

        Args: 
          `z` (float): requested redshift

          `E` (array): requested energies

        Returns:
          `Inu` (array): specific intensity at requested energies

        """ 

        # get interpolating function
        #-----------------------------------------------------------------
        intrpf = self.return_Inu_intrp( z )

        # evaluate at requested energies
        #-----------------------------------------------------------------
        E.units = 'eV'    # just to make sure we have an energy
        lam = self.pc.c / ( E / self.pc.h )
        lam.units = 'cm' 
        log_lam_cm = np.log10( lam.magnitude )
        log_Inu = intrpf( log_lam_cm )
        Inu = 10**log_Inu * self.Inu.units

        return Inu


