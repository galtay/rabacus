""" An isotropic background radiation source module.  
""" 


import numpy as np
from rabacus.utils import utils
from rabacus.constants import physical
from rabacus.constants import units 
from rabacus.hdf5 import h5py_wrap as h5
from rabacus.utils import units_string

from hm12 import HM12Spectrum
from thermal import ThermalSpectrum
from powerlaw import PowerlawSpectrum
from monochromatic import MonochromaticSpectrum
from source import Source


import species


__all__ = ['BackgroundSource']



class BackgroundSource(Source):


    r""" An isotropic background radiation source class.  

    Args: 

      `q_min` (float): minimum photon energy / Rydbergs [dimensionless]

      `q_max` (float): maximum photon energy / Rydbergs [dimensionless]

      `spectrum_type` (str): the spectral shape 
      {``monochromatic``, ``hm12``, ``thermal``, ``powerlaw``, ``user``}
      
    Kwargs:
    
      `Nnu` (int): number of spectral samples, log spaced in energy

      `segmented` (bool): if ``True``, forces spectral samples at H and He 
      ionization thresholds

      `px_fit_type` (str): source to use for photoionization cross section 
      fits {``verner``}
       
      `verbose` (bool): verbose screen output?

      `z` (float): redshift, need if `spectrum_type` = ``hm12``

      `T_eff` (float): effective temperature, need if `spectrum_type` = 
      ``thermal``

      `alpha` (float): powerlaw index, need if `spectrum_type` = ``powerlaw``

      `user_E` (array): energy samples for user defined spectrum. should
      have units of energy. 

      `user_shape` (array): shape of user defined spectrum. should be 
      dimensionless

    Attributes:

      `U` (:class:`~rabacus.constants.units.Units`) 

      `PC` (:class:`~rabacus.constants.physical.PhysicalConstants`)
 
      `H` (:class:`~rabacus.atomic.hydrogen.Hydrogen`)
      
      `He` (:class:`~rabacus.atomic.helium.Helium`)
      
      `PX` (:class:`~rabacus.atomic.photo_xsection.PhotoXsections`)


      `E` (array): Spectrum energy samples

      `lam` (array): Spectrum wavelength samples 

      `nu` (array): Spectrum frequency samples 
      
      `q` (array): Spectrum energy samples / Rydbergs
      
      `Inu` (array): Specific intensity |Inu| 


      `monochromatic` (bool): ``True`` for monochromatic spectra otherwise 
      ``False``

      `grey` (object):  Grey photoionization cross sections and energies.
      Only created if `segmented` = ``True``.

      `log` (object): Convenience object storing the log of other attributes

      `sigma` (object): photoionization cross sections and ratios

      `th` (object): all quantities evaluated at ionization thresholds

      `thin` (object): all optically thin quantities


    Notes:


      The fundamental characterization of the radiation field for this class 
      is the specific intensity |Inu| with units of ``erg/(s cm^2 sr Hz)``
      or units of ``erg/(s cm^2 sr)`` for monochromatic spectra. Integrating 
      over frequency and solid angle yields several other useful 
      characterizations of the radiation field which we summarize below. 
    
      .. math:: F_{u} = 4 \pi \int I_{\nu} d\nu \rightarrow
                \frac{dF_{u}}{d\nu} = 4 \pi I_{\nu}

      .. math:: u = \frac{4 \pi}{c} \int I_{\nu} d\nu \rightarrow
                \frac{du}{d\nu} = \frac{4 \pi}{c} I_{\nu}

      .. math:: F_{n} = 4 \pi \int \frac{I_{\nu}}{h\nu} d\nu \rightarrow
                \frac{dF_{n}}{d\nu} = 4 \pi \frac{I_{\nu}}{h\nu}

      .. math:: n = \frac{4 \pi}{c} \int \frac{I_{\nu}}{h\nu} d\nu \rightarrow
                \frac{dn}{d\nu} = \frac{4 \pi}{c} \frac{I_{\nu}}{h\nu}


      +-----------+--------------------+----------------------+--------+
      | Variable  | Description        | Units                | Symbol |
      +===========+====================+======================+========+
      | Inu       | specific intensity | [erg/(s cm^2 sr Hz)] | |Inu|  |
      +-----------+--------------------+----------------------+--------+
      | thin.Fu   | energy flux        | [erg/(s cm^2)]       | |Fu|   |
      +-----------+--------------------+----------------------+--------+
      | thin.u    | energy density     | [erg/(cm^3)]         | |u|    |
      +-----------+--------------------+----------------------+--------+
      | thin.Fn   | photon flux        | [1/(s cm^2)]         | |Fn|   |
      +-----------+--------------------+----------------------+--------+
      | thin.n    | photon density     | [1/(cm^3)]           | |n|    |
      +-----------+--------------------+----------------------+--------+

      The photoionization and photoheating rates are also integrals over the 
      specific intensity. 

      .. math:: \Gamma_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h \nu} d\nu \rightarrow
                \frac{ d\Gamma_{\rm _{X}} }{d\nu} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h \nu}

      .. math:: H_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h \nu} 
                h(\nu - \nu_{\rm _{X}}) d\nu \rightarrow
                \frac{ dH_{\rm _{X}} }{d\nu} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h \nu} h(\nu-\nu_{\rm _{X}})

      where :math:`\scriptstyle{\rm X}` is one of 
      :math:`\{ {\scriptstyle{\rm HI, \, HeI, \, HeII }} \}`.  
      Note that the :math:`\sigma_{\rm _{X}}` are frequency dependent 
      photoionization cross sections and the :math:`\nu_{\rm _{X}}` are  
      (scalar) frequencies that correspond to the ionization energies for a 
      given species. In many applications, we prefer to use energy as a 
      variable instead of frequency.  The above equations can easily be recast 
      with a change of variable :math:`E=h\nu`. 

      .. math:: \Gamma_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} dE \rightarrow
                \frac{ d\Gamma_{\rm _{X}} }{dE} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E}

      .. math:: H_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} 
                (E - E_{\rm _{X}}) dE \rightarrow
                \frac{ dH_{\rm _{X}} }{dE} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} (E - E_{\rm _{X}})

      where the :math:`E_{\rm _{X}}` are (scalar) ionization energies. 

      +-----------+---------------------+----------------------+--------+
      | Variable  | Description         | Units                | Symbol |
      +===========+=====================+======================+========+
      | thin.H1i  | H1 photo ion. rate  | [1/s]                | |H1i|  |
      +-----------+---------------------+----------------------+--------+
      | thin.He1i | He1 photo ion. rate | [1/s]                | |He1i| |
      +-----------+---------------------+----------------------+--------+
      | thin.He2i | He2 photo ion. rate | [1/s]                | |He2i| |
      +-----------+---------------------+----------------------+--------+
      | thin.H1h  | H1 photo heat rate  | [erg/s]              | |H1h|  |
      +-----------+---------------------+----------------------+--------+
      | thin.He1h | He1 photo heat rate | [erg/s]              | |He1h| |
      +-----------+---------------------+----------------------+--------+
      | thin.He2h | He2 photo heat rate | [erg/s]              | |He2h| |
      +-----------+---------------------+----------------------+--------+

      In addition, the class provides methods for calculating the shielded 
      photoionization/heating rates given an optical depth.  

      .. math:: \Gamma_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} 
                e^{-\tau} dE \rightarrow
                \frac{ d\Gamma_{\rm _{X}} }{dE} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} e^{-\tau}

      .. math:: H_{\rm _{X}} = 4 \pi \int 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} 
                (E - E_{\rm _{X}}) e^{-\tau} dE \rightarrow
                \frac{ dH_{\rm _{X}} }{dE} = 4 \pi 
                \frac{ I_{\nu} \sigma_{\rm _{X}} }{h E} 
                (E - E_{\rm _{X}}) e^{-\tau}

      where :math:`\tau = \sum_{\rm _X} \tau_{\rm _X}`, 
      :math:`\tau_{\rm _X} = N_{\rm _X} \sigma_{\rm _X}`,
      :math:`N_{\rm _X}` are column densities, and
      :math:`\sigma_{\rm _X}` are energy dependent photoionization cross 
      sections. 


    .. note::

      The :math:`\tau_{\rm _X}` are energy dependent optical depths while the 
      parameters to the shielding functions (:func:`shld_H1i`, 
      :func:`shld_H1h`, :func:`shld_He1i`, :func:`shld_He1h`, :func:`shld_He2i`,
      :func:`shld_He2h`) are the :math:`\tau_{\rm _X}` evaluated at the 
      species threhold energy (i.e. scalars). 


    +-----------+---------------------+----------------------+--------+
    | Function  | Description         | Units                | Symbol |
    +===========+=====================+======================+========+
    | shld_H1i  | H1 photo ion. rate  | [1/s]                | |H1i|  |
    +-----------+---------------------+----------------------+--------+
    | shld_He1i | He1 photo ion. rate | [1/s]                | |He1i| |
    +-----------+---------------------+----------------------+--------+
    | shld_He2i | He2 photo ion. rate | [1/s]                | |He2i| |
    +-----------+---------------------+----------------------+--------+
    | shld_H1h  | H1 photo heat rate  | [erg/s]              | |H1h|  |
    +-----------+---------------------+----------------------+--------+
    | shld_He1h | He1 photo heat rate | [erg/s]              | |He1h| |
    +-----------+---------------------+----------------------+--------+
    | shld_He2h | He2 photo heat rate | [erg/s]              | |He2h| |
    +-----------+---------------------+----------------------+--------+

    .. seealso::

      :class:`~rabacus.rad_src.plane.PlaneSource`, 
      :class:`~rabacus.rad_src.point.PointSource`


    """ 


    def __init__( self, 
                  q_min, q_max, spectrum_type, 
                  Nnu = 128, 
                  segmented = True, # puts spectrum points at E_th
                  px_fit_type = 'verner',
                  verbose = False,
                  z = None,         # only for spectrum_type = 'hm12'
                  T_eff = None,     # only for spectrum_type = 'thermal'
                  alpha = None,      # only for spectrum_type = 'powerlaw'
                  user_E = None,      # photon E array for spectrum_type='user'
                  user_shape = None,  # shape for spectrum_type='user'
                  ):   


        # attach specific input background
        #-----------------------------------------------------------------
        self.source_type = 'background' 


        # do generic initialization 
        #-----------------------------------------------------------------
        self.initialize( q_min, q_max, spectrum_type, Nnu, segmented, 
                         px_fit_type, verbose, z, T_eff, alpha, user_E, 
                         user_shape )


        # set spectral shape (from yval)
        #-----------------------------------------------------------------
        if spectrum_type == 'thermal':
            spec = ThermalSpectrum( self, self.T_eff ) 

        elif spectrum_type == 'powerlaw':
            spec = PowerlawSpectrum( self, self.alpha )

        elif spectrum_type == 'hm12':
            spec = HM12Spectrum( self, self.z )
            
        elif spectrum_type == 'monochromatic':
            spec = MonochromaticSpectrum( self )
        

        # set primary spectrum shape from spec.yvals
        #-----------------------------------------------------------------
        num = self.U.erg

        if self.monochromatic:
            den = self.U.s * self.U.cm**2 * self.U.sr 
            units = num / den
            self.Inu = np.ones( self.Nnu ) * spec.yvals * units
        else:
            den = self.U.Hz  * self.U.s * self.U.cm**2 * self.U.sr 
            units = num / den
            if spectrum_type == 'user':
                self.Inu = np.ones( self.Nnu ) * user_shape * units
            else:
                self.Inu = np.ones( self.Nnu ) * spec.yvals * units 
            self.Inu.__doc__ = 'specific intensity' 

            
        self.log.Inu = np.log10( self.Inu.magnitude )


        # set integrals
        #-----------------------------------------------------------------
        self.set_integrands()
        self.thin = OpticallyThinQuantities( self )


        # calculate grey X-sections and energies
        # (note segmented = False if monochromatic = True)
        #-----------------------------------------------------------------
        if self.segmented:
            self.grey = GreyQuantities( self )









    def set_integrands( self ):

        """ Defines tabulated functions of photon energy/frequency which 
        represent integrands for optically thin quantities. All attributes 
        defined in this function are integrands of the fundamental quantity, 
        specific intensity or Inu which has units [erg/(s Hz cm^2 sr)].

        """

        # dont need to set these for monochromatic spectra
        #--------------------------------------------------------
        if self.monochromatic:
            return

        four_pi_sr = 4 * np.pi * self.U.sr


        # first we define differentials wrt frequency.  these are 
        # just for convenience. 
        #--------------------------------------------------------

        self._dFu_over_dnu = four_pi_sr * self.Inu 
        self._dFu_over_dnu.units = 'erg/s/cm^2/Hz'

        self._du_over_dnu = self._dFu_over_dnu / self.PC.c
        self._du_over_dnu.units = 'erg/cm^3/Hz'

        self._dFn_over_dnu = self._dFu_over_dnu / ( self.PC.h * self.nu )
        self._dFn_over_dnu.units = '1/s/cm^2/Hz'

        self._dn_over_dnu = self._dFn_over_dnu / self.PC.c
        self._dn_over_dnu.units = '1/cm^3/Hz'

        #--------------------------------------------------------

        self._dH1i_over_dnu = self._dFn_over_dnu * self.sigma.H1
        self._dH1i_over_dnu.units = '1/s/Hz'
        self._dH1h_over_dnu = self._dH1i_over_dnu * (self.E - self.th.E_H1)
        self._dH1h_over_dnu.units = 'erg/s/Hz'

        self._dHe1i_over_dnu = self._dFn_over_dnu * self.sigma.He1
        self._dHe1i_over_dnu.units = '1/s/Hz'
        self._dHe1h_over_dnu = self._dHe1i_over_dnu * (self.E - self.th.E_He1)
        self._dHe1h_over_dnu.units = 'erg/s/Hz'

        self._dHe2i_over_dnu = self._dFn_over_dnu * self.sigma.He2
        self._dHe2i_over_dnu.units = '1/s/Hz'
        self._dHe2h_over_dnu = self._dHe2i_over_dnu * (self.E - self.th.E_He2)
        self._dHe2h_over_dnu.units = 'erg/s/Hz'



        # now we define differentials wrt energy 
        #--------------------------------------------------------

        self._dFu_over_dE = self._dFu_over_dnu / self.PC.h
        self._dFu_over_dE.units = 'erg/s/cm^2/eV'

        self._du_over_dE = self._dFu_over_dE / self.PC.c
        self._du_over_dE.units = 'erg/cm^3/eV'

        self._dFn_over_dE = self._dFu_over_dnu / ( self.PC.h * self.E )
        self._dFn_over_dE.units = '1/s/cm^2/eV'

        self._dn_over_dE = self._dFn_over_dE / self.PC.c
        self._dn_over_dE.units = '1/cm^3/eV'

        #--------------------------------------------------------

        self._dH1i_over_dE = self._dFn_over_dE * self.sigma.H1
        self._dH1i_over_dE.units = '1/s/eV'        
        self._dH1h_over_dE = self._dH1i_over_dE * (self.E - self.th.E_H1)   
        self._dH1h_over_dE.units = 'erg/s/eV'        

        self._dHe1i_over_dE = self._dFn_over_dE * self.sigma.He1
        self._dHe1i_over_dE.units = '1/s/eV'
        self._dHe1h_over_dE = self._dHe1i_over_dE * (self.E - self.th.E_He1)
        self._dHe1h_over_dE.units = 'erg/s/eV'        

        self._dHe2i_over_dE = self._dFn_over_dE * self.sigma.He2
        self._dHe2i_over_dE.units = '1/s/eV'        
        self._dHe2h_over_dE = self._dHe2i_over_dE * (self.E - self.th.E_He2)
        self._dHe2h_over_dE.units = 'erg/s/eV'





        



# functions that calculate shielded integrals 
#====================================================================

    def return_attenuation( self, tauH1_th, tauHe1_th, tauHe2_th ):
        """ 
        Calculate the attenuation array.  First the energy dependence of tau
        for each absorbing species (H1, He1, He2) is calculated using the input
        scalar optical depths at the ionization thresholds.  Next an energy 
        dependent total tau is calculated and exponentiated to arrive at the 
        attenuation array atten = np.exp(-tau). 

        Args:

          `tauH1_th` (float): H1 optical depth at the H1 ionizing threshold

          `tauHe1_th` (float): He1 optical depth at the He1 ionizing threshold

          `tauHe2_th` (float): He2 optical depth at the He2 ionizing threshold

        Returns:

          `atten` (array): attenuation as a function of energy

        """

        assert isinstance(tauH1_th, float)
        assert isinstance(tauHe1_th, float)
        assert isinstance(tauHe2_th, float)

        tauH1 = tauH1_th * self.sigma.H1_ratio
        tauHe1 = tauHe1_th * self.sigma.He1_ratio
        tauHe2 = tauHe2_th * self.sigma.He2_ratio
        tau = tauH1 + tauHe1 + tauHe2
        atten = np.exp( -tau )
        assert atten.size == self.E.size

        return atten


    def shld_H1i(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HI photoionization rate
        after passing through a column with optical depth tau 
        
        Args:
          `tauH1_th` (float): H1 optical depth at the H1 ionizing threshold
        
          `tauHe1_th` (float): He1 optical depth at the He1 ionizing threshold

          `tauHe2_th` (float): He2 optical depth at the He2 ionizing threshold

        Returns:
          `H1i` (float): attenuated H1 photoionization rate 
           
        """ 
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            H1i = self.thin.H1i * atten
        else:
            H1i = utils.trap( self.E, self._dH1i_over_dE * atten )
        H1i.units = '1/s'
        return H1i

    def shld_H1h(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HI photoheating rate
        after passing through a column with optical depth tau 

        Args:
          `tauH1_th` (float): H1 optical depth at the H1 ionizing threshold
        
          `tauHe1_th` (float): He1 optical depth at the He1 ionizing threshold

          `tauHe2_th` (float): He2 optical depth at the He2 ionizing threshold

        Returns:
          `H1h` (float): attenuated H1 photoheating rate 
           
        """ 
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            H1h = self.thin.H1h * atten
        else:
            H1h = utils.trap( self.E, self._dH1h_over_dE * atten )
        H1h.units = 'erg/s'
        return H1h



    def shld_He1i(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeI photoionization rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1i`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He1i = self.thin.He1i * atten
        else:
            He1i = utils.trap( self.E, self._dHe1i_over_dE * atten )
        He1i.units = '1/s'
        return He1i

    def shld_He1h(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeI photoheating rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1h`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He1h = self.thin.He1h * atten
        else:
            He1h = utils.trap( self.E, self._dHe1h_over_dE * atten )
        He1h.units = 'erg/s'
        return He1h



    def shld_He2i(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeII photoionization rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1i`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He2i = self.thin.He2i * atten
        else:
            He2i = utils.trap( self.E, self._dHe2i_over_dE * atten )
        He2i.units = '1/s'
        return He2i

    def shld_He2h(self, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeII photoheating rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1h`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He2h = self.thin.He2h * atten
        else:
            He2h = utils.trap( self.E, self._dHe2h_over_dE * atten )
        He2h.units = 'erg/s'
        return He2h



# functions that scale and normalize the spectrum 
#====================================================================

    def scale_spectrum( self, fac ): 
        """ Changes normalization of spectrum and recalculates spectral 
        integrands. 

        Args:
          `fac` (float): multiplicative factor
        """ 

        self.Inu *= fac
        self.log_Inu = np.log10( self.Inu.magnitude )

        self.set_integrands()
        self.thin = OpticallyThinQuantities( self )



    def normalize_n( self, n ):
        """ Normalize the spectrum such that the number density of photons
        between `q_min` and `q_max` is `n`. 

        Args:
          `n` (float): target photon number density 
        """ 

        if not hasattr(n,'units'): 
            raise utils.NeedUnitsError, '\n n must have units, e.g. cm^-3 \n'
        else:
            n.units = 'cm^-3'

        fac = n / self.thin.n
        self.scale_spectrum( fac )
        

    def normalize_H1i( self, H1i ):
        """ Normalize the spectrum such that the H1 photoionization rate 
        integral yields `H1i`. 

        Args:
          `H1i` (float): target photoionization rate
        """ 

        if not hasattr(H1i,'units'): 
            raise utils.NeedUnitsError, '\n H1i must have units, e.g. s^-1 \n'
        else:
            H1i.units = 's^-1'

        fac = H1i / self.thin.H1i
        self.scale_spectrum( fac )
        


    def write_to_file( self, fname, top_group='/' ):

        # write out main group
        #-----------------------------------------------------
        tg = top_group
        h5.cg( fname, tg )

        h5.wa( fname, tg, 'q_min', self.q_min )
        h5.wa( fname, tg, 'q_max', self.q_max )
        h5.wa( fname, tg, 'Nnu', self.Nnu )
        h5.wa( fname, tg, 'monochromatic', self.monochromatic )
        
        attrs = {
            'units': units_string(self.q), 
            'description': 'energy samples / Ry'}
        h5.wd( fname, tg, 'q', self.q, attrs )

        attrs = {
            'units': units_string(self.E),
            'description': 'energy samples'}
        h5.wd( fname, tg, 'E', self.E, attrs )

        attrs = {
            'units': units_string(self.nu),
            'description': 'frequency samples'}
        h5.wd( fname, tg, 'nu', self.nu, attrs )

        attrs = {
            'units': units_string(self.lam),
            'description': 'wavelength samples'}
        h5.wd( fname, tg, 'lam', self.lam, attrs )

        attrs = {
            'units': units_string(self.Inu),
            'description': 'specific intensity'}
        h5.wd( fname, tg, 'Inu', self.Inu, attrs )

        # write out optically thin quantities
        #-----------------------------------------------------
        h5.cg( fname, tg+'/thin' )

        attrs = {
            'units': units_string(self.thin.Fn), 
            'description': 'photon flux'}
        h5.wd( fname, tg+'/thin', 'Fn', self.thin.Fn, attrs )

        attrs = {
            'units': units_string(self.thin.Fu), 
            'description': 'energy flux'}
        h5.wd( fname, tg+'/thin', 'Fu', self.thin.Fu, attrs )

        attrs = {
            'units': units_string(self.thin.n), 
            'description': 'photon number density'}
        h5.wd( fname, tg+'/thin', 'n', self.thin.n, attrs )

        attrs = {
            'units': units_string(self.thin.u), 
            'description': 'photon energy density'}
        h5.wd( fname, tg+'/thin', 'u', self.thin.u, attrs )

        attrs = {
            'units': units_string(self.thin.H1i), 
            'description': 'H I photo-ionization rate'}
        h5.wd( fname, tg+'/thin', 'H1i', self.thin.H1i, attrs )

        attrs = {
            'units': units_string(self.thin.He1i), 
            'description': 'He I photo-ionization rate'}
        h5.wd( fname, tg+'/thin', 'He1i', self.thin.He1i, attrs )

        attrs = {
            'units': units_string(self.thin.He2i), 
            'description': 'He II photo-ionization rate'}
        h5.wd( fname, tg+'/thin', 'He2i', self.thin.He2i, attrs )

        attrs = {
            'units': units_string(self.thin.H1h), 
            'description': 'H I photo-heating rate'}
        h5.wd( fname, tg+'/thin', 'H1h', self.thin.H1h, attrs )

        attrs = {
            'units': units_string(self.thin.He1h), 
            'description': 'He I photo-heating rate'}
        h5.wd( fname, tg+'/thin', 'He1h', self.thin.He1h, attrs )

        attrs = {
            'units': units_string(self.thin.He2h), 
            'description': 'He II photo-heating rate'}
        h5.wd( fname, tg+'/thin', 'He2h', self.thin.He2h, attrs )


        # write out photo-ionization cross-sections
        #-----------------------------------------------------
        h5.cg( fname, tg+'/sigma' )

        attrs = {
            'units': units_string(self.sigma.H1), 
            'description': 'H I photo-ionization cross-section'}
        h5.wd( fname, tg+'/sigma', 'H1', self.sigma.H1, attrs )

        attrs = {
            'units': units_string(self.sigma.He1), 
            'description': 'He I photo-ionization cross-section'}
        h5.wd( fname, tg+'/sigma', 'He1', self.sigma.He1, attrs )

        attrs = {
            'units': units_string(self.sigma.He2), 
            'description': 'He II photo-ionization cross-section'}
        h5.wd( fname, tg+'/sigma', 'He2', self.sigma.He2, attrs )


        # write out threshold values
        #-----------------------------------------------------
        h5.cg( fname, tg+'/th' )

        attrs = {
            'units': units_string(self.th.E_H1), 
            'description': 'H I ionization energy'}
        h5.wd( fname, tg+'/th', 'E_H1', self.th.E_H1, attrs )

        attrs = {
            'units': units_string(self.th.E_He1), 
            'description': 'He I ionization energy'}
        h5.wd( fname, tg+'/th', 'E_He1', self.th.E_He1, attrs )

        attrs = {
            'units': units_string(self.th.E_He2), 
            'description': 'He II ionization energy'}
        h5.wd( fname, tg+'/th', 'E_He2', self.th.E_He2, attrs )


        attrs = {
            'units': units_string(self.th.q_H1), 
            'description': 'H I ionization energy / Ry'}
        h5.wd( fname, tg+'/th', 'q_H1', self.th.q_H1, attrs )

        attrs = {
            'units': units_string(self.th.q_He1), 
            'description': 'He I ionization energy / Ry'}
        h5.wd( fname, tg+'/th', 'q_He1', self.th.q_He1, attrs )

        attrs = {
            'units': units_string(self.th.q_He2), 
            'description': 'He II ionization energy / Ry'}
        h5.wd( fname, tg+'/th', 'q_He2', self.th.q_He2, attrs )


        txt = 'photo-ionization cross-section @ ionization energy'
        
        attrs = {
            'units': units_string(self.th.sigma_H1), 
            'description': 'H I ' + txt}
        h5.wd( fname, tg+'/th', 'sigma_H1', self.th.sigma_H1, attrs )

        attrs = {
            'units': units_string(self.th.sigma_He1), 
            'description': 'He I ' + txt}
        h5.wd( fname, tg+'/th', 'sigma_He1', self.th.sigma_He1, attrs )

        attrs = {
            'units': units_string(self.th.sigma_He2), 
            'description': 'He II ' + txt}
        h5.wd( fname, tg+'/th', 'sigma_He2', self.th.sigma_He2, attrs )


        if self.segmented:

            attrs = {
                'description': 'index of energy sample at H I threshold'}
            h5.wd( fname, tg+'/th', 'i_H1', self.th.i_H1, attrs )

            attrs = {
                'description': 'index of energy sample at He I threshold'}
            h5.wd( fname, tg+'/th', 'i_He1', self.th.i_He1, attrs )

            attrs = {
                'description': 'index of energy sample at He II threshold'}
            h5.wd( fname, tg+'/th', 'i_He2', self.th.i_He2, attrs )


        # write out grey quantities if they exist
        #-----------------------------------------------------
        if hasattr( self, 'grey' ):
            h5.cg( fname, tg+'/grey' )

            attrs = {
                'units': units_string(self.grey.E.H1), 
                'description': 'H I grey energy'}
            h5.wd( fname, tg+'/grey', 'E_H1', self.grey.E.H1, attrs )

            attrs = {
                'units': units_string(self.grey.E.He1), 
                'description': 'He I grey energy'}
            h5.wd( fname, tg+'/grey', 'E_He1', self.grey.E.He1, attrs )

            attrs = {
                'units': units_string(self.grey.E.He2), 
                'description': 'He II grey energy'}
            h5.wd( fname, tg+'/grey', 'E_He2', self.grey.E.He2, attrs )

            attrs = {
                'units': units_string(self.grey.sigma.H1), 
                'description': 'H I grey photo-ionization cross-section'}
            h5.wd( fname, tg+'/grey', 'sigma_H1', self.grey.sigma.H1, attrs )

            attrs = {
                'units': units_string(self.grey.sigma.He1), 
                'description': 'He I photo-ionization cross-section'}
            h5.wd( fname, tg+'/grey', 'sigma_He1', self.grey.sigma.He1, attrs )

            attrs = {
                'units': units_string(self.grey.sigma.He2), 
                'description': 'He II photo-ionization cross-section'}
            h5.wd( fname, tg+'/grey', 'sigma_He2', self.grey.sigma.He2, attrs )



# functions for "grey" quantities
#====================================================================

class GreyQuantities:

    """ This class stores grey photoionization cross-sections and 
    energies for each of the absorbing species (H1, He1, He2). Note that 
    this requires a segmented spectrum. """ 


    def __init__( self, rad_src ):

        
        # make sure we have what we need
        #--------------------------------------------
        assert rad_src.segmented == True

        # setup containers
        #--------------------------------------------
        self.sigma = species.AbsorbingSpecies()
        self.E = species.AbsorbingSpecies()


        # integrate between H1 and He1 thresholds
        #--------------------------------------------
        ii = rad_src.th.i_H1
        ff = rad_src.th.i_He1
        xx = rad_src.E[ii:ff]
        yy = rad_src._dH1i_over_dE[ii:ff]
        H1i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE[ii:ff]
        Fn_thin = utils.trap( xx, yy )
        self.sigma.H1 = H1i_thin / Fn_thin
        self.sigma.H1.__doc__ = 'grey H1 photoionization cross section' 

        log_sig = np.log10( rad_src.sigma.H1[ii:ff].magnitude )
        log_gry = np.log10( self.sigma.H1.magnitude )
        log_E = np.log10( rad_src.E[ii:ff].magnitude )

        i_lo = np.argmin( np.abs( log_sig - log_gry ) )
        if log_sig[i_lo] < log_gry:
            i_hi = i_lo + 1
        else:
            i_hi = i_lo
            i_lo = i_lo - 1

        dlogsig = log_sig[i_hi] - log_sig[i_lo]
        frac = (log_gry - log_sig[i_lo]) / dlogsig
        dlogE = log_E[i_hi] - log_E[i_lo]
        grey_E = 10.0**(log_E[i_lo] + frac * dlogE)
        self.E.H1 = grey_E * rad_src.U.eV
        self.E.H1.__doc__ = 'grey H1 energy' 

        # integrate between He1 and He2 thresholds
        #--------------------------------------------
        ii = rad_src.th.i_He1
        ff = rad_src.th.i_He2
        xx = rad_src.E[ii:ff]
        yy = rad_src._dHe1i_over_dE[ii:ff]
        He1i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE[ii:ff]
        Fn_thin = utils.trap( xx, yy )
        self.sigma.He1 = He1i_thin / Fn_thin
        self.sigma.He1.__doc__ = 'grey He1 photoionization cross section' 
        
        log_sig = np.log10( rad_src.sigma.He1[ii:ff].magnitude )
        log_gry = np.log10( self.sigma.He1.magnitude )
        log_E = np.log10( rad_src.E[ii:ff].magnitude )

        i_lo = np.argmin( np.abs( log_sig - log_gry ) )
        if log_sig[i_lo] < log_gry:
            i_hi = i_lo + 1
        else:
            i_hi = i_lo
            i_lo = i_lo - 1
            
        dlogsig = log_sig[i_hi] - log_sig[i_lo]
        frac = (log_gry - log_sig[i_lo]) / dlogsig
        dlogE = log_E[i_hi] - log_E[i_lo]
        grey_E = 10.0**(log_E[i_lo] + frac * dlogE)
        self.E.He1 = grey_E * rad_src.U.eV
        self.E.He1.__doc__ = 'grey He1 energy' 
        
        # integrate above He2 threshold
        #--------------------------------------------
        ii = rad_src.th.i_He2
        xx = rad_src.E[ii:]
        yy = rad_src._dHe2i_over_dE[ii:]
        He2i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE[ii:]
        Fn_thin = utils.trap( xx, yy )
        self.sigma.He2 = He2i_thin / Fn_thin
        self.sigma.He2.__doc__ = 'grey He2 photoionization cross section' 
        
        log_sig = np.log10( rad_src.sigma.He2[ii:].magnitude )
        log_gry = np.log10( self.sigma.He2.magnitude )
        log_E = np.log10( rad_src.E[ii:].magnitude )
        
        i_lo = np.argmin( np.abs( log_sig - log_gry ) )
        if log_sig[i_lo] < log_gry:
            i_hi = i_lo + 1
        else:
            i_hi = i_lo
            i_lo = i_lo - 1
            
        dlogsig = log_sig[i_hi] - log_sig[i_lo]
        frac = (log_gry - log_sig[i_lo]) / dlogsig
        dlogE = log_E[i_hi] - log_E[i_lo]
        grey_E = 10.0**(log_E[i_lo] + frac * dlogE)
        self.E.He2 = grey_E * rad_src.U.eV
        self.E.He2.__doc__ = 'grey He2 energy' 

            
class OpticallyThinQuantities:

    """ This class stores optically thin  characterizations of the radiation 
    field.  """ 
    
    def __init__ ( self, rad_src ):

        self.PC = physical.PhysicalConstants()
        self.U = units.Units()

        # dont need to integrate for monochromatic spectra
        #-----------------------------------------------------------------
        if rad_src.monochromatic:

            four_pi_sr = 4 * np.pi * self.U.sr

            self.Fu = four_pi_sr * rad_src.Inu
            self.u = self.Fu / self.PC.c
            self.Fn = self.Fu / rad_src.E
            self.n = self.Fn / self.PC.c

            self.H1i = self.Fn * rad_src.sigma.H1
            self.H1h = self.H1i * (rad_src.E - rad_src.th.E_H1)
            self.He1i = self.Fn * rad_src.sigma.He1
            self.He1h = self.He1i * (rad_src.E - rad_src.th.E_He1)
            self.He2i = self.Fn * rad_src.sigma.He2
            self.He2h = self.He2i * (rad_src.E - rad_src.th.E_He2)


        # integrate for polychromatic spectra
        #-----------------------------------------------------------------
        else:

            self.Fu = utils.trap( rad_src.E, rad_src._dFu_over_dE )
            self.u = utils.trap( rad_src.E, rad_src._du_over_dE )
            self.Fn = utils.trap( rad_src.E, rad_src._dFn_over_dE )
            self.n = utils.trap( rad_src.E, rad_src._dn_over_dE )

            self.H1i = utils.trap( rad_src.E, rad_src._dH1i_over_dE )
            self.H1h = utils.trap( rad_src.E, rad_src._dH1h_over_dE ) 
            self.He1i = utils.trap(rad_src.E, rad_src._dHe1i_over_dE)
            self.He1h = utils.trap(rad_src.E, rad_src._dHe1h_over_dE) 
            self.He2i = utils.trap(rad_src.E, rad_src._dHe2i_over_dE)
            self.He2h = utils.trap(rad_src.E, rad_src._dHe2h_over_dE)


        # set units and docstrings
        #-----------------------------------------------------------------
        self.Fu.units = 'erg/s/cm^2'        
        self.Fu.__doc__ = 'optically thin energy flux' 

        self.u.units = 'erg/cm^3'
        self.u.__doc__ = 'optically thin energy density' 

        self.Fn.units = '1/s/cm^2'
        self.Fn.__doc__ = 'optically thin photon number flux' 

        self.n.units = '1/cm^3'
        self.n.__doc__ = 'optically thin photon number density' 

        self.H1i.units = '1/s'            
        self.H1i.__doc__ = 'optically thin H1 photoionization rate' 

        self.H1h.units = 'erg/s'
        self.H1h.__doc__ = 'optically thin H1 photoheating rate' 

        self.He1i.units = '1/s'            
        self.He1i.__doc__ = 'optically thin He1 photoionization rate' 

        self.He1h.units = 'erg/s'
        self.He1h.__doc__ = 'optically thin He1 photoheating rate' 

        self.He2i.units = '1/s'            
        self.He2i.__doc__ = 'optically thin He2 photoionization rate' 

        self.He2h.units = 'erg/s'
        self.He2h.__doc__ = 'optically thin He2 photoheating rate' 




        
