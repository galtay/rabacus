""" An plane parallel radiation source class. """ 


import numpy as np
from rabacus.utils import utils
from rabacus.constants import physical

from hm12 import HM12Spectrum
from thermal import ThermalSpectrum
from powerlaw import PowerlawSpectrum
from monochromatic import MonochromaticSpectrum
from source import Source

import species


__all__ = ['PlaneSource']



class PlaneSource(Source):


    r""" A plane parallel radiation source class.  

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

      `dFu_over_dnu` (array): Energy flux density |dFu_dnu| 


      `monochromatic` (bool): ``True`` for monochromatic spectra otherwise 
      ``False``

      `grey` (object):  Grey photoionization cross sections and energies.
      Only created if `segmented` = ``True``.

      `log` (object): Convenience object storing the log of other attributes

      `sigma` (object): photoionization cross sections and ratios

      `th` (object): all quantities evaluated at ionization thresholds

      `thin` (object): all optically thin quantities

    Notes:


    The fundamental characterization of the radiation field for this class is 
    the energy flux density |dFu_dnu| with units of ``erg/(s cm^2 Hz)`` or 
    energy flux :math:`F_{u}` with units of ``erg/(s cm^2)`` for monochromatic 
    spectra. 
    Integrating over frequency yields several other useful characterizations
    of the radiation field which we summarize below. 
    
    .. math:: F_{u} = \int \frac{dF_{u}}{d\nu} d\nu 

    .. math:: u = \frac{1}{c} \int \frac{dF_{u}}{d\nu} d\nu 
              \rightarrow
              \frac{du}{d\nu} = \frac{1}{c} \frac{dF_{u}}{d\nu}

    .. math:: F_{n} = \int \frac{1}{h \nu} \frac{dF_{u}}{d\nu} d\nu 
              \rightarrow
              \frac{dF_{n}}{d\nu} = \frac{1}{h \nu} \frac{dF_{u}}{d\nu}

    .. math:: n = \frac{1}{c} \int \frac{1}{h \nu} \frac{dF_{u}}{d\nu} d\nu 
              \rightarrow 
              \frac{dn}{d\nu} = \frac{1}{c} \frac{1}{h \nu} \frac{dF_{u}}{d\nu}


    +-----------+---------------------+----------------------+-----------+
    | Variable  | Description         | Units                | Symbol    |
    +===========+=====================+======================+===========+
    | dFu_dnu   | energy flux density | [erg/(s cm^2 Hz)]    | |dFu_dnu| |
    +-----------+---------------------+----------------------+-----------+
    | thin.Fu   | energy flux         | [erg/(s cm^2)]       | |Fu|      |
    +-----------+---------------------+----------------------+-----------+
    | thin.u    | energy density      | [erg/(cm^3)]         | |u|       |
    +-----------+---------------------+----------------------+-----------+
    | thin.Fn   | photon flux         | [1/(s cm^2)]         | |Fn|      |
    +-----------+---------------------+----------------------+-----------+
    | thin.n    | photon density      | [1/(cm^3)]           | |n|       |
    +-----------+---------------------+----------------------+-----------+

    The photoionization and photoheating rates are also integrals over the 
    energy flux density. 

    .. math:: \Gamma_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dF_{u}}{d\nu} d\nu 
              \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{d\nu} =  
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dF_{u}}{d\nu}

    .. math:: H_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dF_{u}}{d\nu} 
              h(\nu - \nu_{\rm _{X}}) d\nu \rightarrow
              \frac{ dH_{\rm _{X}} }{d\nu} =  
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dF_{u}}{d\nu} 
              h(\nu - \nu_{\rm _{X}})

    where :math:`\scriptstyle{\rm X}` is one of 
    :math:`\{ {\scriptstyle{\rm HI, \, HeI, \, HeII }} \}`.  
    Note that the :math:`\sigma_{\rm _{X}}` are frequency dependent 
    photoionization cross-sections and the :math:`\nu_{\rm _{X}}` are  
    (scalar) frequencies that correspond to the ionization energies for a given 
    species. In many applications, we prefer to use energy as a variable 
    instead of frequency.  The above equations can easily be recast with a 
    change of variable :math:`E=h\nu`. 

    .. math:: \Gamma_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE} dE 
              \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{dE} =  
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}


    .. math:: H_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}  
              (E - E_{\rm _{X}}) dE \rightarrow
              \frac{ dH_{\rm _{X}} }{dE} =  
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}  
              (E - E_{\rm _{X}})

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

    .. math:: \Gamma_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}  
              e^{-\tau} dE \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{dE} =
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}e^{-\tau}


    .. math:: H_{\rm _{X}} = \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}  
              (E - E_{\rm _{X}}) e^{-\tau} dE \rightarrow
              \frac{ dH_{\rm _{X}} }{dE} =  
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dF_{u}}{dE}  
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

      :class:`~rabacus.rad_src.background.BackgroundSource`, 
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
        self.source_type = 'plane' 


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
            den = self.U.s * self.U.cm**2 
            units = num / den
            self.Fu = np.ones( self.Nnu ) * spec.yvals * units
            self.Fu.__doc__ = 'energy flux'
            self.log.Fu = np.log10( self.Fu.magnitude )
        else:
            den = self.U.Hz  * self.U.s * self.U.cm**2 
            units = num / den
            if spectrum_type == 'user':
                self.dFu_over_dnu = np.ones( self.Nnu ) * user_shape * units
            else:
                self.dFu_over_dnu = np.ones( self.Nnu ) * spec.yvals * units
                if spectrum_type == 'hm12':
                    self.dFu_over_dnu *= 4 * np.pi
            self.dFu_over_dnu.__doc__ = 'energy flux per unit frequency' 
            self.log.dFu_over_dnu = np.log10( self.dFu_over_dnu.magnitude )


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
        energy flux density or dFu_over_dnu which has units [erg/(s cm^2 Hz)].

        """

        # dont need to set these for monochromatic spectra
        #--------------------------------------------------------
        if self.monochromatic:
            return


        # first we define differentials wrt frequency.  these are 
        # just for convenience. 
        #--------------------------------------------------------

        self._du_over_dnu = self.dFu_over_dnu / self.PC.c
        self._du_over_dnu.units = 'erg/cm^3/Hz'

        self._dFn_over_dnu = self.dFu_over_dnu / ( self.PC.h * self.nu )
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

        self._dFu_over_dE = self.dFu_over_dnu / self.PC.h
        self._dFu_over_dE.units = 'erg/s/cm^2/eV'

        self._du_over_dE = self._dFu_over_dE / self.PC.c
        self._du_over_dE.units = 'erg/cm^3/eV'

        self._dFn_over_dE = self.dFu_over_dnu / ( self.PC.h * self.E )
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

        if self.monochromatic:
            self.Fu *= fac
        else:
            self.dFu_over_dnu *= fac

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
        

    def normalize_Fn( self, Fn ):
        """ Normalize the spectrum such that the photon flux between `q_min` 
        and `q_max` is `Fn`. 

        Args:
          `Fn` (float): target photon flux
        """ 

        if not hasattr(Fn,'units'):
            msg = '\n Fn must have units, e.g. s^-1 cm^-2 \n'
            raise utils.NeedUnitsError, msg
        else:
            Fn.units = 's^-1 * cm^-2'

        fac = Fn / self.thin.Fn
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

        # dont need to integrate for monochromatic spectra
        #-----------------------------------------------------------------
        if rad_src.monochromatic:

            self.Fu = rad_src.Fu
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



            

        
