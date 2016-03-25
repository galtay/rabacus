""" A point radiation source class. """ 

import numpy as np
from rabacus.utils import utils
from rabacus.constants import physical

from hm12 import HM12Spectrum
from thermal import ThermalSpectrum
from powerlaw import PowerlawSpectrum
from monochromatic import MonochromaticSpectrum
from source import Source

import species


__all__ = ['PointSource']



class PointSource(Source):

    r""" A point radiation source class. 
    
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

      `dLu_over_dnu` (array): Energy luminosity density |dLu_dnu| 


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
    is the energy luminosity density |dLu_dnu| with units of ``erg/(s Hz)`` 
    or energy luminosity :math:`L_{u}` with units of ``erg/s`` for 
    monochromatic spectra.  All quantities except for the luminosity are 
    dependent on radius.  We summarize them below, 
    
    .. math:: F_{u}(r) = \frac{1}{4 \pi r^2} \int \frac{dL_{u}}{d\nu} d\nu 

    .. math:: u(r) = \frac{1}{c} \frac{1}{4 \pi r^2} 
              \int \frac{dL_{u}}{d\nu} d\nu 

    .. math:: F_{n}(r) = \frac{1}{4 \pi r^2} \int \frac{1}{h \nu}  
              \frac{dL_{u}}{d\nu} d\nu 

    .. math:: n(r) = \frac{1}{c} \frac{1}{4 \pi r^2} 
              \int \frac{1}{h \nu} \frac{dL_{u}}{d\nu} d\nu 


    +-----------+-----------------------+----------------------+-----------+
    | Variable  | Description           | Units                | Symbol    |
    +===========+=======================+======================+===========+
    | dLu_dnu   | energy lum. density   | [erg/(s Hz)]         | |dLu_dnu| |
    +-----------+-----------------------+----------------------+-----------+
    | thin.Fu   | energy flux           | [erg/(s cm^2)]       | |Fu|      |
    +-----------+-----------------------+----------------------+-----------+
    | thin.u    | energy density        | [erg/(cm^3)]         | |u|       |
    +-----------+-----------------------+----------------------+-----------+
    | thin.Fn   | photon flux           | [1/(s cm^2)]         | |Fn|      |
    +-----------+-----------------------+----------------------+-----------+
    | thin.n    | photon density        | [1/(cm^3)]           | |n|       |
    +-----------+-----------------------+----------------------+-----------+


    The photoionization and photoheating rates are also integrals over the 
    energy luminosity density. 

    .. math:: \Gamma_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int  
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dL_{u}}{d\nu} d\nu 
              \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{d\nu} =  
              \frac{1}{4 \pi r^2} 
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dL_{u}}{d\nu}

    .. math:: H_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int 
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dL_{u}}{d\nu} 
              h(\nu - \nu_{\rm _{X}}) d\nu \rightarrow
              \frac{ dH_{\rm _{X}} }{d\nu} = 
              \frac{1}{4 \pi r^2} 
              \frac{ \sigma_{\rm _{X}} }{h \nu} \frac{dL_{u}}{d\nu} 
              h(\nu - \nu_{\rm _{X}})

    where :math:`\scriptstyle{\rm X}` is one of 
    :math:`\{ {\scriptstyle{\rm HI, \, HeI, \, HeII }} \}`.  
    Note that the :math:`\sigma_{\rm _{X}}` are frequency dependent 
    photoionization cross-sections and the :math:`\nu_{\rm _{X}}` are  
    (scalar) frequencies that correspond to the ionization energies for a given 
    species. In many applications, we prefer to use energy as a variable 
    instead of frequency.  The above equations can easily be recast with a 
    change of variable :math:`E=h\nu`. 

    .. math:: \Gamma_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE} dE 
              \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{dE} =  
              \frac{1}{4 \pi r^2} 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE}


    .. math:: H_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE}  
              (E - E_{\rm _{X}}) dE \rightarrow
              \frac{ dH_{\rm _{X}} }{dE} =  
              \frac{1}{4 \pi r^2}
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE}  
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

    .. math:: \Gamma_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE} e^{-\tau} dE 
              \rightarrow
              \frac{ d\Gamma_{\rm _{X}} }{dE} =  
              \frac{1}{4 \pi r^2} 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE} e^{-\tau}


    .. math:: H_{\rm _{X}}(r) = \frac{1}{4 \pi r^2} \int 
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE}  
              (E - E_{\rm _{X}}) e^{-\tau} dE \rightarrow
              \frac{ dH_{\rm _{X}} }{dE} =  
              \frac{1}{4 \pi r^2}
              \frac{ \sigma_{\rm _{X}} }{E} \frac{dL_{u}}{dE}  
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
      :class:`~rabacus.rad_src.plane.PlaneSource`

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
        self.source_type = 'point' 


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
        # the primary variable for point source is dL/dnu [erg/s/Hz]
        #-----------------------------------------------------------------
        num = self.U.erg

        if self.monochromatic:
            den = self.U.s
            units = num / den
            self.Lu = np.ones( self.Nnu ) * spec.yvals * units
            self.Lu.__doc__ = 'energy luminosity'
        else:
            den = self.U.Hz * self.U.s 
            units = num / den
            if spectrum_type == 'user':
                self.dLu_over_dnu = np.ones( self.Nnu ) * user_shape * units
            else:
                self.dLu_over_dnu = np.ones( self.Nnu ) * spec.yvals * units
            self.dLu_over_dnu.__doc__ = 'energy luminosity per unit frequency' 
            self.log.dLu_over_dnu = np.log10( self.dLu_over_dnu.magnitude )


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
        energy luminosity density or dLu_over_dnu which has units [erg/(s Hz)].

        """

        # dont need to set these for monochromatic spectra
        #--------------------------------------------------------
        if self.monochromatic:
            return

        # first we define differentials wrt frequency.  these are 
        # just for convenience. 
        #--------------------------------------------------------
        self._dLn_over_dnu = self.dLu_over_dnu / ( self.PC.h * self.nu )
        self._dLn_over_dnu.units = '1/s/Hz'

        # now we define differentials wrt energy 
        #--------------------------------------------------------
        self._dLu_over_dE = self.dLu_over_dnu / self.PC.h
        self._dLu_over_dE.units = 'erg/s/eV'

        self._dLn_over_dE = self._dLu_over_dE / ( self.PC.h * self.nu )
        self._dLn_over_dE.units = '1/s/eV'


# integrands dependent on radius need functions 
#====================================================================

    def _dFu_over_dnu( self, r ):
        """ Energy flux per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self.dLu_over_dnu / (4.0 * np.pi * r*r)
        dxdx.units = 'erg/s/Hz/cm^2'
        return dxdx

    def _du_over_dnu( self, r ):
        """ Energy density per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFu_over_dnu(r) / self.PC.c
        dxdx.units = 'erg/cm^3/Hz'
        return dxdx

    #-----------------------------------------------------------------
    def _dFn_over_dnu( self, r ):
        """ Photon flux per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dLn_over_dnu / (4.0 * np.pi * r*r)
        dxdx.units = '1/s/Hz/cm^2'
        return dxdx

    def _dn_over_dnu( self, r ):
        """ Photon density per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dnu(r) / self.PC.c
        dxdx.units = '1/cm^3/Hz'
        return dxdx

    #-----------------------------------------------------------------
    def _dFu_over_dE( self, r ):
        """ Energy flux per photon energy at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFu_over_dnu(r) / self.PC.h
        dxdx.units = 'erg/s/cm^2/eV'
        return dxdx

    def _du_over_dE( self, r ):
        """ Energy density per photon energy at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFu_over_dE(r) / self.PC.c
        dxdx.units = 'erg/cm^3/eV'
        return dxdx

    #-----------------------------------------------------------------
    def _dFn_over_dE( self, r ):
        """ Photon flux per photon energy at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFu_over_dnu(r) / ( self.PC.h * self.E )
        dxdx.units = '1/s/cm^2/eV'
        return dxdx

    def _dn_over_dE( self, r ):
        """ Photon density per photon energy at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dE(r) / self.PC.c
        dxdx.units = '1/cm^3/eV'
        return dxdx

    #-----------------------------------------------------------------
    def _dH1i_over_dnu( self, r ):
        """ H1 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dnu(r) * self.sigma.H1
        dxdx.units = '1/s/Hz'
        return dxdx

    def _dH1h_over_dnu( self, r ):
        """ H1 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dH1i_over_dnu(r) * ( self.E - self.th.E_H1 )
        dxdx.units = 'erg/s/Hz'
        return dxdx

    #-----------------------------------------------------------------
    def _dHe1i_over_dnu( self, r ):
        """ He1 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dnu(r) * self.sigma.He1
        dxdx.units = '1/s/Hz'
        return dxdx

    def _dHe1h_over_dnu( self, r ):
        """ He1 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dHe1i_over_dnu(r) * ( self.E - self.th.E_He1 )
        dxdx.units = 'erg/s/Hz'
        return dxdx

    #-----------------------------------------------------------------
    def _dHe2i_over_dnu( self, r ):
        """ He2 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dnu(r) * self.sigma.He2
        dxdx.units = '1/s/Hz'
        return dxdx

    def _dHe2h_over_dnu( self, r ):
        """ He2 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dHe2i_over_dnu(r) * ( self.E - self.th.E_He2 )
        dxdx.units = 'erg/s/Hz'
        return dxdx

    #-----------------------------------------------------------------
    def _dH1i_over_dE( self, r ):
        """ H1 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dE(r) * self.sigma.H1
        dxdx.units = '1/s/eV'
        return dxdx

    def _dH1h_over_dE( self, r ):
        """ H1 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dH1i_over_dE(r) * ( self.E - self.th.E_H1 )
        dxdx.units = 'erg/s/eV'
        return dxdx

    #-----------------------------------------------------------------
    def _dHe1i_over_dE( self, r ):
        """ He1 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dE(r) * self.sigma.He1
        dxdx.units = '1/s/eV'
        return dxdx

    def _dHe1h_over_dE( self, r ):
        """ He1 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dHe1i_over_dE(r) * ( self.E - self.th.E_He1 )
        dxdx.units = 'erg/s/eV'
        return dxdx

    #-----------------------------------------------------------------
    def _dHe2i_over_dE( self, r ):
        """ He2 photoion rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dFn_over_dE(r) * self.sigma.He2
        dxdx.units = '1/s/eV'
        return dxdx

    def _dHe2h_over_dE( self, r ):
        """ He2 photoheat rate per frequency at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        dxdx = self._dHe2i_over_dE(r) * ( self.E - self.th.E_He2 )
        dxdx.units = 'erg/s/eV'
        return dxdx





# functions that calculate shielded integrals dependent on radius
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


    def shld_H1i(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HI photoionization rate
        after passing through a column with optical depth tau 

        Args:
          `r` (float): radial distance from point source

          `tauH1_th` (float): H1 optical depth at the H1 ionizing threshold
        
          `tauHe1_th` (float): He1 optical depth at the He1 ionizing threshold

          `tauHe2_th` (float): He2 optical depth at the He2 ionizing threshold

        Returns:
          `H1i` (float): attenuated H1 photoionization rate 
   
        """ 
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            H1i = self.thin.H1i(r) * atten
        else:
            H1i = utils.trap( self.E, self._dH1i_over_dE(r) * atten )
        H1i.units = '1/s'
        return H1i


    def shld_H1h(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HI photoheating rate
        after passing through a column with optical depth tau 

        Args:
          `r` (float): radial distance from point source

          `tauH1_th` (float): H1 optical depth at the H1 ionizing threshold
        
          `tauHe1_th` (float): He1 optical depth at the He1 ionizing threshold

          `tauHe2_th` (float): He2 optical depth at the He2 ionizing threshold

        Returns:
          `H1h` (float): attenuated H1 photoheating rate 
           
        """ 
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            H1h = self.thin.H1h(r) * atten
        else:
            H1h = utils.trap( self.E, self._dH1h_over_dE(r) * atten )
        H1h.units = 'erg/s'
        return H1h


    def shld_He1i(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeI photoionization rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1i`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He1i = self.thin.He1i(r) * atten
        else:
            He1i = utils.trap( self.E, self._dHe1i_over_dE(r) * atten )
        He1i.units = '1/s'
        return He1i

    def shld_He1h(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeI photoheating rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1h`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He1h = self.thin.He1h(r) * atten
        else:
            He1h = utils.trap( self.E, self._dHe1h_over_dE(r) * atten )
        He1h.units = 'erg/s'
        return He1h


    def shld_He2i(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeII photoionization rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1i`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He2i = self.thin.He2i(r) * atten
        else:
            He2i = utils.trap( self.E, self._dHe2i_over_dE(r) * atten )
        He2i.units = '1/s'
        return He2i

    def shld_He2h(self, r, tauH1_th, tauHe1_th, tauHe2_th):
        """ Integrates spectrum to calculate the HeII photoheating rate
        after passing through a column with optical depth tau. 

        .. seealso::

          :func:`shld_H1h`        
        """
        atten = self.return_attenuation( tauH1_th, tauHe1_th, tauHe2_th ) 
        if self.monochromatic:
            He2h = self.thin.He2h(r) * atten
        else:
            He2h = utils.trap( self.E, self._dHe2h_over_dE(r) * atten )
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
            self.Lu *= fac
        else:
            self.dLu_over_dnu *= fac

        self.set_integrands()
        self.thin = OpticallyThinQuantities( self )


    def normalize_Ln( self, Ln ):
        """ Normalize spectrum such that the photon luminosity is Ln. 

        Args: 
          `Ln` (float): target photon luminosity 
        """ 
        if not hasattr(Ln,'units'): 
            raise utils.NeedUnitsError, '\n n must have units of 1/time \n'
        else:
            Ln.units = '1/s'

        fac = Ln / self.thin.Ln
        self.scale_spectrum( fac )


    def normalize_Lu( self, Lu ):
        """ Normalize spectrum such that the energy luminosity is Lu. 

        Args: 
          `Lu` (float): target energy luminosity 
        """         
        if not hasattr(Lu,'units'): 
            raise utils.NeedUnitsError, '\n n must have units of energy/time \n'
        else:
            Lu.units = 'erg/s'

        fac = Lu / self.thin.Lu
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

        r = 1.0 * rad_src.U.kpc

        # integrate between H1 and He1 thresholds
        #--------------------------------------------
        ii = rad_src.th.i_H1
        ff = rad_src.th.i_He1
        xx = rad_src.E[ii:ff]
        yy = rad_src._dH1i_over_dE(r)[ii:ff]
        H1i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE(r)[ii:ff]
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
        yy = rad_src._dHe1i_over_dE(r)[ii:ff]
        He1i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE(r)[ii:ff]
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
        yy = rad_src._dHe2i_over_dE(r)[ii:]
        He2i_thin = utils.trap( xx, yy )
        yy = rad_src._dFn_over_dE(r)[ii:]
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
        self.rad_src = rad_src

        # dont need to integrate for monochromatic spectra
        #-----------------------------------------------------------------
        if rad_src.monochromatic:

            self.Lu = rad_src.Lu
            self.Ln = self.Lu / rad_src.E

        # integrate for polychromatic spectra
        #-----------------------------------------------------------------
        else:

            self.Lu = utils.trap( rad_src.E, rad_src._dLu_over_dE )
            self.Ln = utils.trap( rad_src.E, rad_src._dLn_over_dE )


        # set units and docstrings
        #-----------------------------------------------------------------
        self.Lu.units = 'erg/s'
        self.Lu.__doc__ = 'optically thin energy luminosity'

        self.Ln.units = '1/s'
        self.Ln.__doc__ = 'optically thin photon luminosity'



# functions that calculate integrals dependent on radius
#====================================================================

    def Fu(self,r):
        """ energy flux at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'        
        if self.rad_src.monochromatic:
            Fu = self.Lu / (4.0 * np.pi * r*r)
        else:
            Fu = utils.trap( self.rad_src.E, self.rad_src._dFu_over_dE(r) )
        return Fu

    def u(self,r):
        """ energy density at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'        
        if self.rad_src.monochromatic:
            u = self.Fu(r) / self.PC.c
        else:
            u = utils.trap( self.rad_src.E, self.rad_src._du_over_dE(r) )
        return u

    #-----------------------------------------------------------------
    def Fn(self,r):
        """ photon flux at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'      
        if self.rad_src.monochromatic:
            Fn = self.Ln / (4.0 * np.pi * r*r)
        else:
            Fn = utils.trap( self.rad_src.E, self.rad_src._dFn_over_dE(r) )
        return Fn

    def n(self,r):
        """ photon density at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            n = self.Fn(r) / self.PC.c
        else:
            n = utils.trap( self.rad_src.E, self.rad_src._dn_over_dE(r) )
        return n

    #-----------------------------------------------------------------
    def H1i(self,r):
        """ H1 photoion rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'        
        if self.rad_src.monochromatic:
            H1i = self.Fn(r) * self.rad_src.sigma.H1
        else:
            H1i = utils.trap( self.rad_src.E, self.rad_src._dH1i_over_dE(r) )
        H1i.units = '1/s'
        return H1i

    def H1h(self,r):
        """ H1 photoheating rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            H1h = self.H1i(r) * (self.rad_src.E - self.rad_src.th.E_H1) 
        else:
            H1h = utils.trap( self.rad_src.E, self.rad_src._dH1h_over_dE(r) )
        H1h.units = 'erg/s'
        return H1h

    #-----------------------------------------------------------------
    def He1i(self,r):
        """ He1 photoion rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            He1i = self.Fn(r) * self.rad_src.sigma.He1
        else:
            He1i = utils.trap( self.rad_src.E, self.rad_src._dHe1i_over_dE(r) )
        He1i.units = '1/s'
        return He1i

    def He1h(self,r):
        """ He1 photoheating rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            He1h = self.He1i(r) * (self.rad_src.E - self.rad_src.th.E_He1)
        else:
            He1h = utils.trap( self.rad_src.E, self.rad_src._dHe1h_over_dE(r) )
        He1h.units = 'erg/s'
        return He1h

    #-----------------------------------------------------------------
    def He2i(self,r):
        """ He2 photoion rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            He2i = self.Fn(r) * self.rad_src.sigma.He2
        else:
            He2i = utils.trap( self.rad_src.E, self.rad_src._dHe2i_over_dE(r) )
        He2i.units = '1/s'
        return He2i

    def He2h(self,r):
        """ He2 photoheating rate at a given radius. """ 
        if not hasattr(r,'units'): 
            raise utils.NeedUnitsError, '\n r must have units \n'
        if self.rad_src.monochromatic:
            He2h = self.He2i(r) * (self.rad_src.E - self.rad_src.th.E_He2)
        else:
            He2h = utils.trap( self.rad_src.E, self.rad_src._dHe2h_over_dE(r) )
        He2h.units = 'erg/s'
        return He2h

    #-----------------------------------------------------------------





