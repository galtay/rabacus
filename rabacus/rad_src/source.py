""" Module for source base class. """ 

import sys
import numpy as np

from rabacus.utils import utils

from rabacus.atomic import hydrogen
from rabacus.atomic import helium
from rabacus.atomic import photo_xsection

from rabacus.constants import physical
from rabacus.constants import units 



__all__ = ['Source']



ValidSpectrumTypes = ['monochromatic', 'thermal', 'powerlaw', 'hm12',
                      'user']



class Source:

    r""" Base class for radiation sources.  Derive specific radiation source 
    classes from this class.  Cannot be directly instantiated.  Note that 
    all derived source classes call :func:`initialize` when they are 
    instanciated.

    .. seealso::


      :class:`~rabacus.rad_src.point.PointSource`,
      :class:`~rabacus.rad_src.plane.PlaneSource`,
      :class:`~rabacus.rad_src.background.BackgroundSource`


    """ 

    def __init__( self ): 
        raise utils.BaseClassError

    def initialize( self, q_min, q_max, spectrum_type, Nnu, segmented,
                    px_fit_type, verbose, z, T_eff, alpha, user_E, 
                    user_shape ):

        """ Perform general initialization.

        Args: 

          `q_min` (float): minimum photon energy / Rydbergs [dimensionless]

          `q_max` (float): maximum photon energy / Rydbergs [dimensionless]

          `spectrum_type` (str): the spectral shape 
          {``monochromatic``, ``hm12``, ``thermal``, ``powerlaw``, ``user``}   

          `Nnu` (int): number of spectral samples, log spaced in energy

          `segmented` (bool): if ``True``, forces spectral samples at H and He 
          ionization thresholds

          `px_fit_type` (str): source to use for photoionization cross section 
          fits {``verner``}
       
          `verbose` (bool): verbose screen output?

          `z` (float): redshift, need if `spectrum_type` = ``hm12``

          `T_eff` (float): effective temperature, need if `spectrum_type` = 
          ``thermal``

          `alpha` (float): powerlaw index, need if `spectrum_type` = 
          ``powerlaw``

          `user_E` (array): energy samples for user defined spectrum. should
          have units of energy. 

          `user_shape` (array): shape of user defined spectrum. should be 
          dimensionless

        """ 

        # attach standard input
        #-------------------------------------------------------- 
        self.q_min = q_min
        self.q_max = q_max
        self.spectrum_type = spectrum_type

        self.Nnu = Nnu
        self.segmented = segmented
        self.px_fit_type = px_fit_type
        self.verbose = verbose


        # attach spectrum specific input
        #-------------------------------------------
        if spectrum_type == 'hm12':
            self.z = z
        elif spectrum_type == 'thermal':
            self.T_eff = T_eff
        elif spectrum_type == 'powerlaw':
            self.alpha = alpha
        elif spectrum_type == 'user':
            self.user_E = user_E
            self.user_shape = user_shape
            self.segmented = False

        if spectrum_type == 'monochromatic':
            self.monochromatic = True
            self.Nnu = 1
            self.segmented = False
        else:
            self.monochromatic = False


        self.validate_spectrum_type()

        # attach physical constants and units 
        #--------------------------------------------------------
        self.U = units.Units()
        self.PC = physical.PhysicalConstants()
        self.PX = photo_xsection.PhotoXsections( fit_type = px_fit_type )

        # attach hydrogen and helium atom classes 
        #--------------------------------------------------------
        self.H = hydrogen.Hydrogen( px_fit_type = self.px_fit_type )
        self.He = helium.Helium( px_fit_type = self.px_fit_type )

        # set quantities related to photo-ionization thresholds
        #--------------------------------------------------------
        self.th = PhotoIonizationThresholds( self )


        # special behaviour for user defined spectra
        #-------------------------------------------
        if spectrum_type == 'user':

            self.E = user_E
            self.E.units = 'eV'
            self.E.__doc__ = 'energy samples' 

            self.q_min = self.E.min() / self.PX.Eth_H1
            self.q_max = self.E.max() / self.PX.Eth_H1
            self.Nnu = self.E.size

            self.q = self.E / self.PX.Eth_H1
            self.q.__doc__ = 'energy samples / Rydbergs'
        
            self.nu = self.E / self.PC.h        
            self.nu.units = 'Hz'
            self.nu.__doc__ = 'frequency samples' 
            
            self.lam = self.PC.c / self.nu
            self.lam.units = 'cm'
            self.lam.__doc__ = 'wavelength samples' 
            

        # set photon energy arrays
        #--------------------------------------------------------        
        else:

            self.set_photon_arrays( self.q_min, self.q_max, self.Nnu, 
                                    self.segmented )

        # set X-sections
        #--------------------------------------------------------
        self.sigma = PhotoIonizationCrossSections( self )

        # store log10 quantities 
        #--------------------------------------------------------
        self.log = LogQuantities( self )


    def set_photon_arrays(self, q_min, q_max, Nnu, segmented):
        """ Creates wavelength, frequency, and energy arrays for the spectrum.

        Args: 

          `q_min` (float): minimum photon energy / Rydbergs [dimensionless]

          `q_max` (float): maximum photon energy / Rydbergs [dimensionless]

          `Nnu` (int): number of spectral samples, log spaced in energy

          `segmented` (bool): if ``True``, forces spectral samples at H and He 
          ionization thresholds

        """  

        # we always start with a uniform spacing in log q
        #-----------------------------------------------------------------
        log_q_min = np.log10(q_min)
        log_q_max = np.log10(q_max)
        log_q = np.linspace( log_q_min, log_q_max, Nnu ) 

        if hasattr( log_q, 'magnitude' ):
            q = 10**log_q.magnitude
        else:
            q = 10**log_q

        q = q * self.U.dimensionless
        log_q = log_q * self.U.dimensionless

        # if a segmented spectrum is requested, we find the entries in the 
        # energy array closest to the photo-ionization thresholds and change
        # them to be exactly equal. 
        #-----------------------------------------------------------------
        if segmented:

            # for a segmented spectrum, the requested energy range must 
            # cover the H1, He1, and He2 ionization thresholds 
            #---------------------------------------------------------
            if q.min() > self.th.q_H1:
                raise ValueError('q_min must be <= 1 for segmented spectra')

            if q.max() < self.th.q_He2:
                th_q_He1 = str(self.th.q_He2.magnitude)
                txt = 'q_max must be >= ' + th_q_He1 + ' for segmented spectra'
                raise ValueError(txt)
            
            # find the entry nearest each threshold, change it, and 
            # store the index in the thresholds class 
            #---------------------------------------------------------
            log_q_H1 = np.log10( self.th.q_H1 )
            i_H1 = np.argmin( np.abs( log_q - log_q_H1 ) )
            q[i_H1] = self.th.q_H1 
            self.th.i_H1 = i_H1

            log_q_He1 = np.log10( self.th.q_He1 )
            i_He1 = np.argmin( np.abs( log_q - log_q_He1 ) )
            q[i_He1] = self.th.q_He1 
            self.th.i_He1 = i_He1
            
            log_q_He2 = np.log10( self.th.q_He2 )
            i_He2 = np.argmin( np.abs( log_q - log_q_He2 ) )
            q[i_He2] = self.th.q_He2 
            self.th.i_He2 = i_He2

 
        # attach the photon arrays to the object
        #-----------------------------------------------------------------    
        self.q = q * self.U.dimensionless        
        self.q.__doc__ = 'energy samples / Rydbergs'

        self.E = self.q * self.PX.Eth_H1
        self.E.units = 'eV'
        self.E.__doc__ = 'energy samples' 
        
        self.nu = self.E / self.PC.h        
        self.nu.units = 'Hz'
        self.nu.__doc__ = 'frequency samples' 
        
        self.lam = self.PC.c / self.nu
        self.lam.units = 'cm'
        self.lam.__doc__ = 'wavelength samples' 
        

        

    def validate_spectrum_type( self ):
        """ Performs check to make sure the input is compatible with the 
        source type. """ 


        if not self.spectrum_type in ValidSpectrumTypes:
            msg = '\nspectrum_type = ' + self.spectrum_type + ' \n' 
            msg += 'spectrum_type must be one of: ' 
            msg += str(ValidSpectrumTypes)
            raise utils.InputError( msg )        

        if self.spectrum_type == 'hm12':
            if self.z == None:
                msg = '\nif source type = hm12, must provide z'
                raise utils.InputError( msg )

        if self.spectrum_type == 'thermal':
            if self.T_eff == None:
                msg = '\nif source type = thermal, must provide T_eff'
                raise utils.InputError( msg )

        if self.spectrum_type == 'powerlaw':
            if self.alpha == None:
                msg = '\nif source type = powerlaw, must provide alpha'
                raise utils.InputError( msg )

        if self.spectrum_type == 'monochromatic':
            if self.q_min != self.q_max:
                msg = '\nsource type = monochromatic, but q_min != q_max'
                raise utils.InputError( msg )

        if self.spectrum_type == 'user':
            if self.user_E is None:
                msg = '\nif source type = user, must provide user_E'
                raise utils.InputError( msg )
            if self.user_shape is None:
                msg = '\nif source type = user, must provide user_shape'
                raise utils.InputError( msg )

            if self.user_E.size != self.user_shape.size:
                msg = '\n E and shape muse be same size for source_type = user'
                raise utils.InputError( msg )



class PhotoIonizationThresholds:

    r""" A class that stores quantities related to the photo-ionization 
    thresholds of H1, He1, and He2. """ 

    def __init__( self, rad_src ):

        self.q_H1  = 1.0 * rad_src.U.dimensionless
        self.q_H1.__doc__ = 'H1 photo-ionization threshold energy / Ry' 

        self.q_He1 = rad_src.PX.Eth_He1 / rad_src.PX.Eth_H1
        self.q_He1.__doc__ = 'He1 photo-ionization threshold energy / Ry' 

        self.q_He2 = rad_src.PX.Eth_He2 / rad_src.PX.Eth_H1
        self.q_He2.__doc__ = 'He2 photo-ionization threshold energy / Ry' 

        self.E_H1 = rad_src.PX.Eth_H1
        self.E_H1.__doc__ = 'H1 photo-ionization threshold energy' 

        self.E_He1 = rad_src.PX.Eth_He1
        self.E_He1.__doc__ = 'He1 photo-ionization threshold energy' 

        self.E_He2 = rad_src.PX.Eth_He2
        self.E_He2.__doc__ = 'He2 photo-ionization threshold energy' 

        self.sigma_H1 = rad_src.PX.sigma_H1th
        self.sigma_H1.__doc__ = \
            'H1 photo-ionization cross section at H1 photo-ionization ' + \
            'threshold energy' 

        self.sigma_He1 = rad_src.PX.sigma_He1th
        self.sigma_He1.__doc__ = \
            'He1 photo-ionization cross section at He1 photo-ionization ' + \
            'threshold energy' 

        self.sigma_He2 = rad_src.PX.sigma_He2th
        self.sigma_He2.__doc__ = \
            'He2 photo-ionization cross section at He2 photo-ionization ' + \
            'threshold energy' 


class PhotoIonizationCrossSections:

    r""" A class that stores photoionization cross-sections at the 
    energy samples defined by a source. For each of the absorbing ions 
    {H1, He1, He2}, this class stores the photo-ionization cross-sections
    and the ratio of those cross-sections with the cross-section at the
    photo-ionization energy thresholds. """

    def __init__( self, rad_src ):

        self.H1 = rad_src.PX.sigma_H1( rad_src.E )
        self.H1.__doc__ = \
            'H1 photo-ionization cross section at energy samples'

        self.H1_ratio = self.H1 / rad_src.th.sigma_H1
        self.H1_ratio.__doc__ = \
            'H1 photo-ionization cross section at energy samples ' + \
            'normalized by the cross section at H1 ionization threshold'

        self.He1 = rad_src.PX.sigma_He1( rad_src.E )
        self.He1.__doc__ = \
            'He1 photo-ionization cross section at energy samples'

        self.He1_ratio = self.He1 / rad_src.th.sigma_He1
        self.He1_ratio.__doc__ = \
            'He1 photo-ionization cross section at energy samples ' + \
            'normalized by the cross section at He1 ionization threshold'

        self.He2 = rad_src.PX.sigma_He2( rad_src.E )
        self.He2.__doc__ = \
            'He2 photo-ionization cross section at energy samples'

        self.He2_ratio = self.He2 / rad_src.th.sigma_He2
        self.He2_ratio.__doc__ = \
            'He2 photo-ionization cross section at energy samples ' + \
            'normalized by the cross section at He2 ionization threshold'



class LogQuantities:

    r""" A class that stores log10 quantities (which must be strictly 
    dimensionless). However we define the units the quantity had before the
    log operation in the variable name. """ 

    def __init__( self, rad_src ):

        # put the dimensionless "unit" into a local variable
        #-----------------------------------------------------------------
        dls = rad_src.U.dimensionless 

        # log of photon arrays
        #-----------------------------------------------------------------
        self.q = np.log10( rad_src.q )
        self.q.__doc__ = 'log10 of energy samples / Ry' 

        rad_src.E.units = 'eV'
        self.E_eV = np.log10( rad_src.E.magnitude ) * dls
        self.E_eV.__doc__ = 'log10 of energy samples / eV'

        rad_src.nu.units = 'Hz'
        self.nu_Hz = np.log10( rad_src.nu.magnitude ) * dls
        self.nu_Hz.__doc__ = 'log10 of frequency samples / Hz' 

        rad_src.lam.units = 'cm' 
        self.lam_cm = np.log10( rad_src.lam.magnitude ) * dls
        self.lam_cm.__doc__ = 'log10 of wavelength samples / cm' 


        # log of q_min / q_max
        #-----------------------------------------------------------------
        self.q_min = np.log10( rad_src.q_min ) * dls
        self.q_min.__doc__ = 'log10 of minimum energy sample / Ry'

        self.q_max = np.log10( rad_src.q_max ) * dls
        self.q_max.__doc__ = 'log10 of maximum energy sample / Ry'

        # log of q_H1, q_He1, and q_He2
        #-----------------------------------------------------------------
        self.q_H1 = np.log10( rad_src.th.q_H1 )
        self.q_H1.__doc__ = 'log10 of H1 photo-ionization threshold / Ry'

        self.q_He1 = np.log10( rad_src.th.q_He1 )
        self.q_He1.__doc__ = 'log10 of He1 photo-ionization threshold / Ry'

        self.q_He2 = np.log10( rad_src.th.q_He2 )
        self.q_He2.__doc__ = 'log10 of He2 photo-ionization threshold / Ry'

        # log of min/max frequency in Hz  
        #-----------------------------------------------------------------
        self.nu_min_Hz = np.log10( rad_src.nu.min().magnitude ) * dls
        self.nu_min_Hz.__doc__ = 'log10 of minimum frequency sample / Hz'

        self.nu_max_Hz = np.log10( rad_src.nu.max().magnitude ) * dls
        self.nu_max_Hz.__doc__ = 'log10 of maximum frequency sample / Hz' 

        # log of min/max wavelength in cm  
        #-----------------------------------------------------------------
        self.lam_min_cm = np.log10( rad_src.lam.min().magnitude ) * dls
        self.lam_min_cm.__doc__ = 'log10 of minimum wavelength sample / cm' 

        self.lam_max_cm = np.log10( rad_src.lam.max().magnitude ) * dls
        self.lam_max_cm.__doc__ = 'log10 of maximum wavelength sample / cm' 
