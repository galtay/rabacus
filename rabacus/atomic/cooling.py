""" A module for loading atomic cooling rates. """ 

import numpy as np
import hui_gnedin_97 
import sd93

from rabacus.utils import utils
from rabacus.constants import units
from rabacus.constants import physical


__all__ = ['CoolingRates']



class CoolingRates:

    r""" A cooling rates class.  

    Provides access to collisional ionization, collisional excitation, 
    recombination, bremsstrahlung, and inverse compton cooling rates as a 
    function of temperature. The user can optionally supply photo-heating rates
    which will become attributes of the returned instance.  These cooling and 
    photo-heating rates determine the evolution of temperature in cosmological 
    gas.  


    .. note:: 

      For compton and bremsstrahlung cooling, only the temperature 
      dependent part of rates are stored.  See the functions 
      :func:`return_cooling` and :func:`return_heating` 
      for the proper multiplicative factors. 


    .. note::

      The case A fractions, `fcA_H2`, `fcA_He2`, and `fcA_He3`, 
      interpolate between case A and case B recombination rates 
      (0.0 = case B, 1.0 = case A). 


    Args:

      `T` (array): Temperatures.
 
      `fcA_H2` (array): Case A fraction for HII.

      `fcA_He2` (array): Case A fraction for HeII.
      
      `fcA_He3` (array): Case A fraction for HeIII.


    Kwargs:

      `H1h` (array): HI photoheating rate.

      `He1h` (array): HeI photoheating rate. 
 
      `He2h` (array): HeII photoheating rate. 

      `fit_name` (string): Source of rate fits {'hg97'}
       
      `add_He2di` (bool): If true, recHe2di is added to recHe2. 


    Attributes:

      `u` (:class:`~rabacus.constants.units.Units`): Units. 

      `cicH1` (array): HI collisional ionization cooling rate.

      `cicHe1` (array): HeI collisional ionization cooling rate.

      `cicHe2` (array): HeII collisional ionization cooling rate.

      `cecH1` (array): HI collisional excitation cooling rate.
    
      `cecHe1` (array): HeI collisional excitation cooling rate.

      `cecHe2` (array): HeII collisional excitation cooling rate.

      `recH2` (array): HII radiative recombination cooling rate. 

      `recHe2` (array): HeII radiative recombination cooling rate. 

      `recHe3` (array): HeIII radiative recombination cooling rate. 
    
      `recHe2di` (array): HeII dielectronic recombination cooling rate. 
    
      `compton` (array): Inverse Compton scattering cooling rate. 

      `bremss` (array): Bremsstrahlung radiation cooling rate. 


    """ 



    def __init__(self, T, fcA_H2, fcA_He2, fcA_He3, 
                 H1h=None, He1h=None, He2h=None, 
                 fit_name='hg97', add_He2di=True):



        # create a list of valid rate names + temperature
        #--------------------------------------------------
        self._rate_names = ['H1h', 'He1h', 'He2h', 
                            'cicH1', 'cicHe1', 'cicHe2',
                            'cecH1', 'cecHe1', 'cecHe2',
                            'recH2', 'recHe2', 'recHe3', 
                            'recHe2di', 
                            'compton', 'bremss', 'T']

        # set optionals once
        # T, H1h, He1h, and He2h can be changed using set
        #----------------------------------------------
        self.fit_name = fit_name
        self.add_He2di = add_He2di
        self.u = units.Units()
        self.pc = physical.PhysicalConstants()

        # set source of rates
        #------------------------------------------
        if fit_name == 'hg97':
            self._rates = hui_gnedin_97.HuiGnedin97()
        elif fit_name == 'enzo13':
            self._rates = enzo_13.Enzo13()
        
        # set initial values
        #----------------------------------------------
        self.set( T, fcA_H2, fcA_He2, fcA_He3, 
                  H1h=H1h, He1h=He1h, He2h=He2h )



    def set( self, T, fcA_H2, fcA_He2, fcA_He3, 
             H1h=None, He1h=None, He2h=None ):

        r""" 
        Updates rates using new temperatures, case A fractions, and 
        photo-heating rates. 

        Args:

          `T` (array): Temperatures.
        
          `fcA_H2` (array): Case A fraction for HII.
        
          `fcA_He2` (array): Case A fraction for HeII.
        
          `fcA_He3` (array): Case A fraction for HeIII.
        
        
        Kwargs:
        
          `H1h` (array): HI photoheating rate.
        
          `He1h` (array): HeI photoheating rate. 
        
          `He2h` (array): HeII photoheating rate. 

        """


        # check input
        #------------------------------------------
        self._clean_input( T, fcA_H2, fcA_He2, fcA_He3, H1h, He1h, He2h )

        # collisional ionization coolling
        #------------------------------------------
        self.cicH1 = self._rates.cicH1(self.T)
        self.cicHe1 = self._rates.cicHe1(self.T)
        self.cicHe2 = self._rates.cicHe2(self.T)

        # collisional excitation coolling
        # (some sources dont provide) 
        #------------------------------------------
        if hasattr( self._rates, 'cecH1' ):
            self.cecH1 = self._rates.cecH1(self.T)
        
        if hasattr( self._rates, 'cecHe1' ):
            self.cecHe1 = self._rates.cecHe1(self.T)

        if hasattr( self._rates, 'cecHe2' ):
            self.cecHe2 = self._rates.cecHe2(self.T)

        # dielectronic recombination coolling
        #------------------------------------------
        self.recHe2di = self._rates.recHe2di(self.T)

        # radiative recombination coolling
        #------------------------------------------
        uu = self._rates.recH2a( 1.0e4 * self.u.K ).units

        recH2a = self._rates.recH2a(self.T)
        recH2b = self._rates.recH2b(self.T)
        log_recH2a = np.log10( recH2a.magnitude )
        log_recH2b = np.log10( recH2b.magnitude )
        log_recH2 = log_recH2a * self.fcA_H2 + \
            log_recH2b * (1.0-self.fcA_H2)
        self.recH2 = 10**log_recH2 * uu
 
        recHe2a = self._rates.recHe2a(self.T)
        recHe2b = self._rates.recHe2b(self.T)
        log_recHe2a = np.log10( recHe2a.magnitude )
        log_recHe2b = np.log10( recHe2b.magnitude )
        log_recHe2 = log_recHe2a * self.fcA_He2 + \
            log_recHe2b * (1.0-self.fcA_He2)
        self.recHe2 = 10**log_recHe2 * uu

        recHe3a = self._rates.recHe3a(self.T)
        recHe3b = self._rates.recHe3b(self.T)
        log_recHe3a = np.log10( recHe3a.magnitude )
        log_recHe3b = np.log10( recHe3b.magnitude )
        log_recHe3 = log_recHe3a * self.fcA_He3 + \
            log_recHe3b * (1.0-self.fcA_He3)
        self.recHe3 = 10**log_recHe3 * uu

        # add dielectronic to radiative?
        #------------------------------------------
        if self.add_He2di:
            self.recHe2 = self.recHe2 + self.recHe2di

        # compton cooling 
        #------------------------------------------
        self.compton = self._rates.compton(self.T)

        # bremsstrahlung cooling 
        #------------------------------------------
        self.bremss = self._rates.bremss(self.T)


    def _clean_input( self, T, fcA_H2, fcA_He2, fcA_He3, H1h, He1h, He2h ):

        # check temperature 
        #------------------------------------------
        if not hasattr( T, 'shape' ):
            raise utils.InputError, '\n T must be an array with units. \n'

        if T.shape == ():
            self.T = np.ones(1) * T
        else:
            self.T = T.copy()

        if not hasattr(self.T,'units'): 
            raise utils.NeedUnitsError, '\n T must have units \n'
        else:
            self.T.units = 'K'

        # check fcA
        #------------------------------------------
        if utils.isiterable( fcA_H2 ):
            assert fcA_H2.shape == self.T.shape
            self.fcA_H2 = fcA_H2.copy()
        else:
            self.fcA_H2 = np.ones(self.T.shape) * fcA_H2

        if utils.isiterable( fcA_He2 ):
            assert fcA_He2.shape == self.T.shape
            self.fcA_He2 = fcA_He2.copy()
        else:
            self.fcA_He2 = np.ones(self.T.shape) * fcA_He2

        if utils.isiterable( fcA_He3 ):
            assert fcA_He3.shape == self.T.shape
            self.fcA_He3 = fcA_He3.copy()
        else:
            self.fcA_He3 = np.ones(self.T.shape) * fcA_He3


        # check HI photo-heating rate
        #------------------------------------------
        if H1h == None:
            self.H1h = np.zeros(1) * self.u.erg / self.u.s
        else:
            if H1h.shape == ():
                self.H1h = np.ones(1) * H1h
            else:
                self.H1h = H1h.copy()

            if not hasattr(self.H1h,'units'): 
                raise utils.NeedUnitsError, '\n H1h must have units \n'
            else:
                self.H1h.units = 'erg/s'

            if self.H1h.shape != self.T.shape:
                raise utils.InputError, '\n T and H1h must have same shape \n'

        # check HeI photo-heating rate
        #------------------------------------------
        if He1h == None:
            self.He1h = np.zeros(1) * self.u.erg / self.u.s
        else:
            if He1h.shape == ():
                self.He1h = np.ones(1) * He1h
            else:
                self.He1h = He1h.copy()

            if not hasattr(self.He1h,'units'): 
                raise utils.NeedUnitsError, '\n He1h must have units \n'
            else:
                self.He1h.units = 'erg/s'

            if self.He1h.shape != self.T.shape:
                raise utils.InputError, '\n T and He1h must have same shape \n'

        # check HeII photo-heating rate
        #------------------------------------------
        if He2h == None:
            self.He2h = np.zeros(1) * self.u.erg / self.u.s
        else:
            if He2h.shape == ():
                self.He2h = np.ones(1) * He2h
            else:
                self.He2h = He2h.copy()

            if not hasattr(self.He2h,'units'): 
                raise utils.NeedUnitsError, '\n He2h must have units \n'
            else:
                self.He2h.units = 'erg/s'

            if self.He2h.shape != self.T.shape:
                raise utils.InputError, '\n T and He2h must have same shape \n'




    def extend_dims( self, N ):
        r""" Takes the one dimensional rate arrays (whose size has been 
        determined by the size of the input temperature array) and extends 
        them to 2-D arrays using copies of the original structures such that 
        the new arrays have dimensions (N, self.T.size). This is useful when
        one wants to use array operations to calculate cooling/heating 
        for all possible combinations of temperature and density.  """ 

        for key in dir( self ):
            if key in self._rate_names:
                if hasattr( self, key ):
                    rate = getattr( self, key )
                    newr = np.array( [ rate for i in xrange(N) ] ) * rate.units
                    setattr( self, key, newr )

    def reduce_dims( self ):
        r""" The inverse of :func:`extend_dims`. """ 

        for key in dir( self ):
            if key in self._rate_names:
                if hasattr( self, key ):
                    rate = getattr( self, key )
                    newr = rate[0,:]
                    setattr( self, key, newr )



    def return_cooling( self, nH, nHe, x, z, Hz=None, sd93_logZ=None ):
        """ Returns a cooling rate in erg/(s cm^3). 

        Args: 

          `nH` (array): Number density of hydrogen. 

          `nHe` (array): Number density of helium.

          `x` (ionization fractions): An object which contains ionization 
          fractions.  For example, instances of the classes in 
          :mod:`~rabacus.solvers.one_zone`.

          `z` (float): Redshift.

        Kwargs:

          `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
          desired.

          `sd93_logZ` (float): set to log metallicity (log Z) to use metal 
          cooling in CIE from Sutherland and Dopita 1993

        """

        if nH.shape == (): nH = np.ones(1) * nH
        if nHe.shape == (): nHe = np.ones(1) * nHe

        c1 = nH.shape == nHe.shape == x.H1.shape == self.T.shape

        if not c1:
            msg = '\n nH, nHe, ionization fractions (x.??) and self.T ' + \
                'must have equal shapes \n'
            print 'nH.shape: ', nH.shape
            print 'nHe.shape: ', nHe.shape
            print 'x.H1.shape: ', x.H1.shape
            print 'self.T.shape: ', self.T.shape
            raise utils.InputError, msg


        ne = nH * x.H2 + nHe * (x.He2 + 2.0 * x.He3)

        cool = \
            self.cecH1  * ne * x.H1  * nH  + \
            self.cecHe2 * ne * x.He2 * nHe + \
            self.cicH1  * ne * x.H1  * nH  + \
            self.cicHe1 * ne * x.He1 * nHe + \
            self.cicHe2 * ne * x.He2 * nHe + \
            self.recH2  * ne * x.H2  * nH  + \
            self.recHe2 * ne * x.He2 * nHe + \
            self.recHe3 * ne * x.He3 * nHe 

        Tcmb = 2.725 * (1.0+z)                 
        cool += self.compton * (1.0+z)**4 * (self.T.magnitude - Tcmb) * ne 

        bfac = nH * x.H2 + nHe * x.He2 + 4.0 * nHe * x.He3
        cool += self.bremss * bfac * ne

        # possibly apply hubble cooling
        #-----------------------------------------------
        if Hz != None:
            if not hasattr(Hz,'units'): 
                raise utils.NeedUnitsError, '\n Hz must have units \n'
            else:
                Hz.units = '1/s'

            adia = 3.0 * Hz * self.pc.kb * self.T * ( nH + nHe + ne )
            cool += adia

        # possibly apply metal cooling
        #-----------------------------------------------
        if sd93_logZ != None:
            sd = sd93.SD93()
            Lambda_N = sd.Zcool( self.T, sd93_logZ )
            nt = nH + nHe
            cool += Lambda_N * ne * nt

        return cool


    def return_heating( self, nH, nHe, x ):
        """ Returns a heating rate in erg/(s cm^3). 

        Args: 

          `nH` (array): Number density of hydrogen. 

          `nHe` (array): Number density of helium.

          `x` (ionization fractions): An object which contains ionization 
          fractions.  For example, instances of the classes in 
          :mod:`~rabacus.solvers.one_zone`

        """

        if nH.shape == (): nH = np.ones(1) * nH
        if nHe.shape == (): nHe = np.ones(1) * nHe

        c1 = nH.shape == nHe.shape == x.H1.shape == self.T.shape

        if not c1:
            msg = '\n nH, nHe, ionization fractions (x.??) and self.T ' + \
                'must have equal shapes \n'
            print 'nH.shape: ', nH.shape
            print 'nHe.shape: ', nHe.shape
            print 'x.H1.shape: ', x.H1.shape
            print 'self.T.shape: ', self.T.shape
            raise utils.InputError, msg
 
        heat = x.H1  * self.H1h  * nH  + \
               x.He1 * self.He1h * nHe + \
               x.He2 * self.He2h * nHe

        return heat
