""" Wrapper to make the fortran calling identical to the python calling. """ 

import rabacus_fc 
import numpy as np


from rabacus.utils import utils
from rabacus.constants import units
from rabacus.constants import physical


__all__ = ['ChemistryRates', 'CoolingRates']



class ChemistryRates:

    r""" A chemistry rates class. 

    Provides access to collisional ionization and recombination rates as a 
    function of temperature.  The user can optionally supply photoionization 
    rates which will become attributes of the returned instance.  These 
    chemistry and photoionization rates determine the evolution of the 
    ionization fractions for primordial cosmological chemistry.         


    Args:

      `T` (array): Temperatures.
 
      `fcA_H2` (array): Case A fraction for HII.

      `fcA_He2` (array): Case A fraction for HeII.
      
      `fcA_He3` (array): Case A fraction for HeIII.

    .. note::

      The case A fractions, `fcA_H2`, `fcA_He2`, and `fcA_He3`, 
      interpolate between case A and case B recombination rates 
      (0.0 = case B, 1.0 = case A). 

    Kwargs:

      `H1i` (array): HI photoionization rate.

      `He1i` (array): HeI photoionization rate. 
 
      `He2i` (array): HeII photoionization rate. 

      `fit_name` (string): Source of rate fits {``hg97``}
       
      `add_He2di` (bool): If true, reHe2di is added to reHe2. 

    Attributes:

      `u` (:class:`~rabacus.constants.units.Units`): Units. 

      `ciH1` (array): HI collisional ionization rate.

      `ciHe1` (array): HeI collisional ionization rate.

      `ciHe2` (array): HeII collisional ionization rate.

      `reH2` (array): HII radiative recombination rate. 

      `reHe2` (array): HeII radiative recombination rate. 

      `reHe3` (array): HeIII radiative recombination rate. 

      `reHe2di` (array): HeII dielectronic recombination rate. 

       
    """

    def __init__( self, T, fcA_H2, fcA_He2, fcA_He3, 
                  H1i=None, He1i=None, He2i=None, 
                  fit_name='hg97', add_He2di=True ):



        # attach units
        #----------------------------------------------
        self.u = units.Units()

        # check input
        #----------------------------------------------
        self._clean_input( T, fcA_H2, fcA_He2, fcA_He3, H1i, He1i, He2i )
        self.fit_name = fit_name
        self.add_He2di = add_He2di

        # choose rates source
        #----------------------------------------------
        if fit_name == 'hg97':
            irate = 1
        elif fit_name == 'enzo13':
            irate = 2

        # choose scalar or vector routines
        #----------------------------------------------
        nn = T.size

        if nn == 1:
            get_kchem = rabacus_fc.chem_cool_rates.get_kchem_s
        else:
            get_kchem = rabacus_fc.chem_cool_rates.get_kchem_v


        (reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2) = get_kchem( 
            self.T, self.fcA_H2, self.fcA_He2, self.fcA_He3, 
            irate, nn )
        

        # attach output 
        #----------------------------------------------
        self.reH2 = reH2 * self.u.cm**3 / self.u.s
        self.reHe2 = reHe2 * self.u.cm**3 / self.u.s
        self.reHe3 = reHe3 * self.u.cm**3 / self.u.s

        self.ciH1 = ciH1 * self.u.cm**3 / self.u.s
        self.ciHe1 = ciHe1 * self.u.cm**3 / self.u.s
        self.ciHe2 = ciHe2 * self.u.cm**3 / self.u.s



    def _clean_input( self, T, fcA_H2, fcA_He2, fcA_He3, H1i, He1i, He2i ):

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

        # check HI photoionization rate
        #------------------------------------------
        if H1i == None:
            self.H1i = np.zeros(self.T.shape) / self.u.s
        else:
            if H1i.shape == ():
                self.H1i = np.ones(1) * H1i
            else:
                self.H1i = H1i.copy()

            if not hasattr(self.H1i,'units'): 
                raise utils.NeedUnitsError, '\n H1i must have units \n'
            else:
                self.H1i.units = '1/s'

            if self.H1i.shape != self.T.shape:
                raise utils.InputError, '\n T and H1i must have same shape \n'


        # check HeI photoionization rate
        #------------------------------------------
        if He1i == None:
           self.He1i = np.zeros(1) / self.u.s
        else:
            if He1i.shape == ():
                self.He1i = np.ones(1) * He1i
            else:
                self.He1i = He1i.copy()

            if not hasattr(self.He1i,'units'): 
                raise utils.NeedUnitsError, '\n He1i must have units \n'
            else:
                self.He1i.units = '1/s'

            if self.He1i.shape != self.T.shape:
                raise utils.InputError, '\n T and He1i must have same shape \n'


        # check HeII photoionization rate
        #------------------------------------------
        if He2i == None:
            self.He2i = np.zeros(1) / self.u.s
        else:
            if He2i.shape == ():
                self.He2i = np.ones(1) * He2i
            else:
                self.He2i = He2i.copy()

            if not hasattr(self.He2i,'units'): 
                raise utils.NeedUnitsError, '\n He2i must have units \n'
            else:
                self.He2i.units = '1/s'

            if self.He2i.shape != self.T.shape:
                raise utils.InputError, '\n T and He2i must have same shape \n'





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

        # attach physical constants and units
        #----------------------------------------------
        self.u = units.Units()
        self.pc = physical.PhysicalConstants()

        # check input
        #----------------------------------------------
        self._clean_input( T, fcA_H2, fcA_He2, fcA_He3, H1h, He1h, He2h )
        self.fit_name = fit_name
        self.add_He2di = add_He2di


        # choose rates source
        #----------------------------------------------
        if fit_name == 'hg97':
            irate = 1

        # choose scalar or vector routines
        #----------------------------------------------
        nn = T.size

        if nn == 1:
            get_kcool = rabacus_fc.chem_cool_rates.get_kcool_s
        else:
            get_kcool = rabacus_fc.chem_cool_rates.get_kcool_v


        (recH2, recHe2, recHe3, cicH1, cicHe1, cicHe2, 
         cecH1, cecHe2, bremss, compton) = get_kcool( 
            self.T, self.fcA_H2, self.fcA_He2, self.fcA_He3, 
            irate, nn )
        

        # attach output 
        #----------------------------------------------
        self.recH2 = recH2 * self.u.erg * self.u.cm**3 / self.u.s
        self.recHe2 = recHe2 * self.u.erg * self.u.cm**3 / self.u.s
        self.recHe3 = recHe3 * self.u.erg * self.u.cm**3 / self.u.s

        self.cicH1 = cicH1 * self.u.erg * self.u.cm**3 / self.u.s
        self.cicHe1 = cicHe1 * self.u.erg * self.u.cm**3 / self.u.s
        self.cicHe2 = cicHe2 * self.u.erg * self.u.cm**3 / self.u.s

        self.cecH1 = cecH1 * self.u.erg * self.u.cm**3 / self.u.s
        self.cecHe2 = cecHe2 * self.u.erg * self.u.cm**3 / self.u.s

        self.bremss = bremss * self.u.erg * self.u.cm**3 / self.u.s
        self.compton = compton * self.u.erg / self.u.s



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



    def return_cooling( self, nH, nHe, x, z, Hz=None ):
        """ Returns a cooling rate in erg/(s cm^3). 

        Args: 

          `nH` (array): Number density of hydrogen. 

          `nHe` (array): Number density of helium.

          `x` (ionization fractions): An object which contains ionization 
          fractions.  For example, instances of the classes in 
          :mod:`~rabacus.f2py.ion_solver`.

          `z` (float): Redshift.

        Kwargs:

          `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
          desired.

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

        if Hz != None:
            if not hasattr(Hz,'units'): 
                raise utils.NeedUnitsError, '\n Hz must have units \n'
            else:
                Hz.units = '1/s'

            adia = 3.0 * Hz * self.pc.kb * self.T * ( nH + nHe + ne )
            cool += adia

        return cool


    def return_heating( self, nH, nHe, x ):
        """ Returns a heating rate in erg/(s cm^3). 

        Args: 

          `nH` (array): Number density of hydrogen. 

          `nHe` (array): Number density of helium.

          `x` (ionization fractions): An object which contains ionization 
          fractions.  For example, instances of the classes in 
          :mod:`~rabacus.f2py.ion_solver`.


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
