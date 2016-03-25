""" A module for loading atomic chemistry rates. """ 

import numpy as np
import hui_gnedin_97 
import badnell_06

from rabacus.utils import utils
from rabacus.constants import units



__all__ = ['ChemistryRates']



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
                  fit_name='hg97', add_He2di=True, 
                  use_badnell=False ):

        
        # create a list of valid rate names + temperature
        #--------------------------------------------------
        self._rate_names = ['H1i', 'He1i', 'He2i', 
                            'ciH1', 'ciHe1', 'ciHe2',
                            'reH2', 'reHe2', 'reHe3', 
                            'reHe2di', 'T']

        # set optionals once
        # T, fcA, H1i, He1i, and He2i can be changed using set
        #----------------------------------------------
        self.fit_name = fit_name
        self.add_He2di = add_He2di
        self.U = units.Units()
        self.use_badnell = use_badnell

        # set source of rates
        #------------------------------------------
        if fit_name == 'hg97':
            self._rates = hui_gnedin_97.HuiGnedin97()
        
        # set initial values
        #----------------------------------------------
        self.set( T, fcA_H2, fcA_He2, fcA_He3, 
                  H1i=H1i, He1i=He1i, He2i=He2i )



    def set( self, T, fcA_H2, fcA_He2, fcA_He3, 
             H1i=None, He1i=None, He2i=None ):

        r""" 
        Updates rates using new temperatures, case A fractions, and 
        photoionization rates. 

        Args:

          `T` (array): Temperatures.
        
          `fcA_H2` (array): Case A fraction for HII.
        
          `fcA_He2` (array): Case A fraction for HeII.
        
          `fcA_He3` (array): Case A fraction for HeIII.
        
        Kwargs:
        
          `H1i` (array): HI photoionization rate.
        
          `He1i` (array): HeI photoionization rate. 
        
          `He2i` (array): HeII photoionization rate. 
        
        """


        # check input
        #------------------------------------------
        self._clean_input( T, fcA_H2, fcA_He2, fcA_He3, H1i, He1i, He2i )

        # collisional ionization 
        #------------------------------------------
        self.ciH1 = self._rates.ciH1(self.T)
        self.ciHe1 = self._rates.ciHe1(self.T)
        self.ciHe2 = self._rates.ciHe2(self.T)

        # dielectronic recombination
        #------------------------------------------
        self.reHe2di = self._rates.reHe2di(self.T)

        # radiative recombination
        #------------------------------------------
        uu = self._rates.reH2a( 1.0e4 * self.U.K ).units

        reH2a = self._rates.reH2a(self.T)
        reH2b = self._rates.reH2b(self.T)
        log_reH2a = np.log10( reH2a.magnitude )
        log_reH2b = np.log10( reH2b.magnitude )
        log_reH2 = log_reH2a * self.fcA_H2 + log_reH2b * (1.0-self.fcA_H2)
        self.reH2 = 10**log_reH2 * uu
 
        reHe2a = self._rates.reHe2a(self.T)
        reHe2b = self._rates.reHe2b(self.T)
        log_reHe2a = np.log10( reHe2a.magnitude )
        log_reHe2b = np.log10( reHe2b.magnitude )
        log_reHe2 = log_reHe2a * self.fcA_He2 + log_reHe2b * (1.0-self.fcA_He2)
        self.reHe2 = 10**log_reHe2 * uu

        reHe3a = self._rates.reHe3a(self.T)
        reHe3b = self._rates.reHe3b(self.T)
        log_reHe3a = np.log10( reHe3a.magnitude )
        log_reHe3b = np.log10( reHe3b.magnitude )
        log_reHe3 = log_reHe3a * self.fcA_He3 + log_reHe3b * (1.0-self.fcA_He3)
        self.reHe3 = 10**log_reHe3 * uu

        # add dielectronic to radiative?
        #------------------------------------------
        if self.add_He2di:
            self.reHe2 = self.reHe2 + self.reHe2di

        # replace recombination rates with those from badnell?
        #------------------------------------------
        if self.use_badnell:

            b06 = badnell_06.Badnell06()
            self.reH2 = b06.reH2a(self.T)
            self.reHe2 = b06.reHe2a(self.T)
            self.reHe2di = b06.reHe2di(self.T)
            self.reHe3 = b06.reHe3a(self.T)
            if self.add_He2di:
                self.reHe2 = self.reHe2 + self.reHe2di



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
            self.H1i = np.zeros(self.T.shape) / self.U.s
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
            self.He1i = np.zeros(1) / self.U.s
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
            self.He2i = np.zeros(1) / self.U.s
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





    def extend_dims( self, N ):
        r""" Takes the one dimensional rate arrays (whose size has been 
        determined by the size of the input temperature array) and extends 
        them to 2-D arrays using copies of the original structures such that 
        the new arrays have dimensions (N, self.T.size). This is useful when
        one wants to use array operations to calculate ionization fractions
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
