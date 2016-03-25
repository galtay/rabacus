""" 
A module that stores rate fits from two Badnell 06 papers
http://adsabs.harvard.edu/abs/2006ApJS..167..334B
http://adsabs.harvard.edu/abs/2006A%26A...447..389B
"""

import numpy as np
from rabacus.constants import physical
from rabacus.constants import units



__all__ = ['Badnell06']



class Badnell06:

    r"""    
    Fits to radiative recombination rates from 
    http://adsabs.harvard.edu/abs/2006ApJS..167..334B

    Fits from dielectronic recombination rates from
    http://adsabs.harvard.edu/abs/2006A%26A...447..389B

    Attributes:

      `u` (:class:`~rabacus.constants.units.Units`)

      `pc` (:class:`~rabacus.constants.physical.PhysicalConstants`)


    
    """

    def __init__(self): 

        self.pc = physical.PhysicalConstants()
        self.u = units.Units()

    
    def _check_input_T(self, T):
        msg = '\n Input variable T must have units of K \n'
        if not hasattr(T,'units'): 
            raise ValueError(msg)
            if not T.units == 'K':
                raise ValueError(msg)

        if T.ndim == 0:
            T = np.array( [T] ) * T.units

        return T





    #-----------------------------------------------------------
    # CHEMISTRY
    #-----------------------------------------------------------

    # recombination rates (caseA) 
    #-----------------------------------------------------------
    def reH2a(self, Tin):
        """ Fit to H2 (Z=1,N=0) recombination rate (caseA) """
        T = self._check_input_T(Tin)
        Tm = T.magnitude
        A = 8.318e-11 
        B = 0.7472
        T0 = 2.965 
        T1 = 7.001e5 
        t1 = np.sqrt(Tm/T0)
        t2 = (1+t1)**(1.0-B)
        t3 = ( 1 + np.sqrt(Tm/T1) )**(1.0+B)
        rate = A / (t1*t2*t3) * self.u.cm**3 / self.u.s
        return rate 

    def reHe2a(self, Tin):
        """ Fit to He2 (Z=2, N=1) recombination rate (caseA) """
        T = self._check_input_T(Tin)
        Tm = T.magnitude
        A = 5.235e-11
        B = 0.6988
        T0 = 7.301 
        T1 = 4.475e6 
        C = 0.0829
        T2 = 1.682e5 
        B = B + C * np.exp(-T2/Tm)
        t1 = np.sqrt(Tm/T0)
        t2 = (1+t1)**(1.0-B)
        t3 = ( 1 + np.sqrt(Tm/T1) )**(1.0+B)
        rate = A / (t1*t2*t3) * self.u.cm**3 / self.u.s
        return rate 


    def reHe3a(self, Tin):
        """ Fit to He3 (Z=2, N=0) recombination rate (caseA) """
        T = self._check_input_T(Tin)
        Tm = T.magnitude
        A = 1.818e-10 
        B = 0.7492
        T0 = 1.017e1 
        T1 = 2.786e6 
        t1 = np.sqrt(Tm/T0)
        t2 = (1+t1)**(1.0-B)
        t3 = ( 1 + np.sqrt(Tm/T1) )**(1.0+B)
        rate = A / (t1*t2*t3) * self.u.cm**3 / self.u.s
        return rate 



    # dielectronic recombination rates 
    #-----------------------------------------------------------
    def reHe2di(self, Tin):
        """ Fit to He2 dielectronic recombination rate """
        T = self._check_input_T(Tin)
        Tm = T.magnitude
        c1 = 5.966e-4 
        c2 = 1.613e-4 
        c3 = -2.223e-5 
        E1 = 4.556e5 
        E2 = 5.552e5 
        E3 = 8.982e5 

        rate = Tm**(-3./2) * \
               ( c1*np.exp(-E1/Tm) + c2*np.exp(-E2/Tm) + c3*np.exp(-E3/Tm) )

        rate = rate * self.u.cm**3 / self.u.s
   
        return rate 







