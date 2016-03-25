""" 
A module that stores rate fits from Hui & Gnedin 97 
http://adsabs.harvard.edu/abs/1997MNRAS.292...27H 
"""

import numpy as np
from rabacus.constants import physical
from rabacus.constants import units



__all__ = ['HuiGnedin97']



class HuiGnedin97:

    r"""    
    Fits to atomic data from 
    http://adsabs.harvard.edu/abs/1997MNRAS.292...27H 

    Ferland3_1e9 are fits to the data from Ferland 1992 
    and are accurate to 2% from 3K - 1e9K
    http://adsabs.harvard.edu/abs/1992ApJ...387...95F

    BurgessSeaton5e3_5e5 are fits to the data from Burgess
    and Seaton 1960 and are accurate to 10% from 5e3K - 5e5K
    http://adsabs.harvard.edu/abs/1960MNRAS.121..471B

    Lotz1e4_1e9 are fits to the data from Lotz 1967 
    and are accurate to 3% from 1e4K to 1e9K
    http://adsabs.harvard.edu/abs/1967ApJS...14..207L

    AP3e4_1e6 are fits to the data from Aldrovandi and Pequignot 1973
    and are accurate to 5% from 3e4K to 1e6K
    http://adsabs.harvard.edu/abs/1973A%26A....25..137A

    Black5e3_5e5 are from Black 1981 with the correction from Cen 1992
    and are accurate to 10% from 5e3 to 5e5
    http://adsabs.harvard.edu/abs/1981MNRAS.197..553B


    Attributes:

      `u` (:class:`~rabacus.constants.units.Units`)

      `pc` (:class:`~rabacus.constants.physical.PhysicalConstants`)

      `T_H1` (float): HI ionization threshold in temperature units. 

      `T_He1` (float): HeI ionization threshold in temperature units. 

      `T_He2` (float): HeII ionization threshold in temperature units. 

    
    """

    def __init__(self): 

        self.pc = physical.PhysicalConstants()
        self.u = units.Units()

        self.T_H1 = 157807. * self.u.K
        self.T_He1 = 285335. * self.u.K
        self.T_He2 = 631515. * self.u.K

    
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


    # collisional ionization rates 
    #-----------------------------------------------------------
    def ciH1(self, Tin):
        """ Fit to H1 collisional ionization rate (Lotz1e4_1e9). """        
        T = self._check_input_T(Tin)        
        lam = 2 * self.T_H1/T
        num = 21.11 * lam**(-1.089)
        den = 1 + (lam/0.354)**(0.874)    
        rate = num / den**(1.101)    
        rate = rate * np.exp(-0.5*lam)    
        rate = rate * T**(-1.5) * self.u.K**(1.5)    
        return rate * self.u.cm**3 / self.u.s

    def ciHe1(self, Tin):
        """ Fit to He1 collisional ionization rate (Lotz1e4_1e9) """        
        T = self._check_input_T(Tin)        
        lam = 2 * self.T_He1/T
        num = 32.38 * lam**(-1.146)
        den = 1 + (lam/0.416)**(0.987)    
        rate = num / den**(1.056)    
        rate = rate * np.exp(-0.5*lam)    
        rate = rate * T**(-1.5) * self.u.K**(1.5)    
        return rate * self.u.cm**3 / self.u.s

    def ciHe2(self, Tin):
        """ Fit to He2 collisional ionization rate (Lotz1e4_1e9) """   
        T = self._check_input_T(Tin)        
        lam = 2 * self.T_He2/T
        num = 19.95 * lam**(-1.089)
        den = 1 + (lam/0.553)**(0.735)    
        rate = num / den**(1.275)    
        rate = rate * np.exp(-0.5*lam)    
        rate = rate * T**(-1.5) * self.u.K**(1.5)    
        return rate * self.u.cm**3 / self.u.s


    # recombination rates (caseA) 
    #-----------------------------------------------------------
    def reH2a(self, Tin):
        """ Fit to H2 recombination rate (caseA) (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_H1/T
        num = 1.269e-13 * lam**(1.503)
        den = 1 + (lam/0.522)**(0.470)
        rate = num / den**(1.923)
        return rate * self.u.cm**3 / self.u.s

    def reHe2a(self, Tin):
        """ Fit to He2 recombination rate (caseA) (BurgessSeaton5e3_5e5) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He1/T
        rate = 3.0e-14 * lam**(0.654)
        rate = rate * self.u.cm**3 / self.u.s
        return rate

    def reHe3a(self, Tin):
        """ Fit to He3 recombination rate (caseA) (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He2/T
        num = 2 * 1.269e-13 * lam**(1.503)
        den = 1 + (lam/0.522)**(0.470)
        rate = num / den**(1.923)
        return rate * self.u.cm**3 / self.u.s


    # recombination rates (caseB) 
    #-----------------------------------------------------------
    def reH2b(self, Tin):
        """ Fit to H2 recombination rate (caseB) (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_H1/T
        num = 2.753e-14 * lam**(1.500)
        den = 1 + (lam/2.740)**(0.407)
        rate = num / den**(2.242)
        return rate * self.u.cm**3 / self.u.s

    def reHe2b(self, Tin):
        """ Fit to He2 recombination rate (caseB) (BurgessSeaton5e3_5e5) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He1/T
        rate = 1.26e-14 * lam**(0.750)
        rate = rate * self.u.cm**3 / self.u.s
        return rate

    def reHe3b(self, Tin):
        """ Fit to He3 recombination rate (caseB) (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He2/T
        num = 2 * 2.753e-14 * lam**(1.500)
        den = 1 + (lam/2.740)**(0.407)
        rate = num / den**(2.242)
        return rate * self.u.cm**3 / self.u.s


    # dielectronic recombination rates 
    #-----------------------------------------------------------
    def reHe2di(self, Tin):
        """ Fit to He2 dielectronic recombination rate (AP3e4_1e6) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He2/T
        rate = 1.90e-3 * np.exp(-0.75 * 0.5 * lam)
        rate = rate * (1.0 + 0.3*np.exp(-0.15 * 0.5 * lam))
        rate = rate * T**(-1.5) * self.u.K**(1.5)    
        return rate * self.u.cm**3 / self.u.s


    #-----------------------------------------------------------
    # COOLING
    #-----------------------------------------------------------


    # collisional ionization cooling rates 
    #-----------------------------------------------------------
    def cicH1(self, Tin):
        """ Fit to H1 collisional ionization cooling rate (Lotz1e4_1e9) """  
        T = self._check_input_T(Tin)        
        ciH1 = self.ciH1(Tin)
        rate = self.pc.kb * self.T_H1 * ciH1
        return rate 

    def cicHe1(self, Tin):
        """ Fit to He1 collisional ionization cooling rate (Lotz1e4_1e9) """   
        T = self._check_input_T(Tin)        
        ciHe1 = self.ciHe1(Tin)
        rate = self.pc.kb * self.T_He1 * ciHe1
        return rate 

    def cicHe2(self, Tin):
        """ Fit to He2 collisional ionization cooling rate (Lotz1e4_1e9) """  
        T = self._check_input_T(Tin)        
        ciHe2 = self.ciHe2(Tin)
        rate = self.pc.kb * self.T_He2 * ciHe2
        return rate 


    # collisional excitation cooling rates 
    #-----------------------------------------------------------
    def cecH1(self, Tin):
        """ Fit to H1 collisional excitation cooling rate (Lotz1e4_1e9) """  
        T = self._check_input_T(Tin)        
        lam = 2 * self.T_H1/T
        num = 7.5e-19 * np.exp( -0.75 * lam / 2 )
        den = 1.0 + np.sqrt(T.magnitude/1.0e5)
        rate = num / den
        return rate * self.u.erg * self.u.cm**3 / self.u.s 



    def cecHe2(self, Tin):
        """ Fit to He2 collisional excitation cooling rate (Lotz1e4_1e9) """   
        T = self._check_input_T(Tin)        
        lam = 2 * self.T_He2/T
        num = 5.54e-17 * (1.0/T.magnitude)**0.397 * np.exp( -0.75 * lam / 2 )
        den = 1.0 + np.sqrt(T.magnitude/1.0e5)
        rate = num / den
        return rate * self.u.erg * self.u.cm**3 / self.u.s 


    # recombination cooling rates (caseA) 
    #-----------------------------------------------------------
    def recH2a(self, Tin):
        """ Fit to H2 recombination cooling rate (caseA) (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_H1/T
        num = 1.778e-29 * lam**(1.965)
        den = 1 + (lam/0.541)**(0.502)
        rate = T * num / den**(2.697)
        return rate * self.u.erg * self.u.cm**3 / self.u.s / self.u.K

    def recHe2a(self, Tin):
        """ Fit to He2 recombination cooling rate (caseA) 
        (BurgessSeaton5e3_5e5) """
        T = self._check_input_T(Tin)
        reHe2a = self.reHe2a(Tin)
        rate = self.pc.kb * T * reHe2a
        return rate

    def recHe3a(self, Tin):
        """ Fit to He3 recombination cooling rate (caseA) 
        (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He2/T
        num = 8 * 1.778e-29 * lam**(1.965)
        den = 1 + (lam/0.541)**(0.502)
        rate = T * num / den**(2.697)
        return rate * self.u.erg * self.u.cm**3 / self.u.s / self.u.K


    # recombination cooling rates (caseB) 
    #-----------------------------------------------------------
    def recH2b(self, Tin):
        """ Fit to H2 recombination cooling rate (caseB) 
        (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_H1/T
        num = 3.435e-30 * lam**(1.970)
        den = 1 + (lam/2.250)**(0.376)
        rate = T * num / den**(3.720)
        return rate * self.u.erg * self.u.cm**3 / self.u.s / self.u.K

    def recHe2b(self, Tin):
        """ Fit to He2 recombination cooling rate (caseB)
        (BurgessSeaton5e3_5e5) """
        T = self._check_input_T(Tin)
        reHe2b = self.reHe2b(Tin)
        rate = self.pc.kb * T * reHe2b
        return rate

    def recHe3b(self, Tin):
        """ Fit to He3 recombination cooling rate (caseB) 
        (Ferland3_1e9) """
        T = self._check_input_T(Tin)
        lam = 2 * self.T_He2/T
        num = 8 * 3.435e-30 * lam**(1.970)
        den = 1 + (lam/2.250)**(0.376)
        rate = T * num / den**(3.720)
        return rate * self.u.erg * self.u.cm**3 / self.u.s / self.u.K


    # dielectronic recombination cooling rates 
    #-----------------------------------------------------------
    def recHe2di(self, Tin):
        """ Fit to He2 dielectronic recombination cooling rate 
        (AP3e4_1e6) """
        T = self._check_input_T(Tin)
        reHe2di = self.reHe2di(Tin)
        rate = 0.75 * self.pc.kb * self.T_He2 * reHe2di
        return rate 



    # Compton cooling rate
    #-----------------------------------------------------------
    def compton(self, Tin ):
        """ Fit to temperature dependent part of compton cooling rate. """
        T = self._check_input_T(Tin)
        Tm = T.magnitude
        rate = np.ones( Tm.shape ) * 5.65e-36
        return rate * self.u.erg / self.u.s 

    # Bremsstrahlung cooling rate
    #-----------------------------------------------------------
    def bremss(self, Tin ):
        """ Fit to temperature dependent part of bremsstrahlung cooling. """ 
        T = self._check_input_T(Tin)
        Tm = T.magnitude

        logT = np.log10( Tm )
        gff = 1.1 + 0.34 * np.exp( -(5.5 - logT)**2.0/3.0 )

        rate = 1.43e-27 * np.sqrt(Tm) * gff
        return rate * self.u.erg * self.u.cm**3 / self.u.s 



#
#
#    # Compton cooling rate
#    #-----------------------------------------------------------
#    def compton(self, Tin, z):
#        """ From Eq. 24.45 in Peebles "Principles of Physical Cosmology"
#        Note: Tcmb should be the temperature of the CMB at the redshift
#        of interest and the returned value should be multiplied by n_e 
#        Tcmb = Tcmb,0 * (1+z) ~ 2.725 * (1+z) K
#        """
#
#        T = self._check_input_T(Tin)
#
#        # calculate temperature of CMB
#        Tcmb0 = 2.725 * T.units
#        Tcmb = Tcmb0 * (1+z)
#
#        # calculate energy density of CMB
#        fac1 = (self.pc.kb * Tcmb)**4 / (self.pc.h*self.pc.c)**3
#        u = 8 * np.pi * fac1 * np.pi**4/15
#
#        # calculate temperature difference term
#        fac2 = (4 * self.pc.sigma_t * self.pc.kb) / (self.pc.m_e * self.pc.c)
#
#        rate = u * fac2 * (T-Tcmb)
#        rate.units = 'erg/s'
#        return rate



