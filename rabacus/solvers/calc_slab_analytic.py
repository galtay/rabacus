""" Provides an analytic solution for a pure hydrogen constant density 
and temperature slab with monochromatic plane parallel radiation incident from 
one side. 
""" 

import numpy as np
from rabacus.atomic import hydrogen
from rabacus.atomic import chemistry
from rabacus.atomic import photo_xsection


__all__ = ['AnalyticSlab']


NUDGE_FAC = 1.0e-4


class AnalyticSlab:

    r"""
    Solves an analytic slab.  The distance variable used in this solution
    is the fully neutral optical depth, tauH.  To convert this to a distance
    provide the hydrogen photoionization cross section at a given frequency.
    

    Args:
      `nH` (float): constant hydrogen number density
      
      `T` (float); constant temperature

      `H1i_0` (float): hydrogen photoionization rate at surface of slab

      `y` (float): quantifies electrons from elements heavier than H,
      ne = nH * (xH2 + y) (throughout slab)

      `fcA` (float) constant case A recombination fraction. 

    Kwargs:
      `Nsamp` (int): number of spatial samples. 

    """
    
    def __init__( self, nH, T, H1i_0, y, fcA, Nsamp=500 ):

        He1i_0 = 1.0e-10 * H1i_0
        He2i_0 = 1.0e-10 * H1i_0

        self.nH = nH
        self.T = T
        self.H1i_0 = H1i_0
        self.y = y
        self.fcA = fcA
        self.Nsamp = Nsamp

        fcA_H2 = fcA
        fcA_He2 = 0.0; fcA_He3 = 0.0

        self.Hatom = hydrogen.Hydrogen()
        self.kchem = chemistry.ChemistryRates(T, 
                                              fcA_H2, 
                                              fcA_He2, 
                                              fcA_He3, 
                                              H1i=H1i_0, 
                                              He1i=He1i_0, 
                                              He2i=He2i_0 ) 

        # calculate the neutral fraction at the surface of the slab
        #-----------------------------------------------------------
        self.xH1_0 = self.Hatom.analytic_soltn_xH1( 
            self.nH, self.y, self.kchem )

        # calculate the neutral fraction in collisional equilibrium
        #-----------------------------------------------------------
        self.xH1_c = self.kchem.reH2 / (self.kchem.reH2 + self.kchem.ciH1 )
        
        # sample the neutral fraction between the two extremes above
        #-----------------------------------------------------------
        self.log_xH1_0 = np.log10( self.xH1_0.magnitude )
        self.log_xH1_c = np.log10( self.xH1_c.magnitude ) 

        # nudge the extremes a bit so we don't get non-sense
        d_log_xH1 = self.log_xH1_c - self.log_xH1_0
        log_xH1_min = self.log_xH1_0 + NUDGE_FAC * d_log_xH1
        log_xH1_max = self.log_xH1_c - NUDGE_FAC * d_log_xH1

        self.log_xH1 = np.linspace( 
            log_xH1_min, log_xH1_max, self.Nsamp )

        self.xH1 = 10**self.log_xH1

        # calculate the optical depth tauH = nH * sigma(E) 
        # corresponding to the neutral fractions. 
        #----------------------------------------------------------
        self.tauH = self.Hatom.analytic_slab_soltn_tau( 
            self.nH, self.kchem, self.xH1 )

        # calculate the photoionization rates corresponding 
        # to the neutral fractions
        #----------------------------------------------------------
        self.H1i = self.Hatom.analytic_soltn_H1i( 
            self.nH, self.y, self.kchem, self.xH1 )


    def set_E( self, E ):
        """ Given a monochromatic photon energy, adds attributes that measure
        distance using tauH. """ 

        if hasattr(E,'units'): 
            E.units = 'eV'
        else:
            raise utils.NeedUnitsError, '\n Input variable E must have units \n'

        self.E = E
        self.PX = photo_xsection.PhotoXsections()
        self.sigma = self.PX.sigma_H1( self.E )
        self.NH = self.tauH / self.sigma
        self.Lslab = self.NH.max() / self.nH
        self.L = self.NH / self.nH
        
        self.L.units = 'kpc' 
