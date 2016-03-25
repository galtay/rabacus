""" Sphere geometry base class. """ 

import rabacus_fc
import numpy as np
from rabacus.utils import utils
from rabacus.constants import units
from rabacus.atomic import chemistry


__all__ = ['SphereBase']


class SphereBase(object):

    """ Base class for spheres. """ 


    def __init__( self ):


        # attach units
        #-----------------------------------------------
        self.U = units.Units()

        # check input 
        #-----------------------------------------------
        assert( self.Edges.size == self.T.size + 1 )
        assert( self.T.shape == self.nH.shape == self.nHe.shape )
        assert( self.T.shape == (self.T.size,) )

        assert( self.rec_meth == 'fixed' or 
                self.rec_meth == 'outward' or
                self.rec_meth == 'radial' or
                self.rec_meth == 'isotropic' )

        if self.find_Teq and self.z==None:
            msg = 'need to supply redshift if finding equilibrium temperature'
            raise ValueError(msg)


        # set units
        #-----------------------------------------------
        self.Edges.units = 'cm'
        self.T.units = 'K'
        self.nH.units = '1.0/cm^3'
        self.nHe.units = '1.0/cm^3'


        # format input
        #-----------------------------------------------
        self.format_for_fortran()



    def format_for_fortran( self ):

        # re-form the input for fortran
        #---------------------------------------------------
        self.E_eV = self.rad_src.E
        self.E_eV.units = 'eV'

        if self.rad_src.source_type == 'point':
            if self.rad_src.monochromatic:
                self.shape = self.rad_src.Lu
            else:
                self.shape = self.rad_src.dLu_over_dnu

        elif self.rad_src.source_type == 'plane':
            if self.rad_src.monochromatic:
                self.shape = self.rad_src.Fu
            else:
                self.shape = self.rad_src.dFu_over_dnu

        elif self.rad_src.source_type == 'background':
            if self.rad_src.monochromatic:
                self.shape = self.rad_src.Inu
            else:
                self.shape = self.rad_src.Inu

        else:
            raise ValueError('source type not recognized')


        if self.rec_meth == 'fixed':
            self.i_rec_meth = 1
        elif self.rec_meth == 'outward':
            self.i_rec_meth = 2
        elif self.rec_meth == 'radial':
            self.i_rec_meth = 3
        elif self.rec_meth == 'isotropic':
            self.i_rec_meth = 4

        if self.rad_src.px_fit_type == 'verner':
            self.i_photo_fit = 1

        if self.atomic_fit_name == 'hg97':
            self.i_rate_fit = 1

        if self.find_Teq:
            self.i_find_Teq = 1
        else:
            self.i_find_Teq = 0

        if self.thin:
            self.i_thin = 1
        else:
            self.i_thin = 0
        
        self.Nl = self.nH.size
        self.Nnu = self.E_eV.size




    def __post__( self ):

        """ Post calculation variables """


        # cgs units 
        #------------------------------------------------        
        self.dr = np.abs( self.Edges[1:] - self.Edges[0:-1] )
        self.dr.units = 'cm' 

        self.r_c = (self.Edges[0:-1] + self.Edges[1:]) * 0.5 
        self.r_c.units = 'cm'


        # electron density 
        #------------------------------------------------        
        self.ne = self.nH * self.xH2 + \
            self.nHe * ( self.xHe2 + 2.0 * self.xHe3 ) 


        # calculate effective photoionization rates
        # (i.e. the photoionization rates that result in 
        # the same ionization fractions if caseA rates 
        # are used).  Note, we assume that the recombination 
        # radiation is peaked near the ionization thresholds 
        # and therefore does not alter the temperature.
        #------------------------------------------------        
        kA = chemistry.ChemistryRates( 
            self.T, fcA_H2=1.0, fcA_He2=1.0, fcA_He3=1.0,
            fit_name=self.atomic_fit_name )        
        
        self.H1i_eff = kA.reH2 * self.ne * self.xH2 / self.xH1 - \
            kA.ciH1 * self.ne

        self.He1i_eff = kA.reHe2 * self.ne * self.xHe2 / self.xHe1 - \
            kA.ciHe1 * self.ne

        self.He2i_eff = kA.reHe3 * self.ne * self.xHe3 / self.xHe2 - \
            kA.ciHe2 * self.ne

        # total column densities
        #------------------------------------------------

        self.NH1_thru = np.sum( self.dr * self.nH * self.xH1 ) * 2
        self.logNH1_thru = np.log10( self.NH1_thru.magnitude ) 

        self.NHe1_thru = np.sum( self.dr * self.nHe * self.xHe1 ) * 2
        self.logNHe1_thru = np.log10( self.NHe1_thru.magnitude ) 

        self.NHe2_thru = np.sum( self.dr * self.nHe * self.xHe2 ) * 2
        self.logNHe2_thru = np.log10( self.NHe2_thru.magnitude ) 
        

        # differential column densities
        #------------------------------------------------
        self.dNH1 = self.dr * self.nH * self.xH1
        self.dNHe1 = self.dr * self.nHe * self.xHe1
        self.dNHe2 = self.dr * self.nHe * self.xHe2

        # optical depths at threshold
        #------------------------------------------------
        self.dtau_H1_th = self.dNH1 * self.rad_src.th.sigma_H1
        self.dtau_He1_th = self.dNHe1 * self.rad_src.th.sigma_He1
        self.dtau_He2_th = self.dNHe2 * self.rad_src.th.sigma_He2

        # column density as a function of impact parameter
        #------------------------------------------------
        (NH, NHe, NH1, NHe1, NHe2) = \
            rabacus_fc.sphere_base.set_column_vs_impact(
            self.xH1, 
            self.xH2,
            self.xHe1,
            self.xHe2,
            self.xHe3,
            self.Edges, 
            self.nH, 
            self.nHe, 
            self.T, 
            self.Nl )
            
        self.NH_b = NH * self.NH1_thru.units
        self.NHe_b = NHe * self.NH1_thru.units
        self.NH1_b = NH1 * self.NH1_thru.units
        self.NHe1_b = NHe1 * self.NH1_thru.units
        self.NHe2_b = NHe2 * self.NH1_thru.units

        self.NH2_b = self.NH_b - self.NH1_b
        self.NHe3_b = self.NHe_b - self.NHe1_b - self.NHe2_b

        # check units
        #------------------------------------------------
        self.r_c.units = 'kpc'

