""" Slab geometry base class. """ 

import numpy as np
from rabacus.utils import utils
from rabacus.constants import units
from rabacus.atomic import chemistry
from rabacus.hdf5 import h5py_wrap as h5
from rabacus.cosmology import jeans


__all__ = ['SlabBase']


class Averaged:
    """ Averaged Quantities. """ 
    pass 

class Central:
    """ Quantities at the center of the slab. """ 
    pass 

class Weighted:
    """ Weighted integrated quantities. """ 
    pass



class SlabBase(object):

    """ Base class for slabs. """ 

    def __init__( self ):

        """ Initialization for all slab classes """
        
        # attach units
        #-----------------------------------------------
        self.U = units.Units()

        # check input 
        #-----------------------------------------------
        assert( self.Edges.size == self.T.size + 1 )
        assert( self.T.shape == self.nH.shape == self.nHe.shape )
        assert( self.T.shape == (self.T.size,) )

        assert( self.rec_meth == 'fixed' or 
                self.rec_meth == 'ray' )

        if self.find_Teq and self.z==None:
            msg = 'need to supply redshift if finding equilibrium temperature'
            raise utils.InputError(msg)


        # set units
        #-----------------------------------------------
        self.Edges.units = 'cm' 
        self.T.units = 'K'
        self.nH.units = '1.0/cm^3'
        self.nHe.units = '1.0/cm^3'


        # format input
        #-----------------------------------------------
        self.format_for_fortran()


    def format_for_fortran(self):

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
        elif self.rec_meth == 'ray':
            self.i_rec_meth = 2

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



    def calc_weighted( self, wa ):
        """ Calculates all integrals weighted by wa. """ 

        wq = Averaged()

        xx = self.z_c

        # <xH1>_w 
        #------------------------------------------------
        yy = wa * self.xH1
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.xH1 = num / den

        # <xH2>_w 
        #------------------------------------------------
        yy = wa * self.xH2
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.xH2 = num / den

        # <xHe1>_w 
        #------------------------------------------------
        yy = wa * self.xHe1
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.xHe1 = num / den

        # <xHe2>_w 
        #------------------------------------------------
        yy = wa * self.xHe2
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.xHe2 = num / den

        # <xHe3>_w 
        #------------------------------------------------
        yy = wa * self.xHe3
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.xHe3 = num / den

        # <T>_w 
        #------------------------------------------------
        yy = wa * self.T
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.T = num / den

        # <mu>_w 
        #------------------------------------------------
        yy = wa * self.mu
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.mu = num / den

        # <H1i>_w 
        #------------------------------------------------
        yy = wa * self.H1i
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1i = num / den

        yy = wa * self.H1i_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1i_src = num / den

        yy = wa * self.H1i_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1i_rec = num / den

        # <He1i>_w 
        #------------------------------------------------
        yy = wa * self.He1i
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1i = num / den

        yy = wa * self.He1i_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1i_src = num / den

        yy = wa * self.He1i_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1i_rec = num / den

        # <He2i>_w 
        #------------------------------------------------
        yy = wa * self.He2i
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2i = num / den

        yy = wa * self.He2i_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2i_src = num / den

        yy = wa * self.He2i_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2i_rec = num / den

        # <H1h>_w 
        #------------------------------------------------
        yy = wa * self.H1h
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1h = num / den

        yy = wa * self.H1h_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1h_src = num / den

        yy = wa * self.H1h_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.H1h_rec = num / den

        # <He1h>_w 
        #------------------------------------------------
        yy = wa * self.He1h
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1h = num / den

        yy = wa * self.He1h_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1h_src = num / den

        yy = wa * self.He1h_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He1h_rec = num / den

        # <He2h>_w 
        #------------------------------------------------
        yy = wa * self.He2h
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2h = num / den

        yy = wa * self.He2h_src
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2h_src = num / den

        yy = wa * self.He2h_rec
        num = utils.trap( xx, yy )
        den = utils.trap( xx, wa )
        wq.He2h_rec = num / den

        return wq




    def __post__( self ):

        """ Post calculation variables """


        # make ionization fractions quantities
        #------------------------------------------------        
        self.xH1 = self.xH1 * self.U.dimensionless
        self.xH2 = self.xH2 * self.U.dimensionless

        self.xHe1 = self.xHe1 * self.U.dimensionless
        self.xHe2 = self.xHe2 * self.U.dimensionless
        self.xHe3 = self.xHe3 * self.U.dimensionless

        # calculate mean molecular weight
        #------------------------------------------------        
        J = jeans.Jeans()
        self.mu = J.mu( self.xH2, self.xHe2, self.xHe3 )

        # create averaged values attribute
        #------------------------------------------------        
        self.avg = Averaged()

        # cgs units 
        #------------------------------------------------        
        self.dz = self.Edges[1:] - self.Edges[0:-1]
        self.dz.units = 'cm' 

        self.z_c = self.Edges[0:-1] + 0.5 * self.dz 
        self.z_c.units = 'cm'


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
        kchemA = chemistry.ChemistryRates( 
            self.T, fcA_H2 = 1.0, fcA_He2 = 1.0, fcA_He3 = 1.0,
            H1i = self.H1i, He1i = self.He1i, He2i = self.He2i, 
            fit_name = self.atomic_fit_name )        
        
        self.H1i_eff = kchemA.reH2 * self.ne * self.xH2 / self.xH1 - \
            kchemA.ciH1 * self.ne

        self.He1i_eff = kchemA.reHe2 * self.ne * self.xHe2 / self.xHe1 - \
            kchemA.ciHe1 * self.ne

        self.He2i_eff = kchemA.reHe3 * self.ne * self.xHe3 / self.xHe2 - \
            kchemA.ciHe2 * self.ne

        # total column densities
        #------------------------------------------------
        self.NH1_thru = np.sum( self.dz * self.nH * self.xH1 )
        self.logNH1_thru = np.log10( self.NH1_thru.magnitude )

        self.NHe1_thru = np.sum( self.dz * self.nHe * self.xHe1 )
        self.logNHe1_thru = np.log10( self.NHe1_thru.magnitude )

        self.NHe2_thru = np.sum( self.dz * self.nHe * self.xHe2 )
        self.logNHe2_thru = np.log10( self.NHe2_thru.magnitude )
        

        # differential column densities
        #------------------------------------------------
        self.dNH1 = self.dz * self.nH * self.xH1
        self.dNHe1 = self.dz * self.nHe * self.xHe1
        self.dNHe2 = self.dz * self.nHe * self.xHe2

        # optical depths at threshold
        #------------------------------------------------
        self.dtau_H1_th = self.dNH1 * self.rad_src.th.sigma_H1
        self.dtau_He1_th = self.dNHe1 * self.rad_src.th.sigma_He1
        self.dtau_He2_th = self.dNHe2 * self.rad_src.th.sigma_He2


        # central values
        #================================================
        #================================================
        Nl2 = self.Nl / 2

        central = Central() 

        central.T = self.T[Nl2]
        central.mu = self.mu[Nl2]
        central.ne = self.ne[Nl2]
        central.nH = self.nH[Nl2]
        central.nHe = self.nHe[Nl2]

        central.xH1 = self.xH1[Nl2]
        central.xH2 = self.xH2[Nl2]
        central.xHe1 = self.xHe1[Nl2]
        central.xHe2 = self.xHe2[Nl2]
        central.xHe3 = self.xHe3[Nl2]

        central.H1i = self.H1i[Nl2]
        central.He1i = self.He1i[Nl2]
        central.He2i = self.He2i[Nl2]

        central.H1i_src = self.H1i_src[Nl2]
        central.He1i_src = self.He1i_src[Nl2]
        central.He2i_src = self.He2i_src[Nl2]

        central.H1i_rec = self.H1i_rec[Nl2]
        central.He1i_rec = self.He1i_rec[Nl2]
        central.He2i_rec = self.He2i_rec[Nl2]

        central.H1h = self.H1h[Nl2]
        central.He1h = self.He1h[Nl2]
        central.He2h = self.He2h[Nl2]

        central.H1h_src = self.H1h_src[Nl2]
        central.He1h_src = self.He1h_src[Nl2]
        central.He2h_src = self.He2h_src[Nl2]

        central.H1h_rec = self.H1h_rec[Nl2]
        central.He1h_rec = self.He1h_rec[Nl2]
        central.He2h_rec = self.He2h_rec[Nl2]


        self.avg.cen = central


        # weighted quantities
        #================================================
        #================================================

        # volume weighted
        self.avg.v_w = self.calc_weighted( self.dz )

        # density weighted 
        self.avg.nH_w = self.calc_weighted( self.nH )
        self.avg.nHe_w = self.calc_weighted( self.nHe )

        # ion weighted 
        self.avg.nH1_w = self.calc_weighted( self.nH * self.xH1 )
        self.avg.nHe1_w = self.calc_weighted( self.nHe * self.xHe1 )
        self.avg.nHe2_w = self.calc_weighted( self.nHe * self.xHe2 )


        # check units
        #------------------------------------------------

        self.z_c.units = 'kpc'











