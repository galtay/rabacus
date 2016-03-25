""" Tests the slab solver against an analytic monochromatic solution. """ 

import sys
import time
import numpy as np
import pylab as plt
import rabacus as ra

import unittest


class TestSlabAnalytic(unittest.TestCase):

    def test_slab_analytic(self):
        """ test slab_analytic """ 

        TOL = 1.0e-3

        # create a polychromatic HM12 source
        #-----------------------------------------------------
        q_min = 1.0; q_max = 4.0e2; z = 3.0
        hm12 = ra.rad_src.BackgroundSource( q_min, q_max, 'hm12', z=z )


        # Calculate an analytic monochromatic solution
        #-----------------------------------------------------
        nH = 1.5e-3 / ra.U.cm**3
        T = np.ones(1) * 1.0e4 * ra.U.K
        H1i = np.ones(1) * hm12.thin.H1i
        y = 0.0

        fcA = 1.0
        ana = ra.solvers.AnalyticSlab( nH, T, H1i, y, fcA ) 
        ana.set_E( hm12.grey.E.H1 )


        # use grey energy to create a monochromatic source
        #-----------------------------------------------------
        q_mono = hm12.grey.E.H1 / hm12.PX.Eth_H1
        mono = ra.rad_src.PlaneSource( q_mono, q_mono, 'monochromatic' )
        mono.normalize_H1i( hm12.thin.H1i )


        # solve a slab
        #-----------------------------------------------------
        Nl = 4096
        nH = np.ones(Nl) * nH
        nHe = nH * 1.0e-10
        Tslab = np.ones(Nl) * T

        Lslab = ana.Lslab * 1.1
        Ledges = np.linspace( 0.0 * ra.U.kpc, Lslab, Nl+1 )

        slab = ra.f2py.SolveSlab( Ledges, Tslab, nH, nHe, mono,
                                  rec_meth='fixed', fixed_fcA=1.0 )
        
        # make interpolating function for slab solution
        #-----------------------------------------------------
        ana.L.units = 'kpc'
        slab.z_c.units = 'kpc'

        fit = ra.scipy.interp1d( slab.z_c, np.log10(slab.xH1) )

        # evaluate fit at analytic solution points
        # assert that the errors between the analytic
        # and numerical solution are less than 1%
        #-----------------------------------------------------
        c1 = ana.L > slab.z_c.min()
        c2 = ana.L < slab.z_c.max()

        indx = np.where( c1 & c2 )
        L = ana.L[indx]
        xH1_ana = ana.xH1[indx]
        xH1_num = 10**fit( L )

        err = np.abs( xH1_num - xH1_ana ) / xH1_ana
        ok = not np.any( err > TOL )
        self.assertTrue( ok )


            


suite = unittest.TestLoader().loadTestsFromTestCase(TestSlabAnalytic)
        

