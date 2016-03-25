""" Test the point radiation source. """ 

import sys
import pylab as plt
import numpy as np
import rabacus as ra

import unittest



class TestSourcePoint(unittest.TestCase):

    def test_creation(self):
        """ create each spectral type """ 
        q_min = 1.0; q_max = 4.0e2; st='hm12'; z=3.0
        hm12 = ra.rad_src.PointSource( q_min, q_max, st, z=z )

        q_min = 1.0; q_max = 5.0; st='thermal'; T_eff=1.0e5 * ra.U.K
        t1e5 = ra.rad_src.PointSource( q_min, q_max, st, T_eff=T_eff )

        q_min = 1.0; q_max = 5.0; st='powerlaw'; alpha=-1.0
        p1p0 = ra.rad_src.PointSource( q_min, q_max, st, alpha=alpha )

        q_min = 5.0; q_max = 5.0; st='monochromatic'
        mono = ra.rad_src.PointSource( q_min, q_max, st )

        self.assertTrue( True )


    def test_functions(self):
        """ call all functions """ 

        tauH1_th = tauHe1_th = tauHe2_th = 1.0e1
        r = 1.0 * ra.U.kpc

        # for polychromatic
        #----------------------------------------------------------
        q_min = 1.0; q_max = 4.0e2; st='hm12'; z=3.0
        hm12 = ra.rad_src.PointSource( q_min, q_max, st, z=z )

        H1i = hm12.shld_H1i( r, tauH1_th, tauHe1_th, tauHe2_th )
        H1h = hm12.shld_H1h( r, tauH1_th, tauHe1_th, tauHe2_th )
        He1i = hm12.shld_He1i( r, tauH1_th, tauHe1_th, tauHe2_th )
        He1h = hm12.shld_He1h( r, tauH1_th, tauHe1_th, tauHe2_th )
        He2i = hm12.shld_He2i( r, tauH1_th, tauHe1_th, tauHe2_th )
        He2h = hm12.shld_He2h( r, tauH1_th, tauHe1_th, tauHe2_th )

        hm12.normalize_Ln( 5.0e48 / ra.U.s )
        hm12.normalize_Lu( 1.0e38 * ra.U.erg / ra.U.s )

        # for monochromatic
        #----------------------------------------------------------
        q_min = 1.0; q_max = 1.0; st='monochromatic'
        mono = ra.rad_src.PointSource( q_min, q_max, st )

        H1i = mono.shld_H1i( r, tauH1_th, tauHe1_th, tauHe2_th )
        H1h = mono.shld_H1h( r, tauH1_th, tauHe1_th, tauHe2_th )
        He1i = mono.shld_He1i( r, tauH1_th, tauHe1_th, tauHe2_th )
        He1h = mono.shld_He1h( r, tauH1_th, tauHe1_th, tauHe2_th )
        He2i = mono.shld_He2i( r, tauH1_th, tauHe1_th, tauHe2_th )
        He2h = mono.shld_He2h( r, tauH1_th, tauHe1_th, tauHe2_th )

        mono.normalize_Ln( 5.0e48 / ra.U.s )
        mono.normalize_Lu( 1.0e38 * ra.U.erg / ra.U.s )

        self.assertTrue( True )


suite = unittest.TestLoader().loadTestsFromTestCase(TestSourcePoint)





