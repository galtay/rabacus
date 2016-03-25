""" Test the python ion_solver module """ 


import sys
import time
import pylab as plt
import numpy as np
import rabacus as ra

import unittest 

TOL = 1.0e-8

H = ra.atomic.hydrogen.Hydrogen()


class TestIonSolver(unittest.TestCase):

    def test_ce_py(self):
        """ test solve_ce_py """ 
        N = 128
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K 
        fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0
        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )
        x_ce = ra.solvers.Solve_CE( kchem )

        nH = 1.0e-3 / ra.U.cm**3
        nHe = nH * 1.0e-1
        y = ( x_ce.He2 + 2 * x_ce.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1(nH, y, kchem)

        err = np.abs( (x_ce.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )


    def test_ce_f2py(self):
        """ test solve_ce_f2py """ 
        N = 128
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K 
        fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0
        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )
        x_ce = ra.f2py.Solve_CE( kchem )

        nH = 1.0e-3 / ra.U.cm**3
        nHe = nH * 1.0e-1
        y = ( x_ce.He2 + 2 * x_ce.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1(nH, y, kchem)

        err = np.abs( (x_ce.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )


    def test_pce_py(self):
        """ test solve_pce_py """ 

        N = 128; Yp = 0.24; Nrho = NT = N; z=3.0

        T = np.linspace( 1.0e4, 1.0e5, NT ) * ra.U.K
        fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0

        hmr = ra.uv_bgnd.HM12_Photorates_Table()
        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)

        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        nH = np.linspace( 1.0e-4, 1.0e-3, Nrho ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)

        x_pce_1D = ra.solvers.Solve_PCE( nH, nHe, kchem )

        y = ( x_pce_1D.He2 + 2 * x_pce_1D.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1( nH, y, kchem )

        err = np.abs( (x_pce_1D.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )        


    def test_pce_f2py(self):
        """ test solve_pce_f2py """ 

        N = 128; Yp = 0.24; Nrho = NT = N; z=3.0

        T = np.linspace( 1.0e4, 1.0e5, NT ) * ra.U.K
        fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0

        hmr = ra.uv_bgnd.HM12_Photorates_Table()
        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)

        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        nH = np.linspace( 1.0e-4, 1.0e-3, Nrho ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)

        x_pce_1D = ra.f2py.Solve_PCE( nH, nHe, kchem )

        y = ( x_pce_1D.He2 + 2 * x_pce_1D.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1( nH, y, kchem )

        err = np.abs( (x_pce_1D.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )        



    def test_pcte_py(self):
        """ test solve_pcte_py """ 

        N = 128; Yp = 0.24; z=3.0
        fcA_H2 = fcA_He2 = fcA_He3 = 0.0

        nH = np.linspace( 1.0e-4, 1.0e-3, N ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K

        hmr = ra.uv_bgnd.HM12_Photorates_Table()

        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)
        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        H1h = np.ones(N) * hmr.H1h(z)
        He1h = np.ones(N) * hmr.He1h(z)
        He2h = np.ones(N) * hmr.He2h(z)
        kcool = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1h=H1h, He1h=He1h, He2h=He2h )

        x_pcte = ra.solvers.Solve_PCTE( nH, nHe, kchem, kcool, z )

        y = ( x_pcte.He2 + 2 * x_pcte.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1( nH, y, x_pcte.kchem )

        err = np.abs( (x_pcte.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )        


    def test_pcte_f2py(self):
        """ test solve_pcte_f2py""" 


        N = 128; Yp = 0.24; z=3.0
        fcA_H2 = fcA_He2 = fcA_He3 = 0.0

        nH = np.linspace( 1.0e-4, 1.0e-3, N ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K

        hmr = ra.uv_bgnd.HM12_Photorates_Table()

        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)
        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        H1h = np.ones(N) * hmr.H1h(z)
        He1h = np.ones(N) * hmr.He1h(z)
        He2h = np.ones(N) * hmr.He2h(z)
        kcool = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1h=H1h, He1h=He1h, He2h=He2h )

        x_pcte = ra.f2py.Solve_PCTE( nH, nHe, kchem, kcool, z )

        y = ( x_pcte.He2 + 2 * x_pcte.He3 ) * nHe / nH
        xH1 = H.analytic_soltn_xH1( nH, y, x_pcte.kchem )

        err = np.abs( (x_pcte.H1 - xH1) / xH1 )
        ok = not np.any( err > TOL )

        self.assertTrue( ok )        




suite = unittest.TestLoader().loadTestsFromTestCase(TestIonSolver)








    



