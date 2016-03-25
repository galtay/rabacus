""" Test the python and fortran atomic rates modules """ 


import sys
import time
import pylab as plt
import numpy as np
import rabacus as ra

import unittest 

TOL = 1.0e-10

class TestAtomicRates(unittest.TestCase):

    def test_chem(self):
        """ test kchem """ 
        N = 128
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K
        fcA_H2 = np.ones(N) * 1.0
        fcA_He2 = np.ones(N) * 1.0
        fcA_He3 = np.ones(N) * 1.0

        kchem_py = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )
        kchem_fc = ra.f2py.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )

        #---------------------------------------------------------
        err = np.abs( (kchem_py.reH2 / kchem_fc.reH2) - 1.0 )
        print '\n reH2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kchem_py.reHe2 / kchem_fc.reHe2) - 1.0 )
        print ' reHe2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kchem_py.reHe3 / kchem_fc.reHe3) - 1.0 )
        print ' reHe3 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        #---------------------------------------------------------
        err = np.abs( (kchem_py.ciH1 / kchem_fc.ciH1) - 1.0 )
        print ' ciH1 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kchem_py.ciHe1 / kchem_fc.ciHe1) - 1.0 )
        print ' ciHe1 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kchem_py.ciHe2 / kchem_fc.ciHe2) - 1.0 )
        print ' ciHe2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )



    def test_cool(self):
        """ test kcool """ 
        N = 128
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K
        fcA_H2 = np.ones(N) * 1.0
        fcA_He2 = np.ones(N) * 1.0
        fcA_He3 = np.ones(N) * 1.0

        kcool_py = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3 )
        kcool_fc = ra.f2py.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3 )

        #---------------------------------------------------------
        err = np.abs( (kcool_py.recH2 / kcool_fc.recH2) - 1.0 )
        print '\n recH2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.recHe2 / kcool_fc.recHe2) - 1.0 )
        print ' recHe2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.recHe3 / kcool_fc.recHe3) - 1.0 )
        print ' recHe3 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        #---------------------------------------------------------
        err = np.abs( (kcool_py.cicH1 / kcool_fc.cicH1) - 1.0 )
        print ' cicH1 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.cicHe1 / kcool_fc.cicHe1) - 1.0 )
        print ' cicHe1 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.cicHe2 / kcool_fc.cicHe2) - 1.0 )
        print ' cicHe2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        #---------------------------------------------------------
        err = np.abs( (kcool_py.cecH1 / kcool_fc.cecH1) - 1.0 )
        print ' cecH1 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.cecHe2 / kcool_fc.cecHe2) - 1.0 )
        print ' cecHe2 max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        #---------------------------------------------------------
        err = np.abs( (kcool_py.bremss / kcool_fc.bremss) - 1.0 )
        print ' bremss max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )

        err = np.abs( (kcool_py.compton / kcool_fc.compton) - 1.0 )
        print ' compton max err: ', err.max()
        ok = not np.any( err > TOL )
        self.assertTrue( ok )



suite = unittest.TestLoader().loadTestsFromTestCase(TestAtomicRates)








    



