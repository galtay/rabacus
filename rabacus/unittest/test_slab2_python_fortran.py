""" 
Tests agreement between the python and fortran versions of the slab2 solver.
""" 

import os
import sys
import time
import numpy as np
import pylab as plt
import rabacus as ra

import unittest


X_TOL = 1.0e-5
T_TOL = 1.0e-5

i_solve_python = True


def make_x_plot( slab, fname ):

    # make output dir for png
    #-----------------------------------------------------
    png_dir = 'rabacus_test_suite_png'
    if not os.path.exists( png_dir ):
        os.system( 'mkdir ' + png_dir )

    # make plots
    #-----------------------------------------------------
    plt.figure()
    plt.plot( slab.z_c, np.log10(slab.xH1), lw=2.0, label='H1' )
    plt.plot( slab.z_c, np.log10(slab.xH2), lw=2.0, label='H2' )
    plt.plot( slab.z_c, np.log10(slab.xHe1), ls='--', lw=2.0, label='He1' )
    plt.plot( slab.z_c, np.log10(slab.xHe2), ls='--', lw=2.0, label='He2' )
    plt.plot( slab.z_c, np.log10(slab.xHe3), ls='--', lw=2.0, label='He3' )
    plt.ylim( -9, 0.3 )
    plt.legend( loc='best' )
    plt.savefig( png_dir + '/' + fname )
    

def make_t_plot( slab, fname ):

    # make output dir for png
    #-----------------------------------------------------
    png_dir = 'rabacus_test_suite_png'
    if not os.path.exists( png_dir ):
        os.system( 'mkdir ' + png_dir )

    # make plots
    #-----------------------------------------------------
    plt.figure()
    plt.plot( slab.z_c, np.log10(slab.T.magnitude), lw=2.0, label='T' )
    plt.ylim( 3.8, 4.7 )
    plt.legend( loc='best' )
    plt.savefig( png_dir + '/' + fname )




# make polychromatic source HM12 source 
#-------------------------------------------------
q_min = 1.0; q_max = 4.0e2; z = 3.0
hm12 = ra.rad_src.PlaneSource( q_min, q_max, 'hm12', z=z )

# create a monochromatic source 
#-------------------------------------------------
q_mono = hm12.grey.E.He2 / hm12.PX.Eth_H1
mono = ra.rad_src.PlaneSource( q_mono, q_mono, 'monochromatic' )
mono.normalize_H1i( hm12.thin.H1i )


# describe slab
#-----------------------------------------------------
Nl = 16
Yp = 0.25
Tslab = np.ones(Nl) * 1.0e4 * ra.U.K
Lslab = 200.0 * ra.U.kpc
Ledges = np.linspace( 0.0 * ra.U.kpc, Lslab, Nl+1 )



class TestSlab2PythonFortran(unittest.TestCase):


    def check_slabs(self, slab_fc, slab_py):

        err = np.abs( 1.0 - slab_fc.T / slab_py.T )
        print 'T max err: ', err.max()
        ok = not np.any( err > T_TOL )
        self.assertTrue( ok )

        err = np.abs( 1.0 - slab_fc.xH1 / slab_py.xH1 )
        print 'xH1 max err: ', err.max()
        ok = not np.any( err > X_TOL )
        self.assertTrue( ok )

        err = np.abs( 1.0 - slab_fc.xH2 / slab_py.xH2 )
        print 'xH2 max err: ', err.max()
        ok = not np.any( err > X_TOL )
        self.assertTrue( ok )

        err = np.abs( 1.0 - slab_fc.xHe1 / slab_py.xHe1 )
        print 'xHe1 max err: ', err.max()
        ok = not np.any( err > X_TOL )
        self.assertTrue( ok )

        err = np.abs( 1.0 - slab_fc.xHe2 / slab_py.xHe2 )
        print 'xHe2 max err: ', err.max()
        ok = not np.any( err > X_TOL )
        self.assertTrue( ok )

        err = np.abs( 1.0 - slab_fc.xHe3 / slab_py.xHe3 )
        print 'xHe3 max err: ', err.max()
        ok = not np.any( err > X_TOL )
        self.assertTrue( ok )




    def test_slab_2_py_fort_mono_fixed_pce(self):
        """ test slab_2_python_fortran mono_fixed_pce """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 1.5e-2 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print '\nsolving fortran slab ...'
        slab_fc = ra.f2py.SolveSlab2( 
            Ledges, Tslab, nH, nHe, mono,
            rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        if i_solve_python:
            print 'solving python slab ...' 
            slab_py = ra.solvers.SolveSlab2( 
                Ledges, Tslab, nH, nHe, mono,
                rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        make_x_plot( slab_fc, 'slab_2_mono_fixed_pce_x.png' )
        if i_solve_python:
            self.check_slabs( slab_fc, slab_py )        



    def test_slab_2_py_fort_poly_fixed_pce(self):
        """ test slab_2_python_fortran_poly_fixed_pce """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 2.0e-3 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print '\nsolving fortran slab ...'
        slab_fc = ra.f2py.SolveSlab2( 
            Ledges, Tslab, nH, nHe, hm12,
            rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        if i_solve_python:
            print 'solving python slab ...' 
            slab_py = ra.solvers.SolveSlab2( 
                Ledges, Tslab, nH, nHe, hm12,
                rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        make_x_plot( slab_fc, 'slab_2_poly_fixed_pce_x.png' )
        if i_solve_python:
            self.check_slabs( slab_fc, slab_py )        



    def test_slab_2_py_fort_mono_fixed_pcte(self):
        """ test slab_2_python_fortran mono_fixed_pcte """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 3.0e-2 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print '\nsolving fortran slab ...'
        slab_fc = ra.f2py.SolveSlab2( 
            Ledges, Tslab, nH, nHe, mono, rec_meth='fixed', 
            fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        if i_solve_python:
            print 'solving python slab ...' 
            slab_py = ra.solvers.SolveSlab2( 
                Ledges, Tslab, nH, nHe, mono, rec_meth='fixed', 
                fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        make_x_plot( slab_fc, 'slab_2_mono_fixed_pcte_x.png' )
        make_t_plot( slab_fc, 'slab_2_mono_fixed_pcte_t.png' )
        if i_solve_python:
            self.check_slabs( slab_fc, slab_py )        



    def test_slab_2_py_fort_poly_fixed_pcte(self):
        """ test slab_2_python_fortran_poly_fixed_pcte """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 2.0e-3 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print '\nsolving fortran slab ...'
        slab_fc = ra.f2py.SolveSlab2( 
            Ledges, Tslab, nH, nHe, hm12, rec_meth='fixed', 
            fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        if i_solve_python:
            print 'solving python slab ...' 
            slab_py = ra.solvers.SolveSlab2( 
                Ledges, Tslab, nH, nHe, hm12, rec_meth='fixed', 
                fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        make_x_plot( slab_fc, 'slab_2_poly_fixed_pcte_x.png' )
        make_t_plot( slab_fc, 'slab_2_poly_fixed_pcte_t.png' )
        if i_solve_python:
            self.check_slabs( slab_fc, slab_py )         


    def test_slab_2_py_fort_poly_thresh_pcte(self):
        """ test slab_2_python_fortran_poly_thresh_pcte """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 2.0e-3 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print '\nsolving fortran slab ...'
        slab_fc = ra.f2py.SolveSlab2( 
            Ledges, Tslab, nH, nHe, hm12, rec_meth='thresh', 
            find_Teq=True, z=3.0, tol=1.0e-8 )

        if i_solve_python:
            print 'solving python slab ...' 
            slab_py = ra.solvers.SolveSlab2( 
                Ledges, Tslab, nH, nHe, hm12, rec_meth='thresh', 
                find_Teq=True, z=3.0, tol=1.0e-8 )

        make_x_plot( slab_fc, 'slab_2_poly_thresh_pcte_x.png' )
        make_t_plot( slab_fc, 'slab_2_poly_thresh_pcte_t.png' )
        if i_solve_python:
            self.check_slabs( slab_fc, slab_py )         





suite = unittest.TestLoader().loadTestsFromTestCase(TestSlab2PythonFortran)

















