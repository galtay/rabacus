""" 
Tests agreement between the python and fortran versions of the iliev
sphere solver.
""" 

import os

import pytest
import unittest
import sys
import time
import numpy as np
import pylab as plt
import rabacus as ra

pytestmark = pytest.mark.slow




# The Python and Fortran solvers are independent implementations. They
# agree to ~1e-7 in most layers, but differences up to ~1e-4 (relative)
# appear at ionization fronts. Ionization fractions far below 1 (e.g.
# xHe3 ~ 1e-26 in fully shielded layers) are compared with an absolute
# floor since their relative differences are round-off noise.
X_RTOL = 1.0e-4
X_ATOL = 1.0e-8
T_RTOL = 1.0e-5

i_solve_python = True



def make_r_plot( sphere, fname ):

    # make output dir for png
    #-----------------------------------------------------
    png_dir = 'rabacus_test_suite_png'
    if not os.path.exists( png_dir ):
        os.system( 'mkdir ' + png_dir )

    # make plots
    #-----------------------------------------------------
    plt.figure()
    plt.plot( sphere.r_c, np.log10(sphere.xH1), lw=2.0, label='H1' )
    plt.plot( sphere.r_c, np.log10(sphere.xH2), lw=2.0, label='H2' )
    plt.plot( sphere.r_c, np.log10(sphere.xHe1), ls='--', lw=2.0, label='He1' )
    plt.plot( sphere.r_c, np.log10(sphere.xHe2), ls='--', lw=2.0, label='He2' )
    plt.plot( sphere.r_c, np.log10(sphere.xHe3), ls='--', lw=2.0, label='He3' )
    plt.ylim( -9, 0.3 )
    plt.legend( loc='best' )
    plt.savefig(  png_dir + '/' + fname )
    

def make_t_plot( sphere, fname ):

    # make output dir for png
    #-----------------------------------------------------
    png_dir = 'rabacus_test_suite_png'
    if not os.path.exists( png_dir ):
        os.system( 'mkdir ' + png_dir )

    # make plot
    #-----------------------------------------------------
    plt.figure()
    plt.plot( sphere.r_c, np.log10(sphere.T.magnitude), lw=2.0, label='T' )
    plt.ylim( 3.8, 4.7 )
    plt.legend( loc='best' )
    plt.savefig( png_dir + '/' + fname )




# make polychromatic source HM12 source 
#-------------------------------------------------
q_min = 1.0; q_max = 4.0e2; z = 3.0
hm12 = ra.rad_src.PointSource( q_min, q_max, 'hm12', z=z )
hm12.normalize_Ln( 2.0e53 / ra.u.s )

# create a monochromatic source 
#-------------------------------------------------
q_mono = hm12.grey.E.He2 / hm12.PX.Eth_H1
mono = ra.rad_src.PointSource( q_mono, q_mono, 'monochromatic' )
mono.normalize_Ln( 2.0e53 / ra.u.s )



# describe sphere
#-----------------------------------------------------
Nl = 16
Yp = 0.25
Tsphere = np.ones(Nl) * 1.0e4 * ra.u.K
Rsphere = 200.0 * ra.u.kpc
Redges = np.linspace( 0.0 * ra.u.kpc, Rsphere.rescale("kpc"), Nl+1 )



class TestIspherePythonFortran(unittest.TestCase):


    def check_spheres(self, sphere_fc, sphere_py):

        np.testing.assert_allclose(
            sphere_fc.T.magnitude, sphere_py.T.magnitude, rtol=T_RTOL )

        for x in ['xH1', 'xH2', 'xHe1', 'xHe2', 'xHe3']:
            np.testing.assert_allclose(
                np.asarray( getattr(sphere_fc, x) ),
                np.asarray( getattr(sphere_py, x) ),
                rtol=X_RTOL, atol=X_ATOL, err_msg=x )


    def test_isphere_py_fort_mono_fixed_pce(self):
        """ test isphere_python_fortran mono_fixed_pce """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 1.5e-2 / ra.u.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print('\nsolving fortran sphere ...')
        sphere_fc = ra.f2py.SphereStromgren( 
            Redges, Tsphere, nH, nHe, mono,
            rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        if i_solve_python:
            print('solving python sphere ...') 
            sphere_py = ra.solvers.IlievSphere( 
                Redges, Tsphere, nH, nHe, mono,
                rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        make_r_plot( sphere_fc, 'isphere_mono_fixed_pce_x.png' )
        if i_solve_python:
            self.check_spheres( sphere_fc, sphere_py )        


    def test_isphere_py_fort_poly_fixed_pce(self):
        """ test isphere_python_fortran_poly_fixed_pce """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 2.0e-3 / ra.u.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print('\nsolving fortran sphere ...')
        sphere_fc = ra.f2py.SphereStromgren( 
            Redges, Tsphere, nH, nHe, hm12,
            rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        if i_solve_python:
            print('solving python sphere ...') 
            sphere_py = ra.solvers.IlievSphere( 
                Redges, Tsphere, nH, nHe, hm12,
                rec_meth='fixed', fixed_fcA=0.0, tol=1.0e-7 )

        make_r_plot( sphere_fc, 'isphere_poly_fixed_pce_x.png' )
        if i_solve_python:
            self.check_spheres( sphere_fc, sphere_py )        



    def test_isphere_py_fort_mono_fixed_pcte(self):
        """ test isphere_python_fortran mono_fixed_pcte """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 3.0e-2 / ra.u.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print('\nsolving fortran sphere ...')
        sphere_fc = ra.f2py.SphereStromgren( 
            Redges, Tsphere, nH, nHe, mono, rec_meth='fixed', 
            fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        if i_solve_python:
            print('solving python sphere ...') 
            sphere_py = ra.solvers.IlievSphere( 
                Redges, Tsphere, nH, nHe, mono, rec_meth='fixed', 
                fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        make_r_plot( sphere_fc, 'isphere_mono_fixed_pcte_x.png' )
        make_t_plot( sphere_fc, 'isphere_mono_fixed_pcte_t.png' )
        if i_solve_python:
            self.check_spheres( sphere_fc, sphere_py )        



    def test_isphere_py_fort_poly_fixed_pcte(self):
        """ test isphere_python_fortran_poly_fixed_pcte """ 

        # set densities 
        #-----------------------------------------------------
        nH = np.ones(Nl) * 2.0e-3 / ra.u.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)

        print('\nsolving fortran sphere ...')
        sphere_fc = ra.f2py.SphereStromgren( 
            Redges, Tsphere, nH, nHe, hm12, rec_meth='fixed', 
            fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        if i_solve_python:
            print('solving python sphere ...') 
            sphere_py = ra.solvers.IlievSphere( 
                Redges, Tsphere, nH, nHe, hm12, rec_meth='fixed', 
                fixed_fcA=0.0, find_Teq=True, z=3.0, tol=1.0e-7 )

        make_r_plot( sphere_fc, 'isphere_poly_fixed_pcte_x.png' )
        make_t_plot( sphere_fc, 'isphere_poly_fixed_pcte_t.png' )
        if i_solve_python:
            self.check_spheres( sphere_fc, sphere_py )         



suite = unittest.TestLoader().loadTestsFromTestCase(TestIspherePythonFortran)


        

