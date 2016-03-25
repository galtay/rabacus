""" Test the HM12 UV background model. 
http://adsabs.harvard.edu/abs/2012ApJ...746..125H
""" 

import unittest
import numpy as np
import rabacus as ra



H1i_z0p00 = 0.228e-13 / ra.U.s
H1h_z0p00 = 0.889e-13 * ra.U.eV / ra.U.s

He2i_z15p10 = 0.192e-22 / ra.U.s
He2h_z15p10 = 0.154e-19 * ra.U.eV / ra.U.s


class TestHM12(unittest.TestCase):

    def test_photorates_table(self):
        """ verify I/O of photorates table """ 

        TOL = 1.0e-6
        prt = ra.uv_bgnd.HM12_Photorates_Table()

        z = 0.00
        err = np.abs( ( prt.H1i(z) - H1i_z0p00 ) / H1i_z0p00 )
        self.assertTrue( err < TOL )
        err = np.abs( ( prt.H1h(z) - H1h_z0p00 ) / H1h_z0p00 )
        self.assertTrue( err < TOL )

        z = 15.10
        err = np.abs( ( prt.He2i(z) - He2i_z15p10 ) / He2i_z15p10 )
        self.assertTrue( err < TOL )
        err = np.abs( ( prt.He2h(z) - He2h_z15p10 ) / He2h_z15p10 )
        self.assertTrue( err < TOL )


    def test_spectrum_integrals(self):
        """ compare integrated values to tabulated values """ 

        # instantiate objects
        #---------------------------------------------------
        prt = ra.uv_bgnd.HM12_Photorates_Table()
        uvt = ra.uv_bgnd.HM12_UVB_Table()
        px = ra.atomic.PhotoXsections()

        # loop over all redshifts in the photorates table
        # and check that integrating produces the value in 
        # the table
        #---------------------------------------------------
        E_min = 13.6 
        four_pi_sr = 4 * np.pi * ra.U.sr

        for z in prt.z:

            if z < 8.75:
                E_max = E_min * 400
                Npts = 2048
                log_E = np.linspace( np.log10(E_min), np.log10(E_max), Npts )
                E = 10**log_E * ra.U.eV
                TOL = 5.0e-2
            elif z >= 8.75: 
                E_max = E_min * 1000
                Npts = 8192
                log_E = np.linspace( np.log10(E_min), np.log10(E_max), Npts )
                E = 10**log_E * ra.U.eV
                TOL = 1.5e-1

            # create spectrum
            #---------------------------------------------------
            Inu = uvt.return_spectrum_E( z, E )
            
            # integrate and compare to table
            #---------------------------------------------------

            # energy flux per unit frequency
            dFu_dnu = four_pi_sr * Inu

            # photon flux per unit energy 
            dFn_dE = dFu_dnu / ( ra.PC.h * E )


            # test photo-ionization rates
            #---------------------------------------
            dH1i_dE = dFn_dE * px.sigma_H1(E)
            H1i = ra.utils.trap( E, dH1i_dE ) 
            err = np.abs( ( H1i - prt.H1i(z) ) / prt.H1i(z) ).simplified
            self.assertTrue( err < TOL )

            dHe1i_dE = dFn_dE * px.sigma_He1(E)
            He1i = ra.utils.trap( E, dHe1i_dE ) 
            err = np.abs( ( He1i - prt.He1i(z) ) / prt.He1i(z) ).simplified
            self.assertTrue( err < TOL )

            dHe2i_dE = dFn_dE * px.sigma_He2(E)
            He2i = ra.utils.trap( E, dHe2i_dE ) 
            err = np.abs( ( He2i - prt.He2i(z) ) / prt.He2i(z) ).simplified
            self.assertTrue( err < TOL )


            # test photo-heating rates 
            #---------------------------------------
            dH1h_dE = dH1i_dE * ( E - px.Eth_H1 )
            H1h = ra.utils.trap( E, dH1h_dE ) 
            err = np.abs( ( H1h - prt.H1h(z) ) / prt.H1h(z) ).simplified
            self.assertTrue( err < TOL )

            dHe1h_dE = dHe1i_dE * ( E - px.Eth_He1 )
            He1h = ra.utils.trap( E, dHe1h_dE ) 
            err = np.abs( ( He1h - prt.He1h(z) ) / prt.He1h(z) ).simplified
            self.assertTrue( err < TOL )

            dHe2h_dE = dHe2i_dE * ( E - px.Eth_He2 )
            He2h = ra.utils.trap( E, dHe2h_dE ) 
            err = np.abs( ( He2h - prt.He2h(z) ) / prt.He2h(z) ).simplified
            self.assertTrue( err < TOL )




suite = unittest.TestLoader().loadTestsFromTestCase(TestHM12)


