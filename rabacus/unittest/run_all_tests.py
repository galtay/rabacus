"""
Runs all tests. 
"""

import unittest
import test_hm12 
import test_atomic_rates
import test_ion_solver
import test_source_point
import test_source_plane
import test_source_background
import test_slab_analytic

import test_slab1_python_fortran
import test_slab2_python_fortran
import test_isphere_python_fortran


def run_all_tests():
    
    tests = unittest.TestSuite( 
        [test_hm12.suite,
         test_atomic_rates.suite,
         test_ion_solver.suite,
         test_source_point.suite,
         test_source_plane.suite,
         test_source_background.suite,
         test_slab_analytic.suite,
         test_slab1_python_fortran.suite,
         test_slab2_python_fortran.suite,
         test_isphere_python_fortran.suite,] 
        )

    unittest.TextTestRunner(verbosity=2).run(tests)
