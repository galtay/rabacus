""" 
runs short time tests.  should complete in less than a few minutes. 
""" 

import unittest
import test_hm12 
import test_atomic_rates
import test_ion_solver
import test_source_point
import test_source_plane
import test_source_background
import test_slab_analytic


__all__ = ['run_test_short']


def run_test_short():

    tests = unittest.TestSuite( 
        [test_hm12.suite,
         test_atomic_rates.suite,
         test_ion_solver.suite,
         test_source_point.suite,
         test_source_plane.suite,
         test_source_background.suite,
         test_slab_analytic.suite] 
        )
    
    unittest.TextTestRunner(verbosity=2).run(tests)



