"""
Runs all tests. 
"""

import unittest
from . import test_hm12 
from . import test_atomic_rates
from . import test_ion_solver
from . import test_source_point
from . import test_source_plane
from . import test_source_background
from . import test_slab_analytic

from . import test_slab1_python_fortran
from . import test_slab2_python_fortran
from . import test_isphere_python_fortran


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
