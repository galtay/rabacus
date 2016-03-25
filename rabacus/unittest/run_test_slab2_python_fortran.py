""" 
runs python fortran consistancy tests for slab2
""" 

import os
import unittest
import test_slab2_python_fortran


__all__ = ['run_test_slab2_python_fortran']


def run_test_slab2_python_fortran():


    # run tests
    #-----------------------------------------------------
    tests = unittest.TestSuite(
        [test_slab2_python_fortran.suite]
        )

    unittest.TextTestRunner(verbosity=2).run(tests)
