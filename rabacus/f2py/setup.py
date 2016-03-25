""" A python script that sets up the f2py shared object (fc.so). """ 

import os
import sys
import numpy as np


# use openmp directives?
#-------------------------------------------------------
use_openmp = True


# choose which fortran compiler to use
#-------------------------------------------------------
#compiler = 'intelem'
compiler = 'gnu95'


# shared object name
#-------------------------------------------------------
so_name = 'rabacus_fc'


# set compiler flags
#-------------------------------------------------------
if compiler == 'gnu95':
    omp_lib = '-lgomp'
    omp_flag = '-fopenmp'
    f90_flags = '"-O3 -fbounds-check -mtune=native"'


elif compiler == 'intelem':
    omp_lib = '-liomp5'
    omp_flag = '-openmp'
    f90_flags = '"-xHost -O3 -ipo -funroll-loops -heap-arrays -mcmodel=medium"'
#    f90_flags = '"-g -O3 -heap-arrays -mcmodel=medium -check bounds"'

if use_openmp:
    f90_flags = f90_flags[:-1] + ' ' + omp_flag + '"'


# call f2py and compile the fortan code 
#======================================================================
cmnd = '\\rm *.mod'
os.system(cmnd)

cmnd = '\\rm *.so'
os.system(cmnd)

cmnd = '\\rm *.pyf'
os.system(cmnd)



file_list = [
    'types.f90', 
    'utils.f90',
    'zhang_jin.f90',
    'zhang_jin_f.f90',
    'slatec.f90',
    'special_functions.f90',
    'legendre_polynomial.f90',
    'm_mrgrnk.f90',  
    'physical_constants.f90', 
    'geometry.f90', 
    'hui_gnedin_97.f90', 
    'hm12.f90', 
    'chem_cool_rates.f90', 
    'ion_solver.f90', 
    'verner_96.f90', 
    'photo_xsections.f90', 
    'spectra.f90', 
    'source_point.f90', 
    'source_plane.f90', 
    'source_background.f90', 
    'sphere_base.f90', 
    'sphere_stromgren.f90', 
    'sphere_bgnd.f90', 
    'slab_base.f90',
    'slab_plane.f90',
    'slab_bgnd.f90',
    ]



# scan fortran files and make signature file
# this is just a convenience so we can see how
# f2py interprets the fortran code
#================================================

cmnd = 'f2py -h ' + so_name + '.pyf '
for f in file_list:
    cmnd = cmnd + f + ' ' 

print cmnd
os.system(cmnd)


# compile fortran files into a shared object
#================================================

cmnd = 'f2py --noopt -m ' + so_name + ' ' 
cmnd = cmnd + ' --fcompiler=' + compiler
cmnd = cmnd + ' -c --f90flags=' + f90_flags + ' '

if use_openmp:
    cmnd = cmnd + omp_lib + ' '

for f in file_list:
    cmnd = cmnd + f + ' ' 

print cmnd
os.system(cmnd)


