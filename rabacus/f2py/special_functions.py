""" Wrapper to special function fortran module. """ 

import rabacus_fc
import numpy as np



__all__ = ['E1xa', 'E1xb', 'E1z', 'Eix', 'Enxa', 'Enxb']


def E1xa(x):     
    v = rabacus_fc.special_functions.e1xa( x )
    return v

def E1xb(x):
    v = rabacus_fc.special_functions.e1xb( x )
    return v

def E1z(z):
    v = rabacus_fc.special_functions.e1z( z )
    return v

def Eix(x):
    v = rabacus_fc.special_functions.eix( x )
    return v

def Enxa( n, x ):
    v = rabacus_fc.special_functions.enxa( n, x )
    return v[n]

def Enxb( n, x ):
    v = rabacus_fc.special_functions.enxb( n, x )
    return v[n]

