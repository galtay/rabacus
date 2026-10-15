""" Wrapper to special function fortran module. """

from . import rabacus_fc



__all__ = ['E1xa', 'E1xb', 'Enxa', 'Enxb']


def E1xa(x):
    v = rabacus_fc.zhang_jin.e1xa( x )
    return v

def E1xb(x):
    v = rabacus_fc.zhang_jin.e1xb( x )
    return v

def Enxa( n, x ):
    v = rabacus_fc.zhang_jin.enxa( n, x )
    return v[n]

def Enxb( n, x ):
    v = rabacus_fc.zhang_jin.enxb( n, x )
    return v[n]
