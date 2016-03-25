import numpy as np
import collections


__all__ = ['NeedUnitsError', 'InputError', 'isiterable', 'trap', 
           'units_string']


class Error(Exception):
    """Base class for all exceptions."""
    pass

class BaseClassError(Error):
    """ Exception raised for instantiating base class. """
    def __init__(self):
        message = ' Spectrum base class cannot be directly instantiated.'
        Exception.__init__( self, message )

class NeedUnitsError(Error):
    """Exception raised for providing a quantity without units. """
    def __init__(self, message ):
        message = message + \
            ' Units can be found in the U attribute of the idlehands class '
        Exception.__init__( self, message )


class InputError(Error):
    """Exception raised for bad function input.  """
    def __init__(self, message ):
        Exception.__init__( self, message )





def isiterable(obj):
    """ Returns `True` if the given object is iterable.  Taken from the 
    Astropy souce code. """

    # Numpy arrays are in collections.Iterable no matter what, but if you
    # attempt to iterate over a 0-d array, it throws a TypeError.
    if isinstance(obj, np.ndarray) and len(obj.shape) == 0:
        return False

    if isinstance(obj, collections.Iterable):
        return True

    try:
        iter(obj)
        return True
    except TypeError:
        return False


def trap( x, y ):
    """ Trapezoidal integration rule that preseves units. """ 

    dx = x[1:] - x[0:-1]
    I = 0.5 * ( y[0:-1] + y[1:] ) * dx
    I = np.sum(I)
    return I


def units_string( q ):
    tmp = str( q.units )
    txt = tmp.split( ' ' )[1:]
    return txt 
