""" Simple convenience functions to access the h5py library.  """

import sys
import h5py
import numpy as np


__all__ = ['ra', 'raa', 'rd', 'wd', 'wa', 'cg']


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class OverwriteError(Error):
    """Exception raised for attempting to overwrite data. """
    def __init__(self, message ):
        Exception.__init__( self, message )
       


# Reading methods
#=========================================================================

def ra( fname, path, name ):
    """ Read Attribute. Return a single attribute called <name> associated
    with a group or dataset located at <path> in file <fname>.
    e.g. ra( fname, '/PartType0/Coordinates', 'h-scale-exponent' ). """

    attr = None
    with h5py.File( fname, 'r' ) as h5f:
        attr = h5f[path].attrs[name]
    return attr


def raa( fname, path ):
    """ Read All Attributes. Return a dictionary of all attributes associated
    with a group or dataset located at <path> in file <fname>.
    e.g. ra( fname, '/PartType0/Coordinates' ). """

    attrs = {}
    with h5py.File( fname, 'r' ) as h5f:
        attr_names = h5f[path].attrs.keys()
        for name in attr_names:
            attrs[name] = h5f[path].attrs[name]
    return attrs


def rd( fname, path, dtype=None ):
    """ Read Data. Return a dataset located at <path> in file <fname> as
    a numpy array.
    e.g. rd( fname, '/PartType0/Coordinates' ). """


    data = None    
    with h5py.File( fname, 'r' ) as h5f:
        ds = h5f[path]
        if dtype == None:
            dtype = ds.dtype
        data = np.zeros( ds.shape, dtype=dtype )
        data = ds.value
    return data



# Writing methods
#=========================================================================

def wd( fname, path, name, buf, attrs=None, overwrite=False ):
    """ Write Data. Write a dataset stored in <buf> to hdf5 file <fname>
    at location <path>/<name>.  Optionally a dictionary of attributes
    can be provided which will be written with the dataset.
    e.g. wd( fname, '/PartType0', 'Coordinates', pos [,attrs] ). """

    with h5py.File( fname, 'a' ) as h5f:

        ds_path = path + '/' + name
        path_exists = ds_path in h5f
       
        # if we are trying to overwrite a dataset w/o setting overwrite=True
        #--------------------------------------------------------------------
        if path_exists and not overwrite:            
            msg = '\n attempting to overwrite a dataset witout setting ' +  \
                   'overwrite=True \n ' + \
                   'file name: ' + fname + ' \n ' + \
                   'path: ' + path + ' \n ' + \
                   'name: ' +  name
            raise OverwriteError( msg )
            return

        # otherwise delete ds if needed and write new one
        #--------------------------------------------------------------------
        else:

            if path_exists:
                del h5f[ds_path]
       
            h5f.create_dataset( ds_path, data=buf )

            if attrs:
                for k,v in attrs.items():        
                    h5f[ds_path].attrs[k] = v



def wa( fname, path, name, buf ):
    """ Write Attribute.  Write a single attribute stored in <buf> called
    <name> associated with a group or dataset located at <path> in file <fname>.
    e.g. wa( fname, '/PartType0/Coordinates', 'h-scale-exponent', -1.0 ). """

    with h5py.File( fname, 'a' ) as h5f:
        h5f[path].attrs[name] = buf


def cg( fname, path ):
    """ Create Group.  Creates a group at <path> in file <fname>.
    e.g. cg( fname, '/PartType0' ). """

    with h5py.File( fname, 'a' ) as h5f:
        h5f.require_group( path )
