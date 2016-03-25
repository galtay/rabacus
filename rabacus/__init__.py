"""
""" 

import os as _os
radir = _os.path.dirname(_os.path.realpath(__file__))
_os.putenv( 'RABACUS_DIR', radir )



from constants.units import Units
from constants.units import CosmoUnits
from constants.physical import PhysicalConstants

pc = PhysicalConstants()
u = Units()

del PhysicalConstants
del Units


from atomic.chemistry import ChemistryRates
from atomic.cooling import CoolingRates
from atomic.photo_xsection import PhotoXsections
from atomic.hydrogen import Hydrogen
from atomic.helium import Helium
from atomic.sd93 import SD93

from cosmology.jeans import Jeans
from cosmology.general import Cosmology
from cosmology.mass_function.mass_function import MassFunction
from cosmology.transfer_functions.bbks86 import TransferBBKS
from cosmology.parameters.planck.load import PlanckParameters

from uv_bgnd import HM12_Photorates_Table
from uv_bgnd import HM12_UVB_Table

from f2py import Solve_CE
from f2py import Solve_PCE
from f2py import Solve_PCTE

from f2py import SlabPln
from f2py import Slab2Pln
from f2py import Slab2Bgnd

from f2py import SphereStromgren
from f2py import SphereBgnd

from rad_src import PointSource
from rad_src import PlaneSource
from rad_src import BackgroundSource

from solvers.calc_slab_analytic import AnalyticSlab

import hdf5 

#import constants 


planck13_parameters = PlanckParameters()
planck13_cosmology = Cosmology( planck13_parameters.values )
bbks86 = TransferBBKS( planck13_cosmology.cpdict )


#import atomic
#import cosmology
#import utils
#import uv_bgnd
#import solvers
#import rad_src

#import f2py 

#import unittest

#Planck13_CP = 
#Planck13_CO = 



#del atomic
#del constants
#del cosmology
#del f2py
#del rad_src
#del scipy
#del utils
#del uv_bgnd
#del solvers
