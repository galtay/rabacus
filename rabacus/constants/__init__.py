r"""

The constants package provides classes which handle physical constants and 
units.  In particular,

  - (:class:`~rabacus.constants.units.Units`): 
    Basic Units
 
  - (:class:`~rabacus.constants.units.CosmoUnits`): 
    Cosmological Units (includes factors of `h` and `a`). 

  - (:class:`~rabacus.constants.physical.PhysicalConstants`): 
    Physical Constants


Units and Physical Constants
-----------------------------

All variables in Rabacus are 
`quantities <https://pythonhosted.org/quantities>`_ .  In other words, they 
have a magnitude and associated units.  Whenever the module is imported an 
instance of both :class:`~rabacus.constants.physical.PhysicalConstants` and
:class:`~rabacus.constants.units.Units` will be attached as ``pc`` and ``u``
respectively.  

>>> import rabacus as ra
>>> d = 100 * ra.u.kpc
>>> d
array(100.0) * kpc

Note that scalars are represented as 0-d arrays.  The units of any quantity 
are stored as an attribute and can be changed by string assignment,

>>> d.units
array(1.0) * kpc
>>> d.units = 'cm'
>>> d
array(3.08568025e+23) * cm

Each quantity has a function called ``rescale`` that will return the quantity 
in a different set of units without modifying its current units. 

>>> d.rescale('mile')
array(1.9173528158056945e+18) * mi
>>> d
array(3.08568025e+23) * cm


One can access just the magnitude of any quantity in its current units,

>>> d.magnitude
array(3.08568025e+23)


Attempting to take the log of a quantity with dimensions will raise 
a ValueError so one should always use a magnitude when calculating logarithms.  

>>> import numpy as np
>>> np.log10( d.magnitude ) 
23.489350920801908

Dividing two length units will produce a hybrid unit, 

>>> r = d / (10 * ra.u.kpc)
>>> r
array(3.0856802499999998e+22) * cm/kpc

that can be reduced using the attribute ``simplified``. 

>>> r.simplified
array(10.0) * dimensionless


Cosmological Units
-----------------------------


By creating an instance of :class:`~rabacus.constants.units.CosmoUnits` one 
can also use the hubble parameter and scale factor as units. However one 
must be careful about accidentally mixing unit systems with different values
for `h` and `a`.  When rabacus is imported, an instance of the cosmology class
is created with parameters reported by the Planck mission.  The included Planck
cosmology is called ``planck13_cosmology``.  

>>> import rabacus as ra
>>> ra.planck13_cosmology.H0
array(100.0) * km*hh/(s*Mpc)

Note that the `h` unit is indicated with a double ``hh``.  Similarly the 
scale factor unit is indicated with a double ``aa``.  This is because ``h`` 
and ``a`` are already units in the quantities package.  The above statement
is an expression of the fact that the Hubble parameter is always 100 in units
of km/s/(Mpc/h).  By requesting the quantity in units of km/s/Mpc we recover
the value specific to the Planck cosmology, 

>>> ra.planck13_cosmology.H0.rescale('km/s/Mpc')
array(67.11182) * km/(s*Mpc)

Suppose we create another hubble parameter, 

>>> h = 0.7
>>> H0 = 100 * h * ra.u.km / ra.u.s / ra.u.Mpc
>>> H0
array(70.0) * km/(s*Mpc)

and then request it in units of km/s/(Mpc/h)

>>> H0.units = 'hh*km/s/Mpc'
>>> H0
array(104.30353401233941) * km*hh/(s*Mpc)

We don't get the expected value of 100 because the hubble parameter 
associated with the set of units in ``ra.u`` is the Planck13 value.  
However, if we create a new instance of 
:class:`~rabacus.constants.units.CosmoUnits` we get the correct
result,

>>> cu = ra.CosmoUnits(h=0.7,a=1.0)
>>> H0 = 100 * cu.h * cu.km / cu.s / cu.Mpc
>>> H0
array(100.0) * km*hh/(s*Mpc)
>>> H0.rescale( 'km/s/Mpc' )
array(70.0) * km/(s*Mpc)

This occurs automatically whenever an instance of the 
:class:`~rabacus.cosmology.general.Cosmology` class is 
created, 

>>> di = {'omegam': 0.3, 'omegal': 0.7, 'omegab': 0.05}
>>> di.update( {'h': 0.7, 'sigma8': 0.8, 'ns':1.0, 'Yp':0.25} )
>>> co = ra.Cosmology(di)
>>> H0 = 100 * co.cu.h * co.cu.km / co.cu.s / co.cu.Mpc
>>> H0
array(100.0) * km*hh/(s*Mpc)
>>> H0.rescale( 'km/s/Mpc' )
array(70.0) * km/(s*Mpc)


.. warning :: 

  Always create a new instance of :class:`~rabacus.constants.units.CosmoUnits` 
  when you want a new value of the hubble parameter.  This happens automatically
  when you create an instance of the 
  :class:`~rabacus.cosmology.general.Cosmology` class. 








"""

#from . units import *
#from . physical import *


