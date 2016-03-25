=====================
Geometric Solvers
=====================

The geometric solvers calculate ionization and/or temperature structure 
in simple spherical and planar geometries.  The systems are discretized 
into layers and the attenuation of radiation is calculated in each.  
Sweeps are made through the layers during which the appropriate single zone 
solver is called in each layer.  The sweeps continue until convergence is 
achieved. 

:class:`~rabacus.f2py.slab_plane.SlabPln`:
  A planar gas distribution with plane parallel radiation incident from 
  one side. 

:class:`~rabacus.f2py.slab_plane.Slab2Pln`:
  A planar gas distribution with plane parallel radiation incident from 
  both sides. 

:class:`~rabacus.f2py.slab_bgnd.Slab2Bgnd`:
  A planar gas distribution with an isotropic radiation field incident from 
  both sides. 

:class:`~rabacus.f2py.sphere_stromgren.SphereStromgren`:
  A spherically symmetric gas distribution with a point source at the center.

:class:`~rabacus.f2py.sphere_bgnd.SphereBgnd`:
  A spherically symmetric gas distribution in a uniform background radiation 
  field. 


Positional Arguments
=====================

All geometric classes share a common set of arguments.  The positional 
(i.e. mandatory) arguments describe the gas distribution. 

Args:

  `Edges` (array): Position of discrete elements. 
        
  `T` (array): Temperature in each discrete element.
        
  `nH` (array): Hydrogen number density in each discrete element.
        
  `nHe` (array): Helium number density in each discrete element.
        
  `rad_src` (:class:`~rabacus.rad_src.point.PointSource`, :class:`~rabacus.rad_src.plane.PlaneSource`, or :class:`~rabacus.rad_src.background.BackgroundSource`): Radiation sources. 


The first argument, `Edges`, describes the edges of each layer. 
For slabs, `Edges` is the distance from one of the illuminated surfaces of the 
slab to each layer edge.  For spheres, `Edges` is the distance from 
the center of the sphere to each spherical layer edge.  For a system 
discretized into `N` layers, `Edges` is an array with `N+1` entries.  Most of 
the time, the first element should be zero and the last the thickness of the 
slab or radius of the sphere, but this is not required.  The next three 
arguments, `T`, `nH`, and `nHe`, are arrays with `N` elements and describe the 
initial temperature and density in each layer.  The final required argument(s) 
should be instances of classes describing sources of radiation. 


Keyword Arguments
=====================

Many keyword (i.e. optional) arguments are shared among the geometric classes 
and can be used to change the default solver behaviour. 


Recombination 
----------------------------------

The default solver behaviour is to use case A recombination rates in each 
discrete element.  These keyword arguments alter this behaviour. 

Kwargs:

  `rec_meth` (string): How to treat recombinations {``fixed``, ``outward``, 
  ``radial``, ``ray``}

  `fixed_fcA` (float): If `rec_meth` = ``fixed``, constant case A fraction. 



If `rec_meth` = ``fixed`` a constant case A fraction will be used in each
discrete element.  The keyword `fixed_fcA` sets this value.  If `fixed_fcA` 
is set to ``1.0`` then case A rates are used throughout and recombination 
radiation is ignored.  If `fixed_fcA` is set equal to ``0.0`` then case B 
rates are used throughout.  This is known as the on-the-spot approximation. 
Setting `fixed_fcA` to a float between ``0.0`` and ``1.0`` will produce 
results bracketed by case A and case B results.  

In spherical geometries, one can also choose to set `rec_meth` to ``outward``, 
``radial``, or ``isotropic``.  If this is done, case A recombination rates are 
used, and recombination radiation is transported through the density field.  
For `rec_meth` = ``outward``, all recombination radiation is transported 
radially outward (see [Ritzerveld05]_).  For `rec_meth` = ``radial``, 
recombination radiation is transported radially inward and outward.  For 
`rec_meth` = ``isotropic``, recombination radiation is transported 
isotropically from each layer.  The last is the most realistic option.   

In planar geometries, one can choose to set `rec_meth` to ``ray``.  If this 
is done, case A rates are used and recombination radiation contributes an 
isotropic flux from each slab layer which is transported through the density 
field. 



Other 
----------------------------------

Kwargs:

  `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve for 
  equilibrium T (i.e. use :class:`~rabacus.f2py.ion_solver.Solve_PCTE` as the
  single zone solver in each element).

  `z` (float): Redshift, only need if `find_Teq` = ``True``

  `tol` (float): tolerance for all convergence tests

  `thin` (bool): if ``True`` radiation is not attenuated by passage through
  absorbing gas. 


If the `find_Teq` keyword is set to ``True``, photo collisional thermal 
equilibrium will be found in each discrete element and the user must provide
a redshift, `z`.  If the `thin` keyword is set to ``True``, it is assumed that
the radiation is not attenuated by passage through the discrete elements and
optically thin values for the photoionization and/or photoheating rates are 
used.  The argument `tol` is a general slider in which lower values produce 
longer run times but more accurate results.  It is the deviation from unity
allowed in the sum of the electron number density over all discrete elements
between consecutive iterations. 



References
==========

.. [Ritzerveld05] http://arxiv.org/abs/astro-ph/0506637
