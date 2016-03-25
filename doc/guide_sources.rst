=====================
Radiation Sources
=====================

Rabacus provides classes to handle various radiation sources.  These objects
can be useful in their own right, but most likely you will be passing them 
as arguments to the geometric solvers.  The following sources are available 
in Rabacus,

 - :class:`~rabacus.rad_src.point.PointSource`

 - :class:`~rabacus.rad_src.plane.PlaneSource`

 - :class:`~rabacus.rad_src.background.BackgroundSource`


.. note:: 

  The most important thing to remember is that source spectra (in general) 
  need to be normalized after being created.  See the section on 
  :ref:`section-normalization` for more details. 

.. note:: 

  All examples assume the numpy and rabacus packages have been imported. ::

    import numpy as np
    import rabacus as ra

.. testsetup:: sources

  import numpy as np
  import rabacus as ra
  q_min = 1.0; q_max = 5.0; Nnu = 8
  T_eff = 1.0e5 * ra.u.K
  src8 = ra.PlaneSource( q_min, q_max, 'thermal', T_eff=T_eff, Nnu=8 )


Arguments
================

All source classes share a common set of arguments that determine
their spectrum.  

  `q_min` (float): minimum photon energy 
  
  `q_max` (float): maximum photon energy
  
  `spectrum_type` (string): {``monochromatic``, ``powerlaw``, ``thermal``, 
  ``hm12``, ``user``}


The first two, `q_min` and `q_max`, are dimensionless floats.  They
determine the minimum and maximum photon energies to be considered.
The are interpreted as multiples of the Rydberg energy (``13.6 eV``).  
The third argument, `spectrum_type`, should be set equal to one of the 
following strings 
{``monochromatic``, ``powerlaw``, ``thermal``, ``hm12``, ``user``}.  All of 
the spectral types except ``monochromatic`` require additional keywords when
used.  For ``powerlaw`` spectra, a slope, `alpha`, must also be passed in.  
For ``thermal`` spectra, an effective temperature, `T_eff`, must be supplied.  
For the Haart and Madau 2012 model ``hm12``, one must pass in a redshift, `z`.
Finally, for the user defined spectrum, one must pass in both an array of 
energy samples, `user_E`, and an array indicating spectral shape, `user_shape`.
We will go into more detail in the examples below.    

Point source monochromatic spectrum
-------------------------------------------------

.. testcode:: sources

   q_min = 2.3; q_max = 2.3
   src = ra.PointSource(q_min, q_max, 'monochromatic')
   print src.source_type + ' ' + src.spectrum_type

.. testoutput:: sources

   point monochromatic

Plane source power law spectrum
-------------------------------------------------

.. testcode:: sources

   q_min = 1.0; q_max = 5.0
   alpha = -2.0
   src = ra.PlaneSource(q_min, q_max, 'powerlaw', alpha=alpha)
   print src.source_type + ' ' + src.spectrum_type

.. testoutput:: sources

   plane powerlaw

Point source thermal spectrum
-------------------------------------------------

.. testcode:: sources

   q_min = 1.0; q_max = 5.0
   T_eff = 1.0e5 * ra.u.K
   src = ra.PointSource(q_min, q_max, 'thermal', T_eff=T_eff)
   print src.source_type + ' ' + src.spectrum_type

.. testoutput:: sources

   point thermal

Background source HM12 spectrum
-------------------------------------------------

.. testcode:: sources

   q_min = 1.0; q_max = 5.0
   z = 3.0
   src = ra.BackgroundSource(q_min, q_max, 'hm12', z=z)
   print src.source_type + ' ' + src.spectrum_type

.. testoutput:: sources

   background hm12

Background source user defined spectrum (broken power law)
--------------------------------------------------------------

Note that `q_min` and `q_max` are ignored for user defined spectra,
but they still must be passed in.

.. testcode:: sources

   q_min = q_max = 1.0
   E = np.linspace(1.0, 5.0, 50) * 13.6 * ra.u.eV
   ibrk = 25
   E0 = E[ibrk]
   alpha1 = -1.0
   alpha2 = -2.0
   shp = np.zeros(E.size)
   shp[:ibrk+1] = (E[:ibrk+1]/E0)**alpha1
   shp[ibrk:] = (E[ibrk:]/E0)**alpha2
   src = ra.BackgroundSource(q_min, q_max, 'user', user_E=E, user_shape=shp)
   print src.source_type + ' ' + src.spectrum_type

.. testoutput:: sources

   background user



Energy Samples
================

When a source is created, a spectrum is defined by uniformly (in log space) 
sampling `Nnu` energies between the variables `q_min` and `q_max`.  The number 
of samples, `Nnu`, takes on a default value but can be passed in as a keyword 
to either increase or decrease spectral resolution.  In general, attributes are
defined to characterize the energy, frequency, and wavelength of photons.  


 +-----------+-----------------------+
 | Attribute | Description           | 
 +===========+=======================+
 | q         | energy samples / Ry   | 
 +-----------+-----------------------+
 | E         | energy samples        | 
 +-----------+-----------------------+
 | lam       | wavelengths           | 
 +-----------+-----------------------+
 | nu        | frequencies           | 
 +-----------+-----------------------+


For example, if we create a source with `Nnu=8`, ::

   q_min = 1.0; q_max = 5.0; Nnu = 8
   T_eff = 1.0e5 * ra.u.K
   src8 = ra.PlaneSource( q_min, q_max, 'thermal', T_eff=T_eff, Nnu=Nnu )

Then it will have a dimensionless fundamental photon energy array `q`, 

.. testcode:: sources

   print src8.q

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 1.          1.25849895  1.58381961  1.80808824  
     2.50848455  3.15692518  4.00147059  5.        ] dimensionless

an array with units of energy `E`, 

.. testcode:: sources

   print src8.E

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 13.6        17.11558573  21.53994668 24.59      
     34.11538992 42.93418242  54.42       68.      ] eV

a wavelength array with units of length `lam`,

.. testcode:: sources

   print src8.lam

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 9.11648497e-06  7.24393530e-06 5.75601219e-06
     5.04205757e-06  3.63425996e-06 2.88777353e-06
     2.27828364e-06  1.82329699e-06] cm

and a frequency array with units of inverse time, `nu`,

.. testcode:: sources

   print src8.nu

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 3.28846544e+15   4.13853031e+15   5.20833605e+15
     5.94583568e+15   8.24906477e+15   1.03814394e+16
     1.31586978e+16   1.64423272e+16] Hz

Note that energy samples exist exactly at the hydrogen and helium ionizing 
thresholds, {13.6, 24.59, 54.42} eV.  By default, uniform spacing is 
sacrificed to achieve this.  If exactly uniform spacing is needed, set the 
keyword `segmented` to ``False``.   

Spectral Shape
================================

Each source class has a fundamental spectral variable based on the 
characteristics of the radiation source.  When a radiation source is created,
the `spectral_type` variable along with the geometric type
(i.e. point, plane, or background) are used to determine the value and
units of the fundamental spectral variable at each energy sample. 
These variables take one form for polychromatic spectra and another
for monochromatic spectra. The variables are listed in the tables below. 

Monochromatic Sources
----------------------

 +------------+-----------------------+--------------+------------------------+
 | Source     | Fundamental Intensity | Attribute    | Units                  | 
 +============+=======================+==============+========================+
 | Point      | luminosity            | Lu           | ``erg/(s)``            | 
 +------------+-----------------------+--------------+------------------------+
 | Plane      | flux                  | Fu           | ``erg/(s cm^2)``       | 
 +------------+-----------------------+--------------+------------------------+
 | Background | specific intensity    | Inu          | ``erg/(s cm^2 sr)``    | 
 +------------+-----------------------+--------------+------------------------+

For example,

.. testcode:: sources

   pt = ra.PointSource( 2.0, 2.0, 'monochromatic' )
   pl = ra.PlaneSource( 2.0, 2.0, 'monochromatic' )
   bg = ra.BackgroundSource( 2.0, 2.0, 'monochromatic' )

   print pt.Lu
   print pl.Fu
   print bg.Inu

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 1.] erg/s
   [ 1.] erg/(cm**2*s)
   [ 1.] erg/(cm**2*s*sr)


.. note:: 

  In general, all spectra have to be normalized after being
  instanciated.  See the section on :ref:`section-normalization`.


Polychromatic Sources
----------------------

 +------------+-----------------------+--------------+------------------------+
 | Source     | Fundamental Intensity | Attribute    | Units                  | 
 +============+=======================+==============+========================+
 | Point      | luminosity density    | dLu_over_dnu | ``erg/(Hz s)``         | 
 +------------+-----------------------+--------------+------------------------+
 | Plane      | flux density          | dFu_over_dnu | ``erg/(Hz s cm^2)``    | 
 +------------+-----------------------+--------------+------------------------+
 | Background | specific intensity    | Inu          | ``erg/(Hz s cm^2 sr)`` | 
 +------------+-----------------------+--------------+------------------------+

For example,

.. testcode:: sources

   T = 1.0e5 * ra.u.K
   pt = ra.PointSource( 1.0, 5.0, 'thermal', T_eff=T, Nnu=8 )
   pl = ra.PlaneSource( 1.0, 5.0, 'thermal', T_eff=T, Nnu=8 )
   bg = ra.BackgroundSource( 1.0, 5.0, 'thermal', T_eff=T, Nnu=8 )

   print pt.dLu_over_dnu
   print pl.dFu_over_dnu
   print bg.Inu

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   [ 0.13632697  0.1662243   0.18637547  0.18957652  
     0.16102251  0.11392481  0.06087262  0.02452717] erg/(s*Hz)

   [ 0.13632697  0.1662243   0.18637547  0.18957652  
     0.16102251  0.11392481  0.06087262  0.02452717] erg/(cm**2*s*Hz)

   [ 0.13632697  0.1662243   0.18637547  0.18957652  
     0.16102251  0.11392481  0.06087262  0.02452717] erg/(cm**2*s*sr*Hz)

.. note:: 

  In general, all spectra have to be normalized after being
  instanciated.  See the section on :ref:`section-normalization`.

These variables, along with the energy sample variables in the previous 
section, are available as top level attributes in the returned object. 



Optically Thin
================

When a source is created, all quantities that can be calculated by
performing integrals over the fundamental intensity and without 
reference to optical depth are stored in a sub object called ``thin``.  
We conceptually split these variables into two types.  The first type
are those that do not refer to any photo-ionization cross-sections. 

 +-----------+-----------------------+--------------------+
 | Attribute | Description           | Units              |
 +===========+=======================+====================+
 | Ln        | photon luminosity     | ``1/s``            | 
 +-----------+-----------------------+--------------------+
 | Lu        | energy luminosity     | ``erg/s``          |
 +-----------+-----------------------+--------------------+
 | Fn        | photon flux           | ``1/(s cm^2)``     |
 +-----------+-----------------------+--------------------+
 | Fu        | energy flux           | ``erg/(s cm^2)``   |
 +-----------+-----------------------+--------------------+
 | n         | photon density        | ``1/cm^3``         |
 +-----------+-----------------------+--------------------+
 | u         | energy density        | ``erg/cm^3``       |
 +-----------+-----------------------+--------------------+
 
Note that the luminosity variables `Ln` and `Lu` will only be present
in point sources as the concept does not translate to plane or
background sources.  We also note that the API documentation for each
geometric source type 
(:class:`~rabacus.rad_src.point.PointSource`, 
:class:`~rabacus.rad_src.plane.PlaneSource`,
:class:`~rabacus.rad_src.background.BackgroundSource`) 
goes into detail about how each quantity above is calculated.  


The second type of attribute stored in the ``thin`` object are
photo-ionization and heating rates.  These are related to photo-ionization 
cross-sections. 

 +-----------+------------------------------+--------------------+
 | Attribute | Description                  | Units              |
 +===========+==============================+====================+
 | H1i       | H I photo-ionization rate    | ``1/s``            |
 +-----------+------------------------------+--------------------+
 | He1i      | He I photo-ionization rate   | ``1/s``            |
 +-----------+------------------------------+--------------------+
 | He2i      | He II photo-ionization rate  | ``1/s``            |
 +-----------+------------------------------+--------------------+
 | H1h       | H I photo-heating rate       | ``erg/s``          |
 +-----------+------------------------------+--------------------+
 | He1h      | He I photo-heating rate      | ``erg/s``          |
 +-----------+------------------------------+--------------------+
 | He2h      | He II photo-heating rate     | ``erg/s``          |
 +-----------+------------------------------+--------------------+


It's important to note that many of these attributes will be variables in
plane and background sources, but functions of distance in point
sources.   

Here we present some examples.  First, the H I photoionization rate at
a distance of 1 kpc from a thermal point source normalized to emit 1.0e50
photons per second,

.. testcode:: sources

   T = 1.0e5 * ra.u.K
   pt = ra.PointSource( 1.0, 5.0, 'thermal', T_eff=T )
   pt.normalize_Ln( 1.0e50/ra.u.s )

   H1i = pt.thin.H1i( 1.0*ra.u.kpc )
   print H1i.round(15)

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   1.387e-12 1/s


Next, the same quantity for a thermal plane source normalized to have
a flux of 1.0e6 photons per second per square centimeter,

.. testcode:: sources

   T = 1.0e5 * ra.u.K
   pl = ra.PlaneSource( 1.0, 5.0, 'thermal', T_eff=T )
   pl.normalize_Fn( 1.0e6/ra.u.s/ra.u.cm**2 )

   H1i = pl.thin.H1i
   print H1i.round(15)

.. testoutput:: sources
   :options: +NORMALIZE_WHITESPACE

   1.659e-12 1/s




.. _section-normalization:

Normalization
================

Normalization of spectra after creation is crucial to obtaining the results
you expect.  All sources except background sources with spectral type 
``hm12`` should be normalized after being created.  Here we will list the 
normalization functions available for each source class. 

Point
------------

- :func:`~rabacus.rad_src.point.PointSource.normalize_Ln`: 
  normalize to fixed photon luminosity

- :func:`~rabacus.rad_src.point.PointSource.normalize_Lu`:  
  normalize to fixed energy luminosity


Plane
------------

- :func:`~rabacus.rad_src.plane.PlaneSource.normalize_n`: 
  normalize to fixed photon number density

- :func:`~rabacus.rad_src.plane.PlaneSource.normalize_Fn`: 
  normalize to fixed photon flux

- :func:`~rabacus.rad_src.plane.PlaneSource.normalize_H1i`: 
  normalize to fixed optically thin H1 photoionization rate


Background
------------

- :func:`~rabacus.rad_src.background.BackgroundSource.normalize_n`: 
  normalize to fixed photon number density

- :func:`~rabacus.rad_src.background.BackgroundSource.normalize_H1i`: 
  normalize to fixed optically thin H1 photoionization rate



User Defined Spectra
=====================

If the pre defined spectral types do not suite your needs, you can define 
your own.  This is done by setting the keyword `spectral_type` to ``user``
and setting the keywords `user_E` and `user_shape`.  The array `user_E` 
should have units of energy and the `user_shape` array should be 
dimensionless.  Leaving the shape array without units allows for flexibility. 
If the `user_shape` array is passed to a background source it will take 
the units of specific intensity ``erg/(Hz s cm^2 sr)``.  Likewise, if the 
array is passed to a plane source it will take the units ``erg/(Hz s cm^2)`` 
and if passed to a point source it will take the units ``erg/(Hz s)``.
Note that for user defined spectra, `q_min` and `q_max` are
ignored.  For example, to create a background source with a flat 
spectrum. :: 







