r"""
A package containing classes which represent sources of radiation.  To create 
a source, one needs to choose a geometry class and a spectral shape.  The 
choices for geometry are:

- (:class:`~rabacus.rad_src.point.PointSource`) 
  
- (:class:`~rabacus.rad_src.plane.PlaneSource`) 
  
- (:class:`~rabacus.rad_src.background.BackgroundSource`) 


When creating a source, instantiate a geometric source class and indicate
the spectral shape by setting the argument `spectrum_type` equal to one of
{``monochromatic``, ``powerlaw``, ``thermal``, ``hm12``}.  Note that the 
spectral shape does not determine the normalization of the spectrum, therefore
classes will be returned with an arbitrary normalization.  The exception
is the Haardt & Madau 2012 model which is returned with the normalization
native to the model. Each geometric class contains methods to (re)normalize 
the spectrum.   


Segmented Spectra 
-----------------
By default, the spectra created for sources will be segmented.  This means that
an exactly uniform spacing of the energy samples in log energy space will be 
sacrificed to guarantee that samples exist at the ionization thresholds of 
neutral hydrogen, neutral helium, and singly ionized helium
(:math:`E_{\rm _{HI}}`, :math:`E_{\rm _{HeI}}`, :math:`E_{\rm _{HeII}}`).  
If one tries to create a segmented spectrum in which the energy samples do not 
cover the range between :math:`E_{\rm _{HI}}` and :math:`E_{\rm _{HeII}}`  
(see the `q_min` and `q_max` variables) an error will occur.

Grey Quantities
---------------

Frequency averaged or "grey" quantities will be calculated for any segmented 
spectrum.  We show an example for the 
:class:`~rabacus.rad_src.background.BackgroundSource` 
class here, but the principle is the same for all source types. The grey 
photoionization cross sections are defined as follows, 

.. math:: 
  \sigma_{\rm _{HI}}^{\rm grey} = 
  \frac{ 
  \int_{\nu_{\rm HI}}^{\nu_{\rm HeI}} 
      I_{\nu} \sigma_{\rm _{HI}} \frac{d\nu}{\nu}
  }{ 
  \int_{\nu_{\rm HI}}^{\nu_{\rm HeI}} 
      I_{\nu} \frac{d\nu}{\nu}
  }, \quad
  \sigma_{\rm _{HeI}}^{\rm grey} = 
  \frac{ 
  \int_{\nu_{\rm HeI}}^{\nu_{\rm HeII}} 
      I_{\nu} \sigma_{\rm _{HeI}} \frac{d\nu}{\nu}
  }{ 
  \int_{\nu_{\rm HeI}}^{\nu_{\rm HeII}} I_{\nu} \frac{d\nu}{\nu}
  }, \quad
  \sigma_{\rm _{HeII}}^{\rm grey} = 
  \frac{ 
  \int_{\nu_{\rm HeII}}^{\nu_{\rm max}} 
      I_{\nu} \sigma_{\rm _{HeII}} \frac{d\nu}{\nu}
  }{ 
  \int_{\nu_{\rm HeII}}^{\nu_{\rm max}} I_{\nu} \frac{d\nu}{\nu}
  }

We can implicitly define a grey frequency and therefore a grey energy,

.. math:: 
  \sigma_{\rm _{HI}}( \nu_{\rm _{HI}}^{\rm grey} ) =
  \sigma_{\rm _{HI}}^{\rm grey}, \quad
  E_{\rm _{HI}}^{\rm grey} = h \nu_{\rm _{HI}}^{\rm grey}

  

where similar relationships hold for the helium ions. 


Examples
--------

The following code will create a point source with a photon luminosity of 
5.0e48 photons per second and a 1.0e5 K blackbody spectrum between 1 and 
5 Rydbergs:: 


   import rabacus as ra
   q_min = 1.0; q_max = 5.0; T_eff = 1.0e5 * ra.U.K
   pt_t1e5 = ra.rad_src.PointSource( q_min, q_max, 'thermal', T_eff=T_eff )
   pt_t1e5.normalize_Ln( 5.0e48 / ra.U.s )

The following code will create a uniform background source with a spectral 
shape given by the model of Haardt and Madau 2012 at redshift 3 between 1 and 
400 Rydbergs::

   import rabacus as ra
   q_min = 1.0; q_max = 4.0e2; z=3.0
   bgnd_hm12 = ra.rad_src.BackgroundSource( q_min, q_max, 'hm12', z=z ) 

"""


from background import *
from plane import *
from point import *

