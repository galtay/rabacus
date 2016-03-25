r"""
Tools for calculating the ionization and temperature structure in primordial 
(i.e. composed of only hydrogen and helium) cosmological gas.  
The single zone solvers in :mod:`~rabacus.solvers.one_zone` form
the core of this package and contain classes for solving three types of 
ionization/temperature equilibrium,

- (:class:`~rabacus.solvers.one_zone.Solve_CE`)
  **collisional equilibrium** in which ionization fractions are found such 
  that collisional ionizations are balanced by recombinations at a fixed 
  temperature. 

- (:class:`~rabacus.solvers.one_zone.Solve_PCE`)
  **photo-collisional equilibrium** in which ionization fractions are found 
  such that photo and collisional ionizations are balanced by recombinations
  at a fixed temperature. 

- (:class:`~rabacus.solvers.one_zone.Solve_PCTE`)
  **photo-collisional-thermal equilibrium** in which ionization fractions are
  found such that photo and collisional ionizations are balanced by 
  recombinations and a temperature is found such that heating balances 
  cooling. 

The other modules in this package make use of 
:mod:`~rabacus.solvers.one_zone` to calculate the ionization 
and/or temperature structure in the case of specific geometries and sources. 

- (:mod:`~rabacus.solvers.calc_slab`)
  An infinite plane gas distribution with plane parallel radiation incident
  from one side. 

- (:mod:`~rabacus.solvers.calc_slab_2`)
  An infinite plane gas distribution with plane parallel radiation incident
  from both sides. 

- (:mod:`~rabacus.solvers.calc_iliev_sphere`)
  A spherically symmetric gas distribution with a point source at the center.
  


Recombination radiation
-----------------------

Recombination radiation is handled in all solutions using the 
"variable on-the-spot" approximation.  In the standard "on-the-spot" 
approximation, the effect of ionizing recombination photons are modeled by 
reducing the recombination rate.  The recombination rate to all atomic levels 
is called the case A rate.  Recombinations to the ground state can cause 
ionizations and so a case B rate is defined that negelects these 
recombinations, 

.. math::
  \alpha_{\rm B} = \alpha_{\rm A} - \alpha_{\rm 1}

We define a variable called the case A fraction, :math:`f_{\rm A}`, for each 
recombining ionic species 
:math:`\{ {\scriptstyle{\rm HII, \, HeII, \, HeIII }} \}`
such that the recombination rate used in the solution is, 

.. math:: 
  \log \alpha = \log \alpha_{\rm A} f_{\rm A} + 
      \log \alpha_{\rm B} (1-f_{\rm A})

The solutions above discretize spheres into thin shells and planes into thin
slabs.  We will refer to both of these as discrete elements.  The main 
variable that determines what :math:`f_{\rm A}` will be in each discrete 
element is `rec_meth`.   
If `rec_meth` = ``fixed`` then a constant case A fraction, 
:math:`f_{\rm A}` = `fixed_fcA`, is used in each element.  If `rec_meth` = 
``thresh``, the case A fraction will vary in each discrete element.  
For each element, the optical depth is probed out to a fixed 
distance, :math:`L`.  The probability of a recombination photon being absorbed 
is then calculated,

.. math:: P = 1 - \exp[ -\tau(\nu_{\rm _X})]

where :math:`\tau(\nu_{\rm _X})` is the optical depth for a given species at
its ionization threshold.  This probability is then used to interpolate 
between case A and case B rates.  The variable 
`thresh_P_A` (:math:`P_{\rm A}`) determines the probability at which the 
transition to case B starts and the variable `thresh_P_B` (:math:`P_{\rm B}`)
determines the probability at which the transition to case B ends.  
The variable `thresh_xmfp` (:math:`C`) determines the distance 
:math:`L` and represents a multiplier to the fully neutral mean free path in a 
discrete element,

.. math::
  \lambda = ( n \sigma_{\rm _X} )^{-1}, \quad L = C \lambda

where :math:`n` is either :math:`n_{\rm _H}` or :math:`n_{\rm _{He}}`.  
If :math:`P < P_{\rm A}` then :math:`f_{\rm A} = 1`, 
if :math:`P > P_{\rm B}` then :math:`f_{\rm A} = 0`,   
otherwise,
 
.. math::
  f_{\rm A} =   (P_{\rm B} - P) / ( P_{\rm B} - P_{\rm A} ) 

"""



from .one_zone import *
from .calc_slab import *
from .calc_slab_2 import *
from .calc_iliev_sphere import *
from .calc_slab_analytic import *


