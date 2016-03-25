"""
Classes to provide atomic data.  In particular,

  - (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
    Chemistry rates 

  - (:class:`~rabacus.atomic.cooling.CoolingRates`): 
    Cooling rates 
 
  - (:class:`~rabacus.atomic.photo_xsection.PhotoXsections`): 
    Photoionization cross sections

.. note::

  Both chemistry and cooling objects take three case A fraction 
  arguments, `fcA_H2`, `fcA_He2`, and `fcA_He3`.  These variables  
  interpolate between case A and case B recombination rates 
  (0.0 = case B, 1.0 = case A). 


Chemistry Examples:
---------------------

    The following example will return chemistry rates at 128 temperatures 
    logarithmically spaced between 10^4 and 10^5 K using case A recombination
    rates. ::

      import numpy as np
      import rabacus as ra

      NT = 128
      T = 10**np.linspace( 4.0, 5.0, NT ) * ra.U.K
      fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0
      kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )

    The object ``kchem`` can now be used in calls to collisional equilibrium
    solvers such as :class:`~rabacus.solvers.one_zone.Solve_CE`. 
    If you are interested in solutions which include non zero photoionization
    rates, such as :class:`~rabacus.solvers.one_zone.Solve_PCE`, then those 
    rates must be passed to the chemistry object when it is created, ::

      H1i = np.ones( NT ) * 1.0e-12 / ra.U.s
      He1i = np.ones( NT ) * 1.0e-13 / ra.U.s
      He2i = np.ones( NT ) * 1.0e-14 / ra.U.s
      kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1i=H1i, He1i=He1i, He2i=He2i )



Cooling Examples:
-------------------

    The cooling rates object behaves analogously to the chemistry rates object.
    The following example will return cooling rates at 128 temperatures 
    logarithmically spaced between 10^4 and 10^5 K using recombination rates 
    half way between case A and case B. :: 

      import numpy as np
      import rabacus as ra

      NT = 128
      T = 10**np.linspace( 4.0, 5.0, NT ) * ra.U.K
      fcA_H2 = 0.5; fcA_He2 = 0.5; fcA_He3 = 0.5
      kcool = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3 )


    If photoheating rates are needed in a particular solution, they should 
    be passed to the cooling object when it is created,  ::

      H1h = np.ones( NT ) * 1.0e-12 * ra.U.eV / ra.U.s
      He1h = np.ones( NT ) * 1.0e-12 * ra.U.eV / ra.U.s
      He2h = np.ones( NT ) * 1.0e-14 * ra.U.eV / ra.U.s
      kcool = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3,
                                      H1h=H1h, He1h=He1h, He2h=He2h )


    Note that the heating and cooling rates per volume for a given density 
    of ions and redshift is calculated using the functions    
    :func:`~rabacus.atomic.cooling.CoolingRates.return_cooling` and 
    :func:`~rabacus.atomic.cooling.CoolingRates.return_heating`.


Photoionization Cross Section Examples:
----------------------------------------

    The photoionization cross section class contains functions to calculate
    the cross section as well as the threshold ionization energies of common
    ions.  The following code will store the photoionization cross sections of
    HI, HeI and HeII at 200 frequencies between 1 and 100 Rydbergs::

      import numpy as np
      import rabacus as ra

      Nnu = 200
      px = ra.atomic.PhotoXsections()
      q_min = 1.0
      q_max = 1.0e2
      log_E_min = np.log10( q_min * px.Eth_H1.magnitude )
      log_E_max = np.log10( q_max * px.Eth_H1.magnitude )
      E = 10**np.linspace( log_E_min, log_E_max, Nnu ) * px.Eth_H1.units
      sigma_H1 = px.sigma_H1( E )
      sigma_He1 = px.sigma_He1( E )
      sigma_He2 = px.sigma_He2( E )

""" 

#from chemistry import *
#from cooling import *
#from photo_xsection import *
#from hydrogen import *
#from helium import *



