""" Wrapper to make the fortran calling identical to the python calling. """ 

import rabacus_fc 
import numpy as np

from rabacus.atomic import chemistry 
from rabacus.atomic import cooling
from rabacus.constants import units



__all__ = ['Solve_CE', 'Solve_PCE', 'Solve_PCTE']



class Solve_CE:

    r""" 
    Calculates and stores collisional ionization equilibrium (i.e. all 
    photoionization rates are equal to zero) given a chemistry rates object. 
    This is an exact analytic solution. 

    Args:
    
      `kchem` (:class:`~rabacus.f2py.chem_cool_rates.ChemistryRates`): 
      Chemistry rates

    Kwargs: 

      `fpc` (bool): If True, perform floating point error correction. 
       

    Attributes:

      `T` (array): temperature
 
      `T_floor` (float): minimum possible input temperature 
 
      `H1` (array): neutral hydrogen fraction, |xH1| = |nH1| / |nH|

      `H2` (array): ionized hydrogen fraction, |xH2| = |nH2| / |nH|

      `He1` (array): neutral helium fraction,  |xHe1| = |nHe1| / |nHe|

      `He2` (array): singly ionized helium fraction,  |xHe2| = |nHe2| / |nHe|

      `He3` (array): doubly ionized helium fraction,  |xHe3| = |nHe3| / |nHe|


    Notes: 

      This class finds solutions (`H1`, `H2`, `He1`, `He2`, `He3`) to the 
      following equations 
      
      .. math:: 
        \frac{dx_{\rm _{HI}}}{dt} &= 
        - C_{\rm _{HI}}  x_{\rm _{HI}}
        + R_{\rm _{HII}}  x_{\rm _{HII}} = 0 
        \\
        \frac{dx_{\rm _{HII}}}{dt} &= 
        C_{\rm _{HI}}  x_{\rm _{HI}}
        - R_{\rm _{HII}}  x_{\rm _{HII}} = 0
        \\
        \frac{dx_{\rm _{HeI}}}{dt} &= 
        - C_{\rm _{HeI}}  x_{\rm _{HeI}}
        + R_{\rm _{HeII}}  x_{\rm _{HeII}} = 0 
        \\
        \frac{dx_{\rm _{HeII}}}{dt} &= 
        C_{\rm _{HeI}}  x_{\rm _{HeI}}
        -( C_{\rm _{HeII}}  + R_{\rm _{HeII}}  ) x_{\rm _{HeII}} +
        R_{\rm _{HeIII}}  x_{\rm _{HeIII}} = 0 
        \\
        \frac{dx_{\rm _{HeIII}}}{dt} &= 
        C_{\rm _{HeII}}  x_{\rm _{HeII}}
        - R_{\rm _{HeIII}}  x_{\rm _{HeIII}} = 0 


      with the following closure relationship

      .. math:: 
        1 &= x_{\rm _{HI}} + x_{\rm _{HII}}  \\
        1 &= x_{\rm _{HeI}} + x_{\rm _{HeII}} + x_{\rm _{HeIII}} 


      where the :math:`C_{\rm _X}` are collisional ionization rates and the 
      :math:`R_{\rm _X}` are recombination rates.  These solutions depend only
      on temperature. 


    Examples:

      The following code will solve for collisional ionization equilibrium at 
      128 temperatures equally spaced between 1.0e4 and 1.0e5 K using case A 
      recombination rates for all species:: 

        import numpy as np
        import rabacus as ra
        N = 128
        T = np.logspace( 4.0, 5.0, N ) * ra.u.K 
        fcA_H2 = 1.0; fcA_He2 = 1.0; fcA_He3 = 1.0
        kchem = ra.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3 )
        x = ra.Solve_CE( kchem )



    .. seealso::

       :class:`Solve_PCE`, :class:`Solve_PCTE`


    """ 

    def __init__(self, kchem):

        # attach units. 
        #----------------------------------------------
        self.u = units.Units()

        # solve
        #---------------------------------------------
        self.set( kchem )



    def set(self, kchem ):

        """ 
        This function is called when a new instance of :class:`Solve_CE` is 
        created.  Updated ionization fractions can be calculated by calling
        it again with a new kchem object. 
        
        Args:

          `kchem` (:class:`~rabacus.f2py.chem_cool_rates.ChemistryRates`): 
          Chemistry rates

        """ 

        # assert 1-D input
        #----------------------------------------------
        if len(kchem.T.shape) != 1:
            msg = "input arrays must be 1-D for fortran routines." 
            raise ValueError(msg)

        # choose scalar or vector routines
        #----------------------------------------------
        nn = kchem.T.size

        if nn == 1:
            solve = rabacus_fc.ion_solver.solve_ce_s
        else:
            solve = rabacus_fc.ion_solver.solve_ce_v

        (xH1, xH2, xHe1, xHe2, xHe3) = solve( 
            kchem.reH2, kchem.reHe2, kchem.reHe3, 
            kchem.ciH1, kchem.ciHe1, kchem.ciHe2, 
            nn )

        # package output
        #----------------------------------------------
        self.H1 = np.ones(nn) * xH1 * self.u.dimensionless
        self.H2 = np.ones(nn) * xH2 * self.u.dimensionless
        self.He1 = np.ones(nn) * xHe1 * self.u.dimensionless
        self.He2 = np.ones(nn) * xHe2 * self.u.dimensionless
        self.He3 = np.ones(nn) * xHe3 * self.u.dimensionless
        self.T = kchem.T

        self.kchem = kchem







class Solve_PCE:

    r""" 
    Calculates and stores photo collisional ionization equilibrium given the
    density of hydrogen and helium plus a chemistry rates object. This is an 
    iterative solution.  The density and temperature arrays must have 
    the same size and the output arrays are 1-D of that size.  

    .. note::

      The user must pass in photoionization rates for each species when 
      creating the chemistry rates object (see examples below). 


    Args:

      `nH_in` (array): number density of hydrogen

      `nHe_in` (array): number density of helium

      `kchem` (:class:`~rabacus.f2py.chem_cool_rates.ChemistryRates`): 
      Chemistry rates

    Kwargs:

      `tol` (float): convergence tolerance,
      abs(1 - `ne_new` / `ne_old`) < `tol`
       

    Attributes: 

      `nH` (array): H number density.  

      `nHe` (array): He number density.  

      `ne` (array): electron number density 

      `T` (array): temperature 

      `T_floor` (float): minimum possible temperature 

      `H1` (array): neutral hydrogen fraction, |xH1| = |nH1| / |nH|

      `H2` (array): ionized hydrogen fraction, |xH2| = |nH2| / |nH|

      `He1` (array): neutral helium fraction,  |xHe1| = |nHe1| / |nHe|

      `He2` (array): singly ionized helium fraction,  |xHe2| = |nHe2| / |nHe|

      `He3` (array): doubly ionized helium fraction,  |xHe3| = |nHe3| / |nHe|

      `Hatom` (:class:`~rabacus.atomic.hydrogen.Hydrogen`): H atom

      `itr` (int): number of iterations to converge

      `ne_min` (array): electron density in collisional equilibrium

      `ne_max` (array): electron density for full ionization

      `xH1_ana` (array): analytic H neutral fraction (neglects electrons from 
      elements other than H)

      `ne_ana` (array): electron number density implied by `xH1_ana`

      `x_ce` (:class:`Solve_CE`): the collisional equilibrium solutions. 


    Notes: 

      This class finds solutions (`H1`, `H2`, `He1`, `He2`, `He3`) to the 
      following equations 

      .. math:: 
        \frac{dx_{\rm _{HI}}}{dt} &= 
        - (\Gamma_{\rm _{HI}} + C_{\rm _{HI}} n_{\rm _e}) x_{\rm _{HI}}
        + R_{\rm _{HII}} n_{\rm _e} x_{\rm _{HII}} = 0 
        \\
        \frac{dx_{\rm _{HII}}}{dt} &= 
        (\Gamma_{\rm _{HI}} + C_{\rm _{HI}} n_{\rm _e}) x_{\rm _{HI}}       
        - R_{\rm _{HII}} n_{\rm _e} x_{\rm _{HII}} = 0
        \\
        \frac{dx_{\rm _{HeI}}}{dt} &= 
        - (\Gamma_{\rm _{HeI}} + C_{\rm _{HeI}} n_{\rm _e}) x_{\rm _{HeI}}
        + R_{\rm _{HeII}} n_{\rm _e} x_{\rm _{HeII}} = 0 
        \\
        \frac{dx_{\rm _{HeII}}}{dt} &= 
        (\Gamma_{\rm _{HeI}} + C_{\rm _{HeI}} n_{\rm _e}) 
        x_{\rm _{HeI}} - \\
        & \quad ( \Gamma_{\rm _{HeII}} + C_{\rm _{HeII}} n_{\rm _e} + 
        R_{\rm _{HeII}} n_{\rm _e} ) 
        x_{\rm _{HeII}} + \\
        & \quad R_{\rm _{HeIII}} n_{\rm _e} x_{\rm _{HeIII}} = 0 
        \\
        \frac{dx_{\rm _{HeIII}}}{dt} &= 
        (\Gamma_{\rm _{HeII}} + C_{\rm _{HeII}} n_{\rm _e}) 
        x_{\rm _{HeII}}
        - R_{\rm _{HeIII}} n_{\rm _e} x_{\rm _{HeIII}} = 0 


      with the following closure relationships

      .. math:: 
        1 &= x_{\rm _{HI}} + x_{\rm _{HII}} 
        \\
        1 &= x_{\rm _{HeI}} + x_{\rm _{HeII}} + x_{\rm _{HeIII}} 
        \\
        n_{\rm e} &= x_{\rm _{HII}} n_{\rm _{H}} + 
        ( x_{\rm _{HeII}} + 2 x_{\rm _{HeIII}} ) n_{\rm _{He}}
        


      where the :math:`\Gamma_{\rm _X}` are photoionization rates, the 
      :math:`C_{\rm _X}` are collisional ionization rates and the 
      :math:`R_{\rm _X}` are recombination rates.  These solutions depend on
      temperature and density. 


    Examples: 
    
      The following code will solve for photo-collisional ionization 
      equilibrium at 128 density-temperature pairs using recombination rates 
      half way between case A and case B for all species.  Note that the 
      photoionization rate arrays are grouped with the chemistry rates and so 
      should always be the same size as the temperature array:: 

        import numpy as np
        import rabacus as ra

        fcA_H2 = 0.5; fcA_He2 = 0.5; fcA_He3 = 0.5

        N = 128; Nrho = NT = N; Yp = 0.24

        nH = 10**np.linspace( -4.0, -3.0, Nrho ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = 10**np.linspace( 4.0, 5.0, NT ) * ra.U.K

        H1i = np.ones( NT ) * 8.22e-13 / ra.U.s
        He1i = np.ones( NT ) * 4.76e-13 / ra.U.s
        He2i = np.ones( NT ) * 3.51e-15 / ra.U.s

        kchem = ra.f2py.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1i=H1i, He1i=He1i, He2i=He2i )

        x = ra.f2py.Solve_PCE( nH, nHe, kchem )


    .. seealso::

      :class:`Solve_CE`, :class:`Solve_PCTE`


    """ 


    def __init__(self, nH, nHe, kchem, tol=1.0e-8):

        # attach units. 
        #----------------------------------------------
        self.u = units.Units()

        # check input
        #----------------------------------------------
        if nH.shape == ():
            nH = np.ones(1) * nH

        if nHe.shape == ():
            nHe = np.ones(1) * nHe

        if not len(nH.shape) == len(nHe.shape) == len(kchem.T.shape) == 1:
            msg = '\n Input arrays must be one dimensional'
            raise ValueError(msg)

        if not hasattr(nH,'units'):
            msg = '\n Input variable nH must have units \n'
            raise ValueError(msg) 
        else:
            nH.units = 'cm**-3'
        
        if not hasattr(nHe,'units'): 
            msg = '\n Input variable nHe must have units \n'
            raise ValueError(msg)
        else:
            nHe.units = 'cm**-3'

        if nH.shape != nHe.shape: 
            msg = '\n nH and nHe must have the same shape \n'
            raise ValueError(msg)



        # choose scalar or vector routines
        #----------------------------------------------
        nn = nH.size

        if nn == 1:
            solve = rabacus_fc.ion_solver.solve_pce_s
        else:
            solve = rabacus_fc.ion_solver.solve_pce_v

        (xH1, xH2, xHe1, xHe2, xHe3) = solve( 
            nH, nHe, 
            kchem.reH2, kchem.reHe2, kchem.reHe3, 
            kchem.ciH1, kchem.ciHe1, kchem.ciHe2, 
            kchem.H1i,  kchem.He1i,  kchem.He2i, 
            tol, nn )

        # package output
        #----------------------------------------------
        self.H1 = np.ones(nn) * xH1 * self.u.dimensionless
        self.H2 = np.ones(nn) * xH2 * self.u.dimensionless
        self.He1 = np.ones(nn) * xHe1 * self.u.dimensionless
        self.He2 = np.ones(nn) * xHe2 * self.u.dimensionless
        self.He3 = np.ones(nn) * xHe3 * self.u.dimensionless
        self.T = kchem.T

        self.nH = nH
        self.nHe = nHe

        self.kchem = kchem



class Solve_PCTE:

    r""" 
    Stores and calculates photo-collisional-thermal ionization and temperature
    equilibrium given the density of hydrogen and helium plus chemistry and 
    cooling rate objects. In other words, a self-consistent equilibrium 
    temperature and ionization state will be found for the input density and 
    redshift. This is an iterative solution.  All input arrays should be 1-D. 
    Note that the temperature arrays in `kchem` and `kcool` will be altered. 
    It is not important which temperatures are used to initialize the objects 
    kchem and kcool.  It IS important that you initialize them with 
    photoionization and photoheating rates.  Because we are searching for one 
    temperature for each density all input arrays must have the same shape.  

    .. note::

      The user must pass in photoionization rates for each species when 
      creating the chemistry rates object and photoheating rates for each 
      species when creating the cooling rates object (see examples below). 


    Args: 

      `nH` (array): number density of hydrogen

      `nHe` (array): number density of helium

      `kchem` (:class:`~rabacus.f2py.chem_cool_rates.ChemistryRates`): 
      chemistry rates

      `kcool` (:class:`~rabacus.f2py.chem_cool_rates.CoolingRates`): 
      cooling rates

      `z` (float): Redshift.  Important for Compton cooling

    Kwargs:

      `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
      desired.

      `tol` (float): convergence tolerance,
      abs(1 - `ne_new` / `ne_old`) < `tol`

    Attributes:

      `nH` (array): H number density.  

      `nHe` (array): He number density.  

      `ne` (array): electron number density 

      `T_floor` (float): minimum possible temperature 

      `H1` (array): neutral hydrogen fraction, |xH1| = |nH1| / |nH|

      `H2` (array): ionized hydrogen fraction, |xH2| = |nH2| / |nH|

      `He1` (array): neutral helium fraction,  |xHe1| = |nHe1| / |nHe|

      `He2` (array): singly ionized helium fraction,  |xHe2| = |nHe2| / |nHe|

      `He3` (array): doubly ionized helium fraction,  |xHe3| = |nHe3| / |nHe|

      `Teq` (array): equilibrium temperature 

      `cool` (array): cooling rate at final state (should be equal to `heat`)

      `heat` (array): heating rate at final state (should be equal to `cool`)
       
      `itr` (int): number of iterations to converge
       
      `Tcmb` (float): CMB temperature at input redshift


    Notes: 

      This class finds solutions (`H1`, `H2`, `He1`, `He2`, `He3`, `Teq`) to 
      the following equations 


      .. math:: 
        \frac{dx_{\rm _{HI}}}{dt} &= 
        - (\Gamma_{\rm _{HI}} + C_{\rm _{HI}} n_{\rm _e}) x_{\rm _{HI}}
        + R_{\rm _{HII}} n_{\rm _e} x_{\rm _{HII}} = 0 
        \\
        \frac{dx_{\rm _{HII}}}{dt} &= 
        (\Gamma_{\rm _{HI}} + C_{\rm _{HI}} n_{\rm _e}) x_{\rm _{HI}}       
        - R_{\rm _{HII}} n_{\rm _e} x_{\rm _{HII}} = 0
        \\
        \frac{dx_{\rm _{HeI}}}{dt} &= 
        - (\Gamma_{\rm _{HeI}} + C_{\rm _{HeI}} n_{\rm _e}) x_{\rm _{HeI}}
        + R_{\rm _{HeII}} n_{\rm _e} x_{\rm _{HeII}} = 0 
        \\
        \frac{dx_{\rm _{HeII}}}{dt} &= 
        (\Gamma_{\rm _{HeI}} + C_{\rm _{HeI}} n_{\rm _e}) 
        x_{\rm _{HeI}} - \\
        & \quad ( \Gamma_{\rm _{HeII}} + C_{\rm _{HeII}} n_{\rm _e} + 
        R_{\rm _{HeII}} n_{\rm _e} ) 
        x_{\rm _{HeII}} + \\
        & \quad R_{\rm _{HeIII}} n_{\rm _e} x_{\rm _{HeIII}} = 0 
        \\
        \frac{dx_{\rm _{HeIII}}}{dt} &= 
        (\Gamma_{\rm _{HeII}} + C_{\rm _{HeII}} n_{\rm _e}) 
        x_{\rm _{HeII}}
        - R_{\rm _{HeIII}} n_{\rm _e} x_{\rm _{HeIII}} = 0 
        \\
        \frac{du}{dt} &= \mathcal{H} - \Lambda_{\rm c} = 0


      with the following closure relationships

      .. math:: 
        1 &= x_{\rm _{HI}} + x_{\rm _{HII}} 
        \\
        1 &= x_{\rm _{HeI}} + x_{\rm _{HeII}} + x_{\rm _{HeIII}} 
        \\
        n_{\rm e} &= x_{\rm _{HII}} n_{\rm _{H}} + 
        ( x_{\rm _{HeII}} + 2 x_{\rm _{HeIII}} ) n_{\rm _{He}}
        \\
        u &= \frac{3}{2} ( n_{\rm _H} + n_{\rm _{He}} + n_{\rm _e} )
        k_{\rm b} T


      In the equations above, the :math:`\Gamma_{\rm _X}` are photoionization 
      rates, the :math:`C_{\rm _X}` are collisional ionization rates, the 
      :math:`R_{\rm _X}` are recombination rates, :math:`u` is the internal 
      energy, :math:`\mathcal{H}` is the heating function, and 
      :math:`\Lambda_{\rm c}` is the cooling function.  These solutions depend 
      on temperature and density. 



    Examples:
    
      The following code will solve for photo-collisional-thermal ionization 
      equilibrium at 128 densities using case B recombination rates for all 
      species.  The photoionization and heating rates are taken from the Haardt
      and Madau 2012 model.  Note that the photoionization rate arrays are 
      grouped with the chemistry rates and the photoheating rate arrays are 
      grouped with the cooling rates:: 

        import numpy as np
        import rabacus as ra

        N = 128; Yp = 0.24; z=3.0
        fcA_H2 = 0.0; fcA_He2 = 0.0; fcA_He3 = 0.0

        nH = 10**np.linspace( -4.0, -3.0, N ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = 10**np.linspace( 4.0, 5.0, N ) * ra.U.K

        hmr = ra.uv_bgnd.HM12_Photorates_Table()

        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)
        kchem = ra.f2py.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1i=H1i, He1i=He1i, He2i=He2i )

        H1h = np.ones(N) * hmr.H1h(z)
        He1h = np.ones(N) * hmr.He1h(z)
        He2h = np.ones(N) * hmr.He2h(z)
        kcool = ra.f2py.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3,
                                      H1h=H1h, He1h=He1h, He2h=He2h )

        x = ra.f2py.Solve_PCTE( nH, nHe, kchem, kcool, z )


    .. seealso::

      :class:`Solve_CE`, :class:`Solve_PCE`


    """ 


    def __init__(self, nH, nHe, kchem, kcool, z, Hz=None, tol=1.0e-8):

        # attach keywords
        #----------------------------------------------
        if Hz == None:
            Hz = 0.0
        else:
            if not hasattr(Hz,'units'): 
                raise ValueError, '\n Hz must have units \n'
            else:
                Hz.units = '1/s'

        self.Hz = Hz
        self.tol = tol

        # attach units. 
        #----------------------------------------------
        self.u = units.Units()

        # assert 1-D input
        #----------------------------------------------
        if len(nH.shape) != 1:
            msg = "input arrays must be 1-D for fortran routines." 
            raise ValueError(msg)

        if len(nHe.shape) != 1:
            msg = "input arrays must be 1-D for fortran routines." 
            raise ValueError(msg)

        if not hasattr(nH,'units'):
            msg = '\n Input variable nH must have units \n'
            raise ValueError(msg) 
        else:
            nH.units = 'cm**-3'
        
        if not hasattr(nHe,'units'): 
            msg = '\n Input variable nHe must have units \n'
            raise ValueError(msg)
        else:
            nHe.units = 'cm**-3'

        if len(kchem.T.shape) != 1:
            msg = "input arrays must be 1-D for fortran routines." 
            raise ValueError(msg)

        if len(kcool.T.shape) != 1:
            msg = "input arrays must be 1-D for fortran routines." 
            raise ValueError(msg)

        # set fortran variables
        #-----------------------------------------------
        if kchem.fit_name == 'hg97':
            irate = 1

        # choose scalar or vector routines
        #----------------------------------------------
        nn = nH.size

        if nn == 1:
            solve = rabacus_fc.ion_solver.solve_pcte_s
        else:
            solve = rabacus_fc.ion_solver.solve_pcte_v

        (xH1, xH2, xHe1, xHe2, xHe3, T) = solve( 
            nH, nHe, 
            kchem.H1i,  kchem.He1i,  kchem.He2i, 
            kcool.H1h,  kcool.He1h,  kcool.He2h,
            z, self.Hz, 
            kchem.fcA_H2, kchem.fcA_He2, kchem.fcA_He3, 
            irate, tol, nn )


        # package output
        #----------------------------------------------
        self.H1 = np.ones(nn) * xH1 * self.u.dimensionless
        self.H2 = np.ones(nn) * xH2 * self.u.dimensionless
        self.He1 = np.ones(nn) * xHe1 * self.u.dimensionless
        self.He2 = np.ones(nn) * xHe2 * self.u.dimensionless
        self.He3 = np.ones(nn) * xHe3 * self.u.dimensionless
        self.Teq = T * self.u.K

        self.nH = nH
        self.nHe = nHe

        # get chem and cooling rates at Teq
        #----------------------------------------------
        kchem_out = chemistry.ChemistryRates( 
            self.Teq,
            kchem.fcA_H2,
            kchem.fcA_He2,
            kchem.fcA_He3,
            H1i = kchem.H1i,
            He1i = kchem.He1i,
            He2i = kchem.He2i, 
            fit_name = kchem.fit_name,
            add_He2di = kchem.add_He2di  
            )
        
        kcool_out = cooling.CoolingRates( 
            self.Teq,
            kcool.fcA_H2,
            kcool.fcA_He2,
            kcool.fcA_He3, 
            H1h = kcool.H1h,
            He1h = kcool.He1h,
            He2h = kcool.He2h, 
            fit_name = kcool.fit_name,
            add_He2di = kcool.add_He2di  
            )

        self.kchem = kchem_out
        self.kcool = kcool_out
