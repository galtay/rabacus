""" Single zone ionization and temperature solvers for primordial gas.  """ 

import sys
import numpy as np

from rabacus.atomic import chemistry 
from rabacus.atomic import cooling
from rabacus.atomic import hydrogen
from rabacus.constants import units
from rabacus.utils import utils

__all__ = ['Solve_CE', 'Solve_PCE', 'Solve_PCTE']



MAX_ITR = 5000
TFLOOR = 7.5e3
TMAX = 1.0e5



class Solve_CE:

    r""" 
    Calculates and stores collisional ionization equilibrium (i.e. all 
    photoionization rates are equal to zero) given a chemistry rates object. 
    This is an exact analytic solution. 

    Args:
    
      `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
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
        x_ce = ra.solvers.Solve_CE( kchem )



    .. seealso::

       :class:`Solve_PCE`, :class:`Solve_PCTE`


    """ 

    def __init__(self, kchem, fpc=True):

        # attach optional arguments
        #---------------------------------------------
        self.fpc = fpc

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

          `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
          Chemistry rates

        """ 

        # attach rates and copy temperature 
        #----------------------------------------------
        self.kchem = kchem
        self.T = kchem.T.copy()

        # solve hydrogen
        #----------------------------------------------
        DH = kchem.reH2 + kchem.ciH1
        self.H1 = kchem.reH2 / DH
        self.H2 = kchem.ciH1 / DH
        # floating point error correction
        if self.fpc:
            s = self.H1 + self.H2
            self.H1 *= 1.0/s
            self.H2 *= 1.0/s

        # solve helium
        #----------------------------------------------
        DHe = kchem.reHe2 * kchem.reHe3 + \
              kchem.reHe3 * kchem.ciHe1 + \
              kchem.ciHe1 * kchem.ciHe2
        self.He1 = kchem.reHe2 * kchem.reHe3 / DHe
        self.He2 = kchem.reHe3 * kchem.ciHe1 / DHe
        self.He3 = kchem.ciHe1 * kchem.ciHe2 / DHe
        # floating point error correction
        if self.fpc:
            s = self.He1 + self.He2 + self.He3
            self.He1 *= 1.0/s
            self.He2 *= 1.0/s
            self.He3 *= 1.0/s


        




class Solve_PCE:

    r""" 
    Stores and calculates photo-collisional ionization equilibrium given the
    density of hydrogen and helium plus a chemistry rates object. This is an 
    iterative solution.  The shape of the output depends on the shape of the 
    input arrays.  There are two important numbers, Nrho and NT.  Nrho is the 
    size of the density arrays `nH_in` and `nHe_in` (both of which must be the 
    same).  NT is the size of the temperature array in `kchem`.  If 
    `two_dim_output` = ``False``, the density and temperature arrays must have 
    the same size and the output arrays are 1-D of that size.  If 
    `two_dim_output` = ``True``, all the main output arrays will have shape 
    (Nrho,NT).  In this case, an ionization state solution will be produced 
    for each possible combination of input density and temperature.  

    .. note::

      The user must pass in photoionization rates for each species when 
      creating the chemistry rates object (see examples below). 


    Args:

      `nH_in` (array): number density of hydrogen

      `nHe_in` (array): number density of helium

      `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): chemistry 
      rates

    Kwargs:

      `fpc` (bool): If ``True``, perform floating point error correction. 

      `verbose` (bool): If ``True``, produce verbose output
       
      `two_dim_output` (bool): If ``True``, produce a solution for every 
      possible density-temperature pair 

      `tol` (float): convergence tolerance,
      abs(1 - `ne_new` / `ne_old`) < `tol`
       

    Attributes: 

      `nH` (array): H number density.  Matches output array shape (i.e. 
      equal to `nH_in` if `two_dim_output` = False.)

      `nHe` (array): He number density.  Matches output array shape. (i.e. 
      equal to `nHe_in` if `two_dim_output` = False.)

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
    
      The 1-D output is the default setting.  The following code will solve for 
      photo-collisional ionization equilibrium at 128 density-temperature 
      pairs using recombination rates half way between case A and case B for 
      all species.  Note that the photoionization rate arrays are grouped with
      the chemistry rates and so should always be the same size as the 
      temperature array:: 

        import numpy as np
        import rabacus as ra

        fcA_H2 = 0.5; fcA_He2 = 0.5; fcA_He3 = 0.5

        N = 128; Nrho = NT = N; Yp = 0.24

        nH = np.linspace( 1.0e-4, 1.0e-3, Nrho ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = np.linspace( 1.0e4, 1.0e5, NT ) * ra.U.K

        H1i = np.ones( NT ) * 8.22e-13 / ra.U.s
        He1i = np.ones( NT ) * 4.76e-13 / ra.U.s
        He2i = np.ones( NT ) * 3.51e-15 / ra.U.s

        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        x_pce_1D = ra.solvers.Solve_PCE( nH, nHe, kchem )

      If the `two_dim_output` keyword is set to True, a solution will be 
      produced for all possible combinations of density and temperature.  In 
      the following example, the output arrays will have shape (Nrho,NT)::

        import numpy as np
        import rabacus as ra

        fcA_H2 = 0.5; fcA_He2 = 0.5; fcA_He3 = 0.5

        Nrho = 128; NT = 32; Yp = 0.24

        nH = np.linspace( 1.0e-4, 1.0e-3, Nrho ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = np.linspace( 1.0e4, 1.0e5, NT ) * ra.U.K

        H1i = np.ones( NT ) * 8.22e-13 / ra.U.s
        He1i = np.ones( NT ) * 4.76e-13 / ra.U.s
        He2i = np.ones( NT ) * 3.51e-15 / ra.U.s

        kchem = ra.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                   H1i=H1i, He1i=He1i, He2i=He2i )

        x_pce_2D = ra.solvers.Solve_PCE( nH, nHe, kchem, 
                                         two_dim_output=True )

    .. seealso::

      :class:`Solve_CE`, :class:`Solve_PCTE`


    """ 

    def __init__( self, nH_in, nHe_in, kchem, fpc=True, verbose=False, 
                  two_dim_output=False, tol=1.0e-10 ):

        # attach optional values
        #----------------------------------------------
        self.fpc = fpc
        self.verbose = verbose
        self.two_dim_output = two_dim_output
        self.tol = tol

        # attach hydrogen instance
        #----------------------------------------------
        self.Hatom = hydrogen.Hydrogen()


        # check temperature floor 
        #----------------------------------------------
        self.T_floor = TFLOOR * kchem.T.units
        if np.any( kchem.T < self.T_floor ):
            msg = 'No temperatures can be below T_floor = ' + \
                str(TFLOOR) + ' K \n' + 'kchem.T.min() = ' + \
                str(kchem.T.min())
            raise utils.InputError, msg

        # solve
        #----------------------------------------------
        self.set( nH_in, nHe_in, kchem )


    def set(self, nH_in, nHe_in, kchem ):

        """ 
        This function is called when a new instance of :class:`Solve_PCE` is 
        created.  Updated ionization fractions can be calculated by calling
        it again with new densities and a new kchem object. 

        Args:

          `nH_in` (array): number density of hydrogen

          `nHe_in` (array): number density of helium

          `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
          chemistry rates

        """ 


        # check input
        #----------------------------------------------
        if nH_in.shape == ():
            nH = np.ones(1) * nH_in
        else:
            nH = nH_in.copy()

        if nHe_in.shape == ():
            nHe = np.ones(1) * nHe_in
        else:
            nHe = nHe_in.copy()

        if not len(nH.shape) == len(nHe.shape) == len(kchem.T.shape) == 1:
            msg = '\n Input arrays must be one dimensional'
            raise utils.InputError, msg

        if not hasattr(nH,'units'):
            msg = '\n Input variable nH must have units \n'
            raise utils.NeedUnitsError, msg

        if not hasattr(nHe,'units'): 
            msg = '\n Input variable nHe must have units \n'
            raise utils.NeedUnitsError, msg

        if nH.shape != nHe.shape: 
            msg = '\n nH and nHe must have the same shape \n'
            raise utils.InputError, msg

        if self.two_dim_output == False and nH.shape != kchem.T.shape:
            msg = '\n For two_dim_output=False, nH, nHe, and kchem.T ' + \
                'must have the same shape \n'
            print 'nH.shape: ', nH.shape
            print 'nHe.shape: ', nHe.shape
            print 'kchem.T.shape: ', kchem.T.shape
            raise utils.InputError, msg


        # attach input
        #------------------------------------------
        self.nH_input = nH.copy()
        self.nHe_input = nHe.copy()
        self.kchem = kchem

        
        # screen report
        #--------------------------------------
        if self.verbose:
            print 'calculating photo/collisional ionization equilibrium' 
            print ' Nrho,NT: ', (nH.size,kchem.T.size)
            print ' two_dim_output: ', self.two_dim_output
            print ' nH min/max: ', nH.min(), nH.max()
            print ' T  min/max: ', kchem.T.min(), kchem.T.max()
            print 


        # create some output arrays
        #--------------------------------------
        if self.two_dim_output:

            self.nH = np.ones( [nH.size, kchem.T.size] ) * nH.units
            self.nHe = np.ones( [nH.size, kchem.T.size] ) * nHe.units
            self.T = np.ones( [nH.size, kchem.T.size] ) * kchem.T.units
            yy = np.zeros( [nH.size, kchem.T.size] ) 

            for i in range(kchem.T.size):
                self.nH[:,i] = nH
                self.nHe[:,i] = nHe

            for i in range(nH.size):
                self.T[i,:] = kchem.T

            self.kchem.extend_dims( nH.size )

            if self.kchem.T.shape != self.nH.shape:
                msg = '\n shapes dont match after kchem.extend \n'
                raise utils.InputError, msg

        else:

            self.nH = nH.copy()
            self.nHe = nHe.copy()
            self.T = kchem.T.copy()
            yy = np.zeros( nH.size ) 


        # calculate coll. ion. eq. 
        #--------------------------------------
        self.x_ce = Solve_CE( self.kchem ) 

        yy = ( self.x_ce.He2 + 2.0 * self.x_ce.He3 ) * self.nHe / self.nH   


        # calculate analytic H with zero heavy element electons
        #--------------------------------------
        self.xH1_ana = \
            self.Hatom.analytic_soltn_xH1( self.nH, yy, self.kchem )

        self.ne_ana = (1.0 - self.xH1_ana) * self.nH + \
            (self.x_ce.He2 + 2.0 * self.x_ce.He3) * self.nHe


        # calculate the minimum and maximum ne
        # (one for every temperature and nH combo)
        # ne_min has dimensions [nH.size, x_ce.kchem.T.size]
        #--------------------------------------
        self.ne_min = self.x_ce.H2 * self.nH + \
            (self.x_ce.He2 + 2.0 * self.x_ce.He3) * self.nHe
        self.ne_max = self.nH + 2.0 * self.nHe

        self.ne_left = self.ne_min.copy()
        self.ne_right = self.ne_max.copy()

        self.x_ce.ne = self.ne_min.copy()

        # choose an initial guess for ne
        #--------------------------------------
        self.ne = self.ne_ana.copy()

        ne_old = self.ne.copy()

        x_ne = ImplicitIonfracs( self.ne, self.nH, self.nHe, self.kchem )
        
        err = np.ones( self.ne.shape ) * 1.0e20
        itr = 0

        # iterate until converged
        #--------------------------------------
        while np.any( err > self.tol ):
            
            # calculate new ionization fractions with ne
            #---------------------------------------------
            x_ne.set( self.ne, self.nH, self.nHe, self.kchem )

            # calculate new ne from ionization fractions
            #---------------------------------------------
            self.ne = x_ne.H2 * self.nH + \
                (x_ne.He2 + 2.0 * x_ne.He3) * self.nHe

            # dummy check that it is between min and max values
            #---------------------------------------------
            if np.any( (self.ne / self.ne_max) > (1.0 + self.tol) ):
                print 'ne > ne_max !!!'
                indx = np.where( self.ne > self.ne_max )
                for i in indx[0]:
                    print 'ne,ne_max: ', self.ne[i], self.ne_max[i]
                sys.exit(1)

            if np.any( (self.ne / self.ne_min) < (1.0 - self.tol) ):
                print 'ne < ne_min !!!'
                indx = np.where( self.ne < self.ne_min )
                for i in indx[0]:
                    print 'ne,ne_min: ', self.ne[i], self.ne_min[i]
                sys.exit(1)

            # calculate new errors and store value in ne_old
            #---------------------------------------------
            ne_ratio = self.ne / ne_old


            if np.any( ne_ratio > 1.0 ):
                indx = np.where( ne_ratio > 1.0 )
                self.ne[indx] = ( ne_old[indx] + self.ne_right[indx] ) * 0.5
                self.ne_left[indx] = ne_old[indx]

            if np.any( ne_ratio < 1.0 ): 
                indx = np.where( ne_ratio < 1.0 )
                self.ne[indx] = ( self.ne_left[indx] + ne_old[indx] ) * 0.5
                self.ne_right[indx] = ne_old[indx]



            err = np.abs( 1 - ne_ratio )


            # check number of iterations
            #---------------------------------------------
            itr += 1
            self.itr = itr

            if itr > MAX_ITR:
                txt = ' Solve_PCE has not converged after ' 
                txt += str(MAX_ITR) + ' iterations'
                print txt

                print 'xH1: ', x_ne.H1
                print 'xH2: ', x_ne.H2
                print 'ne: ', self.ne
                print 'ne_old: ', ne_old
                print 'ne_ratio: ', ne_ratio
                print 'err: ', err
                print 'xH1 ce: ', self.x_ce.H1
                print 'xH2 ce: ', self.x_ce.H2

                sys.exit(1)

            ne_old = self.ne.copy()


        if self.two_dim_output == True:
            self.kchem.reduce_dims()

        self.H1 = x_ne.H1
        self.H2 = x_ne.H2
        self.He1 = x_ne.He1
        self.He2 = x_ne.He2
        self.He3 = x_ne.He3






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

      `nH_in` (array): number density of hydrogen

      `nHe_in` (array): number density of helium

      `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
      chemistry rates

      `kcool` (:class:`~rabacus.atomic.cooling.CoolingRates`): 
      cooling rates

      `z` (float): Redshift.  Important for Compton cooling

    Kwargs:

      `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
      desired.

      `sd93_logZ` (float): set to log metallicity (log Z) to use metal 
      cooling in CIE from Sutherland and Dopita 1993

      `fpc` (bool): If ``True``, perform floating point error correction. 

      `verbose` (bool): If ``True``, produce verbose output

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
      grouped with the chemistry rates and the photo-heating rate arrays are 
      grouped with the cooling rates:: 

        import numpy as np
        import rabacus as ra

        N = 128; Yp = 0.24; z=3.0
        fcA_H2 = fcA_He2 = fcA_He3 = 0.0

        nH = np.linspace( 1.0e-4, 1.0e-3, N ) / ra.U.cm**3
        nHe = nH * 0.25 * Yp / (1-Yp)
        T = np.linspace( 1.0e4, 1.0e5, N ) * ra.U.K

        hmr = ra.uv_bgnd.HM12_Photorates_Table()

        H1i = np.ones(N) * hmr.H1i(z)
        He1i = np.ones(N) * hmr.He1i(z)
        He2i = np.ones(N) * hmr.He2i(z)
        kchem = ra.atomic.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3,
                                          H1i=H1i, He1i=He1i, He2i=He2i )

        H1h = np.ones(N) * hmr.H1h(z)
        He1h = np.ones(N) * hmr.He1h(z)
        He2h = np.ones(N) * hmr.He2h(z)
        kcool = ra.atomic.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3,
                                        H1h=H1h, He1h=He1h, He2h=He2h )

        x_pcte = ra.solvers.Solve_PCTE( nH, nHe, kchem, kcool, z )



    .. seealso::

      :class:`Solve_CE`, :class:`Solve_PCE`


    """ 

    def __init__( self, nH_in, nHe_in, kchem, kcool, z, 
                  Hz=None, sd93_logZ=None, fpc=True, verbose=False, 
                  tol=1.0e-10 ):

        # attach optional arguments
        #----------------------------------------------
        self.Hz = Hz
        self.sd93_logZ = sd93_logZ
        self.fpc = fpc
        self.verbose = verbose        
        self.T_floor = TFLOOR * kchem.T.units
        self.tol = tol

        # solve
        #----------------------------------------------
        self.set( nH_in, nHe_in, kchem, kcool, z )


        
    def set( self, nH_in, nHe_in, kchem, kcool, z ):

        """ 
        This function is called when a new instance of :class:`Solve_PCTE` is 
        created.  Updated solutions can be calculated by calling it again with 
        new densities and new chemistry and cooling objects. 

        Args: 
    
          `nH_in` (array): number density of hydrogen

          `nHe_in` (array): number density of helium

          `kchem` (:class:`~rabacus.atomic.chemistry.ChemistryRates`): 
          chemistry rates

          `kcool` (:class:`~rabacus.atomic.cooling.CoolingRates`): 
          cooling rates

          `z` (float): Redshift.  Important for Compton and Bremsstrahlung 
          cooling
        
        """ 

        # check nH input
        #----------------------------------------------
        if nH_in.shape == ():
            nH = np.ones(1) * nH_in
        else:
            nH = nH_in.copy()

        if not hasattr(nH,'units'): 
            msg = '\n Input variable nH must have units \n'
            raise utils.NeedUnitsError, msg

        # check nHe input
        #----------------------------------------------
        if nHe_in.shape == ():
            nHe = np.ones(1) * nHe_in
        else:
            nHe = nHe_in.copy()

        if not hasattr(nHe,'units'): 
            msg = '\n Input variable nHe must have units \n'
            raise utils.NeedUnitsError, msg

        # check all arrays are 1-D 
        #----------------------------------------------
        if not len(nH.shape) == \
               len(nHe.shape) == \
               len(kchem.T.shape) == \
               len(kchem.H1i.shape) == \
               len(kchem.He1i.shape) == \
               len(kchem.He2i.shape) == \
               len(kcool.T.shape) == \
               len(kcool.H1h.shape) == \
               len(kcool.He1h.shape) == \
               len(kcool.He2h.shape) == 1:
            msg = '\n Input arrays must be one dimensional'
            raise utils.InputError, msg

        # check consistency of sizes
        #----------------------------------------------
        if not nH.shape == nHe.shape == kchem.T.shape == kcool.T.shape == \
                kchem.H1i.shape == kchem.He1i.shape == kchem.He2i.shape == \
                kcool.H1h.shape == kcool.He1h.shape == kcool.He2h.shape:

            msg = '\n nH, nHe, kchem.T, kcool.T, \n' +\
                'kchem.H1i, kchem.He1i, kchem.He2i, \n' +\
                'kcool.H1h, kcool.He1h, kcool.He2h \n' +\
                'must all have the same shape and size \n'

            print 'nH/nHe.shape: ', nH.shape, nHe.shape
            print 'kchem/kcool.T.shape: ', kchem.T.shape, kcool.T.shape
            print 'kchem.H1i/He1i/He2i.shape: ', kchem.H1i.shape, \
                kchem.He1i.shape, kchem.He2i.shape
            print 'kcool.H1h/He1h/He2h.shape: ', kcool.H1h.shape, \
                kcool.He1h.shape, kcool.He2h.shape

            raise utils.InputError, msg

        # make copies of input
        #--------------------------------------
        self.nH = nH.copy()
        self.nHe = nHe.copy()
        self.z = z
        self.Tcmb = 2.725 * (1.0+z) * kchem.T.units


        # bracket the equilibrium temperature.  Teq is bracketed if 
        # dudT at Tmin is positive and dudT and Tmax is negative
        #
        # Note that for some density / rates combinations the equilibrium
        # temperature will be below the temperature floor (and dudT at Tmin
        # will be negative).  For these entries we simply return the 
        # equilibrium ionization state at the temperature floor.  This will 
        # only happen in the case of low temperatures and low photoionization
        # and/or heating rates. 
        #================================================================


        # dudT for Tmin
        #------------------------------------------------------------
        Tmin = np.ones( kchem.T.shape ) * TFLOOR * kchem.T.units

        kchem_Tmin = chemistry.ChemistryRates( 
            Tmin, 
            kchem.fcA_H2,
            kchem.fcA_He2,
            kchem.fcA_He3,
            H1i = kchem.H1i.copy(), 
            He1i = kchem.He1i.copy(), 
            He2i = kchem.He2i.copy(),
            fit_name = kchem.fit_name,
            add_He2di = kchem.add_He2di  
            )
        
        kcool_Tmin = cooling.CoolingRates( 
            Tmin, 
            kcool.fcA_H2,
            kcool.fcA_He2,
            kcool.fcA_He3,
            H1h = kcool.H1h.copy(), 
            He1h = kcool.He1h.copy(), 
            He2h = kcool.He2h.copy(),
            fit_name = kcool.fit_name, 
            add_He2di = kcool.add_He2di  
            ) 
        
        x_eq = Solve_PCE( self.nH, self.nHe, kchem_Tmin )
        heat = kcool_Tmin.return_heating( self.nH, self.nHe, x_eq )
        cool = kcool_Tmin.return_cooling( 
            self.nH, self.nHe, x_eq, z, Hz=self.Hz, 
            sd93_logZ=self.sd93_logZ )
        dudTmin = heat - cool



        # dudT for Tmax
        #------------------------------------------------------------
        Tmax = np.ones( kchem.T.shape ) * TMAX * kchem.T.units

        kchem_Tmax = chemistry.ChemistryRates( 
            Tmax, 
            kchem.fcA_H2,
            kchem.fcA_He2,
            kchem.fcA_He3,
            H1i = kchem.H1i.copy(), 
            He1i = kchem.He1i.copy(), 
            He2i = kchem.He2i.copy(),
            fit_name = kchem.fit_name,
            add_He2di = kchem.add_He2di  
            )

        kcool_Tmax = cooling.CoolingRates( 
            Tmax, 
            kcool.fcA_H2,
            kcool.fcA_He2,
            kcool.fcA_He3,
            H1h = kcool.H1h.copy(), 
            He1h = kcool.He1h.copy(), 
            He2h = kcool.He2h.copy(), 
            fit_name = kcool.fit_name,
            add_He2di = kcool.add_He2di  
            ) 
        
        x_eq = Solve_PCE( self.nH, self.nHe, kchem_Tmax )
        heat = kcool_Tmax.return_heating( self.nH, self.nHe, x_eq )
        cool = kcool_Tmax.return_cooling( 
            self.nH, self.nHe, x_eq, z, Hz=self.Hz, 
            sd93_logZ=self.sd93_logZ )
        dudTmax = heat - cool


        # now we need to identify density-temperature pairs for which 
        # Teq is less than T_floor.  This is the case if dudT evaluated
        # at Tmin is less than zero (i.e. the gas would be cooling 
        # even at the minimum temperature.  For these pairs, we simply
        # fix the temperature at the minimum. 
        #---------------------------------------------------------

        b_floor = dudTmin < 0.0               # boolean
        i_floor = np.where( dudTmin < 0.0 )   # index
        n_floor = len( i_floor[0] )           # count


        # solve for PCE for the entries that will stay at T_floor
        #---------------------------------------------------------
        if np.any( b_floor == True ):

            kchem_f = chemistry.ChemistryRates( 
                Tmin[b_floor], 
                kchem.fcA_H2[b_floor],
                kchem.fcA_He2[b_floor],
                kchem.fcA_He3[b_floor],
                H1i = kchem.H1i[b_floor],
                He1i = kchem.He1i[b_floor],
                He2i = kchem.He2i[b_floor], 
                fit_name = kchem.fit_name,
                add_He2di = kchem.add_He2di  
                )
            
            x_f = Solve_PCE( self.nH[b_floor], self.nHe[b_floor], kchem_f )


        # now iterate to find equilibrium T for the entries that will 
        # not stay at T_floor
        #---------------------------------------------------------
        if np.any( b_floor == False ):

            # create sub-arrays
            #---------------------------------------------------------
            Tmin = Tmin[~b_floor]
            Tmax = Tmax[~b_floor]

            Tleft = Tmin.copy()
            Tright = Tmax.copy()

            dudTmin = dudTmin[~b_floor]
            dudTmax = dudTmax[~b_floor]

            nH = self.nH[~b_floor]
            nHe = self.nHe[~b_floor]

            if  np.any( dudTmin < 0 ) or np.any( dudTmax > 0 ):
                print ' !!!! Teq not bracketed !!!! '
                return
                sys.exit(1)
                
            # choose an initial guess for T
            #--------------------------------------
            logT = ( np.log10( Tmin.magnitude ) + 
                     np.log10( Tmax.magnitude ) ) / 2

            logT = np.ones( Tmin.shape ) * logT
            T = 10**logT * Tmin.units

            # set initial rates 
            #--------------------------------------
            kchem_Teq = chemistry.ChemistryRates( 
                T, 
                kchem.fcA_H2[~b_floor],
                kchem.fcA_He2[~b_floor],
                kchem.fcA_He3[~b_floor],
                H1i = kchem.H1i[~b_floor],
                He1i = kchem.He1i[~b_floor],
                He2i = kchem.He2i[~b_floor], 
                fit_name = kchem.fit_name,
                add_He2di = kchem.add_He2di  
                )
            
            kcool_Teq = cooling.CoolingRates( 
                T, 
                kcool.fcA_H2[~b_floor],
                kcool.fcA_He2[~b_floor],
                kcool.fcA_He3[~b_floor],
                H1h = kcool.H1h[~b_floor],
                He1h = kcool.He1h[~b_floor],
                He2h = kcool.He2h[~b_floor], 
                fit_name = kcool.fit_name,
                add_He2di = kcool.add_He2di  
                )


            err = np.ones( T.shape ) * 1.0e20
            itr = 0
            while np.any( err > self.tol ):

                # calculate dudT with new T
                #---------------------------------------------
                x_eq = Solve_PCE( nH, nHe, kchem_Teq )
                heat = kcool_Teq.return_heating( nH, nHe, x_eq )
                cool = kcool_Teq.return_cooling( 
                    nH, nHe, x_eq, z, Hz=self.Hz, 
                    sd93_logZ=self.sd93_logZ )
                dudT = heat - cool

                # if dudT is negative Teq is between Tmin and T
                # if dudT is positive Teq is between T and Tmax
                #---------------------------------------------
                Told = T.copy()

                indx_neg = np.where( dudT < 0 )[0]
                if len(indx_neg) > 0:
                    Tright[indx_neg] = T[indx_neg]
                    T[indx_neg] = ( Tleft[indx_neg] + T[indx_neg] ) / 2

                indx_pos = np.where( dudT >= 0)[0]
                if len(indx_pos) > 0:
                    Tleft[indx_pos] = T[indx_pos]
                    T[indx_pos] = ( T[indx_pos] + Tright[indx_pos] ) / 2

                # check err against tol
                #---------------------------------------------
                Tratio = T / Told
                err = np.abs( 1 - Tratio )


                # recalculate rates
                #---------------------------------------------
                kchem_Teq.set( T,
                               kchem.fcA_H2[~b_floor],
                               kchem.fcA_He2[~b_floor],
                               kchem.fcA_He3[~b_floor], 
                               H1i = kchem.H1i.copy()[~b_floor], 
                               He1i = kchem.He1i.copy()[~b_floor], 
                               He2i = kchem.He2i.copy()[~b_floor] )
                
                kcool_Teq.set( T, 
                               kcool.fcA_H2[~b_floor],
                               kcool.fcA_He2[~b_floor],
                               kcool.fcA_He3[~b_floor], 
                               H1h = kcool.H1h.copy()[~b_floor], 
                               He1h = kcool.He1h.copy()[~b_floor], 
                               He2h = kcool.He2h.copy()[~b_floor] )
                

                # check number of iterations
                #---------------------------------------------
                itr += 1
                self.itr = itr

                if itr > MAX_ITR:
                    txt = ' have not converged after ' 
                    txt += str(MAX_ITR) + ' iterations'
                    print txt
                    
                    print x_eq.H1, x_eq.H2
                    print x_eq.He1, x_eq.He2, x_eq.He3
                    print T, Told
                    print Tratio
                    print err
                    
                    
                    sys.exit(1)


        # create storage for output
        #-----------------------------------------------
        self.H1 = np.ones( self.nH.size )
        self.H2 = np.ones( self.nH.size )
        self.He1 = np.ones( self.nH.size )
        self.He2 = np.ones( self.nH.size )
        self.He3 = np.ones( self.nH.size )
        self.ne = np.ones( self.nH.size ) * self.nH.units
        self.Teq = np.ones( self.nH.size ) * kchem.T.units

        # assign values 
        #-----------------------------------------------
        if np.any( b_floor == False ):
            self.H1[~b_floor] = x_eq.H1
            self.H2[~b_floor] = x_eq.H2
            self.He1[~b_floor] = x_eq.He1
            self.He2[~b_floor] = x_eq.He2
            self.He3[~b_floor] = x_eq.He3
            self.ne[~b_floor] = x_eq.ne
            self.Teq[~b_floor] = T
            
        if np.any( b_floor == True ):
            self.H1[b_floor] = x_f.H1
            self.H2[b_floor] = x_f.H2
            self.He1[b_floor] = x_f.He1
            self.He2[b_floor] = x_f.He2
            self.He3[b_floor] = x_f.He3
            self.ne[b_floor] = x_f.ne
            self.Teq[b_floor] = x_f.T


        # create output chemistry and cool rates that match Teq
        #--------------------------------------------------------
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

        x_eq = Solve_PCE( self.nH, self.nHe, kchem_out )
        heat = kcool_out.return_heating( self.nH, self.nHe, x_eq )
        cool = kcool_out.return_cooling( 
            self.nH, self.nHe, x_eq, z, Hz=self.Hz, 
            sd93_logZ=self.sd93_logZ )
        dudT = heat - cool

        self.kchem = kchem_out
        self.kcool = kcool_out 
        self.heat = heat
        self.cool = cool
        self.dudT = dudT


class ImplicitIonfracs:

    """ 

    Stores and calculates ne dependent ionization equilibrium for six species 
    models. Input is, 
    
      ne:    number density of electrons
      nH:    number density of hydrogen
      nHe:   number density of helium 
      k:     an instance of ChemistryRates

    """ 


    def __init__(self, ne, nH, nHe, kchem, fpc=True):

        # calculate
        #----------------------------------------------
        self.set( ne, nH, nHe, kchem, fpc )


    def set(self, ne, nH, nHe, kchem, fpc=True):

        # check input
        #----------------------------------------------
        if not hasattr(ne,'units'): 
            raise NeedUnitsError, '\n Input variable ne must have units \n'

        if not hasattr(nH,'units'): 
            raise NeedUnitsError, '\n Input variable nH must have units \n'

        if not hasattr(nHe,'units'): 
            raise NeedUnitsError, '\n Input variable nHe must have units \n'

        msg = '\n ne, nH, nHe, and rates arrays must have the same shape \n'
        if not ne.shape == nH.shape == nHe.shape == kchem.T.shape: 
            print ne.shape, nH.shape, nHe.shape, kchem.T.shape
            raise InputError, msg


        # attach input
        #----------------------------------------------
        self.ne = ne
        self.nH = nH
        self.nHe = nHe
        self.kchem = kchem

        # solve hydrogen
        #----------------------------------------------
        R2 = kchem.reH2 * ne
        CG1 = kchem.ciH1 * ne + kchem.H1i

        DH = R2 + CG1
        self.H1 = R2 / DH
        self.H2 = CG1 / DH
        # floating point error correction
        if fpc:
            s = self.H1 + self.H2
            self.H1 *= 1.0/s
            self.H2 *= 1.0/s

        # solve helium
        #----------------------------------------------
        R2 = kchem.reHe2 * ne
        R3 = kchem.reHe3 * ne
        CG1 = kchem.ciHe1 * ne + kchem.He1i
        CG2 = kchem.ciHe2 * ne + kchem.He2i

        DHe = (R2 * R3) + (R3 * CG1) + (CG1 * CG2)
        self.He1 = R2  * R3  / DHe
        self.He2 = R3  * CG1 / DHe
        self.He3 = CG1 * CG2 / DHe
        # floating point error correction
        if fpc:
            s = self.He1 + self.He2 + self.He3
            self.He1 *= 1.0/s
            self.He2 *= 1.0/s
            self.He3 *= 1.0/s






