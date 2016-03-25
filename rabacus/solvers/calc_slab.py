""" Solves a slab with plane parallel radiation incident from one side. """

import numpy as np

import one_zone as ozn
from rabacus.atomic import chemistry
from rabacus.atomic import cooling
from rabacus.constants import units
from rabacus.utils import utils



__all__ = ['Slab']


class Slab:

    """ 

    Solves a slab geometry with plane parallel radiation incident from one 
    side.  The slab is divided into `Nl` layers. 

    Args: 
      `Edges` (array): Distance to layer edges 

      `T` (array): Temperature in each layer 

      `nH` (array): hydrogen number density in each layer 

      `nHe` (array): helium number density in each layer 

      `rad_src` (:class:`~rabacus.rad_src.plane.PlaneSource`): Plane source

    .. note::
      The arrays `T`, `nH`, and `nHe` must all be the same size.  This size
      determines the number of shells, `Nl`.  `Edges` determines the positions
      of the layer edges and must have `Nl+1` entries. 

    Kwargs:

      `rec_meth` (string): How to treat recombinations {``fixed``, ``thresh``}

      `fixed_fcA` (float): If `rec_meth` = ``fixed``, constant caseA fraction 

      `thresh_P_A` (float): If `rec_meth` = ``thresh``, this is the probability 
      of absorption, P = 1 - exp(-tau), at which the transition to caseB rates 
      begins. 

      `thresh_P_B` (float): If `rec_meth` = ``thresh``, this is the probability 
      of absorption, P = 1 - exp(-tau), at which the transition to caseB ends.

      `thresh_xmfp` (float): Determines the distance to probe when calculating
      optical depth for the caseA to caseB transition.  Specifies a multiple
      of the mean free path for fully neutral gas.  In equation form, 
      L = thresh_xmfp / (n * sigma)
      
      `atomic_fit_name` (string): Source for atomic rate fits {``hg97``}

      `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve for 
      equilibrium T 

      `z` (float): Redshift, only need if `find_Teq` = ``True``

      `verbose` (bool): Verbose output? 

      `tol` (float): tolerance for all convergence tests

      `thin` (bool): if ``True`` only solves optically thin
     
    Attributes:

       `U` (:class:`~rabacus.constants.units.Units`)

       `z_c` (array): distance from surface of slab to center of layer

       `dz` (array): thickness of layer

       `NH_c` (array): H column density from surface of slab to center of layer 

       `dNH` (array): H column density through layer

       `NH1_c` (array): HI column density from surface of slab to center of 
       layer

       `dNH1` (array): HI column density through layer

       `H1i` (array): HI photo-ionization rate  

       `H1h` (array): HI photo-heating rate  

       `xH1` (array): H neutral fraction nHI / nH
       
       `xH2` (array): H ionized fraction nHII / nH

       `ne` (array): electron number density 

       `fcA_H2` (array): HII case A fraction 
 
       `cool` (array): cooling rate [erg / (cm^3 K)]

       `heat` (array): heating rate [erg / (cm^3 K)]

       `heatH1` (array): contribution to `heat` from H1 photoheating

       `dtauH1_th` (array): HI optical depth at H1 ionizing threshold
       through layer 

       `tauH1_th_lo` (array): HI optical depth below this layer

       `tauH1_th_hi` (array): HI optical depth above this layer 

       `NH1_thru` (float): HI column density through entire slab
 
       `Nl` (int): Number of layers

       `itr` (int): Number of iterations to converge


    .. note::

      For many of the attributes above, there are analagous versions for 
      helium.  We also note that the source of the atomic rates fit is stored
      in the variable `atomic_fit_name`, but the source of the photoionization
      cross section fits are stored in the point source object in the variable
      `px_fit_type`. 


    Examples:    

      The following code will solve a monochromatic source incident 
      onto a slab with fixed recombination rates and temperatues. In this
      case we use a Haardt and Madau 2012 spectrum to calculate a grey photon 
      energy which is then used to create a monochromatic plane source :: 

        import numpy as np
        import rabacus as ra
        
        # create Haardt and Madau 2012 source
        #---------------------------------------------------------------
        q_min = 1.0; q_max = 4.0e2; z = 3.0
        hm12 = ra.rad_src.PlaneSource( q_min, q_max, 'hm12', z=z )
 
        # get grey photon energy and create a monochromatic plane source
        #---------------------------------------------------------------
        q_mono = hm12.grey.E.He2 / hm12.PX.Eth_H1
        mono = ra.rad_src.PlaneSource( q_mono, q_mono, 'monochromatic' )
        mono.normalize_H1i( hm12.thin.H1i )

        # describe slab
        #---------------------------------------------------------------
        Nl = 512; Yp = 0.24
        nH = np.ones(Nl) * 1.5e-2 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)
        Tslab = np.ones(Nl) * 1.0e4 * ra.U.K
        Lslab = 200.0 * ra.U.kpc
        Edges = np.linspace( 0.0 * ra.U.kpc, Lslab, Nl+1 )

        # solve slab
        #---------------------------------------------------------------
        slab = ra.solvers.Slab( Edges, Tslab, nH, nHe, mono,
                                rec_meth='fixed', fixed_fcA=1.0 )

    """ 

    def __init__(self, 
                 Edges, 
                 T, 
                 nH, 
                 nHe, 
                 rad_src,
#
                 rec_meth = "fixed", 
                 fixed_fcA = 1.0, 
                 thresh_P_A = 0.01, 
                 thresh_P_B = 0.99,
                 thresh_xmfp = 3.0e1, 
                 atomic_fit_name = "hg97",    
                 find_Teq = False,          
                 z = None,                 
                 verbose = False,
                 tol = 1.0e-10,
                 thin = False,
                 ):


        # attach units
        #-----------------------------------------------
        self.U = units.Units()

        # check input 
        #-----------------------------------------------
        assert( Edges.size == T.size + 1 )
        assert( T.shape == nH.shape == nHe.shape )
        assert( T.shape == (T.size,) )

        assert( rec_meth == 'fixed' or rec_meth == 'thresh' )

        if find_Teq and z==None:
            msg = 'need to supply redshift if finding equilibrium temperature'
            raise utils.InputError(msg)

        if rad_src.source_type != 'plane':
            msg = 'source type needs to be plane'
            raise utils.InputError(msg)


        # set units
        #-----------------------------------------------
        Edges.units = 'cm' 
        T.units = 'K'
        nH.units = '1.0/cm^3'
        nHe.units = '1.0/cm^3'


        # attach input
        #-----------------------------------------------
        self.Edges = Edges.copy()
        self.T = T.copy()
        self.nH = nH.copy()
        self.nHe = nHe.copy()
        self.rad_src = rad_src
        self.rec_meth = rec_meth 

        if self.rec_meth == 'fixed':
            self.fixed_fcA = fixed_fcA

        elif self.rec_meth == 'thresh':
            assert thresh_xmfp > 0.0
            assert thresh_P_B > thresh_P_A
            self.thresh_xmfp = thresh_xmfp
            self.thresh_P_A = thresh_P_A
            self.thresh_P_B = thresh_P_B
            self.thresh_dPAB = thresh_P_B - thresh_P_A
            self.thresh_tau_A = -np.log( 1.0 - thresh_P_A )
            self.thresh_tau_B = -np.log( 1.0 - thresh_P_B )

        self.atomic_fit_name = atomic_fit_name
        self.find_Teq = find_Teq
        self.verbose = verbose
        self.tol = tol
        self.thin = thin

        if find_Teq:
            self.z = z


        # initialize slab
        #-----------------------------------------------
        self.init_slab()
        self.set_optically_thin()

        if self.thin:
            return

        # solve slab (sweep until convergence)
        #-----------------------------------------------------------
        conv_old = np.sum( self.ne )
        not_converged = True
        self.itr = 0

        while not_converged:
            self.sweep_slab()
            conv_new = np.sum( self.ne )
            if np.abs( conv_new/conv_old - 1.0 ) < self.tol:
                not_converged = False
            conv_old = conv_new
            self.itr += 1


        # finalize slab
        #-----------------------------------------------------------
        self.finalize_slab()







    def init_slab( self ):
        """ Initialize slab values.  """

        if self.verbose:
            print 'begin initialization' 

        # instantiate atomic rates (optically thin)
        #-----------------------------------------------
        if self.rec_meth == 'fixed':
            fcA_H2 = self.fixed_fcA
            fcA_He2 = self.fixed_fcA
            fcA_He3 = self.fixed_fcA
        elif self.rec_meth == 'thresh':
            fcA_H2 = 1.0
            fcA_He2 = 1.0
            fcA_He3 = 1.0

        kchem = chemistry.ChemistryRates( 
            self.T[0], fcA_H2, fcA_He2, fcA_He3, 
            H1i = self.rad_src.thin.H1i, 
            He1i = self.rad_src.thin.He1i, 
            He2i = self.rad_src.thin.He2i,
            fit_name = self.atomic_fit_name
            )

        kcool = cooling.CoolingRates(
            self.T[0], fcA_H2, fcA_He2, fcA_He3, 
            H1h = self.rad_src.thin.H1h, 
            He1h = self.rad_src.thin.He1h, 
            He2h = self.rad_src.thin.He2h, 
            fit_name = self.atomic_fit_name
            ) 

        self.kchem = kchem
        self.kcool = kcool 

        if self.verbose:
            print '  created kchem and kcool' 

        # setup arrays
        #-----------------------------------------------
        self.Nl = self.nH.size        
        Nl = self.Nl

        self.dz = self.Edges[1:] - self.Edges[0:-1]
        self.z_c = self.Edges[0:-1] + 0.5 * self.dz 
                                                     
        self.xH1 = np.zeros(Nl)
        self.dNH1 = np.zeros(Nl) / self.U.cm**2
        self.dtauH1_th = np.zeros(Nl)
        self.tauH1_th_lo = np.zeros(Nl) # H1 optical depth below layer
        self.tauH1_th_hi = np.zeros(Nl) # H1 optical depth above layer 
        self.H1i = np.zeros(Nl) / self.U.s
        self.H1h = np.zeros(Nl) * self.U.eV / self.U.s
        self.fcA_H2 = np.zeros(Nl)

        self.xHe1 = np.zeros(Nl)
        self.dNHe1 = np.zeros(Nl) / self.U.cm**2
        self.dtauHe1_th = np.zeros(Nl)
        self.tauHe1_th_lo = np.zeros(Nl) # He1 optical depth below layer
        self.tauHe1_th_hi = np.zeros(Nl) # He1 optical depth above layer
        self.He1i = np.zeros(Nl) / self.U.s
        self.He1h = np.zeros(Nl) * self.U.eV / self.U.s
        self.fcA_He2 = np.zeros(Nl)

        self.xHe2 = np.zeros(Nl)
        self.dNHe2 = np.zeros(Nl) / self.U.cm**2
        self.dtauHe2_th = np.zeros(Nl)
        self.tauHe2_th_lo = np.zeros(Nl) # He2 optical depth below layer
        self.tauHe2_th_hi = np.zeros(Nl) # He2 optical depth above layer
        self.He2i = np.zeros(Nl) / self.U.s
        self.He2h = np.zeros(Nl) * self.U.eV / self.U.s
        self.fcA_He3 = np.zeros(Nl)

        self.cool = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heat = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatH1 = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatHe1 = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatHe2 = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)

        self.xH2 = np.zeros(Nl)
        self.xHe3 = np.zeros(Nl)
        self.ne = np.zeros(Nl) / self.U.cm**3

        if self.verbose:
            print '  created simple arrays' 

        # calc NH and NHe between S and center of each layer
        #----------------------------------------
        self.dNH = self.dz * self.nH 
        self.dNHe = self.dz * self.nHe

        self.NH_c = np.zeros(Nl) / self.U.cm**2
        self.NHe_c = np.zeros(Nl) / self.U.cm**2

        for i in xrange(Nl):

            self.NH_c[i] = np.sum( self.dNH[0:i] ) + 0.5 * self.dNH[i]
            self.NHe_c[i] = np.sum( self.dNHe[0:i] ) + 0.5 * self.dNHe[i]

        if self.verbose:
            print '  created NH_c and NHe_c' 

        # calc mean free path when fully neutral for each shell
        # also calculate the bounds for each shell
        #--------------------------------------------------------
        if self.rec_meth == 'thresh':

            L = self.Edges[-1]
            
            self.mfp_H1 = 1.0 / ( self.nH * self.rad_src.th.sigma_H1 )
            self.L_H1 = self.mfp_H1 * self.thresh_xmfp
            if np.any( self.L_H1 > L ):
                indx = np.where( self.L_H1 > L )
                self.L_H1[indx] = L

            self.mfp_He1 = 1.0 / ( self.nHe * self.rad_src.th.sigma_He1 )
            self.L_He1 = self.mfp_He1 * self.thresh_xmfp
            if np.any( self.L_He1 > L ):
                indx = np.where( self.L_He1 > L )
                self.L_He1[indx] = L

            self.mfp_He2 = 1.0 / ( self.nHe * self.rad_src.th.sigma_He2 )
            self.L_He2 = self.mfp_He2 * self.thresh_xmfp
            if np.any( self.L_He2 > L ):
                indx = np.where( self.L_He2 > L )
                self.L_He2[indx] = L


            self.ii_H1 = np.zeros(Nl, dtype=int)
            self.ff_H1 = np.zeros(Nl, dtype=int)
            self.ii_He1 = np.zeros(Nl, dtype=int)
            self.ff_He1 = np.zeros(Nl, dtype=int)
            self.ii_He2 = np.zeros(Nl, dtype=int)
            self.ff_He2 = np.zeros(Nl, dtype=int)

            for i in xrange(Nl):

                self.ii_H1[i],  self.ff_H1[i]  = \
                    self.return_bounds( i, self.L_H1[i] )
  
                self.ii_He1[i], self.ff_He1[i] = \
                    self.return_bounds( i, self.L_He1[i] )

                self.ii_He2[i], self.ff_He2[i] = \
                    self.return_bounds( i, self.L_He2[i] )

            if self.verbose:
                print '  created index arrays for fcA' 


        if self.verbose:
            print 'initialization complete' 



    def return_bounds( self, i, Ltarget ):
        """ Given a shell number and a distance, return the indices such that
        self.r[i] - self.r[ii] > Ltarget and self.r[ff] - self.r[i] > 
        Ltarget. """ 

        L = 0.0 * self.U.cm
        k = 1
        while L < Ltarget:
            ii = i-k
            if ii < 0:
                ii = 0
                L = 1.5 * Ltarget
            else:
                L = L + self.dz[ii]
                k += 1

        L = 0.0 * self.U.cm
        k = 1
        while L < Ltarget:
            ff = i+k
            if ff >= self.Nl:
                ff = self.Nl-1
                L = 1.5 * Ltarget
            else:
                L = L + self.dz[ff]
                k += 1

        return ii,ff



    def set_optically_thin( self ): 

        """ Set optically thin ionzation state """ 

        # initialize a chemistry object
        #----------------------------------------
        if self.rec_meth == 'fixed':
            self.fcA_H2 = np.ones(self.Nl) * self.fixed_fcA
            self.fcA_He2 = np.ones(self.Nl) * self.fixed_fcA
            self.fcA_He3 = np.ones(self.Nl) * self.fixed_fcA
                
        else:
            self.fcA_H2 = np.ones(self.Nl) 
            self.fcA_He2 = np.ones(self.Nl) 
            self.fcA_He3 = np.ones(self.Nl) 


        kchem = chemistry.ChemistryRates( 
            self.T, self.fcA_H2, self.fcA_He2, self.fcA_He3,
            H1i = np.ones(self.Nl) * self.rad_src.thin.H1i, 
            He1i = np.ones(self.Nl) * self.rad_src.thin.He1i, 
            He2i = np.ones(self.Nl) * self.rad_src.thin.He2i,
            fit_name = self.atomic_fit_name, 
            )

        kcool = cooling.CoolingRates( 
            self.T, self.fcA_H2, self.fcA_He2, self.fcA_He3,
            H1h = np.ones(self.Nl) * self.rad_src.thin.H1h, 
            He1h = np.ones(self.Nl) * self.rad_src.thin.He1h, 
            He2h = np.ones(self.Nl) * self.rad_src.thin.He2h,
            fit_name = self.atomic_fit_name, 
            )

                
        if self.find_Teq:
            x = ozn.Solve_PCTE( self.nH, self.nHe, kchem, kcool, self.z, 
                                self.tol )
            self.T = x.Teq
        else:
            x = ozn.Solve_PCE( self.nH, self.nHe, kchem, self.tol )
        
        self.xH1 = x.H1
        self.xH2 = x.H2
        self.xHe1 = x.He1
        self.xHe2 = x.He2
        self.xHe3 = x.He3
        
        # summarize
        #----------------------------------------
        self.dNH1 = self.dNH * self.xH1
        self.dNHe1 = self.dNHe * self.xHe1
        self.dNHe2 = self.dNHe * self.xHe2
        
        self.dtauH1_th = self.dNH1 * self.rad_src.th.sigma_H1 
        self.dtauHe1_th = self.dNHe1 * self.rad_src.th.sigma_He1 
        self.dtauHe2_th = self.dNHe2 * self.rad_src.th.sigma_He2 
        
        self.ne = self.xH2 * self.nH + \
            ( self.xHe2 + 2.0 * self.xHe3 ) * self.nHe 



    def set_taus( self, i ):

        """ Calculate optical depth at the ionization thresholds of each
        species above and below this layer.  Does not include a contribution 
        from the layer itself. 

        Args: 
          `i` (int): layer index
        """ 

        self.tauH1_th_lo[i] = np.sum( self.dtauH1_th[0:i] )
        self.tauH1_th_hi[i] = np.sum( self.dtauH1_th[i+1:self.Nl] )
        
        self.tauHe1_th_lo[i] = np.sum( self.dtauHe1_th[0:i] )
        self.tauHe1_th_hi[i] = np.sum( self.dtauHe1_th[i+1:self.Nl] )
        
        self.tauHe2_th_lo[i] = np.sum( self.dtauHe2_th[0:i] )
        self.tauHe2_th_hi[i] = np.sum( self.dtauHe2_th[i+1:self.Nl] )


    def return_fcA( self, tau_th ):
        """ Given the optical depth at nu_th for any given species, returns
        the case A fraction. """ 

        if tau_th < self.thresh_tau_A:
            fcA = 1.0
        elif tau_th > self.thresh_tau_B:
            fcA = 0.0
        else:
            P = 1.0 - np.exp( -tau_th )
            delta = P - self.thresh_P_A
            fcB = delta / self.thresh_dPAB
            fcA = 1.0 - fcB

        return fcA


    def set_fcA( self, i ):
        """ Calculate case A fraction for each species in this layer.
        
        Args: 
          `i` (int): layer index
        """ 
        if self.rec_meth == 'fixed':

            self.fcA_H2[i] = self.fixed_fcA
            self.fcA_He2[i] = self.fixed_fcA
            self.fcA_He3[i] = self.fixed_fcA

        elif self.rec_meth == 'thresh':

            ii = self.ii_H1[i]; ff = self.ff_H1[i]
            tau_lo = np.sum( self.dtauH1_th[ii:i] )
            fcA_lo = self.return_fcA( tau_lo )
            tau_hi = np.sum( self.dtauH1_th[i+1:ff+1] )
            fcA_hi = self.return_fcA( tau_hi )
            self.fcA_H2[i] = (fcA_lo + fcA_hi) * 0.5

            ii = self.ii_He1[i]; ff = self.ff_He1[i]
            tau_lo = np.sum( self.dtauHe1_th[ii:i] )
            fcA_lo = self.return_fcA( tau_lo )
            tau_hi = np.sum( self.dtauHe1_th[i+1:ff+1] )
            fcA_hi = self.return_fcA( tau_hi )
            self.fcA_He2[i] = (fcA_lo + fcA_hi) * 0.5

            ii = self.ii_He2[i]; ff = self.ff_He2[i]
            tau_lo = np.sum( self.dtauHe2_th[ii:i] )
            fcA_lo = self.return_fcA( tau_lo )
            tau_hi = np.sum( self.dtauHe2_th[i+1:ff+1] )
            fcA_hi = self.return_fcA( tau_hi )
            self.fcA_He3[i] = (fcA_lo + fcA_hi) * 0.5


    def set_photoion_rates( self, i ):
        """ set photoionization rates for this layer 

        Args: 
          `i` (int): layer index
          """ 

        tauH1_th = np.float( self.tauH1_th_lo[i] + 0.5 * self.dtauH1_th[i] )
        tauHe1_th = np.float( self.tauHe1_th_lo[i] + 0.5 * self.dtauHe1_th[i] )
        tauHe2_th = np.float( self.tauHe2_th_lo[i] + 0.5 * self.dtauHe2_th[i] )
        
        self.H1i[i] = self.rad_src.shld_H1i( tauH1_th, tauHe1_th, tauHe2_th )
        self.He1i[i] = self.rad_src.shld_He1i( tauH1_th, tauHe1_th, tauHe2_th ) 
        self.He2i[i] = self.rad_src.shld_He2i( tauH1_th, tauHe1_th, tauHe2_th )


    def set_photoheat_rates( self, i ):
        """ set photoheating rates for this layer

        Args: 
          `i` (int): layer index
          """ 

        tauH1_th = np.float( self.tauH1_th_lo[i] + 0.5 * self.dtauH1_th[i] )
        tauHe1_th = np.float( self.tauHe1_th_lo[i] + 0.5 * self.dtauHe1_th[i] )
        tauHe2_th = np.float( self.tauHe2_th_lo[i] + 0.5 * self.dtauHe2_th[i] )

        self.H1h[i] = self.rad_src.shld_H1h( tauH1_th, tauHe1_th, tauHe2_th )  
        self.He1h[i] = self.rad_src.shld_He1h( tauH1_th, tauHe1_th, tauHe2_th )
        self.He2h[i] = self.rad_src.shld_He2h( tauH1_th, tauHe1_th, tauHe2_th )




    def sweep_slab(self):
        """ Performs one sweep through slab. """ 


        for i in xrange(self.Nl):

            if self.verbose:
                print 'i = ', i

            # calc tauXX above and below this layer
            # (constant during iterations)
            #----------------------------------------
            self.set_taus( i )

            # calculate fcA (average over directions)
            # (constant during iterations)
            #----------------------------------------
            self.set_fcA( i )

            # iterate until we have convergence in the optical
            # depth through this layer
            #--------------------------------------------------
            conv_old = self.dtauH1_th[i] + \
                self.dtauHe1_th[i] + self.dtauHe2_th[i]
            not_converged = True
            itr = 0

            while not_converged:

                # calculate photoion / chemistry rates
                #----------------------------------------
                self.set_photoion_rates( i )
        
                self.kchem.set( self.T[i], 
                                self.fcA_H2[i], 
                                self.fcA_He2[i], 
                                self.fcA_He3[i], 
                                H1i = self.H1i[i], 
                                He1i = self.He1i[i], 
                                He2i = self.He2i[i] )
                
                # calculate photoheat / cooling rates
                #----------------------------------------
                self.set_photoheat_rates( i )
        
                self.kcool.set( self.T[i], 
                                self.fcA_H2[i], 
                                self.fcA_He2[i], 
                                self.fcA_He3[i], 
                                H1h = self.H1h[i], 
                                He1h = self.He1h[i], 
                                He2h = self.He2h[i] )

                # if we are finding the equilibrium temperature 
                # we need to call Solve_PCTE
                #----------------------------------------
                if self.find_Teq:

                    x = ozn.Solve_PCTE( np.ones(1) * self.nH[i], 
                                        np.ones(1) * self.nHe[i], 
                                        self.kchem, 
                                        self.kcool, 
                                        self.z, 
                                        self.tol )

                    self.T[i] = x.Teq

                # otherwise we call Solve_PCE
                #----------------------------------------
                else:

                    x = ozn.Solve_PCE( np.ones(1) * self.nH[i], 
                                       np.ones(1) * self.nHe[i], 
                                       self.kchem, self.tol )
                    

                # set the ionization fractions in the object
                #----------------------------------------            
                self.xH1[i] = x.H1
                self.xH2[i] = x.H2
                self.xHe1[i] = x.He1
                self.xHe2[i] = x.He2
                self.xHe3[i] = x.He3
                
                # calculate heating rates in layer
                #----------------------------------------
                self.heatH1[i] = self.H1h[i] * self.nH[i] * self.xH1[i] 
                self.heatHe1[i] = self.He1h[i] * self.nHe[i] * self.xHe1[i] 
                self.heatHe2[i] = self.He2h[i] * self.nHe[i] * self.xHe2[i] 
                self.heat[i] = self.heatH1[i] + \
                    self.heatHe1[i] + self.heatHe2[i]

                # calculate electron density in layer
                #----------------------------------------
                self.ne[i] = self.xH2[i] * self.nH[i] + \
                    ( self.xHe2[i] + 2 * self.xHe3[i] ) * self.nHe[i] 

                # calculate tauXX through layer
                #----------------------------------------
                self.dNH1[i] = self.dNH[i] * self.xH1[i]
                self.dNHe1[i] = self.dNHe[i] * self.xHe1[i]
                self.dNHe2[i] = self.dNHe[i] * self.xHe2[i]
            
                self.dtauH1_th[i] = self.dNH1[i] * self.rad_src.th.sigma_H1 
                self.dtauHe1_th[i] = self.dNHe1[i] * self.rad_src.th.sigma_He1 
                self.dtauHe2_th[i] = self.dNHe2[i] * self.rad_src.th.sigma_He2 

                # check convergence
                #----------------------------------------
                conv_new = self.dtauH1_th[i] + \
                    self.dtauHe1_th[i] + self.dtauHe2_th[i]
                if np.abs( conv_new/conv_old - 1.0 ) < self.tol:
                    not_converged = False
                conv_old = conv_new
                itr += 1
                



    def finalize_slab( self ):
        """ Calculate some final values using fully solved slab. """ 


        # calculate column densities up to each layer
        #-----------------------------------------------------------
        self.NH1_c = np.zeros(self.Nl) / self.U.cm**2
        self.NHe1_c = np.zeros(self.Nl) / self.U.cm**2
        self.NHe2_c = np.zeros(self.Nl) / self.U.cm**2

        for i in xrange(self.Nl):

            self.NH1_c[i] = np.sum( self.dNH1[0:i] ) + 0.5 * self.dNH1[i]
            self.NHe1_c[i] = np.sum( self.dNHe1[0:i] ) + 0.5 * self.dNHe1[i]
            self.NHe2_c[i] = np.sum( self.dNHe2[0:i] ) + 0.5 * self.dNHe2[i]

        # calculate column densities through whole slab
        #-----------------------------------------------------------
        self.NH1_thru = np.sum( self.dNH * self.xH1 ) 
        self.logNH1_thru = np.log10( self.NH1_thru.magnitude )

        self.NHe1_thru = np.sum( self.dNHe * self.xHe1 ) 
        self.logNHe1_thru = np.log10( self.NHe1_thru.magnitude )

        self.NHe2_thru = np.sum( self.dNHe * self.xHe2 ) 
        self.logNHe2_thru = np.log10( self.NHe2_thru.magnitude )

        # set preferred units
        #-----------------------------------------------------------
        self.z_c.units = 'kpc'

        # delete things we don't need 
        #-----------------------------------------------------------
        del(self.kchem)
        del(self.kcool)

        # calculate effective photoionization rates
        # (i.e. the photoionization rates that result in 
        # the same ionization fractions if caseA rates 
        # are used).  Note, we assume that the recombination 
        # radiation is peaked near the ionization thresholds 
        # and therefore does not alter the temperature.
        #------------------------------------------------        
        
        kchemA = chemistry.ChemistryRates( 
            self.T, fcA_H2 = 1.0, fcA_He2 = 1.0, fcA_He3 = 1.0,
            H1i = self.H1i, He1i = self.He1i, He2i = self.He2i, 
            fit_name = self.atomic_fit_name )
        
        self.H1i_eff = kchemA.reH2 * self.ne * self.xH2 / self.xH1 - \
            kchemA.ciH1 * self.ne

        self.He1i_eff = kchemA.reHe2 * self.ne * self.xHe2 / self.xHe1 - \
            kchemA.ciHe1 * self.ne

        self.He2i_eff = kchemA.reHe3 * self.ne * self.xHe3 / self.xHe2 - \
            kchemA.ciHe2 * self.ne


