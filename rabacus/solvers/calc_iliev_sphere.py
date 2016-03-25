""" Solves a sphere with a point source at the center.  Canonical examples
are Tests 1 and 2 from http://arxiv.org/abs/astro-ph/0603199. """

import numpy as np
import one_zone as ozn
from rabacus.atomic import chemistry
from rabacus.atomic import cooling
from rabacus.constants import units
from rabacus.utils import utils



__all__ = ['IlievSphere']


class IlievSphere:

    """ 

    Stores and calculates equilibrium ionization (and optionally temperature) 
    structure in a spherically symmetric geometry with a point source at the 
    center.  The sphere is divided into `Nl` shells. We will refer to the 
    center of the sphere as "C" and the radius of the sphere as "R". 


    Args:
      `Edges` (array): Radius of all shell edges
        
      `T` (array): Temperature in each shell
        
      `nH` (array): Hydrogen number density in each shell
        
      `nHe` (array): Helium number density in each shell 
        
      `rad_src` (:class:`~rabacus.rad_src.point.PointSource`): Point source

    .. note::
      The arrays `T`, `nH`, and `nHe` must all be the same size.  This size
      determines the number of shells, `Nl`.  `Edges` determines the positions
      of the shell edges and must have `Nl+1` entries. 


    Kwargs:

      `rec_meth` (string): How to treat recombinations 
      {"fixed"}

      `fixed_fcA` (float): If `rec_meth` = "fixed", constant caseA fraction 
      
      `atomic_fit_name` (string): Source for atomic rate fits {"hg97"}

      `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve for 
      equilibrium T 

      `z` (float): Redshift, only need if `find_Teq` = ``True``

      `verbose` (bool): Verbose output? 

      `tol` (float): tolerance for all convergence tests

      `thin` (bool): if ``True`` only solves optically thin
     
    Attributes:

       `U` (:class:`~rabacus.constants.units.Units`)

       `r_c` (array): distance from C to center of shell

       `dr` (array): radial thickness of shell

       `NH_c` (array): H column density from C to center of shell 

       `dNH` (array): H column density through shell

       `NH1_c` (array): HI column density from C to center of shell

       `dNH1` (array): HI column density through shell

       `H1i` (array): HI photo-ionization rate  

       `H1h` (array): HI photo-heating rate  

       `xH1` (array): H neutral fraction nHI / nH
       
       `xH2` (array): H ionized fraction nHII / nH

       `ne` (array): electron number density 

       `fcA_H2` (array): HII case A fraction 
 
       `cool` (array): cooling rate [erg / (cm^3 K)]

       `heat` (array): heating rate [erg / (cm^3 K)]

       `heatH1` (array): contribution to `heat` from H1 photo-heating

       `dtauH1_th` (array): HI optical depth at H1 ionizing threshold
       through shell 

       `tauH1_th_lo` (array): HI optical depth below this shell

       `tauH1_th_hi` (array): HI optical depth above this shell 

       `NH1_thru` (float): HI column density from C to R
 
       `Nl` (int): Number of shells

       `itr` (int): Number of iterations to converge

    .. note::

      For many of the attributes above, there are analagous versions for 
      helium.  We also note that the source of the atomic rates fit is stored
      in the variable `atomic_fit_name`, but the source of the photoionization
      cross section fits are stored in the point source object in the variable
      `px_fit_type`. 


    Examples:    

      The following code will solve a monochromatic point source inside 
      a spherical distribution with fixed recombination rates and temperatues. 
      In this case we use a Haardt and Madau 2012 spectrum to calculate a grey 
      photon energy which is then used to create a monochromatic plane 
      source :: 

        import numpy as np
        import rabacus as ra
        
        # create Haardt and Madau 2012 source
        #---------------------------------------------------------------
        q_min = 1.0; q_max = 4.0e2; z = 3.0
        hm12 = ra.rad_src.PointSource( q_min, q_max, 'hm12', z=z )
 
        # get grey photon energy and create a monochromatic plane source
        #---------------------------------------------------------------
        q_mono = hm12.grey.E.He2 / hm12.PX.Eth_H1
        mono = ra.rad_src.PointSource( q_mono, q_mono, 'monochromatic' )
        mono.normalize_Ln( 5.0e48 / ra.U.s )

        # describe slab
        #---------------------------------------------------------------
        Nl = 512; Yp = 0.24
        nH = np.ones(Nl) * 1.0e-3 / ra.U.cm**3
        nHe = nH * Yp * 0.25 / (1.0-Yp)
        Tslab = np.ones(Nl) * 1.0e4 * ra.U.K
        Rsphere = 6.6 * ra.U.kpc
        Edges = np.linspace( 0.0 * ra.U.kpc, Rsphere, Nl+1 )

        # solve slab
        #---------------------------------------------------------------
        sphere = ra.solvers.IlievSphere( Edges, Tslab, nH, nHe, mono,
                                         rec_meth='fixed', fixed_fcA=0.0 )

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

        assert( rec_meth == 'fixed' )

        if find_Teq and z==None:
            msg = 'need to supply redshift if finding equilibrium temperature'
            raise utils.InputError(msg)

        if rad_src.source_type != 'point':
            msg = 'source type needs to be point'
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

        self.atomic_fit_name = atomic_fit_name
        self.find_Teq = find_Teq
        self.verbose = verbose
        self.tol = tol
        self.thin = thin

        if find_Teq:
            self.z = z


        # initialize sphere
        #-----------------------------------------------
        self.init_sphere()
        self.set_optically_thin()

        if self.thin:
            return

        # solve sphere (sweep until convergence)
        #-----------------------------------------------------------
        conv_old = np.sum( self.ne )
        not_converged = True
        self.itr = 0

        while not_converged:
            self.sweep_sphere()
            conv_new = np.sum( self.ne )
            if np.abs( conv_new/conv_old - 1.0 ) < self.tol:
                not_converged = False
            conv_old = conv_new
            self.itr += 1


        # finalize slab
        #-----------------------------------------------------------
        self.finalize_sphere()







    def init_sphere( self ):
        """ Initialize sphere values.  """

        if self.verbose:
            print 'begin initialization' 

        # set arbitrary radius
        #-----------------------------------------------
        r = 1.0 * self.U.kpc

        # instantiate atomic rates (optically thin)
        #-----------------------------------------------
        fcA_H2 = self.fixed_fcA
        fcA_He2 = self.fixed_fcA
        fcA_He3 = self.fixed_fcA

        kchem = chemistry.ChemistryRates( 
            self.T[0], fcA_H2, fcA_He2, fcA_He3, 
            H1i = self.rad_src.thin.H1i(r), 
            He1i = self.rad_src.thin.He1i(r), 
            He2i = self.rad_src.thin.He2i(r),
            fit_name = self.atomic_fit_name 
            )

        kcool = cooling.CoolingRates(
            self.T[0], fcA_H2, fcA_He2, fcA_He3, 
            H1h = self.rad_src.thin.H1h(r), 
            He1h = self.rad_src.thin.He1h(r), 
            He2h = self.rad_src.thin.He2h(r), 
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

        self.dr = self.Edges[1:] - self.Edges[0:-1]
        self.r_c = self.Edges[0:-1] + 0.5 * self.dr 
                                                     
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

        self.cool    = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heat    = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatH1  = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatHe1 = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)
        self.heatHe2 = np.zeros(Nl) * self.U.eV / (self.U.s * self.U.cm**3)

        self.xH2 = np.zeros(Nl)
        self.xHe3 = np.zeros(Nl)
        self.ne = np.zeros(Nl) / self.U.cm**3

        if self.verbose:
            print '  created simple arrays' 

        # calc NH and NHe between C and center of each layer
        #----------------------------------------
        self.dNH  = self.dr * self.nH 
        self.dNHe = self.dr * self.nHe

        self.NH_c = np.zeros(Nl) / self.U.cm**2
        self.NHe_c = np.zeros(Nl) / self.U.cm**2

        for i in xrange(Nl):

            self.NH_c[i] = np.sum( self.dNH[0:i] ) + 0.5 * self.dNH[i]
            self.NHe_c[i] = np.sum( self.dNHe[0:i] ) + 0.5 * self.dNHe[i]


        if self.verbose:
            print 'initialization complete' 




    def set_optically_thin( self ): 

        """ Set optically thin ionzation state at input temperature.  
        Adjustments will be made during the sweeps. """ 

        # initialize a chemistry object
        #----------------------------------------
        fcA_H2 = self.fixed_fcA
        fcA_He2 = self.fixed_fcA
        fcA_He3 = self.fixed_fcA


        kchem = chemistry.ChemistryRates( 
            self.T[0], fcA_H2, fcA_He2, fcA_He3,
            H1i = self.rad_src.thin.H1i( self.r_c[0] ), 
            He1i = self.rad_src.thin.He1i( self.r_c[0] ),
            He2i = self.rad_src.thin.He2i( self.r_c[0] ),
            fit_name = self.atomic_fit_name, 
            )

        kcool = cooling.CoolingRates( 
            self.T[0], fcA_H2, fcA_He2, fcA_He3,
            H1h = self.rad_src.thin.H1h( self.r_c[0] ), 
            He1h = self.rad_src.thin.He1h( self.r_c[0] ),
            He2h = self.rad_src.thin.He2h( self.r_c[0] ),
            fit_name = self.atomic_fit_name, 
            )


        # loop through layers and set x
        #----------------------------------------
        for i in range(self.Nl):

            H1i = self.rad_src.thin.H1i( self.r_c[i] )
            He1i = self.rad_src.thin.He1i( self.r_c[i] )
            He2i = self.rad_src.thin.He2i( self.r_c[i] )

            kchem.set( self.T[i], fcA_H2, fcA_He2, fcA_He3,
                       H1i=H1i, He1i=He1i, He2i=He2i )

            H1h = self.rad_src.thin.H1h( self.r_c[i] )
            He1h = self.rad_src.thin.He1h( self.r_c[i] )
            He2h = self.rad_src.thin.He2h( self.r_c[i] )

            kcool.set( self.T[i], fcA_H2, fcA_He2, fcA_He3,
                       H1h=H1h, He1h=He1h, He2h=He2h )

            if self.find_Teq:
                x = ozn.Solve_PCTE( self.nH[i], self.nHe[i], 
                                    kchem, kcool, self.z, self.tol )
                self.T[i] = x.Teq
            else:
                x = ozn.Solve_PCE( self.nH[i], self.nHe[i], 
                                   kchem, self.tol )

        
            self.xH1[i] = x.H1
            self.xH2[i] = x.H2
            self.xHe1[i] = x.He1
            self.xHe2[i] = x.He2
            self.xHe3[i] = x.He3

        
        # summarize
        #----------------------------------------
        self.dNH1 = self.dNH * self.xH1
        self.dNHe1 = self.dNHe * self.xHe1
        self.dNHe2 = self.dNHe * self.xHe2
        
        self.dtauH1_th = self.dNH1 * self.rad_src.th.sigma_H1
        self.dtauHe1_th = self.dNHe1 * self.rad_src.th.sigma_He1
        self.dtauHe2_th = self.dNHe2 * self.rad_src.th.sigma_He2
        
        self.ne = self.xH2 * self.nH + \
            ( self.xHe2 + 2 * self.xHe3 ) * self.nHe 



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



    def set_photoion_rates( self, i ):
        """ set photoionization rates for this layer 

        Args: 
          `i` (int): layer index
          """ 

        tauH1_th = np.float( self.tauH1_th_lo[i] + 0.5 * self.dtauH1_th[i] )
        tauHe1_th = np.float( self.tauHe1_th_lo[i] + 0.5 * self.dtauHe1_th[i] )
        tauHe2_th = np.float( self.tauHe2_th_lo[i] + 0.5 * self.dtauHe2_th[i] )
        
        self.H1i[i] = self.rad_src.shld_H1i( self.r_c[i],
                                             tauH1_th, tauHe1_th, tauHe2_th )

        self.He1i[i] = self.rad_src.shld_He1i( self.r_c[i],
                                               tauH1_th, tauHe1_th, tauHe2_th ) 

        self.He2i[i] = self.rad_src.shld_He2i( self.r_c[i],
                                               tauH1_th, tauHe1_th, tauHe2_th )


    def set_photoheat_rates( self, i ):
        """ set photoheating rates for this layer 

        Args: 
          `i` (int): layer index
          """ 

        tauH1_th = np.float( self.tauH1_th_lo[i] + 0.5 * self.dtauH1_th[i] )
        tauHe1_th = np.float( self.tauHe1_th_lo[i] + 0.5 * self.dtauHe1_th[i] )
        tauHe2_th = np.float( self.tauHe2_th_lo[i] + 0.5 * self.dtauHe2_th[i] )
        
        self.H1h[i] = self.rad_src.shld_H1h( self.r_c[i],
                                             tauH1_th, tauHe1_th, tauHe2_th )

        self.He1h[i] = self.rad_src.shld_He1h( self.r_c[i],
                                               tauH1_th, tauHe1_th, tauHe2_th ) 

        self.He2h[i] = self.rad_src.shld_He2h( self.r_c[i],
                                               tauH1_th, tauHe1_th, tauHe2_th )



    def sweep_sphere(self):
        """ Performs one sweep through sphere. """ 

        for i in xrange(self.Nl):

            if self.verbose:
                print 'i = ', i

            # calc tauXX above and below this layer
            # (constant during iterations)
            #----------------------------------------
            self.set_taus( i )

            # set fcA 
            #----------------------------------------
            self.fcA_H2 = np.ones(self.Nl) * self.fixed_fcA
            self.fcA_He2 = np.ones(self.Nl) * self.fixed_fcA
            self.fcA_He3 = np.ones(self.Nl) * self.fixed_fcA
            
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
                # we need to call IonTemp_Eq
                #----------------------------------------
                if self.find_Teq:

                    x = ozn.Solve_PCTE( np.ones(1) * self.nH[i], 
                                        np.ones(1) * self.nHe[i], 
                                        self.kchem, 
                                        self.kcool, 
                                        self.z,
                                        self.tol )

                    self.T[i] = x.Teq

                # otherwise we call Ion_Eq
                #----------------------------------------
                else:

                    x = ozn.Solve_PCE( np.ones(1) * self.nH[i], 
                                       np.ones(1) * self.nHe[i], 
                                       self.kchem, self.tol )


                # set the ionization fraction in the object
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


    def finalize_sphere( self ):
        """ Calculate some final values using fully solved sphere. """ 


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
        self.r_c.units = 'kpc'

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




