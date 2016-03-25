""" Solves a sphere with a point source at the center.  Canonical examples
are Tests 1 and 2 from http://arxiv.org/abs/astro-ph/0603199. """

import rabacus_fc
import numpy as np
from rabacus.atomic import chemistry
from rabacus.constants import units
from sphere_base import SphereBase

__all__ = ['SphereStromgren']


class SphereStromgren(SphereBase):

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
      determines the number of layers, `Nl`.  `Edges` determines the positions
      of the shells that bound the layers and must have `Nl+1` entries. 


    Kwargs:

      `rec_meth` (string): How to treat recombinations 
      {"fixed", "outward", "radial", "isotropic"}

      `fixed_fcA` (float): If `rec_meth` = "fixed", constant caseA fraction 
      
      `atomic_fit_name` (string): Source for atomic rate fits {"hg97"}

      `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve for 
      equilibrium T 

      `z` (float): Redshift, only need if `find_Teq` = ``True``

      `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
      desired.

      `em_H1_fac` (float):  multiplicative factor for H1 recomb emission

      `em_He1_fac` (float):  multiplicative factor for He1 recomb emission

      `em_He2_fac` (float):  multiplicative factor for He2 recomb emission

      `verbose` (bool): Verbose output? 

      `tol` (float): tolerance for all convergence tests

      `thin` (bool): if ``True`` only solves optically thin

      `Nmu` (int): number of polar angle bins

     
    Attributes:

       `U` (:class:`rabacus.constants.units.Units`)

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
                 Hz = None,
                 em_H1_fac=1.0,
                 em_He1_fac=1.0,
                 em_He2_fac=1.0,
                 verbose = False,           
                 tol = 1.0e-4,
                 thin = False,
                 Nmu = 32,
                 ):


        # check input
        #-----------------------------------------------
        if rad_src.source_type != 'point':
            msg = 'source type needs to be point'
            raise ValueError(msg)


        # attach input
        #-----------------------------------------------
        self.Edges = Edges.copy()
        self.T = T.copy()
        self.nH = nH.copy()
        self.nHe = nHe.copy()
        self.rad_src = rad_src

        self.rec_meth = rec_meth 
        self.fixed_fcA = fixed_fcA
        self.atomic_fit_name = atomic_fit_name
        self.find_Teq = find_Teq
        self.verbose = verbose
        self.tol = tol
        self.thin = thin
        self.Nmu = Nmu

        self.em_H1_fac = em_H1_fac
        self.em_He1_fac = em_He1_fac
        self.em_He2_fac = em_He2_fac

        if find_Teq:
            self.z = z
        else:
            self.z = 0.0


        # call base class init
        #-----------------------------------------------
        super(SphereStromgren,self).__init__()

        if Hz == None:
            self.Hz = 0.0 / self.U.s
        else:
            self.Hz = Hz
            self.Hz.units = '1/s'

        
        # call solver
        #-----------------------------------------------
        self.call_fortran_solver()


        # do post-calculation
        #-----------------------------------------------
        super(SphereStromgren,self).__post__()





            

    def call_fortran_solver( self ):
                
        (xH1, xH2, xHe1, xHe2, xHe3, 
         H1i_src, He1i_src, He2i_src, 
         H1i_rec, He1i_rec, He2i_rec, 
         H1h_src, He1h_src, He2h_src, 
         H1h_rec, He1h_rec, He2h_rec) = \
         rabacus_fc.sphere_stromgren.sphere_stromgren_solve( 
             self.Edges, 
             self.nH, 
             self.nHe, 
             self.T, 
             self.Nmu,
             self.E_eV, 
             self.shape, 
             self.i_rec_meth, 
             self.fixed_fcA, 
             self.i_photo_fit, 
             self.i_rate_fit, 
             self.i_find_Teq,
             self.i_thin,
             self.em_H1_fac,
             self.em_He1_fac,
             self.em_He2_fac,
             self.z, 
             self.Hz,
             self.tol,
             self.Nl, 
             self.Nnu )
         
        # attach output with units. 
        #-----------------------------------------------
        self.xH1 = xH1
        self.xH2 = xH2
        self.xHe1 = xHe1
        self.xHe2 = xHe2
        self.xHe3 = xHe3


        self.H1i_src = H1i_src / self.U.s
        self.He1i_src = He1i_src / self.U.s
        self.He2i_src = He2i_src / self.U.s

        self.H1i_rec = H1i_rec / self.U.s
        self.He1i_rec = He1i_rec / self.U.s
        self.He2i_rec = He2i_rec / self.U.s

        self.H1i = self.H1i_src + self.H1i_rec
        self.He1i = self.He1i_src + self.He1i_rec
        self.He2i = self.He2i_src + self.He2i_rec


        self.H1h_src = H1h_src * self.U.erg / self.U.s
        self.He1h_src = He1h_src * self.U.erg / self.U.s
        self.He2h_src = He2h_src * self.U.erg / self.U.s

        self.H1h_rec = H1h_rec * self.U.erg / self.U.s
        self.He1h_rec = He1h_rec * self.U.erg / self.U.s
        self.He2h_rec = He2h_rec * self.U.erg / self.U.s

        self.H1h = self.H1h_src + self.H1h_rec
        self.He1h = self.He1h_src + self.He1h_rec
        self.He2h = self.He2h_src + self.He2h_rec


