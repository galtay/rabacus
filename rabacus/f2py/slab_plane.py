""" Wrapper to make the fortran calling identical to the python calling. """ 

import rabacus_fc
import numpy as np
from slab_base import SlabBase
from rabacus.constants import units
from rabacus.utils import utils
from rabacus.atomic import chemistry
from rabacus.atomic import cooling



__all__ = ['SlabPln', 'Slab2Pln']



class SlabPln(SlabBase):

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

      `rec_meth` (string): How to treat recombinations {"fixed", "ray"}

      `fixed_fcA` (float): If `rec_meth` = "fixed", constant caseA fraction 
      
      `atomic_fit_name` (string): Source for atomic rate fits {"hg97", "enzo13"}

      `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve for 
      equilibrium T 

      `z` (float): Redshift, only need if `find_Teq` = ``True``

      `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
      desired.

      `verbose` (bool): Verbose output? 

      `tol` (float): tolerance for all convergence tests

      `thin` (bool): if ``True`` only solves optically thin


     
    Attributes:

       `U` (:class:`~rabacus.constants.units.Units`)

       `z_c` (array): distance from C to center of layer

       `dz` (array): radial thickness of layer

       `NH_c` (array): H column density from C to center of layer 

       `dNH` (array): H column density through layer

       `NH1_c` (array): HI column density from C to center of layer

       `dNH1` (array): HI column density through layer

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
       through layer 

       `tauH1_th_lo` (array): HI optical depth below this layer

       `tauH1_th_hi` (array): HI optical depth above this layer 

       `NH1_thru` (float): HI column density from C to R
 
       `Nl` (int): Number of layers

       `itr` (int): Number of iterations to converge



      .. note::

        For many of the attributes above, there are analagous versions for 
        helium.  We also note that the source of the atomic rates fit is stored
        in the variable `atomic_fit_name`, but the source of the 
        photoionization cross section fits are stored in the point source 
        object in the variable `px_fit_type`. 
 
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
                 verbose = False,
                 tol = 1.0e-4,
                 thin = False,
                 ):


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

        if find_Teq:
            self.z = z
        else:
            self.z = 0


        # call base class init
        #-----------------------------------------------
        super(SlabPln, self).__init__()

        if Hz == None:
            self.Hz = 0.0 / self.U.s
        else:
            self.Hz = Hz
            self.Hz.units = '1/s'

        # call the fortran solver routine
        #-----------------------------------------------
        self.i_2side = 0

        (xH1, xH2, xHe1, xHe2, xHe3, 
         H1i_src, He1i_src, He2i_src, 
         H1i_rec, He1i_rec, He2i_rec, 
         H1h_src, He1h_src, He2h_src, 
         H1h_rec, He1h_rec, He2h_rec) = \
         rabacus_fc.slab_plane.slab_plane_solve( self.Edges, 
                                                 self.nH, 
                                                 self.nHe, 
                                                 self.T, 
                                                 self.E_eV, 
                                                 self.shape,
                                                 self.i_2side,
                                                 self.i_rec_meth,
                                                 self.fixed_fcA, 
                                                 self.i_photo_fit,
                                                 self.i_rate_fit,
                                                 self.i_find_Teq,
                                                 self.i_thin, 
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



        # do post-calculation
        #-----------------------------------------------
        super(SlabPln,self).__post__()



class Slab2Pln(SlabBase):

    """ 

    Solves a slab geometry with plane parallel radiation incident from both 
    sides.  The slab is divided into `Nl` layers. 

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

      `rec_meth` (string): How to treat recombinations {"fixed", "ray"}

      `fixed_fcA` (float): If `rec_meth` = "fixed", constant caseA fraction 
      
      `atomic_fit_name` (string): Source for atomic rate fits {"hg97"}

      `find_Teq` (bool): If ``False``, use fixed input T, if ``True`` solve 
      for equilibrium T 

      `z` (float): Redshift, only need if `find_Teq` = ``True``

      `Hz` (float): Hubble parameter at `z` if Hubble cooling is 
      desired.

      `verbose` (bool): Verbose output? 

      `tol` (float): tolerance for all convergence tests

      `thin` (bool): if ``True`` only solves optically thin

      `i_rec_ems_prf` (int): if ``0`` recombination radiation is emitted 
      isotropically (i.e. more attenuation), if ``1`` recombination radiation
      is emitted perpendicular to the layers. 
     
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
                 verbose = False,
                 tol = 1.0e-4,
                 thin = False,
                 i_rec_ems_prf = 0
                 ):



        # initialize 
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
        self.i_rec_ems_prf = i_rec_ems_prf

        if find_Teq:
            self.z = z
        else:
            self.z = 0

        # call base class init
        #-----------------------------------------------
        super(Slab2Pln,self).__init__()

        if Hz == None:
            self.Hz = 0.0 / self.U.s
        else:
            self.Hz = Hz
            self.Hz.units = '1/s'

        # check for even number of layers
        #-----------------------------------------------
        if self.nH.size % 2 != 0:
            msg = 'the number of layers must be even'
            raise utils.InputError(msg)



        # call the fortran solver routine
        #-----------------------------------------------
        self.i_2side = 1

        (xH1, xH2, xHe1, xHe2, xHe3, 
         H1i_src, He1i_src, He2i_src, 
         H1i_rec, He1i_rec, He2i_rec, 
         H1h_src, He1h_src, He2h_src,
         H1h_rec, He1h_rec, He2h_rec) = \
         rabacus_fc.slab_plane.slab_plane_solve( self.Edges, 
                                                 self.nH, 
                                                 self.nHe, 
                                                 self.T, 
                                                 self.E_eV, 
                                                 self.shape,
                                                 self.i_2side,
                                                 self.i_rec_meth,
                                                 self.fixed_fcA, 
                                                 self.i_photo_fit,
                                                 self.i_rate_fit,
                                                 self.i_find_Teq,
                                                 self.i_thin, 
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


        # do post-calculation
        #-----------------------------------------------
        super(Slab2Pln,self).__post__()



