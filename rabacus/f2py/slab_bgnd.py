""" Wrapper to make python calling easy. """ 

import rabacus_fc
import numpy as np
from slab_base import SlabBase
from rabacus.constants import units
from rabacus.utils import utils
from rabacus.atomic import chemistry
from rabacus.atomic import cooling
from rabacus.hdf5 import h5py_wrap as h5
from rabacus.utils import units_string

__all__ = ['Slab2Bgnd']



class Slab2Bgnd(SlabBase):

    """ 

    Solves a slab geometry in a uniform background.  Radiation is implicitly 
    incident from both sides of the slab.  The slab is divided into `Nl` 
    layers. 

    Args: 
      `Edges` (array): Distance to layer edges 

      `T` (array): Temperature in each layer 

      `nH` (array): hydrogen number density in each layer 

      `nHe` (array): helium number density in each layer 

      `rad_src` (:class:`~rabacus.rad_src.background.BackgroundSource`): 
      Background source


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

        if z != None:
            self.z = z
        else:
            self.z = 0.0


        # call base class init
        #-----------------------------------------------
        super(Slab2Bgnd,self).__init__()

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


        # call solver 
        #-----------------------------------------------
        self.call_fortran_solver()


        # do post-calculation
        #-----------------------------------------------
        super(Slab2Bgnd,self).__post__()



    def call_fortran_solver( self ):

        (xH1, xH2, xHe1, xHe2, xHe3, 
         H1i_src, He1i_src, He2i_src, 
         H1i_rec, He1i_rec, He2i_rec, 
         H1h_src, He1h_src, He2h_src,
         H1h_rec, He1h_rec, He2h_rec) = \
         rabacus_fc.slab_bgnd.slab_bgnd_solve( 
             self.Edges, 
             self.nH, 
             self.nHe, 
             self.T, 
             self.E_eV, 
             self.shape,
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







    def write_to_file_generic( 
            self, fname, top_group, grp_name, obj,
            write_temperature=True, write_number_density=True,
            write_ionization_fractions=True, write_photoion_rates=True,
            write_photoheat_rates=True, write_mu=True ):

        """ 
        writes a lot of attributes to file fname in hdf5 group, 

          fname/top_group/grp_name

        the attributes are taken from obj. 
        """ 

        # write layers
        #=====================================================
        tg = top_group
        gg = tg+'/'+grp_name
        h5.cg( fname, gg )


        if write_temperature:

            # temperatures
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.T), 
                'description': 'temperature'}
            h5.wd( fname, gg, 'T', obj.T, attrs )


        if write_mu:

            # mean molecular weight
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.mu), 
                'description': 'mean molecular weight'}
            h5.wd( fname, gg, 'mu', obj.mu, attrs )


        if write_number_density:

            # number density
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.nH), 
                'description': 'hydrogen number density'}
            h5.wd( fname, gg, 'nH', obj.nH, attrs )
            
            attrs = {
                'units': units_string(obj.nHe), 
                'description': 'helium number density'}
            h5.wd( fname, gg, 'nHe', obj.nHe, attrs )
            
            attrs = {
                'units': units_string(obj.ne), 
                'description': 'electron number density'}
            h5.wd( fname, gg, 'ne', obj.ne, attrs )


        if write_ionization_fractions:

            # ionization fractions
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.xH1), 
                'description': 'H I ionization fraction = nHI / nH'}
            h5.wd( fname, gg, 'xH1', obj.xH1, attrs )

            attrs = {
                'units': units_string(obj.xH2), 
                'description': 'H II ionization fraction = nHII / nH'}
            h5.wd( fname, gg, 'xH2', obj.xH2, attrs )

            attrs = {
                'units': units_string(obj.xHe1), 
                'description': 'He I ionization fraction = nHeI / nHe'}
            h5.wd( fname, gg, 'xHe1', obj.xHe1, attrs )
            
            attrs = {
                'units': units_string(obj.xHe2), 
                'description': 'He II ionization fraction = nHeII / nHe'}
            h5.wd( fname, gg, 'xHe2', obj.xHe2, attrs )
            
            attrs = {
                'units': units_string(obj.xHe3), 
                'description': 'He III ionization fraction = nHeIII / nHe'}
            h5.wd( fname, gg, 'xHe3', obj.xHe3, attrs )
    

        if write_photoion_rates:

            # total photo-ionization rates
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.H1i), 
                'description': 'total H I photo-ionization rate'}
            h5.wd( fname, gg, 'H1i', obj.H1i, attrs )

            attrs = {
                'units': units_string(obj.He1i), 
                'description': 'total He I photo-ionization rate'}
            h5.wd( fname, gg, 'He1i', obj.He1i, attrs )
            
            attrs = {
                'units': units_string(obj.He2i), 
                'description': 'total He II photo-ionization rate'}
            h5.wd( fname, gg, 'He2i', obj.He2i, attrs )


            # background source photo-ionization rates
            #-----------------------------------------------------
            txt = 'due to background source' 

            attrs = {
                'units': units_string(obj.H1i_src), 
                'description': 'H I photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'H1i_src', obj.H1i_src, attrs )

            attrs = {
                'units': units_string(obj.He1i_src), 
                'description': 'He I photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'He1i_src', obj.He1i_src, attrs )

            attrs = {
                'units': units_string(obj.He2i_src), 
                'description': 'He II photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'He2i_src', obj.He2i_src, attrs )


            # recombination photo-ionization rates
            #-----------------------------------------------------
            txt = 'due to recombination photons' 

            attrs = {
                'units': units_string(obj.H1i_rec), 
                'description': 'H I photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'H1i_rec', obj.H1i_rec, attrs )

            attrs = {
                'units': units_string(obj.He1i_rec), 
                'description': 'He I photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'He1i_rec', obj.He1i_rec, attrs )

            attrs = {
                'units': units_string(obj.He2i_rec), 
                'description': 'He II photo-ionization rate ' + txt}
            h5.wd( fname, gg, 'He2i_rec', obj.He2i_rec, attrs )


        if write_photoheat_rates:

            # total photo-heating rates
            #-----------------------------------------------------
            attrs = {
                'units': units_string(obj.H1h), 
                'description': 'total H I photo-heating rate'}
            h5.wd( fname, gg, 'H1h', obj.H1h, attrs )
            
            attrs = {
                'units': units_string(obj.He1h), 
                'description': 'total He I photo-heating rate'}
            h5.wd( fname, gg, 'He1h', obj.He1h, attrs )
            
            attrs = {
                'units': units_string(obj.He2h), 
                'description': 'total He II photo-heating rate'}
            h5.wd( fname, gg, 'He2h', obj.He2h, attrs )


            # background source photo-heating rates
            #-----------------------------------------------------
            txt = 'due to background source' 

            attrs = {
                'units': units_string(obj.H1h_src), 
                'description': 'H I photo-heating rate ' + txt}
            h5.wd( fname, gg, 'H1h_src', obj.H1h_src, attrs )
            
            attrs = {
                'units': units_string(obj.He1h_src), 
                'description': 'He I photo-heating rate ' + txt}
            h5.wd( fname, gg, 'He1h_src', obj.He1h_src, attrs )
            
            attrs = {
                'units': units_string(obj.He2h_src), 
                'description': 'He II photo-heating rate ' + txt}
            h5.wd( fname, gg, 'He2h_src', obj.He2h_src, attrs )
            

            # recombination photo-heating rates
            #-----------------------------------------------------
            txt = 'due to recombination photons' 
            
            attrs = {
                'units': units_string(obj.H1h_rec), 
                'description': 'H I photo-heating rate ' + txt}
            h5.wd( fname, gg, 'H1h_rec', obj.H1h_rec, attrs )
            
            attrs = {
                'units': units_string(obj.He1h_rec), 
                'description': 'He I photo-heating rate ' + txt}
            h5.wd( fname, gg, 'He1h_rec', obj.He1h_rec, attrs )
            
            attrs = {
                'units': units_string(obj.He2h_rec), 
                'description': 'He II photo-heating rate ' + txt}
            h5.wd( fname, gg, 'He2h_rec', obj.He2h_rec, attrs )
            




    def write_to_file( self, fname, top_group='/', write_layers=True,
                       write_centrals=True, write_averages=True, 
                       write_source=True ):

        """ Writes a background slab out to an hdf5 file. """ 

        # write out to file
        #=====================================================
        tg = top_group
        h5.cg( fname, tg )

        # write out source information
        #-----------------------------------------------------
        if write_source:
            self.rad_src.write_to_file( fname, top_group=tg+'/Src' )

        # write config
        #-----------------------------------------------------
        h5.cg( fname, tg+'/Config' )

        h5.wa( fname, tg+'/Config', 'tol', self.tol )
        h5.wa( fname, tg+'/Config', 'z', self.z )
        h5.wa( fname, tg+'/Config', 'Nl', self.Nl )
        h5.wa( fname, tg+'/Config', 'rec_meth', self.rec_meth )
        h5.wa( fname, tg+'/Config', 'atomic_fit_name', self.atomic_fit_name )
        h5.wa( fname, tg+'/Config', 'fixed_fcA', self.fixed_fcA )

        # write total column densities
        #-----------------------------------------------------
        attrs = {
            'units': units_string(self.NH1_thru), 
            'description': 'H I column density = sum(dz * nH * xH1)'}
        h5.wd( fname, tg, 'NH1', self.NH1_thru, attrs )

        attrs = {
            'units': units_string(self.NHe1_thru), 
            'description': 'He I column density = sum(dz * nHe * xHe1)'}
        h5.wd( fname, tg, 'NHe1', self.NHe1_thru, attrs )

        attrs = {
            'units': units_string(self.NHe2_thru), 
            'description': 'He II column density = sum(dz * nHe * xHe2)'}
        h5.wd( fname, tg, 'NHe2', self.NHe2_thru, attrs )


        # write slab thickness
        #-----------------------------------------------------
        Lslab = self.Edges[-1] - self.Edges[0]
        Lslab.units = 'kpc'

        attrs = {
            'units': units_string(Lslab), 
            'description': 'Thickness of slab.'}
        h5.wd( fname, tg, 'Lslab', Lslab, attrs )

        # write layers
        #-----------------------------------------------------
        if write_layers:

            grp_name = 'Layers'
            obj = self
            self.write_to_file_generic( fname, top_group, grp_name, obj )
        
            attrs = {
                'units': units_string(self.Edges), 
                'description': 'position of slab edges'}
            h5.wd( fname, tg+'/Layers', 'Edges', self.Edges, attrs )

            attrs = {
                'units': units_string(self.z_c), 
                'description': 'central position of layer'}
            h5.wd( fname, tg+'/Layers', 'z_c', self.z_c, attrs )
            
            attrs = {
                'units': units_string(self.dz), 
                'description': 'thickness of each layer'}
            h5.wd( fname, tg+'/Layers', 'dz', self.dz, attrs )

            attrs = {
                'units': units_string(self.H1i_eff), 
                'description': 
                ('H I effective photo-ionization rate.'  
                 '(for fixed_fcA, the photo-ion rate that results in the'
                 'same xHI if case A rates are used)')}
            h5.wd( fname, tg+'/Layers', 'H1i_eff', self.H1i_eff, attrs )

            attrs = {
                'units': units_string(self.He1i_eff), 
                'description': 
                ('He I effective photo-ionization rate.'  
                 '(for fixed_fcA, the photo-ion rate that results in the'
                 'same xHeI if case A rates are used)')}
            h5.wd( fname, tg+'/Layers', 'He1i_eff', self.He1i_eff, attrs )

            attrs = {
                'units': units_string(self.He2i_eff), 
                'description': 
                ('He II effective photo-ionization rate.'  
                 '(for fixed_fcA, the photo-ion rate that results in the'
                 'same xHeII if case A rates are used)')}
            h5.wd( fname, tg+'/Layers', 'He2i_eff', self.He2i_eff, attrs )



        # write central values
        #-----------------------------------------------------
        if write_centrals:

            grp_name = 'Central'
            obj = self.avg.cen
            self.write_to_file_generic( fname, top_group, grp_name, obj )

        # write averages
        #-----------------------------------------------------
        if write_averages:

            h5.cg( fname, tg+'/Average' )            

            grp_name = 'Average/v_w'
            obj = self.avg.v_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )

            grp_name = 'Average/nH_w'
            obj = self.avg.nH_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )

            grp_name = 'Average/nH1_w'
            obj = self.avg.nH1_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )

            grp_name = 'Average/nHe_w'
            obj = self.avg.nHe_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )

            grp_name = 'Average/nHe1_w'
            obj = self.avg.nHe1_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )

            grp_name = 'Average/nHe2_w'
            obj = self.avg.nHe2_w
            self.write_to_file_generic( fname, top_group, grp_name, obj,
                                        write_number_density=False )
