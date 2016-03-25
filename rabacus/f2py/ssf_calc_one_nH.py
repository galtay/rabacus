""" Calculates a slab solution for one density and one redshift """ 

import fc
import numpy as np
import idlehands as ih

class UserDefinedType:
    pass 

# Fits to use for photoion-xsection and atomic rates
#----------------------------------------------------------------
i_photo_fit = 1  # verner 96
i_rate_fit = 1   # hui gnedin 97

fcA = 0.0
i_caseA = 0


# Create PhotoXsections object.
#----------------------------------------------------------------
PX = ih.atomic.PhotoXsections(fit='verner')


# Cosmo Objects
#----------------------------------------------------------------
Planck = ih.cosmology.parameters.planck.load.PlanckParameters()
Cosmo = ih.cosmology.Cosmology( Planck.values )
Yp = Planck.values['yhe']
fg = Planck.values['omegab'] / Planck.values['omegam']
Jeans = ih.cosmology.Jeans( Yp, fg )


# Number of frequency bins, number of layers, spectral range
#--------------------------------------------------------
Nnu = 1000
Nl = 5000

q_min = 1.0
q_max = 1.0e2


# baseline spectra 
#--------------------------------------------------------
z = 3.02
hm12_bgnd = ih.source_class.BackgroundSource( 
    q_min, q_max, 'hm12', Npts=Nnu, z=z )
hm12_pln = ih.source_class.PlaneSource( 
    q_min, q_max, 'hm12', z=z )
hm12_tab = ih.uv_bgnd.HM12_Photorates_Table()


# make polychromatic spectra
#-----------------------------------------
E_min = q_min * PX.Eth_H1
E_max = q_max * PX.Eth_H1

hm12_fc = ih.f2py.wrap.SourcePlane_HM12( E_min, E_max, Nnu, z )
hm12_fc.normalize_n( hm12_bgnd.n_thin )

PR = UserDefinedType()

PR.H1i = hm12_fc.H1i_thin
PR.He1i = hm12_fc.He1i_thin
PR.He2i = hm12_fc.He2i_thin

PR.H1h = hm12_fc.H1h_thin
PR.He1h = hm12_fc.He1h_thin
PR.He2h = hm12_fc.He2h_thin


# set up slab variables
#-----------------------------------------
Delta = 4.0e3
nH_z = Cosmo.nH_critz( z )
nH = Delta * nH_z
nHe = nH * Yp / ( 4.0 * (1.0-Yp) )

print 'slab: '
print 'Delta: ', Delta
print 'nH: ', nH
print 'nHe: ', nHe
print 


pcte = ih.f2py.wrap.SolvePCTE( nH, nHe, z, fcA, i_rate_fit, PR )

print 'one zone:'
print 'xH1: ', pcte.xH1
print 'xH2: ', pcte.xH2
print 'xHe1: ', pcte.xHe1
print 'xHe2: ', pcte.xHe2
print 'xHe3: ', pcte.xHe3
print 'T: ', pcte.T
print 




# solve slab
#==============================================================
slab = UserDefinedType()

slab.i_2side = 1
slab.i_caseA = 0
slab.i_photo_fit = i_photo_fit
slab.i_rate_fit = i_rate_fit
slab.i_find_Teq = 1




# find initial Jeans scales
#----------------------------------------
mu = Jeans.mu( pcte.xH2, pcte.xHe2, pcte.xHe3 )
LJ = Jeans.L( nH, pcte.T, mu ) 
LJ_old = LJ

print 'Jeans: '
print 'mu: ', mu
print 'LJ: ', LJ
print 

TOL = 1.0e-4
err = 1.0e20
itr = 0



while ( err > TOL ):

    # create slab
    #----------------------------------------
    LJ.units = 'cm'
    slab.Ledges = np.linspace( 0.0 * ih.U.kpc, LJ, Nl+1 )
    slab.T = np.ones( Nl ) * pcte.T
    slab.nH = np.ones( Nl ) * nH
    slab.nHe = np.ones( Nl ) * nHe
    slab.Lcen = np.linspace( 0.0 * ih.U.kpc, LJ, Nl )

    # solve slab
    #----------------------------------------
    ih.f2py.wrap.SlabPlane( slab, hm12_fc, z )

    Tcen = ( slab.T[Nl/2-1] + slab.T[Nl/2] ) * 0.5
    xH2cen = ( slab.xH2[Nl/2-1] + slab.xH2[Nl/2] ) * 0.5
    xHe2cen = ( slab.xHe2[Nl/2-1] + slab.xHe2[Nl/2] ) * 0.5
    xHe3cen = ( slab.xHe3[Nl/2-1] + slab.xHe3[Nl/2] ) * 0.5

    print 'slab (min/max): '
    print 'xH1: ', slab.xH1.min(), slab.xH1.max()
    print 'xH2: ', slab.xH2.min(), slab.xH2.max()
    print 'xHe1: ', slab.xHe1.min(), slab.xHe1.max()
    print 'xHe2: ', slab.xHe2.min(), slab.xHe2.max()
    print 'xHe3: ', slab.xHe3.min(), slab.xHe3.max()
    print 'T: ', slab.T.min(), slab.T.max()
    print

    yy = slab.nH * slab.xH1 * slab.T
    num = fc.utils.utils_integrate_trap( slab.Lcen, yy )
    
    yy = slab.nH * slab.xH1
    den = fc.utils.utils_integrate_trap( slab.Lcen, yy )
    
    TH1 = num/den
    
    print 'Centered (min/max): '
    print 'TH1: ', TH1
    print 'Tcen: ', Tcen
    print 'xH2cen: ', xH2cen
    print 'xHe2cen: ', xHe2cen
    print 'xHe3cen: ', xHe3cen
    print
    
    mu = Jeans.mu( xH2cen, xHe2cen, xHe3cen )
    LJ = Jeans.L( nH, Tcen, mu ) 
    
    LJ.units = 'kpc' 
    LJ_old.units = 'kpc'

    print 'Jeans: '
    print 'itr: ', itr
    print 'mu: ', mu
    print 'LJ: ', LJ
    print 'LJ_old: ', LJ_old
    print 

    err = np.abs( 1.0 - LJ / LJ_old )
    itr += 1
       
    LJ_old = LJ
