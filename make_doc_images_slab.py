import time
import numpy as np
import pylab as plt
import rabacus as ra


def plot_slab_analytic( s, ana, fname ):
    """ put plot of analytic vs. rabacus solution into fname """ 

    plt.figure( figsize=(7,7) )
    s.Edges.units = 'kpc'
    s.z_c.units = 'kpc'
    xx = s.z_c 
    L = s.Edges[-1]

    ana.L.units = 'kpc'

    plt.plot( xx, np.log10( s.xH1 ), 
              color='red', ls='-', label = r'$x_{\rm HI}$' )
    plt.plot( xx, np.log10( s.xH2 ), 
              color='red', ls='--', label = r'$x_{\rm HII}$' )
        
    plt.scatter( ana.L, np.log10( ana.xH1 ), s=20, 
                 color='black', marker='o', label='analytic' )

    plt.xlim( -L/20, L+L/20 )
    plt.xlabel( 'z_c [kpc]' )

    plt.ylim( -4.0, 0.2 )
    plt.ylabel( 'log 10 ( x )' )

    plt.grid()
    plt.legend(loc='best', ncol=2)
    plt.tight_layout()
    plt.savefig( 'doc/img/x_' + fname )


def plot_slab_x( s, fname, plot_H=True, plot_He=True ):
    """ put plot of ionization fractions from slab `s` into fname """ 

    plt.figure( figsize=(7,7) )
    s.Edges.units = 'kpc'
    s.z_c.units = 'kpc'
    xx = s.z_c 
    L = s.Edges[-1]

    if plot_He:
        plt.plot( xx, np.log10( s.xHe1 ), 
                  color='green', ls='-', label = r'$x_{\rm HeI}$' )
        plt.plot( xx, np.log10( s.xHe2 ), 
                  color='green', ls='--', label = r'$x_{\rm HeII}$' )
        plt.plot( xx, np.log10( s.xHe3 ), 
                  color='green', ls=':', label = r'$x_{\rm HeIII}$' )

    if plot_H:
        plt.plot( xx, np.log10( s.xH1 ), 
                  color='red', ls='-', label = r'$x_{\rm HI}$' )
        plt.plot( xx, np.log10( s.xH2 ), 
                  color='red', ls='--', label = r'$x_{\rm HII}$' )
        

    plt.xlim( -L/20, L+L/20 )
    plt.xlabel( 'z_c [kpc]' )

    plt.ylim( -4.0, 0.2 )
    plt.ylabel( 'log 10 ( x )' )

    plt.grid()
    plt.legend(loc='best', ncol=2)
    plt.tight_layout()
    plt.savefig( 'doc/img/x_' + fname )


def plot_slab_T( s, fname ):
    """ put plot of temperature from slab `s` into fname """ 

    plt.figure( figsize=(7,7) )
    s.Edges.units = 'kpc'
    s.z_c.units = 'kpc'
    xx = s.z_c 
    L = s.Edges[-1]

    plt.plot( xx, s.T, 
              color='black', ls='-', label = r'$T$' )

    plt.xlim( -L/20, L+L/20 )
    plt.xlabel( 'z_c [kpc]' )

    plt.ylim( 8.0e3, 2.6e4 )
    plt.ylabel( 'T [K]' )

    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig( 'doc/img/T_' + fname )





# setup slab
#=================================================================

Nl = 512
T = np.ones(Nl) * 1.0e4 * ra.u.K

Lslab = 1.5e2 * ra.u.kpc
Edges = np.linspace( 0.0 * ra.u.kpc, Lslab, Nl+1 )
Yp = 0.24
nH = np.ones(Nl) * 2.2e-3 / ra.u.cm**3
nHe = nH * 0.25 * Yp / ( 1.0 - Yp )
nHe_null = np.ones(Nl) * 1.0e-15 / ra.u.cm**3


# setup sources 
#=================================================================
z = 3.0

q_min = 1.0
q_max = 4.0e2
src_hm12 = ra.PlaneSource( q_min, q_max, 'hm12', z=z )

q_mono = src_hm12.grey.E.H1 / src_hm12.th.E_H1
q_min = q_mono
q_max = q_mono
src_mono = ra.PlaneSource( q_min, q_max, 'monochromatic' )
src_mono.normalize_H1i( src_hm12.thin.H1i )



# solve slabs
#=================================================================


slabs = {}

# Slab 1 Tests
#=================================================================


# optically thin case B mono fix T
#-----------------------------------------------------------------
key = 'thin_caseB_mono_fixT'
t1 = time.time()
slabs[key] = ra.SlabPln( Edges, T, nH, nHe, src_mono, 
                           fixed_fcA=0.0, thin=True )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png' )



# radiative transfer case B mono fix T
#-----------------------------------------------------------------
ana = ra.AnalyticSlab( 
    nH[0], T[0], src_hm12.thin.H1i, y=0.0, fcA=0.0 )
ana.set_E( src_hm12.grey.E.H1 )

key = 'rt_caseB_mono_fixT'
t1 = time.time()
slabs[key] = ra.SlabPln( Edges, T, nH, nHe, src_mono, 
                           fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_analytic( slabs[key], ana, 'slab_analytic.png', )


# add polychromatic source. 
#-----------------------------------------------------------------
key = 'rt_caseB_hm12_fixT'
t1 = time.time()
slabs[key] = ra.SlabPln( Edges, T, nH, nHe, src_hm12, 
                           fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )


# radiative transfer case B hm12 evo T
#-----------------------------------------------------------------
key = 'rt_caseB_hm12_evoT'
t1 = time.time()
slabs[key] = ra.SlabPln( Edges, T, nH, nHe, src_hm12, 
                           fixed_fcA=0.0, find_Teq=True, z=z )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )
plot_slab_T( slabs[key], 'slab_' + key + '.png', )


# radiative transfer ray hm12 evo T
#-----------------------------------------------------------------
key = 'rt_ray_hm12_evoT'
t1 = time.time()
slabs[key] = ra.SlabPln( Edges, T, nH, nHe, src_hm12, 
                           rec_meth='ray', find_Teq=True, z=z )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )
plot_slab_T( slabs[key], 'slab_' + key + '.png', )



# Slab 2 Tests
#=================================================================


# radiative transfer case B mono fix T
#-----------------------------------------------------------------
key = '2_rt_caseB_mono_fixT'
t1 = time.time()
slabs[key] = ra.Slab2Pln( Edges, T, nH, nHe, src_mono, 
                            fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )


# add polychromatic source. 
#-----------------------------------------------------------------
key = '2_rt_caseB_hm12_fixT'
t1 = time.time()
slabs[key] = ra.Slab2Pln( Edges, T, nH, nHe, src_hm12, 
                            fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )


# radiative transfer case B hm12 evo T
#-----------------------------------------------------------------
key = '2_rt_caseB_hm12_evoT'
t1 = time.time()
slabs[key] = ra.Slab2Pln( Edges, T, nH, nHe, src_hm12, 
                            fixed_fcA=0.0, find_Teq=True, z=z )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )
plot_slab_T( slabs[key], 'slab_' + key + '.png', )


# radiative transfer ray hm12 evo T
#-----------------------------------------------------------------
key = '2_rt_ray_hm12_evoT'
t1 = time.time()
slabs[key] = ra.Slab2Pln( Edges, T, nH, nHe, src_hm12, 
                            rec_meth='ray', find_Teq=True, z=z )
t2 = time.time()
print key + ' ' + str(t2-t1)
plot_slab_x( slabs[key], 'slab_' + key + '.png', )
plot_slab_T( slabs[key], 'slab_' + key + '.png', )
