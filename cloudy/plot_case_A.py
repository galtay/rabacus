import time 
import numpy as np
import rabacus as ra



# setup Stromgren sphere
#=================================================================

Nl = 512
T = np.ones(Nl) * 1.0e4 * ra.U.K

Rsphere = 6.6 * ra.U.kpc
Edges = np.linspace( 0.0 * ra.U.kpc, Rsphere, Nl+1 )
nH = np.ones(Nl) * 1.0e-3 / ra.U.cm**3
nHe = np.ones(Nl) * 8.7e-5 / ra.U.cm**3
nHe_null = np.ones(Nl) * 1.0e-15 / ra.U.cm**3

Ln = 5.0e48 / ra.U.s  # set photon luminosity
Nnu = 128

q_mono = 1.2
q_min = q_mono
q_max = q_mono
src_mono = ra.rad_src.PointSource( q_min, q_max, 'monochromatic' )
src_mono.normalize_Ln( Ln )



#=================================================================
# Iliev 06 tests
#=================================================================


# radiative transfer case A mono fix T
#-----------------------------------------------------------------
t1 = time.time()
s = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=1.0 )
t2 = time.time()
print 'time: ', str(t2-t1)


plt.figure()
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red', ls='-' )
fname = 'raicevic_case_A.ovr'
dat = np.loadtxt( fname )
plt.plot( dat[:,0], dat[:,6], color='blue', ls='--' )

plt.show()

