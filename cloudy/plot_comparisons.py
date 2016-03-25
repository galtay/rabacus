import time 
import numpy as np
import rabacus as ra
import scipy


# setup Stromgren sphere
#=================================================================

Rsphere = 6.6 * ra.U.kpc

Ln = 5.0e48 / ra.U.s  # set photon luminosity
Nnu = 128

q_mono = 1.1
q_min = q_mono
q_max = q_mono
src_mono = ra.rad_src.PointSource( q_min, q_max, 'monochromatic' )
src_mono.normalize_Ln( Ln )

q_min = 1.0
q_max = 10.0
T_eff = 1.0e5 * ra.U.K
src_thrm = ra.rad_src.PointSource( q_min, q_max, 'thermal', T_eff=T_eff,
                                   Nnu=Nnu)
src_thrm.normalize_Ln( Ln )


spheres = {}


def read_cloudy( fbase ):

    fname = fbase + '.out'
    f = open( fname, 'r' )
    lines = f.readlines()
    f.close()
    line = lines[-2]
    nzone = int( line.partition(':')[2].partition( 'zones' )[0] )
    
    fname = fbase + '.ovr'
    dat_all = np.loadtxt( fname )
    dat = dat_all[-nzone:,:]
    return dat


def set_sphere( dat ):

    Nl = dat[:,0].size 
    T = np.ones(Nl) * 1.0e4 * ra.U.K
    nH = np.ones(Nl) * 1.0e-3 / ra.U.cm**3
    nHe = np.ones(Nl) * 8.7e-5 / ra.U.cm**3
    nHe_null = np.ones(Nl) * 1.0e-15 / ra.U.cm**3

    r_c = dat[:,0] * ra.U.cm
    Edges = np.zeros( Nl+1 ) * ra.U.kpc
    Edges[0] = 0.0 * ra.U.cm
    Edges[-1] = Rsphere
    for i in range(1,Nl):
        Edges[i] = ( r_c[i-1] + r_c[i] ) * 0.5
    return (Edges, T, nH, nHe, nHe_null) 

#=================================================================
# Iliev 06 tests
#=================================================================


# radiative transfer case B mono fix T
#-----------------------------------------------------------------
key = 'iliev_test_1'
dat = read_cloudy( 'iliev_test_1' )
(Edges, T, nH, nHe, nHe_null) = set_sphere( dat )
t1 = time.time()
spheres[key] = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)

plt.figure()
plt.plot( dat[:,0], dat[:,6], color='blue' )
s = spheres['iliev_test_1']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red' )
plt.show()




# radiative transfer case B thermal
#-----------------------------------------------------------------
key = 'iliev_test_2'
dat = read_cloudy( 'iliev_test_2' )
(Edges, T, nH, nHe, nHe_null) = set_sphere( dat )
t1 = time.time()
spheres[key] = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_thrm, find_Teq=True, z=0.0, fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)

plt.figure()
plt.plot( dat[:,0], dat[:,6], color='blue' )
s = spheres['iliev_test_2']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red' )
plt.show()

plt.figure()
plt.plot( dat[:,0], dat[:,1], color='blue' )
s = spheres['iliev_test_2']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.T.magnitude ), color='red' )
plt.show()





#=================================================================
# Raiveciv 14 Tests
#=================================================================

# radiative transfer case A mono fix T
#-----------------------------------------------------------------
key = 'raicevic_case_A'
dat = read_cloudy( 'raicevic_case_A' )
(Edges, T, nH, nHe, nHe_null) = set_sphere( dat )
t1 = time.time()
spheres[key] = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=1.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)

# radiative transfer case B mono fix T
#-----------------------------------------------------------------
key = 'raicevic_case_B'
dat = read_cloudy( 'raicevic_case_B' )
(Edges, T, nH, nHe, nHe_null) = set_sphere( dat )
t1 = time.time()
spheres[key] = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=0.0 )
t2 = time.time()
print key + ' ' + str(t2-t1)

# radiative transfer outward mono fix T
#-----------------------------------------------------------------
key = 'raicevic_outward'
dat = read_cloudy( 'raicevic_outward' )
(Edges, T, nH, nHe, nHe_null) = set_sphere( dat )
t1 = time.time()

Nl = dat[:,0].size 
T = np.ones(Nl) * 1.0e4 * ra.U.K
nH = np.ones(Nl) * 1.0e-3 / ra.U.cm**3
nHe = np.ones(Nl) * 8.7e-5 / ra.U.cm**3
nHe_null = np.ones(Nl) * 1.0e-15 / ra.U.cm**3
Edges = np.linspace( 0.0 * ra.U.kpc, Rsphere, Nl+1 )

spheres[key] = ra.f2py.StromgrenSphere( 
    Edges, T, nH, nHe_null, src_mono, rec_meth='outward' )
t2 = time.time()
print key + ' ' + str(t2-t1)



plt.figure()
s = spheres['raicevic_case_A']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red', ls='-' )
dat = read_cloudy( 'raicevic_case_A' )
plt.plot( dat[:,0], dat[:,6], color='blue', ls='-' )

s = spheres['raicevic_case_B']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red', ls=':' )
dat = read_cloudy( 'raicevic_case_B' )
plt.plot( dat[:,0], dat[:,6], color='blue', ls=':' )

s = spheres['raicevic_outward']
s.r_c.units = 'cm'
plt.plot( s.r_c, np.log10( s.xH1 ), color='red', ls='--' )
dat = read_cloudy( 'raicevic_outward' )
plt.plot( dat[:,0], dat[:,6], color='blue', ls='--' )

plt.xlim( 0.5e22, 2.1e22 )
plt.ylim( -2.5, 0.1 )



plt.show()

sys.exit(1)
