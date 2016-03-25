
import os
import numpy as np
import rabacus as ra


z = 3.0
cloudy_exe = '/home/galtay/Downloads/Software/c13.03/source/cloudy.exe'



def write_spectrum_file( z, NH1, NHe1, NHe2, fac=1.0 ):

    # generate HM12 spectrum
    #---------------------------------------------
    Nnu = 100
    q_min = 1.0e-2
    q_max = 1.0e6
    uvb = ra.BackgroundSource( q_min, q_max, 'hm12', z=z, Nnu=Nnu )

    tau = uvb.sigma.H1 * NH1 + uvb.sigma.He1 * NHe1 + uvb.sigma.He2 * NHe2

    interpolate_command = []

    uvb.E.units = 'Ry_inf'
    logE = np.log10( uvb.E.magnitude )
    logInu = np.log10( uvb.Inu.magnitude * fac * np.exp(-tau.magnitude) ) 

    for i in range(Nnu):
        interpolate_command.append( 
            ' (' + str(logE[i]) + ' ' + str(logInu[i]) + ') ' )

    # note f(nu) takes 4 pi Jnu in [erg/(s Hz cm^2)]
    # save continuum outputs in 4 pi nu Jnu
    # also note f(nu) interprets E_norm as log if <= 0 and linear if > 0

    f = open( 'spectrum.ini', 'w' )
    f.write( 'interpolate\n' )
    for i in range(Nnu):
        if i%1 == 0:
            f.write('continue ')
            f.write( interpolate_command[i] )
            f.write('\n')
    logInu_norm = np.log10( 4*np.pi*10**logInu[Nnu/2] )
    logE_norm = logE[Nnu/2]
    if logE_norm > 0.0:
        f.write('f(nu) = ' + str(logInu_norm) + ' at ' + 
                str(10**logE_norm) + '\n')
    else:
        f.write('f(nu) = ' + str(logInu_norm) + ' at ' + 
                str(logE_norm) + '\n')

    f.close()








# this includes the 11 elements from the OWLS cooling model
# H, He, C, N, O, Ne, Si, Mg, S, Ca, Fe
#-----------------------------------------------------------
abundance_lines_solar = [
    'abundances GASS10',
#    'element hydrogen off',
#    'element helium off', 
    'element lithium off',
    'element beryllium off',
    'element boron  off',
    'element fluorine  off',
    'element phosphor off',
    'element chlorine off',
    'element potassium off',
    'element scandium  off',
    'element titanium off',
    'element vanadium off',
    'element chromium off',
    'element manganese off',
    'element cobalt off',
    'element copper off',
    'element zinc  off',
#    'element carbon off',
#    'element nitrogen off',
#    'element oxygen off',
#    'element neon off',
    'element sodium off',
#    'element magnesium off',
    'element aluminium off',
#    'element silicon off',
#    'element sulphur off',
    'element argon off',
#    'element calcium off',
#    'element iron off',
    'element nickel off',
]


abundance_lines_primordial = [
    'abundances GASS10',
#    'element hydrogen off',
#    'element helium off', 
    'element lithium off',
    'element beryllium off',
    'element boron  off',
    'element fluorine  off',
    'element phosphor off',
    'element chlorine off',
    'element potassium off',
    'element scandium  off',
    'element titanium off',
    'element vanadium off',
    'element chromium off',
    'element manganese off',
    'element cobalt off',
    'element copper off',
    'element zinc  off',
    'element carbon off',
    'element nitrogen off',
    'element oxygen off',
    'element neon off',
    'element sodium off',
    'element magnesium off',
    'element aluminium off',
    'element silicon off',
    'element sulphur off',
    'element argon off',
    'element calcium off',
    'element iron off',
    'element nickel off',
]






log_T_lo = 3.50
log_T_hi = 7.50
d_log_T = 0.25
NT = int( np.rint( (log_T_hi - log_T_lo)/d_log_T ) + 1 ) 


cloudy_lines = [
    'init "spectrum.ini"',
    'stop zone 1',
    'iterate to convergence',
    'no molecules',
#    'no grain physics',
    'hden -2',
    'set dr 0', 
    'set WeakHeatCool 1.0e-3',
    'constant temperature, T=4 vary',
    'grid from '+str(log_T_lo) + \
    ' to ' + str(log_T_hi) + \
    ' in ' + str(d_log_T) + ' dex steps',
    'CMB redshift ' + str(z),
    'save overview ".ovr" last',
    'save continuum ".con" last',
    'save cooling each ".cool" last', 
    'save heating ".heat" last',
#    'save gammaa ".gamma" last',
    'save grain heating ".grain" last', 
    'save grid ".grd" no clobber',
    'print last',
]

# table read file = "spectrum.txt"
# coronal equilibrium

# grid command
# chemical composition
# abundances are set by number relative to hydrogen n(X)/n(H)
# continuum specified between 10 m and 100 MeV !!! 
# use CMB command
# background gives a simple approx. of UV background (includes CMB)
# title

# continuum output units are Ryd
# continuum output units are 4 pi nu Jnu [erg cm^-2 s^-1] 
# 1.001e-8 Ryd to 7.354e6 Ryd

# save cooling [each]
# save heating

# save gammas?

def write_input_file( fname, abundance_commands, cloudy_commands ):
    f = open( fname, 'w' )
    for cmnd in abundance_commands:
        f.write( cmnd + '\n' )
    for cmnd in cloudy_commands:
        f.write( cmnd + '\n' )
    f.close()




# collisional equilibrium
#-----------------------------------------------------------

z = 3.0
NH1 = 0.0 / ra.u.cm**2
NHe1 = 0.0 / ra.u.cm**2
NHe2 = 0.0 / ra.u.cm**2
write_spectrum_file( z, NH1, NHe1, NHe2, fac=1.0e-20 )

f_in = 'primordial_ce.in'
write_input_file( f_in, abundance_lines_primordial, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )

f_in = 'solar_ce.in'
write_input_file( f_in, abundance_lines_solar, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )


# photo collisional equilibrium
#-----------------------------------------------------------

z = 3.0
NH1 = 0.0 / ra.u.cm**2
NHe1 = 0.0 / ra.u.cm**2
NHe2 = 0.0 / ra.u.cm**2
write_spectrum_file( z, NH1, NHe1, NHe2 )

f_in = 'primordial_pce.in'
write_input_file( f_in, abundance_lines_primordial, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )

f_in = 'solar_pce.in'
write_input_file( f_in, abundance_lines_solar, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )




NH1 = 1.0e18 / ra.u.cm**2
NHe1 = 1.0e17 / ra.u.cm**2
NHe2 = 1.0e17 / ra.u.cm**2
write_spectrum_file( z, NH1, NHe1, NHe2 )

f_in = 'primordial_pce_18.in'
write_input_file( f_in, abundance_lines_primordial, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )

f_in = 'solar_pce_18.in'
write_input_file( f_in, abundance_lines_solar, cloudy_lines )
cmnd = cloudy_exe + ' ' + f_in 
os.system( cmnd )



