import sys
import numpy as np
import pylab as plt



#
# read overview file
#----------------------------------------
def read_overview_file( fbase, NT ):

    tmprtr = np.zeros( (NT) )
    xH1 = np.zeros( (NT) )
    xH2 = np.zeros( (NT) )
    xHe1 = np.zeros( (NT) )
    xHe2 = np.zeros( (NT) )
    xHe3 = np.zeros( (NT) )

    fname = fbase+'.ovr'
    f = open( fname, 'r' )
    ovr_header = f.readline().split('\t')
    ovr_lines = f.readlines()
    f.close()

    for i,line in enumerate(ovr_lines):
        ovr_lines[i] = ovr_lines[i].split('\t')
    
    for i,line in enumerate(ovr_lines):
        #print line
        if len(line) != 1:
            #print '**', i/2, line[1], line[2]
            ii = i/2
            tmprtr[ii] = line[1]
            xH1[ii] = line[6]
            xH2[ii] = line[7]
            xHe1[ii] = line[8]
            xHe2[ii] = line[9]
            xHe3[ii] = line[10]


            
    return 10**tmprtr, 10**xH1, 10**xH2, 10**xHe1, 10**xHe2, 10**xHe3


#
# read continuum file
#----------------------------------------
def read_continuum_file( fbase ):

    fname = fbase+'.con'
    f = open( fname, 'r' )
    dum = f.readline()
    break_check = ''
    E_list = []
    I_list = []
    while break_check != '#':
        line = f.readline()
        break_check=line[0]
        line_bits = line.split()
        if break_check != '#':
            E_list.append( np.float(line_bits[0]) )
            I_list.append( np.float(line_bits[1]) )
        #print line.split()

    E_Ry = np.array( E_list )
    I_nu = np.array( I_list )

    return E_Ry, I_nu

#
# read grid file
#----------------------------------------
def read_grid_file( fbase ):

    fname = fbase+'.grd'
    f = open( fname, 'r' )
    grd_header = f.readline().split('\t')
    grd_lines = f.readlines()
    f.close()

    NT = len(grd_lines)
    grd_keys = []
    consT = {}
    for line in grd_lines:
        split_line = line.split('\t')
        key = 'grid'+split_line[0]
        grd_keys.append( key )
        consT[key] = split_line[6]

    return NT


#
# read cooling file
#----------------------------------------
def read_cooling_file( fbase, NT ):

    tmprtr = np.zeros( (NT) )
    ctot = np.zeros( (NT) )
    cH = np.zeros( (NT) )
    cHe = np.zeros( (NT) )
    ccomp = np.zeros( (NT) )
    ceeff = np.zeros( (NT) )
    cFFcm = np.zeros( (NT) )
#    chvFB = np.zeros( (NT) )
#    cH2p = np.zeros( (NT) )
#    cHDro = np.zeros( (NT) )
#    cH2ln = np.zeros( (NT) )
#    cHdfb = np.zeros( (NT) )
#    cCTsC = np.zeros( (NT) )
#    cH2cX = np.zeros( (NT) )
#    cdust = np.zeros( (NT) )
#    cmolecule = np.zeros( (NT) )

    fname = fbase+'.cool'
    f = open( fname, 'r' )
    cool_header = f.readline().split('\t')
    cool_lines = f.readlines()
    f.close()

    for i,line in enumerate(cool_lines):
        cool_lines[i] = cool_lines[i].split('\t')
    
    for i,line in enumerate(cool_lines):
        #print line
        if len(line) != 1:
            #print '**', i/2, line[1], line[2]
            ii = i/2
            tmprtr[ii] = line[1]
            ctot[ii] = line[2]
            cH[ii] = line[3]
            cHe[ii] = line[4]

            ccomp[ii] = line[44]
            ceeff[ii] = line[43]
            cFFcm[ii] = line[41]
            #chvFB[ii] = line[42]
            #cH2p[ii] = line[40]
            #cHDro[ii] = line[39]
            #cH2ln[ii] = line[38]
            #cHdfb[ii] = line[37]
            #cCTsC[ii] = line[36]
            #cH2cX[ii] = line[35]
            #cdust[ii] = line[34]
            #cmolecule[ii] = line[33]
        
        # this line describes the grid number
        else:
            key = line[0].split('--')[1].strip()
            #print 'key = ', key



    return tmprtr, ctot, cH, cHe, ccomp, ceeff, cFFcm




#
# read heating file
#----------------------------------------
def read_heating_file( fbase, NT ):

    tmprtr = np.zeros( (NT) )
    htot = np.zeros( (NT) )
    ctot = np.zeros( (NT) )
    hH1 = np.zeros( (NT) )
    hHe1 = np.zeros( (NT) )
    hHe2 = np.zeros( (NT) )

    fname = fbase+'.heat'
    f = open( fname, 'r' )
    heat_header = f.readline().split('\t')
    heat_lines = f.readlines()
    f.close()

    for i,line in enumerate(heat_lines):
        heat_lines[i] = heat_lines[i].split('\t')
    
    for i,line in enumerate(heat_lines):
        #print line
        if len(line) != 1:
            #print '**', i/2, line[1], line[2]
            ii = i/2
            tmprtr[ii] = line[1]
            htot[ii] = line[2]
            ctot[ii] = line[3]

            for i,entry in enumerate(line):

                if entry == 'H  1':
                    hH1[ii] = np.float(line[i+1]) * htot[ii]
                if entry == 'He 1':
                    hHe1[ii] = np.float(line[i+1]) * htot[ii]
                if entry == 'He 2':
                    hHe2[ii] = np.float(line[i+1]) * htot[ii]

                #hH1[ii] = np.float(line[5]) * htot[ii]
                #hHe1[ii] = np.float(line[7]) * htot[ii]
                #hHe2[ii] = np.float(line[9]) * htot[ii]

    return tmprtr, htot, ctot, hH1, hHe1, hHe2



nH = 10**(-2)

fig = plt.figure( figsize=(20,10) )

ax1 = plt.subplot2grid( (2,4), (0, 0) )   
ax2 = plt.subplot2grid( (2,4), (1, 0) )   

ax3 = plt.subplot2grid( (2,4), (0, 1) )
ax4 = plt.subplot2grid( (2,4), (1, 1) )

ax5 = plt.subplot2grid( (2,4), (0, 2) )
ax6 = plt.subplot2grid( (2,4), (1, 2) )

ax7 = plt.subplot2grid( (2,4), (0, 3) )
ax8 = plt.subplot2grid( (2,4), (1, 3) )



# primordial gas cooling
#----------------------------------------------------------
fbase = 'primordial_pce'
NT = read_grid_file( fbase )
tmprtr, ctot, cH, cHe, ccomp, ceeff, cFFcm = read_cooling_file( fbase, NT )
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax1.loglog( tmprtr, htot / nH**2, 
            color='red', lw=3.0, ls='--', label=r'$\mathcal{H}$' )
ax1.loglog( tmprtr, ctot / nH**2, 
            color='blue', lw=3.0, ls='-', label=r'$\Lambda$' )
ax1.loglog( tmprtr, np.abs(ctot-htot) / nH**2, 
            color='black', ls='--', lw=3.0, 
            label=r'$|\Lambda - \mathcal{H}|$' )


fbase = 'primordial_pce_18'
NT = read_grid_file( fbase )
tmprtr, ctot, cH, cHe, ccomp, ceeff, cFFcm = read_cooling_file( fbase, NT )
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax1.loglog( tmprtr, htot / nH**2, color='red', ls='--', lw=1.0 )
ax1.loglog( tmprtr, ctot / nH**2, color='blue', lw=1.0, ls='-' )
ax1.loglog( tmprtr, np.abs(ctot-htot) / nH**2, color='black', lw=1.0, ls='--' )


ax1.set_xlim( 1.0e4/2, 1.0e5 )
ax1.set_ylim( 1.0e-25, 1.0e-21 )
#ax1.set_xlabel( 'T [K]', fontsize=20 )
ax1.set_xticklabels( [] )
ax1.set_ylabel( r'$\Lambda / n_{\rm H}^2$', fontsize=20 )


# solar gas cooling
#----------------------------------------------------------
fbase = 'solar_pce'
NT = read_grid_file( fbase )
tmprtr, ctot, cH, cHe, ccomp, ceeff, cFFcm = read_cooling_file( fbase, NT )
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax2.loglog( tmprtr, htot / nH**2, color='red', lw=3.0, ls='--' )
ax2.loglog( tmprtr, ctot / nH**2, color='blue', lw=3.0, ls='-' )
ax2.loglog( tmprtr, np.abs(ctot-htot) / nH**2, color='black', ls='--', lw=3.0 )



fbase = 'solar_pce_18'
NT = read_grid_file( fbase )
tmprtr, ctot, cH, cHe, ccomp, ceeff, cFFcm = read_cooling_file( fbase, NT )
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax2.loglog( tmprtr, htot / nH**2, color='red', ls='--', lw=1.0 )
ax2.loglog( tmprtr, ctot / nH**2, color='blue', lw=1.0, ls='-' )
ax2.loglog( tmprtr, np.abs(ctot-htot) / nH**2, color='black', lw=1.0, ls='--' )


ax2.set_xlim( 1.0e4/2, 1.0e5 )
ax2.set_ylim( 1.0e-25, 1.0e-21 )
ax2.set_xlabel( 'T [K]', fontsize=20 )
ax2.set_ylim( 1.0e-25, 1.0e-21 )
ax2.set_ylabel( r'$\Lambda / n_{\rm H}^2$', fontsize=20 )



# ionization fractions primordial
#----------------------------------------------------------
fbase = 'primordial_pce'
tmprtr, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file( fbase, NT )
ax3.loglog( tmprtr, xH1, lw=3.0,  color='red', ls='-' )   
ax3.loglog( tmprtr, xH2, lw=3.0,  color='red', ls='--' ) 
ax3.loglog( tmprtr, xHe1, lw=3.0, color='blue', ls='-' )  
ax3.loglog( tmprtr, xHe2, lw=3.0, color='blue', ls='--' ) 
ax3.loglog( tmprtr, xHe3, lw=3.0, color='blue', ls=':' )

fbase = 'primordial_pce_18'
tmprtr, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file( fbase, NT )
ax3.loglog( tmprtr, xH1, lw=1.0, color='red', ls='-' )   
ax3.loglog( tmprtr, xH2, lw=1.0, color='red', ls='--' ) 
ax3.loglog( tmprtr, xHe1, lw=1.0,color='blue', ls='-' )  
ax3.loglog( tmprtr, xHe2, lw=1.0,color='blue', ls='--' )  
ax3.loglog( tmprtr, xHe3, lw=1.0,color='blue', ls=':' ) 

ax3.set_xlim( 1.0e4/2, 1.0e5 )
ax3.set_ylim( 1.0e-5, 2.0 )
ax3.set_xticklabels( [] )
ax3.set_ylabel( 'x' )



# ionization fractions primordial
#----------------------------------------------------------
fbase = 'solar_pce'
tmprtr, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file( fbase, NT )
ax4.loglog( tmprtr, xH1, lw=3.0,  color='red', ls='-', label='HI' )   
ax4.loglog( tmprtr, xH2, lw=3.0,  color='red', ls='--', label='HII' ) 
ax4.loglog( tmprtr, xHe1, lw=3.0, color='blue', ls='-', label='HeI' )  
ax4.loglog( tmprtr, xHe2, lw=3.0, color='blue', ls='--', label='HeII' ) 
ax4.loglog( tmprtr, xHe3, lw=3.0, color='blue', ls=':', label='HeIII' )

fbase = 'solar_pce_18'
tmprtr, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file( fbase, NT )
ax4.loglog( tmprtr, xH1, lw=1.0, color='red', ls='-' )   
ax4.loglog( tmprtr, xH2, lw=1.0, color='red', ls='--' ) 
ax4.loglog( tmprtr, xHe1, lw=1.0,color='blue', ls='-' )  
ax4.loglog( tmprtr, xHe2, lw=1.0,color='blue', ls='--' )  
ax4.loglog( tmprtr, xHe3, lw=1.0,color='blue', ls=':' ) 

ax4.set_xlim( 1.0e4/2, 1.0e5 )
ax4.set_ylim( 1.0e-5, 2.0 )
ax4.set_xlabel( 'T [K]', fontsize=20 )
ax4.set_ylabel( 'x' )


# spectrum - continuum
#----------------------------------------------------------
fbase = 'primordial_pce'
E_Ry, I_nu = read_continuum_file( fbase )
ax5.loglog( E_Ry, I_nu, color='black', lw=3.0, ls='-' )
ax5.set_xlim( 1.0e-2, 1.0e2 )
ax5.set_ylim( 1.0e-7, 1.0e-2 )

fbase = 'primordial_pce_18'
E_Ry, I_nu = read_continuum_file( fbase )
ax5.loglog( E_Ry, I_nu, color='black', lw=1.0, ls='-' )
ax5.set_xlim( 2.0e-1, 5.0e1 )
ax5.set_ylim( 1.0e-8, 2.0e-3 )

ax5.set_xticklabels( [] )
ax5.set_ylabel( r'$4 \pi \, \nu \, J_{\nu}$', fontsize=20 ) 


# spectrum - continuum
#----------------------------------------------------------
fbase = 'solar_pce'
E_Ry, I_nu = read_continuum_file( fbase )
ax6.loglog( E_Ry, I_nu, color='black', lw=3.0, ls='-' )
ax6.set_xlim( 1.0e-2, 1.0e2 )
ax6.set_ylim( 1.0e-7, 1.0e-2 )

fbase = 'solar_pce_18'
E_Ry, I_nu = read_continuum_file( fbase )
ax6.loglog( E_Ry, I_nu, color='black', lw=1.0, ls='-' )
ax6.set_xlim( 2.0e-1, 5.0e1 )
ax6.set_ylim( 1.0e-8, 2.0e-3 )

ax6.set_xlabel( 'E [Ry]' )
ax6.set_ylabel( r'$4 \pi \, \nu \, J_{\nu}$', fontsize=20 ) 




# heating breakdown
#----------------------------------------------------------
fbase = 'primordial_pce'
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax7.loglog( tmprtr, htot / nH**2, 
            color='black', lw=3.0, ls='-', label=r'$\mathcal{H}$' )
ax7.loglog( tmprtr, hH1 / nH**2, 
            color='red', lw=3.0, ls='-', label=r'${\rm HI}$' )
ax7.loglog( tmprtr, hHe1 / nH**2, 
            color='blue', lw=3.0, ls='-', label=r'${\rm HeI}$' )
ax7.loglog( tmprtr, hHe2 / nH**2, 
            color='blue', lw=3.0, ls='--', label=r'${\rm HeII}$' )




fbase = 'primordial_pce_18'
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax7.loglog( tmprtr, htot / nH**2, 
            color='black', lw=1.0, ls='-', label=r'$\mathcal{H}$' )
ax7.loglog( tmprtr, hH1 / nH**2, 
            color='red', lw=1.0, ls='-', label=r'${\rm HI}$' )
ax7.loglog( tmprtr, hHe1 / nH**2, 
            color='blue', lw=1.0, ls='-', label=r'${\rm HeI}$' )
ax7.loglog( tmprtr, hHe2 / nH**2, 
            color='blue', lw=1.0, ls='--', label=r'${\rm HeII}$' )

ax7.set_xlim( 1.0e4/2, 1.0e5 )
ax7.set_xticklabels( [] )
ax7.set_ylim( 1.0e-27, 5.0e-23 )





fbase = 'solar_pce'
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax8.loglog( tmprtr, htot / nH**2, 
            color='black', lw=3.0, ls='-', label=r'$\mathcal{H}$' )
ax8.loglog( tmprtr, hH1 / nH**2, 
            color='red', lw=3.0, ls='-', label=r'${\rm HI}$' )
ax8.loglog( tmprtr, hHe1 / nH**2, 
            color='blue', lw=3.0, ls='-', label=r'${\rm HeI}$' )
ax8.loglog( tmprtr, hHe2 / nH**2, 
            color='blue', lw=3.0, ls='--', label=r'${\rm HeII}$' )


fbase = 'solar_pce_18'
tmprtr, htot, ctot, hH1, hHe1, hHe2 = read_heating_file( fbase, NT )
ax8.loglog( tmprtr, htot / nH**2, 
            color='black', lw=1.0, ls='-', label=r'$\mathcal{H}$' )
ax8.loglog( tmprtr, hH1 / nH**2, 
            color='red', lw=1.0, ls='-', label=r'${\rm HI}$' )
ax8.loglog( tmprtr, hHe1 / nH**2, 
            color='blue', lw=1.0, ls='-', label=r'${\rm HeI}$' )
ax8.loglog( tmprtr, hHe2 / nH**2, 
            color='blue', lw=1.0, ls='--', label=r'${\rm HeII}$' )



ax8.set_xlim( 1.0e4/2, 1.0e5 )
ax8.set_xlabel( 'T [K]', fontsize=20 )
ax8.set_ylim( 1.0e-26, 5.0e-23 )





xtxt = 0.1
ytxt = 0.85

ax1.text( xtxt, ytxt, 'primordial - HM12', horizontalalignment='left', 
          fontsize=20, verticalalignment='center', transform = ax1.transAxes )

ax2.text( xtxt, ytxt, 'solar - HM12', horizontalalignment='left', 
          fontsize=20, verticalalignment='center', transform = ax2.transAxes )


ax1.legend(ncol=3, fontsize=12)
ax4.legend(ncol=2, fontsize=12, loc='lower left')
plt.tight_layout()
