import rabacus as ra
import numpy as np
import pylab as plt

con_fname = 'test.con'
dat = np.loadtxt( con_fname, usecols=(0,1) )
cldy_Ry = dat[:,0] * ra.u.Ry_inf
cldy_nu = cldy_Ry / ra.pc.h
cldy_nu.units = 'Hz'
cldy_4pi_nu_Jnu = dat[:,1] * ra.u.erg / (ra.u.s * ra.u.cm**2)
cldy_Jnu = cldy_4pi_nu_Jnu / (4.0*np.pi*ra.u.sr * cldy_nu)

z = 3.0
hm12 = ra.HM12_UVB_Table()
indx = np.where( cldy_Ry > hm12.E[-1] )
ra_Ry = cldy_Ry[indx]

ra_Jnu = hm12.return_spectrum_E( z, ra_Ry )



NH1 = 1.0e20 / ra.u.cm**2
sigmaH1 = uvb.PX.sigma_H1( ra_Ry )
tauH1 = NH1 * sigmaH1

NHe1 = 1.0e16 / ra.u.cm**2
sigmaHe1 = uvb.PX.sigma_He1( ra_Ry )
tauHe1 = NHe1 * sigmaHe1




tau = tauH1 + tauHe1

ra_Ry.units = 'Ry_inf'

#plt.loglog( cldy_Ry, cldy_4pi_nu_Jnu )
plt.loglog( cldy_Ry, cldy_Jnu, color='blue' )
plt.loglog( ra_Ry, ra_Jnu, ls='--', color='red' )
plt.loglog( ra_Ry, ra_Jnu*np.exp(-tau), ls='--', color='green' )



