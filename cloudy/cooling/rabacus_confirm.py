import rabacus as ra
import pylab as plt
import numpy as np

z = 3.0
Nnu = 100
q_min = 1.0e-2
q_max = 1.0e6
uvb = ra.BackgroundSource( q_min, q_max, 'hm12', z=z, Nnu=Nnu )

NT=100

T = np.logspace( 4.0, 5.0, NT ) * ra.u.K

nH = np.ones( NT ) * 1.0e-2 / ra.u.cm**3
nHe = nH * 10**(-1.0701)

H1i = np.ones(T.size) * uvb.thin.H1i
He1i = np.ones(T.size) * uvb.thin.He1i
He2i = np.ones(T.size) * uvb.thin.He2i

H1h = np.ones(T.size) * uvb.thin.H1h
He1h = np.ones(T.size) * uvb.thin.He1h
He2h = np.ones(T.size) * uvb.thin.He2h

fcA_H2 = 1.0
fcA_He2 = 1.0
fcA_He3 = 1.0

kchem = ra.ChemistryRates( T, fcA_H2, fcA_He2, fcA_He3, 
                           H1i=H1i, He1i=He1i, He2i=He2i )

kcool = ra.CoolingRates( T, fcA_H2, fcA_He2, fcA_He3, 
                         H1h=H1h, He1h=He1h, He2h=He2h )



x_pce = ra.Solve_PCE( nH, nHe, kchem )

heat_H1 = H1h * nH * x_pce.H1
heat_He1 = He1h * nHe * x_pce.He1
heat_He2 = He2h * nHe * x_pce.He2  
heat = heat_H1 + heat_He1 + heat_He2


plt.loglog( T, heat, color='black', lw=2.0 )
plt.loglog( T, heat_H1 )
plt.loglog( T, heat_He1 )
plt.loglog( T, heat_He2 )
