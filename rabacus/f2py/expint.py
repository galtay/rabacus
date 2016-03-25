import numpy as np
import rabacus as ra

G_E1 = np.array( [3.0667, -18.3308, 61.4953, -130.3433, 161.0925, -104.1769, 27.2146] )

def expint_pfast( x ):
    r = x / (1.0 + x)
    e1 = 0.0
    for i in range(7):
        e1 = e1 + G_E1[i] * r**i
    return e1


def create_E1_lookup( N, x_lo=1.0e-6, x_hi=10.0 ):

    log_x_lo = np.log10( x_lo )
    log_x_hi = np.log10( x_hi )

    log_x = np.linspace( log_x_lo, log_x_hi, N )
    x = 10**log_x
    E1_table = np.ones( N )
    log_E1_table = np.ones( N )
    for i in range(N):
        E1_table[i] = ra.f2py.rabacus_fc.slab_base.expint( 1, x[i] )
    log_E1_table = np.log10( E1_table )
    return (log_x, log_E1_table)
