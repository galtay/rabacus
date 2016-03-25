import numpy as np
import rabacus as ra


def read_cool( fname ):

    """ Return a dictionary with cooling rates from Cloudy 13 assuming
    the command:  'save cooling each ".cool"' was given. """ 

    f = open( fname, 'r' )
    h = f.readline().split('\t')
    f.close()
    dat = np.loadtxt( fname )

    cool = {}
    for i,k in enumerate( h ):
        if i == 0:
            pass
        elif i == 1:
            cool['Temp'] = dat[:,i] * ra.u.K
        elif i == 2:
            cool['Ctot'] = dat[:,i] * ra.u.erg / ra.u.cm**3 / ra.u.s
        else:
            cool[k.strip()] = dat[:,i] * ra.u.erg / ra.u.cm**3 / ra.u.s

    return cool



