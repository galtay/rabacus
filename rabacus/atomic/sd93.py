""" 
A module that stores rate fits from Sutherland and Dopita 93
http://adsabs.harvard.edu//abs/1993ApJS...88..253S
"""

import os
import numpy as np
from rabacus.constants import physical
from rabacus.constants import units
from scipy.interpolate import interp1d



__all__ = ['SD93']



class SD93:

    r""" 
    Fits to atomic data from 
    http://adsabs.harvard.edu//abs/1993ApJS...88..253S
    Provides acess to Lambda_N from Table 6
    """


    def __init__(self):

        """ Read and clean data. """ 

        # create physical constants and units 
        #--------------------------------------------------------
        self.pc = physical.PhysicalConstants()
        self.u = units.Units()

        # read data
        #-----------------------------------------------------
        local = os.path.dirname(os.path.realpath(__file__))
        fname = local + '/sd93/zcool_sd93.dat'        
        dat = np.loadtxt( fname )
        self._dat = dat
        self._logZ = np.array( [-3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5] )
        self._NZ = self._logZ.size

        # get temperatures
        #-----------------------------------------------------
        self._logT = dat[:,0]
        self._T = 10**self._logT * self.u.K

        # get Lambda.  this is the raw data from the file
        # column 1 is temperature 2-9 are metallicity starting 
        # with nil and then moving through the _logZ array above
        #-----------------------------------------------------
        lamu = self.u.cm**3 * self.u.erg / self.u.s
        self._logLambda = np.array( dat[:,1:] ) 
        self._Lambda = 10**self._logLambda * lamu

        # make copies of zero metallicity column
        #--------------------------------------------------------------
        Z0 = np.array( [self._dat[:,1] for i in range(self._NZ+1)] )
        Z0 = Z0.transpose()
        Zf = dat[:,1:]


        # subtract the nill values from the finite metallicity values
        #--------------------------------------------------------------
        self._Z0 = Z0
        self._Zf = Zf
        self._Lambda_Zonly = 10**self._Zf - 10**self._Z0

        # set minimum value to make log safe and add units
        #--------------------------------------------------------------
        indx = np.where( self._Lambda_Zonly <= 0.0 )
        self._Lambda_Zonly[indx] = 1.0e-50
        self._logLambda_Zonly = np.log10( self._Lambda_Zonly )
        self._Lambda_Zonly = self._Lambda_Zonly * lamu



    def _check_input_T(self, T):
        msg = '\n Input variable T must have units of K \n'
        if not hasattr(T,'units'): 
            raise ValueError(msg)
            if not T.units == 'K':
                raise ValueError(msg)

        if T.ndim == 0:
            T = np.array( [T] ) * T.units

        return T


    def Zcool( self, Tin, logZ ):
        """ Interpolate metal cooling function at a fixed metallicity and 
        many temperatures. 

        Args:
          `Tin` (array like): temperatures
        
          `logZ` (float): log metallicity 
        """ 


        # set iZ and Z_frac
        #------------------------------
        if logZ <= self._logZ[0]:
            iZ = 0
            Z_frac = 0.0
        elif logZ >= self._logZ[-1]:
            iZ = self._NZ - 2
            Z_frac = 1.0
        else:
            iZ = 0
            for i in range(self._NZ):
                if logZ < self._logZ[i]:
                    iZ = i-1
                    break
            dZ = self._logZ[iZ+1] - self._logZ[iZ]
            Z_frac = ( logZ - self._logZ[iZ] ) / dZ


        # interpolate cooling between two metallicities
        #------------------------------------------------
        Lambda_lo = self._Lambda_Zonly[:,iZ]
        Lambda_hi = self._Lambda_Zonly[:,iZ+1]
        Lambda = Lambda_lo + Z_frac * ( Lambda_hi - Lambda_lo )
        logLambda = np.log10( Lambda.magnitude )

        # create interpolating function lam(T)
        #------------------------------------------------
        interpf = interp1d( self._logT, logLambda )


        # clean and clip input temperature
        #------------------------------
        T = self._check_input_T(Tin)        
        logT = np.log10( T.magnitude )
        logT = np.clip( logT, self._logT[0], self._logT[-1] )

        log_lam = interpf( logT )
        lam = 10**log_lam
        return lam * self._Lambda_Zonly.units


    def cool( self, Tin, logZ ):
        """ Interpolate cooling function at a fixed metallicity and 
        many temperatures. 

        Args:
          `Tin` (array like): temperatures
        
          `logZ` (float): log metallicity 
        """ 

        # set iZ and Z_frac
        #------------------------------
        if logZ <= self._logZ[0]:
            iZ = 0
            Z_frac = 0.0
        elif logZ >= self._logZ[-1]:
            iZ = self._NZ - 2
            Z_frac = 1.0
        else:
            iZ = 0
            for i in range(self._NZ):
                if logZ < self._logZ[i]:
                    iZ = i-1
                    break
            dZ = self._logZ[iZ+1] - self._logZ[iZ]
            Z_frac = ( logZ - self._logZ[iZ] ) / dZ


        # interpolate cooling between two metallicities
        #------------------------------------------------
        Lambda_lo = self._Lambda[:,iZ]
        Lambda_hi = self._Lambda[:,iZ+1]
        Lambda = Lambda_lo + Z_frac * ( Lambda_hi - Lambda_lo )
        logLambda = np.log10( Lambda.magnitude )

        # create interpolating function lam(T)
        #------------------------------------------------
        interpf = interp1d( self._logT, logLambda )


        # clean and clip input temperature
        #------------------------------
        T = self._check_input_T(Tin)        
        logT = np.log10( T.magnitude )
        logT = np.clip( logT, self._logT[0], self._logT[-1] )

        log_lam = interpf( logT )
        lam = 10**log_lam
        return lam * self._Lambda.units

