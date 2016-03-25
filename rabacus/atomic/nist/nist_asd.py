""" 
Atomic Data from the NIST Atomic Spectra Database (NIST ASD) 
http://physics.nist.gov/asd
"""

import numpy as np
import os.path

from rabacus.constants import physical
from rabacus.constants import units
from rabacus.utils import utils



__all__ = ['NIST_ASD']



class NIST_ASD:

    """ 
    .. note::
      Reference as ...
 
      Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2011). 
      NIST Atomic Spectra Database (ver. 5.0), [Online].
      Atomic Data from the NIST Atomic Spectra Database. 
      Available at: http://physics.nist.gov/asd [2013, July 29]. 
      National Institute of Standards and Technology, Gaithersburg, MD. 
      
    Attributes:

      `U` (:class:`~rabacus.constants.units.Units`) 

      `PC` (:class:`~rabacus.constants.physical.PhysicalConstants`) 

    """

    def __init__( self ):

        """ Reads data file and stores energies. """ 

        # create physical constants and units 
        #--------------------------------------------------------
        self.PC = physical.PhysicalConstants()
        self.U = units.Units()

        # set local directory
        #-----------------------------------------------------
        local = os.path.dirname(os.path.realpath(__file__))
        self._local = local
    
        # read and store ionization energies
        #-----------------------------------------------------
        IonizationEnergies = {}

        remove_syms = ["(", ")", "[", "]"]

        fname = local + '/ion_enrg.dat' 
        self._fname = fname

        f = open( fname, 'r' )

        for line in f.readlines():
            if not line.startswith( '#' ):
                words = line.split()
                Z = int( words[0] )
                Esym = words[1]
                Irom = words[2]
                N = Z - int( words[3] ) + 1
                Shells = words[4]
                Eraw = words[5]
                
                for sym in remove_syms:
                    Eraw = Eraw.replace( sym, "" )
                Eraw = float(Eraw)

                #print Z, Esym, Irom, Ne, Shells, Eraw * self.U.Ry_inf

                IonizationEnergies[(Z,N)] = Eraw * self.U.Ry_inf

        f.close()

        self._IonizationEnergies = IonizationEnergies 


    def set_ionization_energy_unit( self, unit ):
        """ Sets the default unit for the return value of 
        :func:`get_ionization_energy`

        Args:
          `unit` (string): a string that represents a unit, e.g. "eV"

        """

        for k,v in self._IonizationEnergies.items():
            v.units = unit
            self._IonizationEnergies[k] = v

    def get_ionization_energy(self, Z, N):
        """ Returns threshold ionization energy for ions defined by Z and N.

        Args: 

          `Z` (int): atomic number (number of protons)        

          `N` (int): electron number (number of electrons)

        Returns:

          `Eth` (float): ionization energy

        """

        assert isinstance( Z, int )
        assert isinstance( N, int )

        E = self._IonizationEnergies[(Z,N)]
        return E

        
