""" Generic access to all photoionization cross section routines. Ions are 
labeled with Z (atomic number) and N (electron number).  For convenience, we 
provide named routines for H1, He1, and He2.  All access to photoionization
cross sections should be through this module. """ 

from rabacus.utils import utils
from verner.photox import verner_photox


__all__ = ['PhotoXsections']



class PhotoXsections:
    """ General photoionization cross sections class. 

    Kwargs:

      `fit_type` (string): fits come from Verner {``verner``}
      

    Attributes:

      `Eth_H1` (float): H1 ionization threshold energy

      `Eth_He1` (float): He1 ionization threshold energy

      `Eth_He2` (float): He2 ionization threshold energy

    """ 

    

    def __init__(self, fit_type='verner' ):


        self.fit_type = fit_type
        self.VALID_FITS = ['verner']

        if self.fit_type == 'verner':

            self.v96 = verner_photox.PhotoXsections_Verner96()
            self.Eth_H1 = self.v96.return_Eth( 1,1 )
            self.Eth_He1 = self.v96.return_Eth( 2,2 )
            self.Eth_He2 = self.v96.return_Eth( 2,1 )

        else:

            msg = 'PhotoX fit_type not recognized\n' + \
                  'fit_type: ' + fit_type + '\n' + \
                  'valid fit types: ' + str(self.VALID_FITS) + '\n'
            raise utils.InputError( msg )


        self.sigma_H1th = self.sigma_th( 1, 1 )
        self.sigma_He1th = self.sigma_th( 2, 2 )
        self.sigma_He2th = self.sigma_th( 2, 1 )




    def Eth( self, Z, N ):
        """ Returns threshold ionization energy for ions defined by Z and N.

        Args:

          `Z` (int): atomic number (number of protons)        

          `N` (int): electron number (number of electrons)

        Returns:

          `Eth` (float): ionization energy

        """

        assert isinstance( Z, int )
        assert isinstance( N, int )

        if self.fit_type == 'verner':
            Eth = self.v96.return_Eth( Z, N )
        return Eth


    def sigma( self, Z, N, E ):
        """ Returns a photo-ionization cross section for an ion defined by 
        `Z` and `N` at energies `E`.

        Args:

          `Z` (int): atomic number (number of protons)        

          `N` (int): electron number (number of electrons)

          `E` (array): calculate cross section at these energies

        Returns:

          `sigma` (array): photoionization cross sections

        """
    
        assert isinstance( Z, int )
        assert isinstance( N, int )

        if self.fit_type == 'verner':
            sigma = self.v96.return_fit( Z=Z, N=N, E=E )
        return sigma


    def sigma_th( self, Z, N ):
        """ Returns photo-ionization cross sections for ions defined by Z and N
        at their ionization thresholds.

        Args:

          `Z` (int): atomic number (number of protons)        

          `N` (int): electron number (number of electrons)

        Returns:

          `sigma` (float): photoionization cross sections

        """
        
        assert isinstance( Z, int )
        assert isinstance( N, int )

        Eth = self.Eth( Z, N )
        sigma = self.sigma( Z, N, Eth )
        return sigma



    def sigma_H1( self, E ):
        """ Convenience function that calls :func:`sigma` with `Z` and `N` set
        for neutral hydrogen. 

        Args:

          `E` (array): calculate cross section at these energies

        Returns:

          `sigma` (array): photoionization cross sections

        """

        if self.fit_type == 'verner':
            sigma = self.v96.return_fit( Z=1, N=1, E=E )
        return sigma


    def sigma_He1( self, E ):
        """ Convenience function that calls :func:`sigma` with `Z` and `N` set
        for neutral helium. 

        Args:

          `E` (array): calculate cross section at these energies

        Returns:

          `sigma` (array): photoionization cross sections

        """

        if self.fit_type == 'verner':
            sigma = self.v96.return_fit( Z=2, N=2, E=E )
        return sigma


    def sigma_He2( self, E ):
        """ Convenience function that calls :func:`sigma` with `Z` and `N` set
        for singly ionized helium. 

        Args:

          `E` (array): calculate cross section at these energies

        Returns:

          `sigma` (array): photoionization cross sections


        """
 
        if self.fit_type == 'verner':
            sigma = self.v96.return_fit( Z=2, N=1, E=E )
        return sigma


