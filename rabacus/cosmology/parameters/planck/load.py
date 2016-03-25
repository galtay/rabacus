""" I/O routines to read in Plank data """ 

import os as _os
import string as _string
import numpy as _np

__all__ = ['PlanckParameters']

class _Error(Exception):
    """Base class for exceptions in this module."""
    pass


class _InvalidDataSets(_Error):
    """Exception raised for requesting an invalid group of data sets. """ 
    def __init__(self, message ):
        Exception.__init__( self, message )


class PlanckParameters:
    r""" Load Planck cosmological parameters.  

    Args:
      `data_sets` (str): A string indicating which data sets to include. 

    Possible values for `data_sets` include, 

      - ``base_planck_lowl``
      - ``base_planck_lowl_lowLike``
      - ``base_planck_lowl_post_lensing``
      - ``base_planck_lowl_lowLike_highL``
      - ``base_planck_lowl_lowLike_post_lensing``
      - ``base_planck_lowl_lowLike_highL_post_lensing``    
  

    Attributes:
      values (dict): cosmological parameter values

      latex (dict): latex representation of parameter name
      
       
    Notes:
 
      Description of data sets. 

      ``base``: 6 parameter LCDM model

      ``planck``: high-l planck temperature data 

      ``lowl``: low-l planck temperature data

      ``lowLike``: low-l WMAP9 polarization data (WP) 

      ``lensing``: planck lensing power spectrum

      ``highL``: high-l data from ACT and SPT
        

      The default data sets ``base_plank_lowl`` will return the values in 
      column 1 of table 2 in http://arxiv.org/abs/1303.5076

""" 


    def __init__( self, data_sets='base_planck_lowl' ):

        self.ValidDataSets = ['base_planck_lowl',
                              'base_planck_lowl_lowLike',
                              'base_planck_lowl_post_lensing',
                              'base_planck_lowl_lowLike_highL',
                              'base_planck_lowl_lowLike_post_lensing',
                              'base_planck_lowl_lowLike_highL_post_lensing']


        if not data_sets in self.ValidDataSets:
            msg = '\n'
            msg += 'received data_sets = ' + data_sets + '\n'
            msg += '\ndata_sets must be one of ... \n'
            for ds in self.ValidDataSets:
                msg += '  ' + ds + '\n'
            raise _InvalidDataSets(msg)

        self.data_sets = data_sets
        fname_in = data_sets + '.minimum'
        local = _os.path.dirname(_os.path.realpath(__file__))
        fname = local + '/' + fname_in
    
        self.values = {}
        self.latex = {}
        f = open( fname, 'r' )
        for line in f:
            pieces = line.split()
            if len(pieces) >= 4:
                if pieces[0].isdigit():
                    self.values[pieces[2]] = _np.float64(pieces[1])
                    self.latex[pieces[2]]  = _string.join(pieces[3:])
        f.close()

        self.values['h'] = self.values['H0'] / 100.0
        self.latex['h'] = 'h'

        self.values['omegab'] = self.values['omegabh2'] / self.values['h']**2
        self.latex['omegab'] = '\Omega_b'

        self.values['omegac'] = self.values['omegach2'] / self.values['h']**2
        self.latex['omegac'] = '\Omega_c'

        self.values['omegar'] = self.values['omegam'] / (1.0+self.values['zeq'])
        self.latex['omegar'] = '\Omega_r'

        # copy yheused into yhe and make alias Yp
        self.values['yhe'] = self.values['yheused']
        self.latex['yhe'] = self.latex['yheused']

        self.values['Yp'] = self.values['yhe']
        self.latex['Yp'] = r'Y_p'


    def __str__(self):
        print 'PlanckParameters Instance'
        print '  data sets: ', self.data_sets
        print '  OmegaM:  %5.3f' % self.values['omegam']
        print '  OmegaB:  %5.3f' % self.values['omegab']
        print '  OmegaC:  %5.3f' % self.values['omegac']
        print '  OmegaL:  %5.3f' % self.values['omegal']
        print '  OmegaR:  %5.3f' % self.values['omegar']
        print '  h:       %5.3f' % self.values['h']
        print '  sigma_8: %5.3f' % self.values['sigma8']
        print '  n_s:     %5.3f' % self.values['ns']
        print '  Y_P:     %5.3f' % self.values['yhe']
        print 
        return ''



if __name__ == "__main__":

    """ Get the values presented in Tables 2 and 5 of 
    http://arxiv.org/abs/1303.5076 """

    data_sets = 'base_planck_lowl'
    pp_base = PlanckParameters( data_sets=data_sets ) 

    data_sets = 'base_planck_lowl_post_lensing'
    pp_lens = PlanckParameters( data_sets=data_sets )

    data_sets = 'base_planck_lowl_lowLike'
    pp_WP = PlanckParameters( data_sets=data_sets )

    data_sets = 'base_planck_lowl_lowLike_highL'
    pp_WP_highL = PlanckParameters( data_sets=data_sets )

    data_sets = 'base_planck_lowl_lowLike_highL_post_lensing'
    pp_WP_highL_lens = PlanckParameters( data_sets=data_sets )

    print
    print 'Showing values from Table 2 http://arxiv.org/abs/1303.5076 \n'

    print '           [OmegaL, OmegaM, OmegaB, OmegaC, H0,   sigma8, ns]:'
    fs = '%.5f'

    pb = pp_base.values
    print 'Base      :', fs % pb['omegal'], fs % pb['omegam'], \
        fs % pb['omegab'], fs % pb['omegac'], '%.2f' % pb['H0'], \
        fs % pb['sigma8'], fs % pb['ns']

    p = pp_lens.values
    print 'Base+lens :', fs % p['omegal'], fs % p['omegam'], \
        fs % p['omegab'], fs % p['omegac'], '%.2f' % p['H0'], \
        fs % p['sigma8'], fs % p['ns']

    p = pp_WP.values
    print 'Base+WP   :', fs % p['omegal'], fs % p['omegam'], \
        fs % p['omegab'], fs % p['omegac'], '%.2f' % p['H0'], \
        fs % p['sigma8'], fs % p['ns']

    print
    print 'Showing percent difference from base \n'

    p = pp_lens.values
    print 'Base+lens :', \
        fs % ( 100.0*(p['omegal'] - pb['omegal']) / pb['omegal']), \
        fs % ( 100.0*(p['omegam'] - pb['omegam']) / pb['omegam']), \
        fs % ( 100.0*(p['omegab'] - pb['omegab']) / pb['omegab']), \
        fs % ( 100.0*(p['omegac'] - pb['omegac']) / pb['omegac']), \
        fs % ( 100.0*(p['H0'] - pb['H0']) / pb['H0']), \
        fs % ( 100.0*(p['sigma8'] - pb['sigma8']) / pb['sigma8']), \
        fs % ( 100.0*(p['ns'] - pb['ns']) / pb['ns'])

    p = pp_WP.values
    print 'Base+WP   :', \
        fs % ( 100.0*(p['omegal'] - pb['omegal']) / pb['omegal']), \
        fs % ( 100.0*(p['omegam'] - pb['omegam']) / pb['omegam']), \
        fs % ( 100.0*(p['omegab'] - pb['omegab']) / pb['omegab']), \
        fs % ( 100.0*(p['omegac'] - pb['omegac']) / pb['omegac']), \
        fs % ( 100.0*(p['H0'] - pb['H0']) / pb['H0']), \
        fs % ( 100.0*(p['sigma8'] - pb['sigma8']) / pb['sigma8']), \
        fs % ( 100.0*(p['ns'] - pb['ns']) / pb['ns'])



    print
    print
    print 'Showing values from Table 5 http://arxiv.org/abs/1303.5076 \n'

    pb = pp_base.values
    print 'Base            :', fs % pb['omegal'], fs % pb['omegam'], \
        fs % pb['omegab'], fs % pb['omegac'], '%.2f' % pb['H0'], \
        fs % pb['sigma8'], fs % pb['ns']

    p = pp_WP_highL.values
    print 'Base+WP+hiL     :', fs % p['omegal'], fs % p['omegam'], \
        fs % p['omegab'], fs % p['omegac'], '%.2f' % p['H0'], \
        fs % p['sigma8'], fs % p['ns']

    p = pp_WP_highL_lens.values
    print 'Base+WP+hiL+lens:', fs % p['omegal'], fs % p['omegam'], \
        fs % p['omegab'], fs % p['omegac'], '%.2f' % p['H0'], \
        fs % p['sigma8'], fs % p['ns']

    print
    print 'Showing percent difference from base \n'

    p = pp_WP_highL.values
    print 'Base+WP+hiL     :', \
        fs % ( 100.0*(p['omegal'] - pb['omegal']) / pb['omegal']), \
        fs % ( 100.0*(p['omegam'] - pb['omegam']) / pb['omegam']), \
        fs % ( 100.0*(p['omegab'] - pb['omegab']) / pb['omegab']), \
        fs % ( 100.0*(p['omegac'] - pb['omegac']) / pb['omegac']), \
        fs % ( 100.0*(p['H0'] - pb['H0']) / pb['H0']), \
        fs % ( 100.0*(p['sigma8'] - pb['sigma8']) / pb['sigma8']), \
        fs % ( 100.0*(p['ns'] - pb['ns']) / pb['ns'])

    p = pp_WP_highL_lens.values
    print 'Base+WP+hiL+lens:', \
        fs % ( 100.0*(p['omegal'] - pb['omegal']) / pb['omegal']), \
        fs % ( 100.0*(p['omegam'] - pb['omegam']) / pb['omegam']), \
        fs % ( 100.0*(p['omegab'] - pb['omegab']) / pb['omegab']), \
        fs % ( 100.0*(p['omegac'] - pb['omegac']) / pb['omegac']), \
        fs % ( 100.0*(p['H0'] - pb['H0']) / pb['H0']), \
        fs % ( 100.0*(p['sigma8'] - pb['sigma8']) / pb['sigma8']), \
        fs % ( 100.0*(p['ns'] - pb['ns']) / pb['ns'])


    print






