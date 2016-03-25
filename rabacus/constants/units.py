""" 
A wrapper to the quantities package
(https://github.com/python-quantities/python-quantities) 
providing units.  
""" 

import quantities as pq
import physical

__all__ = ['Units', 'CosmoUnits']

pc = physical.PhysicalConstants()

class Units:
    """ Standard units. 

    Attributes:

      `dimensionless`: no units

      `angstrom`: angstrom

      `nm`: nanometer

      `um`: micrometer

      `cm`: centimeter

      `m`: meter

      `km`: kilometer

      `pc`: parsec

      `AU`: atronomical unit

      `kpc`: kiloparsec

      `Mpc`: megaparsec

      `s`: second

      `yr`: year

      `Myr`: megayear

      `Gyr`: gigayear

      `g`: gram

      `kg`: kilogram

      `Msun`: solar mass

      `erg`: erg

      `eV`: electron volt

      `J`: joule

      `Ry_inf`: rydberg energy unit (infinite nucleus mass compared to electron)

      `Ry_H`: rydberg energy unit (reduced mass for H electron)

      `K`: kelvin

      `Hz`: Hertz

      `sr`: steradian
      
    """

    

    def __init__(self):

        # dimensionless
        #------------------------------------------------
        self.dimensionless = pq.dimensionless

        # length
        #------------------------------------------------
        self.angstrom = pq.angstrom

        self.nm = pq.nm
        self.um = pq.um
        self.cm = pq.cm
        self.m = pq.m
        self.km = pq.km

        self.pc = pq.pc 

        self.AU  = pq.UnitQuantity( 'astronomical unit',
                                     pc.AU,
                                     symbol='AU' )

        self.kpc = pq.UnitQuantity( 'kiloparsec',
                                     1.0e3 * self.pc,
                                     symbol='kpc' )
        
        self.Mpc = pq.UnitQuantity( 'megaparsec',
                                     1.0e6 * self.pc,
                                     symbol='Mpc' )

        # time
        #------------------------------------------------
        self.s = pq.s
        self.yr = pq.yr

        self.Myr = pq.UnitQuantity( 'megayear',
                                     1.0e6 * self.yr,
                                     symbol='Myr' )

        self.Gyr = pq.UnitQuantity( 'gigayear',
                                     1.0e9 * self.yr,
                                     symbol='Gyr' )

        # mass
        #------------------------------------------------
        self.g = pq.g
        self.kg = pq.kg

        self.Msun = pq.UnitQuantity( 'solar mass',
                                      pc.GMsun / pc.G,
                                      symbol='Msun' )


        # energy
        #------------------------------------------------
        self.erg = pq.erg
        self.eV = pq.eV
        self.J = pq.J

        self.Ry_inf = pq.UnitQuantity( 'Rydberg_infinity',
                                        pc.Ry_inf,
                                        symbol='Ry_inf' )

        m_H = pc.m_p + pc.m_e
        self.Ry_H = pq.UnitQuantity( 'Rydberg_H',
                                     pc.Ry_inf * (1.0 - pc.m_e/m_H),
                                     symbol='Ry_H' )


        # temperature
        #------------------------------------------------
        self.K = pq.K


        # frequency
        #------------------------------------------------
        self.Hz = pq.Hz


        # solid angle 
        #------------------------------------------------
        self.sr = pq.sr





class CosmoUnits(Units):    
    """ Provides cosmological units. 

    Args: 

      `h` (float): hubble parameter, H_0 = 100 `h` km / s / Mpc

      `a` (float): scale factor for co-moving conversions

    """


    def __init__(self, h, a):

        # define new units
        #----------------------------------------------------
        self.h = pq.UnitQuantity( 'dimensionless hubble parameter',
                                  h,
                                  symbol='hh' )

        self.a = pq.UnitQuantity( 'scale factor a = 1/(1+z)',
                                  a,
                                  symbol='aa' )

        Units.__init__(self)

#        self.H0 = 100.0 * self.h * self.km / self.s / self.Mpc








