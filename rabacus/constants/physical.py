""" 
Physical constants from NIST  
(http://physics.nist.gov/cuu/Constants/index.html). 
We also make use of constants defined by the IAU-2009 Report 
(pdf included in the doc/ref directory). 
""" 

import quantities as pq


__all__ = ['PhysicalConstants']


class PhysicalConstants:
    """ Physical constants.  

    Attributes:
      `m_p`: proton mass. 

      `m_n`: neutraon mass. 

      `m_e`: electron mass. 

      `c`: speed of light in vacuum. 

      `h`: Planck constant. 

      `kb`: Boltzmann constant. 

      `Rbohr`: Bohr radius.

      `sigma_t`: Thomson cross-section.

      `alpha`: fine-structure constant.

      `G`: Newtonian constant of gravitation. 

      `Ry_inf`: Rydberg constant times hc.

      `sfkb`: Stefan-Boltzmann constant.

      `GMsun`: solar mass parameter.

      `Msun`: solar mass. 

      `AU`: IAU defined astronomical unit. 
      

    """

    def __init__(self):



        #
        # we convert the CODATA values from mks to cgs units 
        #
        
        # values from CODATA
        #-------------------------------------------------------------
        self.m_p = 1.672621777e-27 * pq.kg
        self.m_p.units = pq.g

        self.m_n = 1.674927351e-27 * pq.kg
        self.m_n.units = pq.g

        self.m_e = 9.10938291e-31 * pq.kg
        self.m_e.units = pq.g

        self.c = 2.99792458e8 * pq.m / pq.s
        self.c.units = pq.cm / pq.s

        self.h = 6.62606957e-34 * pq.J * pq.s
        self.h.units = pq.erg * pq.s

        self.kb = 1.3806488e-23 * pq.J / pq.K
        self.kb.units = pq.erg / pq.K

        self.Rbohr = 0.52917721092e-10 * pq.m
        self.Rbohr.units = pq.cm

        self.sigma_t = 0.6652458734e-28 * pq.m**2
        self.sigma_t.units =   pq.cm**2

        self.alpha = 7.2973525698e-3 * pq.dimensionless

        self.G = 6.67384e-11 * pq.m**3 * pq.kg**(-1) * pq.s**(-2)
        self.G.units = pq.cm**3 * pq.g**(-1) * pq.s**(-2)
        
        self.Ry_inf = 13.60569253e0 * pq.eV
        
        self.sfkb = 5.670373e-8 * pq.W * pq.m**(-2) * pq.K**(-4)
        self.sfkb.units = pq.erg / pq.s * pq.cm**(-2) * pq.K**(-4)


        # values from IAU-2009
        #-------------------------------------------------------------

        # Solar mass parameter = G * M_sun
        self.GMsun = 1.327124421e20 * pq.m**3 * pq.s**(-2)
        self.GMsun.units = pq.cm**3 * pq.s**(-2)

        # Solar mass = GM_sun / G
        self.Msun = self.GMsun / self.G
        self.Msun.units = pq.g

        # IAU definition of au
        self.AU = 1.49597870700e11 * pq.m
        self.AU.units = pq.cm
 




