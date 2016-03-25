import time
import numpy as np
import pylab as plt
import rabacus as ra




def plot_sphere_x( s, fname, style, plot_H=True, plot_He=True ):
  """ put plot of ionization fractions from sphere `s` into fname """

  plt.figure( figsize=(7,7) )

  s.r_c.units = 'kpc'
  s.Edges.units = 'kpc'
  if style == 'I':
      xx = s.r_c / s.Edges[-1]
  else:
      xx = s.r_c

  if plot_He:
      plt.plot( xx, np.log10( s.xHe1 ),
                color='green', ls='-', label = r'$x_{\rm HeI}$' )
      plt.plot( xx, np.log10( s.xHe2 ),
                color='green', ls='--', label = r'$x_{\rm HeII}$' )
      plt.plot( xx, np.log10( s.xHe3 ),
                color='green', ls=':', label = r'$x_{\rm HeIII}$' )

  if plot_H:
      plt.plot( xx, np.log10( s.xH1 ),
                color='red', ls='-', label = r'$x_{\rm HI}$' )
      plt.plot( xx, np.log10( s.xH2 ),
                color='red', ls='--', label = r'$x_{\rm HII}$' )

  if style == 'I':
      plt.xlim( -0.05, 1.05 )
      plt.xlabel( 'r_c / Rsphere' )
  else:
      plt.xlim( -0.50, 15.2 )
      plt.xlabel( 'r_c [kpc]' )

  if style == 'I':
      plt.ylim( -6.0, 0.5 )
  else:
      plt.ylim( -3.1, 0.1 )
  plt.ylabel( 'log 10 ( x )' )

  plt.grid()
  plt.legend(loc='lower center', ncol=2)
  plt.tight_layout()
  plt.savefig( 'doc/img/x_' + fname )


def plot_sphere_T( s, fname, style ):
  """ put plot of temperature from sphere `s` into fname """

  plt.figure( figsize=(7,7) )

  s.r_c.units = 'kpc'
  s.Edges.units = 'kpc'
  if style == 'I':
      xx = s.r_c / s.Edges[-1]
  else:
      xx = s.r_c

  plt.plot( xx, s.T,
            color='black', ls='-', label = r'$T$' )

  if style == 'I':
      plt.xlim( -0.05, 1.05 )
      plt.xlabel( 'r_c / Rsphere' )
  else:
      plt.xlim( -0.50, 15.2 )
      plt.xlabel( 'r_c [kpc]' )

  plt.ylim( 0.0, 5.0e4 )
  plt.ylabel( 'T [K]' )

  plt.grid()
  plt.legend(loc='best')
  plt.tight_layout()
  plt.savefig( 'doc/img/T_' + fname )




Ln = 5.0e48 / ra.u.s  # set photon luminosity

q_mono = 1.0
q_min = q_mono
q_max = q_mono
src_mono = ra.PointSource( q_min, q_max, 'monochromatic' )
src_mono.normalize_Ln( Ln )

q_min = 1.0
q_max = 10.0
T_eff = 1.0e5 * ra.u.K
src_thrm = ra.PointSource( q_min, q_max, 'thermal', T_eff=T_eff )
src_thrm.normalize_Ln( Ln )

q_min = 1.0
q_max = 10.0
alpha = -1.0
src_pwr1 = ra.PointSource( q_min, q_max, 'powerlaw', alpha=alpha )
src_pwr1.normalize_Ln( Ln )




spheres = {}




Nl = 512
T = np.ones(Nl) * 1.0e4 * ra.u.K

Rsphere = 6.6 * ra.u.kpc
Edges = np.linspace( 0.0 * ra.u.kpc, Rsphere, Nl+1 )
nH = np.ones(Nl) * 1.0e-3 / ra.u.cm**3
nHe = np.ones(Nl) * 8.7e-5 / ra.u.cm**3
nHe_null = np.ones(Nl) * 1.0e-15 / ra.u.cm**3




key = 'thin_caseB_mono_fixT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=0.0, thin=True )

plot_sphere_x(
    spheres[key], 'strm_sphere_' + key + '.png', 'I', plot_He=False )




key = 'rt_caseB_mono_fixT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=0.0 )

plot_sphere_x(
    spheres[key], 'strm_sphere_' + key + '.png', 'I', plot_He=False )




key = 'rt_caseB_thrm_evoT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_thrm, fixed_fcA=0.0, find_Teq=True, z=0.0 )

plot_sphere_x(
    spheres[key], 'strm_sphere_' + key + '.png', 'I', plot_He=False )
plot_sphere_T(
    spheres[key], 'strm_sphere_' + key + '.png', 'I' )




Rsphere = 15.0 * ra.u.kpc
Edges = np.linspace( 0.0 * ra.u.kpc, Rsphere, Nl+1 )




key = 'rt_caseA_thrm_fixT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe, src_thrm, fixed_fcA=1.0 )

plot_sphere_x( spheres[key], 'strm_sphere_' + key + '.png', 'F' )



key = 'rt_caseB_thrm_fixT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe, src_thrm, fixed_fcA=0.0 )

plot_sphere_x( spheres[key], 'strm_sphere_' + key + '.png', 'F' )



key = 'rt_caseB_pwr1_evoT'

spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe, src_pwr1, fixed_fcA=0.0, find_Teq=True, z=0.0 )

plot_sphere_x( spheres[key], 'strm_sphere_' + key + '.png', 'F' )
plot_sphere_T( spheres[key], 'strm_sphere_' + key + '.png', 'F' )





Rsphere = 6.6 * ra.u.kpc
Edges = np.linspace( 0.0 * ra.u.kpc, Rsphere, Nl+1 )


key = 'rt_caseA_mono_fixT'
spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, fixed_fcA=1.0 )


key = 'rt_outward_mono_fixT'
spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, rec_meth='outward' )


key = 'rt_radial_mono_fixT'
spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, rec_meth='radial' )


key = 'rt_isotropic_mono_fixT'
spheres[key] = ra.SphereStromgren(
    Edges, T, nH, nHe_null, src_mono, rec_meth='isotropic' )



plt.figure( figsize=(7,7) )

key = 'rt_caseA_mono_fixT'
s = spheres[key]
s.r_c.units = 'kpc'
s.Edges.units = 'kpc'
xx = s.r_c / s.Edges[-1]

key = 'rt_caseA_mono_fixT'
s = spheres[key]
plt.plot( xx, np.log10( s.xH1 ),
          color='blue', ls=':', label = r'case A' )

key = 'rt_caseB_mono_fixT'
s = spheres[key]
plt.plot( xx, np.log10( s.xH1 ),
          color='black', ls='--', label = r'case B' )

key = 'rt_outward_mono_fixT'
s = spheres[key]
plt.plot( xx, np.log10( s.xH1 ),
          color='cyan', ls='-', label = r'outward' )

key = 'rt_radial_mono_fixT'
s = spheres[key]
plt.plot( xx, np.log10( s.xH1 ),
          color='green', ls='-', label = r'radial' )

key = 'rt_isotropic_mono_fixT'
s = spheres[key]
plt.plot( xx, np.log10( s.xH1 ),
          color='red', ls='-', label = r'isotropic' )

plt.xlim( -0.05, 1.05 )
plt.xlabel( 'r_c / Rsphere' )

plt.ylim( -6.0, 0.5 )
plt.ylabel( 'log 10 ( x )' )

plt.grid()
plt.legend(loc='best', ncol=2)
plt.tight_layout()
plt.savefig( 'doc/img/x_raicevic.png' )








plt.figure( figsize=(7,7) )

key = 'rt_isotropic_mono_fixT'
s = spheres[key]
s.r_c.units = 'kpc'
s.Edges.units = 'kpc'
xx = s.r_c / (5.4*ra.u.kpc) # s.Edges[-1]

geo = 4.0 * np.pi * s.r_c**2
I0 = s.H1i_src[0] * geo[0]

# source I/I0
#------------------------------------
Is = s.H1i_src * geo
yy = Is / I0
plt.plot( xx, yy,
          color='red', ls='-', label='source' )

# recomb I/I0
#------------------------------------
Id = s.H1i_rec * geo
yy = Id / I0
plt.plot( xx, yy,
          color='red', ls='--', label='diffuse' )

# ratio diffuse/source
#------------------------------------
yy = Id / Is
plt.plot( xx, yy,
          color='red', ls=':', label='ratio' )

plt.ylabel( r'$I/I_0$', fontsize=25 )
plt.ylim( 0.0, 1.05 )

plt.xlabel( 'r_c / R_strom' )
plt.xlim( -0.05, 1.05 )

plt.grid()
plt.legend(loc='best', ncol=1)
plt.tight_layout()
plt.savefig( 'doc/img/I_raicevic_isotropic.png' )





plt.figure( figsize=(7,7) )

key = 'rt_radial_mono_fixT'
s = spheres[key]
s.r_c.units = 'kpc'
s.Edges.units = 'kpc'
xx = s.r_c / (5.4*ra.u.kpc) # s.Edges[-1]

geo = 4.0 * np.pi * s.r_c**2
I0 = s.H1i_src[0] * geo[0]

# source I/I0
#------------------------------------
Is = s.H1i_src * geo
yy = Is / I0
plt.plot( xx, yy,
          color='green', ls='-', label='source' )

# recomb I/I0
#------------------------------------
Id = s.H1i_rec * geo
yy = Id / I0
plt.plot( xx, yy,
          color='green', ls='--', label='diffuse' )

# ratio diffuse/source
#------------------------------------
yy = Id / Is
plt.plot( xx, yy,
          color='green', ls=':', label='ratio' )

plt.ylabel( r'$I/I_0$', fontsize=25 )
plt.ylim( 0.0, 1.05 )

plt.xlabel( 'r_c / R_strom' )
plt.xlim( -0.05, 1.05 )

plt.grid()
plt.legend(loc='best', ncol=1)
plt.tight_layout()
plt.savefig( 'doc/img/I_raicevic_radial.png' )
