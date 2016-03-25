!
! Routines for point sources of radiation.  For these types of sources, 
! the componenets of a spectrum are interpreted in the following way ... 
!
!    E_eV: energies sampled [eV]
!    shape: dLu/dnu [erg/(s Hz)] (for polychromatic point sources)
!    shape: Lu [erg/s] (for monochromatic point sources)
!
! This module supplies functions and subroutines for calculating 
! various integrals over the spectrum. 
!
!
! 
! Luminosity: 
! ----------------------
!   Lu_thin:  energy luminosity [erg/s]
!   Ln_thin:  photon luminosity [1/s]
!
! Flux: 
! ----------------------
!   Fu_thin:  energy flux at a given radius [erg/s/cm^2] 
!   Fn_thin:  photon flux at a given radius [1/s/cm^2] 
!
! Density: 
! ----------------------
!   u_thin:  energy density at a given radius [erg/cm^3] 
!   n_thin:  photon density at a given radius [1/cm^3] 
!
! Optically Thin Photoionization/heating: 
! ----------------------
!   H1i/He1i/He2i_thin: H1/He1/He2 opt. thin photoionization rate at R [1/s]
!   H1h/He1h/He2h_thin: H1/He1/He2 opt. thin photoheating rate at R [erg/s]
!
! Shielded Photoionization/heating: 
! ----------------------
!   H1i/He1i/He2i_shld: H1/He1/He2 shielded photoionization rate at R [1/s]
!   H1h/He1h/He2h_shld: H1/He1/He2 shielded photoheating rate at R [erg/s]
!
! Normalize: 
! ----------------------
!   Ln: normalize the spectrum shape such that it returns desired Ln_thin
!
!------------------------------------------------------------------------

module source_point
  use types
  use physical_constants, only: h_eV_s, h_erg_s, c_cm_s
  use physical_constants, only: eV_2_erg, erg_2_eV
  use physical_constants, only: pi
  use utils, only: utils_integrate_trap
  implicit none


contains


!====================================================================
!====================================================================
!
! Luminosity integrals
!
!====================================================================
!====================================================================

  !
  ! Returns energy luminosity Lu [erg/s].
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !
  !  Polychromatic:
  !    Lu = int (dLu/dnu) / (h) * dE * eV_2_erg [erg/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_Lu_thin( E_eV, dLu_over_dnu, nn ) result( Lu )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    integer(int32), intent(in) :: nn
    real(real64) :: Lu
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Lu = dLu_over_dnu(0)
    else
       yy = dLu_over_dnu / h_erg_s * eV_2_erg
       Lu = utils_integrate_trap( E_eV, yy, nn ) 
    end if
  end function pt_src_return_Lu_thin


  !
  ! Returns photon luminosity Ln [1/s].
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !    Ln = Lu / E * erg_2_eV [1/s] 
  !
  !  Polychromatic:
  !    Ln = int (dLu/dnu) / (E h) dE [1/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_Ln_thin( E_eV, dLu_over_dnu, nn ) result( Ln )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    integer(int32), intent(in) :: nn
    real(real64) :: Ln
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then 
       Ln = dLu_over_dnu(0) / E_eV(0) * erg_2_eV
    else
       yy = dLu_over_dnu / (E_eV * h_erg_s)
       Ln = utils_integrate_trap( E_eV, yy, nn ) 
    end if
  end function pt_src_return_Ln_thin


!====================================================================
!====================================================================
!
! Flux integrals
!
!====================================================================
!====================================================================


  !
  ! Returns energy flux Fu [erg/s/cm^2] at given radius [cm].
  !
  !   Fu = Lu / (4 pi r^2) [erg/(s cm^2)]
  !
  !------------------------------------------------------------------------
  function pt_src_return_Fu_thin( E_eV, dLu_over_dnu, radius, nn ) result( Fu )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: Fu
    real(real64) :: Lu
    Lu = pt_src_return_Lu_thin( E_eV, dLu_over_dnu, nn )
    Fu = Lu / (4.0d0 * pi * radius*radius)
  end function pt_src_return_Fu_thin


  !
  ! Returns photon flux Fn [1/s/cm^2] at given radius [cm].
  !
  !   Fn = Ln / (4 pi r^2) [1/(s cm^2)]
  !
  !------------------------------------------------------------------------
  function pt_src_return_Fn_thin( E_eV, dLu_over_dnu, radius, nn ) result( Fn )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: Fn
    real(real64) :: Ln
    Ln = pt_src_return_Ln_thin( E_eV, dLu_over_dnu, nn )
    Fn = Ln / (4.0d0 * pi * radius*radius)
  end function pt_src_return_Fn_thin


!====================================================================
!====================================================================
!
! Volume density integrals
!
!====================================================================
!====================================================================


  !
  ! Returns energy density u [erg/cm^3] at given radius [cm].
  !
  !    u = Fu / c [erg/cm^3]
  !
  !------------------------------------------------------------------------
  function pt_src_return_u_thin( E_eV, dLu_over_dnu, radius, nn ) result( u )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: u
    real(real64) :: Fu
    Fu = pt_src_return_Fu_thin( E_eV, dLu_over_dnu, radius, nn )
    u = Fu / c_cm_s
  end function pt_src_return_u_thin

  !
  ! Returns photon density n [1/cm^3] at given radius [cm].
  !
  !    n = Fn / c [1/cm^3]
  !
  !------------------------------------------------------------------------
  function pt_src_return_n_thin( E_eV, dLu_over_dnu, radius, nn ) result( n )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: n
    real(real64) :: Fn
    Fn = pt_src_return_Fn_thin( E_eV, dLu_over_dnu, radius, nn )
    n = Fn / c_cm_s
  end function pt_src_return_n_thin


!====================================================================
!====================================================================
!
! Optically thin photoionization/heating rates 
!
!====================================================================
!====================================================================

  !
  ! Returns H1/He1/He2 photoionization rate (i) [1/s] at given radius [cm] 
  ! The rate returned is determined by which sigma is passed in. 
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !    i = Fn * sigma [1/s]
  !
  !  Polychromatic:
  !     dFn/dE = (dLu/dnu) / (E h 4 pi r^2) [1/(eV s cm^2]
  !     i = int dFn/dE * sigma * dE [1/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_ionrate_thin( E_eV, dLu_over_dnu, sigma, &
       radius, nn ) result( ionrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: ionrate
    real(real64) :: Fn
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Fn = pt_src_return_Fn_thin( E_eV, dLu_over_dnu, radius, nn )
       ionrate = Fn * sigma(0)
    else
       yy = dLu_over_dnu / (E_eV * h_erg_s * 4.0d0 * pi * radius*radius) 
       yy = yy * sigma
       ionrate = utils_integrate_trap( E_eV, yy, nn )
    end if
  end function pt_src_return_ionrate_thin


  !
  ! Returns H1/He1/He2 photoheating rate (h) [erg/s] at given radius [cm] 
  ! The rate returned is determined by which sigma and E_th is passed in. 
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !    h = Fn * sigma * (E - Eth) * eV_2_erg [erg/s]
  !
  !  Polychromatic:
  !      dFn/dE = (dLu/dnu) / (E h 4 pi r^2) [1/(eV s cm^2]
  !       di/dE = dFn/dE * sigma [1/eV s]
  !           h = int di/dE * (E-Eth) * eV_2_erg dE [erg/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_heatrate_thin( E_eV, dLu_over_dnu, sigma, E_th, &
       radius, nn ) result( heatrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: E_th
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: heatrate
    real(real64) :: Fn
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Fn = pt_src_return_Fn_thin( E_eV, dLu_over_dnu, radius, nn )
       heatrate = Fn * sigma(0) * ( E_eV(0) - E_th ) * eV_2_erg 
    else
       yy = dLu_over_dnu / (E_eV * h_erg_s * 4.0d0 * pi * radius*radius)
       yy = yy * sigma
       yy = yy * (E_eV - E_th) * eV_2_erg
       heatrate = utils_integrate_trap( E_eV, yy, nn )
    end if
  end function pt_src_return_heatrate_thin




!====================================================================
!====================================================================
!
! Shielded photoionization/heating rates 
!
!====================================================================
!====================================================================


  !
  ! Calculates attenuation array = exp(-tau)
  !
  !------------------------------------------------------------------------
  function pt_src_set_attenuation( tauH1, tauHe1, tauHe2, nn ) result( atten )
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: atten
    real(real64), dimension(0:nn-1) :: tau  

    tau = tauH1 + tauHe1 + tauHe2
    atten = exp( -tau )

  end function pt_src_set_attenuation



  !
  ! Returns H1/He1/He2 photoionization rate [1/s] at given radius [cm] 
  ! The rate calculated is determined by the sigma passed in. 
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !    ionrate = Fn * sigma * exp(-tau) [1/s]
  !
  !  Polychromatic:
  !    dFn/dE = (dLu/dnu) / (E h 4 pi r^2) [1/(eV s cm^2]
  !    ionrate = int dFn/dE * sigma * exp(-tau) dE [1/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_ionrate_shld( E_eV, dLu_over_dnu, sigma, &
       tauH1, tauHe1, tauHe2, radius, nn ) result( ionrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: ionrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten

    atten = pt_src_set_attenuation( tauH1, tauHe1, tauHe2, nn )

    if ( nn == 1 ) then
       ionrate = pt_src_return_ionrate_thin( E_eV, dLu_over_dnu, sigma, &
                                             radius, nn )
       ionrate = ionrate * atten(0)
    else
       yy = dLu_over_dnu / (E_eV * h_erg_s * 4.0d0 * pi * radius*radius) 
       yy = yy * sigma * atten
       ionrate = utils_integrate_trap( E_eV, yy, nn )
    end if
  end function pt_src_return_ionrate_shld


  !
  ! Returns H1/He1/He2 photoioheating rate [erg/s] at given radius [cm] 
  ! The rate calculated is determined by the sigma and E_th passed in. 
  !
  !  Monochromatic:
  !    dLu_over_dnu is actually Lu [erg/s] 
  !    heatrate = Fn * sigma * exp(-tau) * (E-Eth) [erg/s]
  !
  !  Polychromatic:
  !    dFn/dE = (dLu/dnu) / (E h 4 pi r^2) [1/(eV s cm^2]
  !    heatrate = int dFn/dE * sigma * exp(-tau) * (E-Eth) dE [1/s]
  !
  !------------------------------------------------------------------------
  function pt_src_return_heatrate_shld( E_eV, dLu_over_dnu, sigma, &
       E_th, tauH1, tauHe1, tauHe2, radius, nn ) result( heatrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: dLu_over_dnu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: E_th
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    real(real64), intent(in) :: radius
    integer(int32), intent(in) :: nn
    real(real64) :: heatrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten

    atten = pt_src_set_attenuation( tauH1, tauHe1, tauHe2, nn )

    if ( nn == 1 ) then
       heatrate = pt_src_return_heatrate_thin( E_eV, dLu_over_dnu, sigma, &
                                               E_th, radius, nn )
       heatrate = heatrate * atten(0)
    else
       yy = dLu_over_dnu / (E_eV * h_erg_s * 4.0d0 * pi * radius*radius) 
       yy = yy * sigma * atten
       yy = yy * (E_eV - E_th) * eV_2_erg
       heatrate = utils_integrate_trap( E_eV, yy, nn )
    end if
  end function pt_src_return_heatrate_shld


!====================================================================
!====================================================================
!
! Normalize spectrum
!
!====================================================================
!====================================================================


  !
  ! Normalize shape such that it will produce a photon luminosity of Ln_goal. 
  !------------------------------------------------------------------------
  subroutine pt_src_normalize_shape_Ln( E_eV, shape, Ln_goal, nn ) 
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(inout) :: shape
    real(real64), intent(in) :: Ln_goal
    integer(int32), intent(in) :: nn
    real(real64) :: Ln

    Ln = pt_src_return_Ln_thin( E_eV, shape, nn )
    shape = ( Ln_goal / Ln ) * shape 

  end subroutine pt_src_normalize_shape_Ln

  






end module source_point
