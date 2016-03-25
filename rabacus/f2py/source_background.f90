!
! Routines for isotropic background sources of radiation.  For these sources, 
! the componenets of a spectrum are interpreted in the following way ... 
!
!    E_eV: energies sampled [eV]
!    shape: Inu [erg/(s Hz cm^2 sr)] (for polychromatic point sources)
!    shape: Inu [erg/(s cm^2 sr)] (for monochromatic point sources)
!
! This module supplies functions and subroutines for calculating 
! various integrals over the spectrum. 
!
!
! 
!
! Flux: 
! ----------------------
!   Fu_thin:  energy flux [erg/s/cm^2] 
!   Fn_thin:  photon flux [1/s/cm^2] 
!
! Density: 
! ----------------------
!   u_thin:  energy density [erg/cm^3] 
!   n_thin:  photon density [1/cm^3] 
!
! Optically Thin Photoionization/heating: 
! ----------------------
!   H1i/He1i/He2i_thin: H1/He1/He2 opt. thin photoionization rate [1/s]
!   H1h/He1h/He2h_thin: H1/He1/He2 opt. thin photoheating rate [erg/s]
!
! Shielded Photoionization/heating: 
! ----------------------
!   H1i/He1i/He2i_shld: H1/He1/He2 shielded photoionization rate [1/s]
!   H1h/He1h/He2h_shld: H1/He1/He2 shielded photoheating rate [erg/s]
!
! Normalize: 
! ----------------------
!   Fn: normalize the spectrum shape such that it returns desired Fn_thin
!    n: normalize the spectrum shape such that it returns desired n_thin
!------------------------------------------------------------------------

module source_background
  use types
  use physical_constants, only: h_eV_s, h_erg_s, c_cm_s
  use physical_constants, only: eV_2_erg, erg_2_eV
  use physical_constants, only: pi
  use utils, only: utils_integrate_trap
!  use special_functions, only: E2
  use zhang_jin_f, only: e2xa
  implicit none


  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: four_pi = 4.0d0 * pi


contains



!====================================================================
!====================================================================
!
! Flux integrals
!
!====================================================================
!====================================================================


  !
  ! Returns energy flux Fu [erg/s/cm^2].
  !
  !  Monochromatic:
  !    Fu = 4 pi Inu [erg/(s cm^2)] 
  !
  !  Polychromatic:
  !    Fu = int (4 pi Inu) / (h) * dE * eV_2_erg [erg/s/cm^2]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_Fu_thin( E_eV, Inu, nn ) result( Fu )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    integer(int32), intent(in) :: nn
    real(real64) :: Fu
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Fu = four_pi * Inu(0)
    else
       yy = four_pi * Inu / h_erg_s * eV_2_erg
       Fu = utils_integrate_trap( E_eV, yy, nn ) 
    end if

  end function bgnd_src_return_Fu_thin

  !
  ! Returns photon flux Fn [1/s].
  !
  !  Monochromatic:
  !    Inu is actually Fu [erg/s] 
  !    Fn = Fu / E * erg_2_eV [1/(s cm^2)]
  !
  !  Polychromatic:
  !    Fn = int (dFu/dnu) / (E h) * dE [1/(s cm^2)]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_Fn_thin( E_eV, Inu, nn ) result( Fn )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    integer(int32), intent(in) :: nn
    real(real64) :: Fn
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then 
       Fn = four_pi * Inu(0) / E_eV(0) * erg_2_eV
    else
       yy = four_pi * Inu / (E_eV * h_erg_s)
       Fn = utils_integrate_trap( E_eV, yy, nn ) 
    end if

  end function bgnd_src_return_Fn_thin


!====================================================================
!====================================================================
!
! Volume density integrals
!
!====================================================================
!====================================================================


  !
  ! Returns energy density u [erg/cm^3] 
  !
  !    u = Fu / c [erg/cm^3]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_u_thin( E_eV, Inu, nn ) result( u )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    integer(int32), intent(in) :: nn
    real(real64) :: u
    real(real64) :: Fu
    Fu = bgnd_src_return_Fu_thin( E_eV, Inu, nn )
    u = Fu / c_cm_s

  end function bgnd_src_return_u_thin

  !
  ! Returns photon density n [1/cm^3] 
  !
  !    n = Fn / c [1/cm^3]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_n_thin( E_eV, Inu, nn ) result( n )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    integer(int32), intent(in) :: nn
    real(real64) :: n
    real(real64) :: Fn
    Fn = bgnd_src_return_Fn_thin( E_eV, Inu, nn )
    n = Fn / c_cm_s

  end function bgnd_src_return_n_thin


!====================================================================
!====================================================================
!
! Optically thin photoionization/heating rates 
!
!====================================================================
!====================================================================

  !
  ! Returns H1/He1/He2 photoionization rate (i) [1/s] 
  ! The rate returned is determined by which sigma is passed in. 
  !
  !  Monochromatic:
  !    Inu is actually Fu [erg/(s cm^2)] 
  !    i = Fn * sigma [1/s]
  !
  !  Polychromatic:
  !     dFn/dE = (dFu/dnu) / (E h) [1/(eV s cm^2]
  !     i = int dFn/dE * sigma * dE [1/s]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_ionrate_thin( E_eV, Inu, sigma, nn ) &
       result( ionrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    integer(int32), intent(in) :: nn
    real(real64) :: ionrate
    real(real64) :: Fn
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Fn = bgnd_src_return_Fn_thin( E_eV, Inu, nn )
       ionrate = Fn * sigma(0)
    else
       yy = four_pi * Inu / (E_eV * h_erg_s)
       yy = yy * sigma
       ionrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_ionrate_thin


  !
  ! Returns H1/He1/He2 photoheating rate (h) [erg/s] 
  ! The rate returned is determined by which sigma and E_th is passed in. 
  !
  !  Monochromatic:
  !    Inu is actually Fu [erg/(s cm^2)] 
  !    h = Fn * sigma * (E - Eth) * eV_2_erg [erg/s]
  !
  !  Polychromatic:
  !      dFn/dE = (dFu/dnu) / (E h) [1/(eV s cm^2]
  !       di/dE = dFn/dE * sigma [1/eV s]
  !           h = int di/dE * (E-Eth) * eV_2_erg dE [erg/s]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_heatrate_thin( E_eV, Inu, sigma, E_th, &
       nn ) result( heatrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: E_th
    integer(int32), intent(in) :: nn
    real(real64) :: heatrate
    real(real64) :: Fn
    real(real64), dimension(0:nn-1) :: yy
    if ( nn == 1 ) then
       Fn = bgnd_src_return_Fn_thin( E_eV, Inu, nn )
       heatrate = Fn * sigma(0) * ( E_eV(0) - E_th ) * eV_2_erg 
    else
       yy = four_pi * Inu / (E_eV * h_erg_s)
       yy = yy * sigma
       yy = yy * (E_eV - E_th) * eV_2_erg
       heatrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_heatrate_thin




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
  function bgnd_src_set_attenuation( tauH1, tauHe1, tauHe2, nn ) result( atten )
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: atten
    real(real64), dimension(0:nn-1) :: tau  

    tau = tauH1 + tauHe1 + tauHe2
    atten = exp( -tau )

  end function bgnd_src_set_attenuation



  !
  ! Returns H1/He1/He2 photoionization rate [1/s] 
  ! The rate calculated is determined by the sigma passed in. 
  !
  !  Monochromatic:
  !    Inu is actually Fu [erg/s] 
  !    ionrate = Fn * sigma * exp(-tau) [1/s]
  !
  !  Polychromatic:
  !    dFn/dE = (dFu/dnu) / (E h) [1/(eV s cm^2]
  !    ionrate = int dFn/dE * sigma * exp(-tau) dE [1/s]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_ionrate_shld( E_eV, Inu, sigma, &
       tauH1, tauHe1, tauHe2, nn ) result( ionrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64) :: ionrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten

    atten = bgnd_src_set_attenuation( tauH1, tauHe1, tauHe2, nn )

    if ( nn == 1 ) then
       ionrate = bgnd_src_return_ionrate_thin( E_eV, Inu, sigma, nn )
       ionrate = ionrate * atten(0)
    else
       yy = four_pi * Inu / (E_eV * h_erg_s) 
       yy = yy * sigma * atten
       ionrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_ionrate_shld


  !
  ! Returns H1/He1/He2 photoioheating rate [erg/s] 
  ! The rate calculated is determined by the sigma and E_th passed in. 
  !
  !  Monochromatic:
  !    Inu is actually Fu [erg/s] 
  !    heatrate = Fn * sigma * exp(-tau) * (E-Eth) * eV_2_erg [erg/s]
  !
  !  Polychromatic:
  !    dFn/dE = (dLu/dnu) / (E h) [1/(eV s cm^2]
  !    heatrate = int dFn/dE * sigma * exp(-tau) * (E-Eth) dE [1/s]
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_heatrate_shld( E_eV, Inu, sigma, &
       E_th, tauH1, tauHe1, tauHe2, nn ) result( heatrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: E_th
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64) :: heatrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten

    atten = bgnd_src_set_attenuation( tauH1, tauHe1, tauHe2, nn )

    if ( nn == 1 ) then
       heatrate = bgnd_src_return_heatrate_thin( E_eV, Inu, sigma, &
                                               E_th, nn )
       heatrate = heatrate * atten(0)
    else
       yy = four_pi * Inu / (E_eV * h_erg_s) 
       yy = yy * sigma * atten
       yy = yy * (E_eV - E_th) * eV_2_erg
       heatrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_heatrate_shld



  !
  ! Returns H1/He1/He2 photoionization rate [1/s] 
  ! The rate calculated is determined by the sigma passed in. 
  ! This is for planar gas distributions with uniform backgrounds
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_ionrate_shld_ei( E_eV, Inu, sigma, &
       tauH1, tauHe1, tauHe2, nn ) result( ionrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64) :: ionrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten
    real(real64), dimension(0:nn-1) :: tau  
    integer(int32) :: i

    tau = tauH1 + tauHe1 + tauHe2
    if ( sum(tau) == zero ) then
       ionrate = bgnd_src_return_ionrate_thin( E_eV, Inu, sigma, nn )
       return
    end if

    do i = 0, nn-1
!       atten(i) = tau(i) * Ei(-tau(i)) + exp( -tau(i) )
       atten(i) = e2xa(tau(i)) 
    end do

    if ( nn == 1 ) then
       ionrate = bgnd_src_return_ionrate_thin( E_eV, Inu, sigma, nn )
       ionrate = ionrate * atten(0)
    else
       yy = four_pi * Inu / (E_eV * h_erg_s) 
       yy = yy * sigma * atten
       ionrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_ionrate_shld_ei


  !
  ! Returns H1/He1/He2 photoioheating rate [erg/s] 
  ! The rate calculated is determined by the sigma and E_th passed in. 
  ! This is for planar gas ditribtions with uniform backgrounds
  !
  !------------------------------------------------------------------------
  function bgnd_src_return_heatrate_shld_ei( E_eV, Inu, sigma, &
       E_th, tauH1, tauHe1, tauHe2, nn ) result( heatrate )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(in) :: Inu
    real(real64), dimension(0:nn-1), intent(in) :: sigma
    real(real64), intent(in) :: E_th
    real(real64), dimension(0:nn-1), intent(in) :: tauH1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe1
    real(real64), dimension(0:nn-1), intent(in) :: tauHe2
    integer(int32), intent(in) :: nn
    real(real64) :: heatrate
    real(real64), dimension(0:nn-1) :: yy
    real(real64), dimension(0:nn-1) :: atten
    real(real64), dimension(0:nn-1) :: tau  
    integer(int32) :: i

    tau = tauH1 + tauHe1 + tauHe2
    if ( sum(tau) == zero ) then
       heatrate = bgnd_src_return_heatrate_thin( E_eV, Inu, sigma, E_th, nn )
       return
    end if

    do i = 0, nn-1
!       atten(i) = tau(i) * Ei(-tau(i)) + exp( -tau(i) )
       atten(i) = e2xa(tau(i)) 
    end do

    if ( nn == 1 ) then
       heatrate = bgnd_src_return_heatrate_thin( E_eV, Inu, sigma, &
                                               E_th, nn )
       heatrate = heatrate * atten(0)
    else
       yy = four_pi * Inu / (E_eV * h_erg_s) 
       yy = yy * sigma * atten
       yy = yy * (E_eV - E_th) * eV_2_erg
       heatrate = utils_integrate_trap( E_eV, yy, nn )
    end if

  end function bgnd_src_return_heatrate_shld_ei



!====================================================================
!====================================================================
!
! Normalize spectrum
!
!====================================================================
!====================================================================


  !
  ! Normalize shape such that it will produce a photon flux of Fn_goal. 
  !------------------------------------------------------------------------
  subroutine bgnd_src_normalize_shape_Fn( E_eV, shape, Fn_goal, nn ) 
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(inout) :: shape
    real(real64), intent(in) :: Fn_goal
    integer(int32), intent(in) :: nn
    real(real64) :: Fn

    Fn = bgnd_src_return_Fn_thin( E_eV, shape, nn )
    shape = ( Fn_goal / Fn ) * shape 

  end subroutine bgnd_src_normalize_shape_Fn


  !
  ! Normalize shape such that it will produce a photon density of n_goal. 
  !------------------------------------------------------------------------
  subroutine bgnd_src_normalize_shape_n( E_eV, shape, n_goal, nn ) 
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), dimension(0:nn-1), intent(inout) :: shape
    real(real64), intent(in) :: n_goal
    integer(int32), intent(in) :: nn
    real(real64) :: n

    n = bgnd_src_return_n_thin( E_eV, shape, nn )
    shape = ( n_goal / n ) * shape 

  end subroutine bgnd_src_normalize_shape_n
  



end module source_background
