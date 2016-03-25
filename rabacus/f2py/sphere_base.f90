! 
! base routines for spherical geometries
!
!---------------------------------------------------------

module sphere_base

  use types
  use physical_constants, only: pi
  use physical_constants, only: eV_2_erg, erg_2_eV
  use geometry, only: m_for_ray, ds_segment
  use ion_solver, only: solve_pce, solve_pcte
  use chem_cool_rates, only: get_kchem

  use photo_xsections, only: return_E_H1_th, return_sigma_H1
  use photo_xsections, only: return_E_He1_th, return_sigma_He1
  use photo_xsections, only: return_E_He2_th, return_sigma_He2

  use legendre_polynomial, only: p_quadrature_rule

  implicit none


  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: quarter = 0.25d0
  real(real64), parameter :: half = 0.5d0
  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: two = 2.0d0
  real(real64), parameter :: three = 3.0d0
  real(real64), parameter :: four = 4.0d0





contains


  !======================================================================
  ! set_recomb_photons_outward
  !
  ! Calculates the photoionization and photoheating rates in each layer
  ! due to recombination photons emitted from every other layer. 
  ! The outward only approximation is used in this routine. 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  ! ---------------------------------------------------------------
  !   Edges: radius of the Nl+1 shells [cm] 
  !   nH: hydrogen number density in each layer [cm^-3]
  !   nHe: helium number density in each layer [cm^-3]
  !   Tmprtr: temperature in each layer [K]
  ! ---------------------------------------------------------------
  !   sigma_H1_th: H1 photoion xsection at H1 threshold
  !   sigma_He1_th: He1 photoion xsection at He1 threshold
  !   sigma_He2_th: He2 photoion xsection at He2 threshold
  ! ---------------------------------------------------------------
  !   E_H1_th: H1 ionization threshold [eV]
  !   E_He1_th: He1 ionization threshold [eV]
  !   E_He2_th: He2 ionization threshold [eV]
  ! ---------------------------------------------------------------
  !   i_photo_fit: photo_xsection fits {1=verner96}
  !   i_rate_fit: atomic rate fits {1=hg97}
  !   Nl: number of layers
  !
  ! Output: 
  !   H1i: H1 photoion rate due to recombination photons
  !   He1i: He1 photoion rate due to recombination photons
  !   He2i: He2 photoion rate due to recombination photons
  !   H1h: H1 photoheat rate due to recombination photons
  !   He1h: He1 photoheat rate due to recombination photons
  !   He2h: He2 photoheat rate due to recombination photons
  !
  !======================================================================
  subroutine set_recomb_photons_outward( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_photo_fit, i_rate_fit, &
       H1i, He1i, He2i, H1h, He1h, He2h, Nl )

    ! Arguments
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1), intent(inout) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(inout) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl), intent(inout) :: Edges
    real(real64), dimension(0:Nl-1), intent(inout) ::  nH, nHe
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr
    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th
    integer(int32), intent(in) :: i_photo_fit
    integer(int32), intent(in) :: i_rate_fit
    real(real64), dimension(0:Nl-1), intent(out) :: H1i, He1i, He2i
    real(real64), dimension(0:Nl-1), intent(out) :: H1h, He1h, He2h
    integer(int32), intent(in) :: Nl

    ! Local
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1) :: dr, r_c
    real(real64), dimension(0:Nl-1) :: dNH, dNHe
    real(real64), dimension(0:Nl-1) :: dtau_H1_th
    real(real64), dimension(0:Nl-1) :: dtau_He1_th
    real(real64), dimension(0:Nl-1) :: dtau_He2_th

    integer(int32), parameter :: Nnu = 1
    real(real64), dimension(0:Nnu-1) :: E_eV
    real(real64) :: sigma_H1_at_H1_th
    real(real64) :: sigma_H1_at_He1_th
    real(real64) :: sigma_H1_at_He2_th
    real(real64) :: sigma_He1_at_He1_th
    real(real64) :: sigma_He1_at_He2_th
    real(real64) :: sigma_He2_at_He2_th

    real(real64), dimension(0:Nl-1) :: fcA_H2, fcA_He2, fcA_He3
    real(real64), dimension(0:Nl-1) :: reH2_A, reHe2_A, reHe3_A
    real(real64), dimension(0:Nl-1) :: reH2_B, reHe2_B, reHe3_B
    real(real64), dimension(0:Nl-1) :: reH2_1, reHe2_1, reHe3_1
    real(real64), dimension(0:Nl-1) :: ciH1, ciHe1, ciHe2

    real(real64), dimension(0:Nl-1) :: ne, area, vol
    real(real64), dimension(0:Nl-1) :: em_H2, em_He2, em_He3
    real(real64), dimension(0:Nl-1) :: F_H2, F_He2, F_He3
    real(real64), dimension(0:Nl-1) :: F0_H2, F0_He2, F0_He3

    integer(int32) :: isl, irl
    real(real64) :: geo_fac, tau
    real(real64) :: atten_H1, atten_He1, atten_He2

    real(real64) :: tau_H1_th, sigma_H1, sigma_H1_ra
    real(real64) :: tau_He1_th, sigma_He1, sigma_He1_ra
    real(real64) :: tau_He2_th, sigma_He2, sigma_He2_ra

    logical :: reversed
    real(real64), dimension(0:Nl-1) :: dum

    ! geometry.  
    ! here we make sure we have monotonically decreasing Edges
    ! and indexes starting at 0
    !---------------------------------------------

    reversed = .false.
    if ( monotonic_increasing(  Edges, Nl+1 ) ) then

       reversed = .true. 
       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )

    end if

    call get_dr_and_r_c(  Edges, dr, r_c, Nl )

    dNH = dr * nH
    dNHe = dr * nHe

    dtau_H1_th = dNH * xH1 * sigma_H1_th 
    dtau_He1_th = dNHe * xHe1 * sigma_He1_th
    dtau_He2_th = dNHe * xHe2 * sigma_He2_th

    area = four * pi * r_c*r_c
    vol = abs( four / three * pi * ( Edges(0:Nl-1)**3 - Edges(1:Nl)**3 ) )

    ! calculate cross-sections
    !------------------------------------------------------------
    E_eV = E_He2_th
    sigma_H1_at_He2_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He2_th = sum( return_sigma_He1( E_eV, i_photo_fit, 1 ) )
    sigma_He2_at_He2_th = sigma_He2_th

    E_eV = E_He1_th
    sigma_H1_at_He1_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He1_th = sigma_He1_th

    sigma_H1_at_H1_th = sigma_H1_th

    ! get case A rates
    !------------------------------------------------
    fcA_H2 = one
    fcA_He2 = one
    fcA_He3 = one
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_A, reHe2_A, reHe3_A, &
         ciH1, ciHe1, ciHe2, Nl )

    ! get case B rates
    !------------------------------------------------
    fcA_H2 = zero
    fcA_He2 = zero
    fcA_He3 = zero
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_B, reHe2_B, reHe3_B, &
         ciH1, ciHe1, ciHe2, Nl )
    
    ! calculate recomb. to ground state rates 
    !------------------------------------------------
    reH2_1 = reH2_A - reH2_B
    reHe2_1 = reHe2_A - reHe2_B
    reHe3_1 = reHe3_A - reHe3_B


    ! calculate emission coefficients
    ! photons / ( s cm^3 sr )
    !------------------------------------------------------------
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    em_H2 = ( reH2_1 * nH * xH2 * ne ) / (four * pi)
    em_He2 = ( reHe2_1 * nHe * xHe2 * ne ) / (four * pi)
    em_He3 = ( reHe3_1 * nHe * xHe3 * ne ) / (four * pi)


    ! calculate flux from each recombining layer 
    !------------------------------------------------------------
    F0_H2 = em_H2 * four * pi * vol / area
    F0_He2 = em_He2 * four * pi * vol / area
    F0_He3 = em_He3 * four * pi * vol / area


    ! calculate flux in each solve layer
    !============================================================

    F_H2 = zero    
    F_He2 = zero
    F_He3 = zero

    over_solve_layers: do isl = 0, Nl-1

       ! for each solve layer we build up the contribution from 
       ! recombination radiation in layers at smaller radii. 

       tau_H1_th = zero
       tau_He1_th = zero
       tau_He2_th = zero

       over_recomb_layers: do irl = isl, Nl-1

          ! increment optical depths 
          !-----------------------------------------------------------
          if ( irl == isl ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl == isl+1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl > isl+1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl-1) * half
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl-1) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl-1) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          end if

          ! geometric factor
          !-----------------------------------------------------------
          geo_fac = ( r_c(irl) / r_c(isl) ) * ( r_c(irl) / r_c(isl) )

          ! for H2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_H1_th 
          sigma_H1_ra = sigma_H1 / sigma_H1_th
          tau = tau_H1_th * sigma_H1_ra
          F_H2(isl) = F_H2(isl) +  F0_H2(irl) * geo_fac * exp(-tau) 
          
          ! for He2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He1_th
          sigma_He1 = sigma_He1_at_He1_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th           
          sigma_He1_ra = sigma_He1 / sigma_He1_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra

          F_He2(isl) = F_He2(isl) + F0_He2(irl) * geo_fac * exp(-tau)

          ! for He3 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He2_th
          sigma_He1 = sigma_He1_at_He2_th
          sigma_He2 = sigma_He2_at_He2_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th 
          sigma_He1_ra = sigma_He1 / sigma_He1_th
          sigma_He2_ra = sigma_He2 / sigma_He2_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra + &
                tau_He2_th * sigma_He2_ra

          F_He3(isl) = F_He3(isl) + F0_He3(irl) * geo_Fac * exp(-tau)
          
       end do over_recomb_layers
       
    end do over_solve_layers


    ! photo-ionization rates from fluxes
    !----------------------------------------------
    H1i = F_H2 * sigma_H1_at_H1_th + &
          F_He2 * sigma_H1_at_He1_th + &
          F_He3 * sigma_H1_at_He2_th

    He1i = F_He2 * sigma_He1_at_He1_th + &
           F_He3 * sigma_He1_at_He2_th

    He2i = F_He3 * sigma_He2_at_He2_th


    ! photo-heating rates from fluxes
    !----------------------------------------------
    H1h = F_He2 * ( E_He1_th - E_H1_th ) * sigma_H1_at_He1_th + &
          F_He3 * ( E_He2_th - E_H1_th ) * sigma_H1_at_He2_th
    H1h = H1h * eV_2_erg

    He1h = F_He2 * sigma_He1_at_He1_th + &
           F_He3 * sigma_He1_at_He2_th
    He1h = He1h * eV_2_erg

    He2h = zero
    He2h = He2h * eV_2_erg



    if ( reversed ) then 
       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )
    end if

  end subroutine set_recomb_photons_outward



  !======================================================================
  ! set_recomb_photons_radial
  !
  ! Calculates the photoionization and photoheating rates in each layer
  ! due to recombination photons emitted from every other layer. 
  ! The radial only approximation is used in this routine.  It is 
  ! convenient to have Edges be monotonically increasing and begin
  ! indexing layers with 1 in this routine. 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  ! ---------------------------------------------------------------
  !   Edges: radius of all Nl+1 shells  [cm] 
  !   nH: hydrogen number density in each layer [cm^-3]
  !   nHe: helium number density in each layer [cm^-3]
  !   Tmprtr: temperature in each layer [K]
  ! ---------------------------------------------------------------
  !   sigma_H1_th: H1 photoion xsection at H1 threshold
  !   sigma_He1_th: He1 photoion xsection at He1 threshold
  !   sigma_He2_th: He2 photoion xsection at He2 threshold
  ! ---------------------------------------------------------------
  !   E_H1_th: H1 ionization threshold [eV]
  !   E_He1_th: He1 ionization threshold [eV]
  !   E_He2_th: He2 ionization threshold [eV]
  ! ---------------------------------------------------------------
  !   i_photo_fit: photo_xsection fits {1=verner96}
  !   i_rate_fit: atomic rate fits {1=hg97}
  !   Nl: number of layers
  !
  ! Output: 
  !   H1i: H1 photoion rate due to recombination photons
  !   He1i: He1 photoion rate due to recombination photons
  !   He2i: He2 photoion rate due to recombination photons
  !   H1h: H1 photoheat rate due to recombination photons
  !   He1h: He1 photoheat rate due to recombination photons
  !   He2h: He2 photoheat rate due to recombination photons
  !
  !======================================================================
  subroutine set_recomb_photons_radial( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_photo_fit, i_rate_fit, &
       H1i, He1i, He2i, H1h, He1h, He2h, Nl )

    ! Arguments
    !---------------------------------------------------------------
    real(real64), dimension(1:Nl), intent(inout) :: xH1, xH2
    real(real64), dimension(1:Nl), intent(inout) :: xHe1, xHe2, xHe3
    real(real64), dimension(1:Nl+1), intent(inout) :: Edges
    real(real64), dimension(1:Nl), intent(inout) ::  nH, nHe
    real(real64), dimension(1:Nl), intent(inout) :: Tmprtr
    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th
    integer(int32), intent(in) :: i_photo_fit
    integer(int32), intent(in) :: i_rate_fit
    real(real64), dimension(1:Nl), intent(out) :: H1i, He1i, He2i
    real(real64), dimension(1:Nl), intent(out) :: H1h, He1h, He2h
    integer(int32), intent(in) :: Nl



    ! Local
    !---------------------------------------------------------------
    real(real64), dimension(1:Nl) :: dr, r_c
    real(real64), dimension(1:Nl) :: dNH, dNHe
    real(real64), dimension(1:Nl) :: dtau_H1_th
    real(real64), dimension(1:Nl) :: dtau_He1_th
    real(real64), dimension(1:Nl) :: dtau_He2_th

    integer(int32), parameter :: Nnu = 1
    real(real64), dimension(0:Nnu-1) :: E_eV
    real(real64) :: sigma_H1_at_H1_th
    real(real64) :: sigma_H1_at_He1_th
    real(real64) :: sigma_H1_at_He2_th
    real(real64) :: sigma_He1_at_He1_th
    real(real64) :: sigma_He1_at_He2_th
    real(real64) :: sigma_He2_at_He2_th

    real(real64), dimension(1:Nl) :: fcA_H2, fcA_He2, fcA_He3
    real(real64), dimension(1:Nl) :: reH2_A, reHe2_A, reHe3_A
    real(real64), dimension(1:Nl) :: reH2_B, reHe2_B, reHe3_B
    real(real64), dimension(1:Nl) :: reH2_1, reHe2_1, reHe3_1
    real(real64), dimension(1:Nl) :: ciH1, ciHe1, ciHe2

    real(real64), dimension(1:Nl) :: ne, area, vol
    real(real64), dimension(1:Nl) :: em_H2, em_He2, em_He3
    real(real64), dimension(1:Nl) :: F_H2, F_He2, F_He3
    real(real64), dimension(1:Nl) :: F0_H2, F0_He2, F0_He3

    integer(int32) :: isl, irl, il
    real(real64) :: geo_fac, tau
    real(real64) :: atten_H1, atten_He1, atten_He2

    real(real64) :: tau_H1_th_lo, tau_H1_th_hi, tau_H1_th
    real(real64) :: tau_He1_th_lo, tau_He1_th_hi, tau_He1_th
    real(real64) :: tau_He2_th_lo, tau_He2_th_hi, tau_He2_th

    real(real64) :: sigma_H1, sigma_H1_ra
    real(real64) :: sigma_He1, sigma_He1_ra
    real(real64) :: sigma_He2, sigma_He2_ra

    real(real64), dimension(1:Nl) :: dum
    logical :: reversed


    ! geometry (need monotonic increasing for this one)
    !---------------------------------------------
    reversed = .false.
    if ( .not. monotonic_increasing( Edges, Nl+1 ) ) then

       reversed = .true.

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )
       
    end if

    if ( .not. monotonic_increasing( Edges, Nl+1 ) ) then
       write(*,*) 'cant make Edges monotonic increasing' 
       stop
    end if

    call get_dr_and_r_c(  Edges, dr, r_c, Nl )

    dNH = dr * nH
    dNHe = dr * nHe

    dtau_H1_th = dNH * xH1 * sigma_H1_th 
    dtau_He1_th = dNHe * xHe1 * sigma_He1_th
    dtau_He2_th = dNHe * xHe2 * sigma_He2_th

    area = four * pi * r_c*r_c
    vol = abs( four / three * pi * ( Edges(1:Nl)**3 - Edges(2:Nl+1)**3 ) )

    ! calculate cross-sections
    !------------------------------------------------------------
    E_eV = E_He2_th
    sigma_H1_at_He2_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He2_th = sum( return_sigma_He1( E_eV, i_photo_fit, 1 ) )
    sigma_He2_at_He2_th = sigma_He2_th

    E_eV = E_He1_th
    sigma_H1_at_He1_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He1_th = sigma_He1_th

    sigma_H1_at_H1_th = sigma_H1_th

    ! get case A rates
    !------------------------------------------------
    fcA_H2 = one
    fcA_He2 = one
    fcA_He3 = one
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_A, reHe2_A, reHe3_A, &
         ciH1, ciHe1, ciHe2, Nl )

    ! get case B rates
    !------------------------------------------------
    fcA_H2 = zero
    fcA_He2 = zero
    fcA_He3 = zero
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_B, reHe2_B, reHe3_B, &
         ciH1, ciHe1, ciHe2, Nl )
    
    ! calculate recomb. to ground state rates 
    !------------------------------------------------
    reH2_1 = reH2_A - reH2_B
    reHe2_1 = reHe2_A - reHe2_B
    reHe3_1 = reHe3_A - reHe3_B


    ! calculate emission coefficients
    ! photons / ( s cm^3 sr )
    !------------------------------------------------------------
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    em_H2 = ( reH2_1 * nH * xH2 * ne ) / (four * pi)
    em_He2 = ( reHe2_1 * nHe * xHe2 * ne ) / (four * pi)
    em_He3 = ( reHe3_1 * nHe * xHe3 * ne ) / (four * pi)


    ! calculate flux from each recombining layer 
    !------------------------------------------------------------
    F0_H2 = em_H2 * four * pi * vol / area
    F0_He2 = em_He2 * four * pi * vol / area
    F0_He3 = em_He3 * four * pi * vol / area


    ! calculate flux in each solve layer 
    !============================================================

    F_H2 = zero    
    F_He2 = zero
    F_He3 = zero

    over_solve_layers: do isl = 1, Nl

       ! for each solve layer we build up the contribution from 
       ! recombination radiation in layers at smaller and larger
       ! radii.  cycle over il = 0 b/c it doesn't exist 

       tau_H1_th = zero
       tau_He1_th = zero
       tau_He2_th = zero

       over_lower_recomb_layers: do il = isl, -Nl, -1

          if ( il == 0 ) cycle

          irl = abs(il)

          ! increment optical depths 
          !-----------------------------------------------------------
          if ( irl == isl ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl == isl-1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl < isl-1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl+1) * half
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl+1) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl+1) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          end if

          ! geometric factor
          !-----------------------------------------------------------
          geo_fac = ( r_c(irl) / r_c(isl) ) * ( r_c(irl) / r_c(isl) )

          ! for H2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_H1_th 
          sigma_H1_ra = sigma_H1 / sigma_H1_th
          tau = tau_H1_th * sigma_H1_ra
          F_H2(isl) = F_H2(isl) + half * F0_H2(irl) * geo_fac * exp(-tau) 
          
          ! for He2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He1_th
          sigma_He1 = sigma_He1_at_He1_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th           
          sigma_He1_ra = sigma_He1 / sigma_He1_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra

          F_He2(isl) = F_He2(isl) + half * F0_He2(irl) * geo_fac * exp(-tau)

          ! for He3 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He2_th
          sigma_He1 = sigma_He1_at_He2_th
          sigma_He2 = sigma_He2_at_He2_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th 
          sigma_He1_ra = sigma_He1 / sigma_He1_th
          sigma_He2_ra = sigma_He2 / sigma_He2_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra + &
                tau_He2_th * sigma_He2_ra

          F_He3(isl) = F_He3(isl) + half * F0_He3(irl) * geo_Fac * exp(-tau)
          
       end do over_lower_recomb_layers



       tau_H1_th = zero
       tau_He1_th = zero
       tau_He2_th = zero

       over_upper_recomb_layers: do il = isl, Nl

          irl = abs(il)

          ! increment optical depths 
          !-----------------------------------------------------------
          if ( irl == isl ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl == isl+1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
          else if ( irl > isl+1 ) then
             tau_H1_th = tau_H1_th + dtau_H1_th(irl) * half
             tau_H1_th = tau_H1_th + dtau_H1_th(irl-1) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl) * half
             tau_He1_th = tau_He1_th + dtau_He1_th(irl-1) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl) * half
             tau_He2_th = tau_He2_th + dtau_He2_th(irl-1) * half
          end if

          ! geometric factor
          !-----------------------------------------------------------
          geo_fac = ( r_c(irl) / r_c(isl) ) * ( r_c(irl) / r_c(isl) )

          ! for H2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_H1_th 
          sigma_H1_ra = sigma_H1 / sigma_H1_th
          tau = tau_H1_th * sigma_H1_ra
          F_H2(isl) = F_H2(isl) + half * F0_H2(irl) * geo_fac * exp(-tau) 
          
          ! for He2 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He1_th
          sigma_He1 = sigma_He1_at_He1_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th           
          sigma_He1_ra = sigma_He1 / sigma_He1_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra

          F_He2(isl) = F_He2(isl) + half * F0_He2(irl) * geo_fac * exp(-tau)

          ! for He3 recombination photons
          !--------------------------------------------------------
          sigma_H1 = sigma_H1_at_He2_th
          sigma_He1 = sigma_He1_at_He2_th
          sigma_He2 = sigma_He2_at_He2_th

          sigma_H1_ra = sigma_H1 / sigma_H1_th 
          sigma_He1_ra = sigma_He1 / sigma_He1_th
          sigma_He2_ra = sigma_He2 / sigma_He2_th

          tau = tau_H1_th * sigma_H1_ra + &
                tau_He1_th * sigma_He1_ra + &
                tau_He2_th * sigma_He2_ra

          F_He3(isl) = F_He3(isl) + half * F0_He3(irl) * geo_Fac * exp(-tau)
          
       end do over_upper_recomb_layers
       


    end do over_solve_layers


    ! photo-ionization rates from fluxes
    !----------------------------------------------
    H1i = F_H2 * sigma_H1_at_H1_th + &
          F_He2 * sigma_H1_at_He1_th + &
          F_He3 * sigma_H1_at_He2_th

    He1i = F_He2 * sigma_He1_at_He1_th + &
           F_He3 * sigma_He1_at_He2_th

    He2i = F_He3 * sigma_He2_at_He2_th


    ! photo-heating rates from fluxes
    !----------------------------------------------
    H1h = F_He2 * ( E_He1_th - E_H1_th ) * sigma_H1_at_He1_th + &
          F_He3 * ( E_He2_th - E_H1_th ) * sigma_H1_at_He2_th
    H1h = H1h * eV_2_erg

    He1h = F_He2 * sigma_He1_at_He1_th + &
           F_He3 * sigma_He1_at_He2_th
    He1h = He1h * eV_2_erg

    He2h = zero
    He2h = He2h * eV_2_erg


    if ( reversed ) then

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )
       
    end if



  end subroutine set_recomb_photons_radial




  !======================================================================
  ! set_recomb_photons_isotropic
  !
  ! Calculates the photoionization and photoheating rates in each layer
  ! due to recombination photons emitted from every other layer. 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  ! ---------------------------------------------------------------
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH: hydrogen number density in each layer [cm^-3]
  !   nHe: helium number density in each layer [cm^-3]
  !   Tmprtr: temperature in each layer [K]
  !   Nmu: number of polar angle samples
  ! ---------------------------------------------------------------
  !   sigma_H1_th: H1 photoion xsection at H1 threshold
  !   sigma_He1_th: He1 photoion xsection at He1 threshold
  !   sigma_He2_th: He2 photoion xsection at He2 threshold
  ! ---------------------------------------------------------------
  !   E_H1_th: H1 ionization threshold [eV]
  !   E_He1_th: He1 ionization threshold [eV]
  !   E_He2_th: He2 ionization threshold [eV]
  ! ---------------------------------------------------------------
  !   i_photo_fit: photo_xsection fits {1=verner96}
  !   i_rate_fit: atomic rate fits {1=hg97}
  ! ---------------------------------------------------------------
  !   em_H1_fac: multiplicative factor for H1 recomb emission
  !   em_He1_fac: multiplicative factor for He1 recomb emission
  !   em_He2_fac: multiplicative factor for He2 recomb emission
  ! ---------------------------------------------------------------
  !   Nl: number of layers
  !
  ! Output: 
  !   H1i: H1 photoion rate due to recombination photons
  !   He1i: He1 photoion rate due to recombination photons
  !   He2i: He2 photoion rate due to recombination photons
  !   H1h: H1 photoheat rate due to recombination photons
  !   He1h: He1 photoheat rate due to recombination photons
  !   He2h: He2 photoheat rate due to recombination photons
  !
  !======================================================================
  subroutine set_recomb_photons_isotropic( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, Nmu, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_photo_fit, i_rate_fit, &
       em_H1_fac, em_He1_fac, em_He2_fac, &
       H1i, He1i, He2i, H1h, He1h, He2h, Nl )

    ! Arguments
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1), intent(inout) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(inout) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl), intent(inout) :: Edges
    real(real64), dimension(0:Nl-1), intent(inout) ::  nH, nHe
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr
    integer(int32), intent(in) :: Nmu
    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th
    integer(int32), intent(in) :: i_photo_fit
    integer(int32), intent(in) :: i_rate_fit
    real(real64), intent(in) :: em_H1_fac, em_He1_fac, em_He2_fac
    real(real64), dimension(0:Nl-1), intent(out) :: H1i, He1i, He2i
    real(real64), dimension(0:Nl-1), intent(out) :: H1h, He1h, He2h
    integer(int32), intent(in) :: Nl

    ! Local
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1) :: dk_at_EH1_from_H1
    real(real64), dimension(0:Nl-1) :: dk_at_EHe1_from_H1
    real(real64), dimension(0:Nl-1) :: dk_at_EHe1_from_He1
    real(real64), dimension(0:Nl-1) :: dk_at_EHe2_from_H1
    real(real64), dimension(0:Nl-1) :: dk_at_EHe2_from_He1
    real(real64), dimension(0:Nl-1) :: dk_at_EHe2_from_He2

    real(real64), dimension(0:Nl-1) :: dk_at_EH1, dk_at_EHe1, dk_at_EHe2

    integer(int32), parameter :: Nnu = 1
    real(real64), dimension(0:Nnu-1) :: E_eV
    real(real64) :: sigma_H1_at_H1_th
    real(real64) :: sigma_H1_at_He1_th
    real(real64) :: sigma_H1_at_He2_th
    real(real64) :: sigma_He1_at_He1_th
    real(real64) :: sigma_He1_at_He2_th
    real(real64) :: sigma_He2_at_He2_th

    real(real64), dimension(0:Nl-1) :: fcA_H2, fcA_He2, fcA_He3
    real(real64), dimension(0:Nl-1) :: reH2_A, reHe2_A, reHe3_A
    real(real64), dimension(0:Nl-1) :: reH2_B, reHe2_B, reHe3_B
    real(real64), dimension(0:Nl-1) :: reH2_1, reHe2_1, reHe3_1
    real(real64), dimension(0:Nl-1) :: ciH1, ciHe1, ciHe2
    real(real64), dimension(0:Nl-1) :: ne

    real(real64), dimension(0:Nl-1) :: em_EH1, em_EHe1, em_EHe2
    real(real64), dimension(0:Nl-1) :: J_EH1, J_EHe1, J_EHe2

    real(real64) :: mu
    real(real64), dimension(0:Nmu-1) :: xi, wi
    real(real64), dimension(0:Nmu-1) :: I_EH1, I_EHe1, I_EHe2

    real(real64) :: tau_EH1, tau_EHe1, tau_EHe2
    integer(int32) :: jnl, m, j, imu, il, j0
    real(real64) :: ds, dI, dmu

    logical :: reversed
    real(real64), dimension(0:Nl-1) :: dum

    ! geometry (need monotonic decreasing for this one)
    !---------------------------------------------
    reversed = .false.
    if ( monotonic_increasing( Edges, Nl+1 ) ) then

       reversed = .true.

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )
       
    end if

    if ( .not. monotonic_decreasing( Edges, Nl+1 ) ) then
       write(*,*) 'cant make Edges monotonic decreasing' 
       stop
    end if


    ! calculate cross-sections
    !------------------------------------------------------------
    E_eV = E_He2_th
    sigma_H1_at_He2_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He2_th = sum( return_sigma_He1( E_eV, i_photo_fit, 1 ) )
    sigma_He2_at_He2_th = sigma_He2_th

    E_eV = E_He1_th
    sigma_H1_at_He1_th = sum( return_sigma_H1( E_eV, i_photo_fit, 1 ) )
    sigma_He1_at_He1_th = sigma_He1_th

    sigma_H1_at_H1_th = sigma_H1_th


    ! partial opacity term at energy X from species Y
    !------------------------------------------------------------
    dk_at_EH1_from_H1 = nH * xH1 * sigma_H1_at_H1_th
    dk_at_EHe1_from_H1 = nH * xH1 * sigma_H1_at_He1_th
    dk_at_EHe2_from_H1 = nH * xH1 * sigma_H1_at_He2_th

    dk_at_EHe1_from_He1 = nHe * xHe1 * sigma_He1_at_He1_th
    dk_at_EHe2_from_He1 = nHe * xHe1 * sigma_He1_at_He2_th

    dk_at_EHe2_from_He2 = nHe * xHe2 * sigma_He2_at_He2_th


    dk_at_EH1 = dk_at_EH1_from_H1
    dk_at_EHe1 = dk_at_EHe1_from_H1 + dk_at_EHe1_from_He1
    dk_at_EHe2 = dk_at_EHe2_from_H1 + dk_at_EHe2_from_He1 + dk_at_EHe2_from_He2

    ! get case A rates
    !------------------------------------------------
    fcA_H2 = one
    fcA_He2 = one
    fcA_He3 = one
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_A, reHe2_A, reHe3_A, &
         ciH1, ciHe1, ciHe2, Nl )

    ! get case B rates
    !------------------------------------------------
    fcA_H2 = zero
    fcA_He2 = zero
    fcA_He3 = zero
    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2_B, reHe2_B, reHe3_B, &
         ciH1, ciHe1, ciHe2, Nl )
    
    ! calculate recomb. to ground state rates 
    !------------------------------------------------
    reH2_1 = (reH2_A - reH2_B) 
    reHe2_1 = (reHe2_A - reHe2_B)
    reHe3_1 = (reHe3_A - reHe3_B)


    ! calculate emission coefficients.  em_X is for emission at
    ! energy = ionization threshold of X
    ! photons / ( s cm^3 sr )
    !------------------------------------------------------------
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    em_EH1 = ( reH2_1 * nH * xH2 * ne ) / (four * pi)
    em_EHe1 = ( reHe2_1 * nHe * xHe2 * ne ) / (four * pi) 
    em_EHe2 = ( reHe3_1 * nHe * xHe3 * ne ) / (four * pi) 

    em_EH1 = em_EH1 * em_H1_fac 
    em_EHe1 = em_EHe1 * em_He1_fac 
    em_EHe2 = em_EHe2 * em_He2_fac 


    ! calculate xi and wi for Gauss-Legendre Quadrature
    !============================================================

    call p_quadrature_rule( Nmu, xi, wi )


    ! loop over layers and calculate specific intensity as a function
    ! of direction (mu)
    !============================================================

    J_EH1 = zero
    J_EHe1 = zero
    J_EHe2 = zero


    ! begin parallel region
    !============================================================

    !$omp parallel private( il, imu, mu, m, j0, j, jnl, ds, dI, &
    !$omp& I_EH1, I_EHe1, I_EHe2, tau_EH1, tau_EHe1, tau_EHe2 )

    !$omp do
    over_layers_1: do il = 0, Nl-1

       I_EH1 = zero
       I_EHe1 = zero
       I_EHe2 = zero
       
       over_mu_1: do imu = 0, Nmu-1

          mu = xi(imu)
          m = m_for_ray( Edges, il, mu, Nl )

          tau_EH1 = zero
          tau_EHe1 = zero
          tau_EHe2 = zero

          if ( mu <= zero ) then
             j0 = il
          else
             j0 = 2*m - il
          end if

          over_segements: do j = j0, 2*m
             jnl = m-abs(m-j)
             ds = ds_segment( Edges, il, mu, j, Nl )
             if ( j == j0 ) ds = ds * half
             
             ! tau at different energies
             !----------------------------------------
             tau_EH1 = tau_EH1 + ds * dk_at_EH1(jnl)
             tau_EHe1 = tau_EHe1 + ds * dk_at_EHe1(jnl)
             tau_EHe2 = tau_EHe2 + ds * dk_at_EHe2(jnl)
             
             dI = em_EH1(jnl) * ds * exp( -tau_EH1 )
             I_EH1(imu) = I_EH1(imu) + dI
             
             dI = em_EHe1(jnl) * ds * exp(-tau_EHe1 )
             I_EHe1(imu) = I_EHe1(imu) + dI
             
             dI = em_EHe2(jnl) * ds * exp(-tau_EHe2 )
             I_EHe2(imu) = I_EHe2(imu) + dI
             
          end do over_segements


       end do over_mu_1


       J_EH1(il) = sum( I_EH1 * wi )
       J_EHe1(il) = sum( I_EHe1 * wi )
       J_EHe2(il) = sum( I_EHe2 * wi )


    end do over_layers_1
    !$omp end do

    !$omp end parallel 

    J_EH1 = J_EH1 * half 
    J_EHe1 = J_EHe1 * half 
    J_EHe2 = J_EHe2 * half 

        
    ! photo-ionization rates from angle averaged specific intensity
    !----------------------------------------------
    H1i = ( J_EH1 * sigma_H1_at_H1_th + &
            J_EHe1 * sigma_H1_at_He1_th + &
            J_EHe2 * sigma_H1_at_He2_th ) * four * pi

    He1i = ( J_EHe1 * sigma_He1_at_He1_th + &
             J_EHe2 * sigma_He1_at_He2_th ) * four * pi

    He2i = ( J_EHe2 * sigma_He2_at_He2_th ) * four * pi


    ! photo-heating rates from angle averaged specific intensity
    !----------------------------------------------
    H1h = ( J_EHe1 * ( E_He1_th - E_H1_th ) * sigma_H1_at_He1_th + &
            J_EHe2 * ( E_He2_th - E_H1_th ) * sigma_H1_at_He2_th ) * four * pi
    H1h = H1h * eV_2_erg

    He1h = ( J_EHe1 * sigma_He1_at_He1_th + &
             J_EHe2 * sigma_He1_at_He2_th ) * four * pi
    He1h = He1h * eV_2_erg

    He2h = zero
    He2h = He2h * eV_2_erg

    if ( reversed ) then

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i, He1i, He2i, &
            dum, dum, dum, &
            H1h, He1h, He2h, &
            dum, dum, dum, Nl  )
       
    end if



  end subroutine set_recomb_photons_isotropic




  !======================================================================
  ! set_column_vs_impact
  !
  ! Calculates the column density as a function of impact parameter. 
  ! One ray for each layer. 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  ! ---------------------------------------------------------------
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH_cc: hydrogen number density in each shell [cm^-3]
  !   nHe_cc: helium number density in each shell [cm^-3]
  ! ---------------------------------------------------------------
  !   Nl: number of layers
  !
  ! Output: 
  !   NH: H column
  !   NHe: He column
  !   NH1: H1 column
  !   NHe1: He1 column
  !   NHe2: He2 column
  !
  !======================================================================
  subroutine set_column_vs_impact( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH_cc, nHe_cc, Tmprtr, &
       NH, NHe, NH1, NHe1, NHe2, Nl )

    ! Arguments
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1), intent(inout) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(inout) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl), intent(inout) :: Edges
    real(real64), dimension(0:Nl-1), intent(inout) ::  nH_cc, nHe_cc
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr
    real(real64), dimension(0:Nl-1), intent(out) :: NH, NHe
    real(real64), dimension(0:Nl-1), intent(out) :: NH1, NHe1, NHe2

    integer(int32), intent(in) :: Nl

    ! Local
    !---------------------------------------------------------------
    integer(int32) :: jnl, m, j, imu, il, j0
    real(real64) :: ds, mu

    logical :: reversed
    real(real64), dimension(0:Nl-1) :: dum

    ! geometry (need monotonic decreasing for this one)
    !---------------------------------------------
    reversed = .false.
    if ( monotonic_increasing( Edges, Nl+1 ) ) then

       reversed = .true.

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH_cc, nHe_cc, Tmprtr, &
            dum, dum, dum, &
            dum, dum, dum, &
            dum, dum, dum, &
            dum, dum, dum, Nl  )
       
    end if

    if ( .not. monotonic_decreasing( Edges, Nl+1 ) ) then
       write(*,*) 'cant make Edges monotonic decreasing' 
       stop
    end if


    NH = zero
    NHe = zero
    NH1 = zero
    NHe1 = zero
    NHe2 = zero

    over_layers: do il = 0, Nl-1
       
       mu = zero
       m = m_for_ray( Edges, il, mu, Nl )

       over_segements: do j = 0, 2*m

          jnl = m-abs(m-j)
          ds = ds_segment( Edges, il, mu, j, Nl )

          NH(il) = NH(il) + ds * nH_cc(jnl)
          NHe(il) = NHe(il) + ds * nHe_cc(jnl)

          NH1(il) = NH1(il) + ds * nH_cc(jnl) * xH1(jnl)
          NHe1(il) = NHe1(il) + ds * nHe_cc(jnl) * xHe1(jnl)
          NHe2(il) = NHe2(il) + ds * nHe_cc(jnl) * xHe2(jnl)
             
       end do over_segements

    end do over_layers
        

    if ( reversed ) then

       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH_cc, nHe_cc, Tmprtr, &
            NH, NHe, NH1, &
            NHe1, NHe2, dum, &
            dum, dum, dum, &
            dum, dum, dum, Nl  )
       
    end if



  end subroutine set_column_vs_impact





  !======================================================================
  ! set_optical_depths
  !
  ! calculate optical depth above and below a layer
  !
  ! Input:
  !   inl: layer index
  !   dtau_H1_th: H1 optical depth through each layer
  !   dtau_He1_th: He1 optical depth through each layer
  !   dtau_He2_th: He2 optical depth through each layer
  !   Nl: number of layers
  !
  ! Output:
  !   tau_H1_th_lo: H1 optical depth below layer
  !   tau_He1_th_lo: He1 optical depth below layer
  !   tau_He2_th_lo: He2 optical depth below layer
  !   tau_H1_th_hi: H1 optical depth above layer
  !   tau_He1_th_hi: He1 optical depth above layer
  !   tau_He2_th_hi: He2 optical depth above layer
  !
  !======================================================================
  subroutine set_optical_depths( Edges, inl, &
       dtau_H1_th, dtau_He1_th, dtau_He2_th, &
       tau_H1_th_lo, tau_He1_th_lo, tau_He2_th_lo, &
       tau_H1_th_hi, tau_He1_th_hi, tau_He2_th_hi, Nl )

    real(real64), dimension(0:Nl), intent(in) :: Edges
    integer(int32), intent(in) :: inl
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_H1_th
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_He1_th
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_He2_th
    real(real64), intent(out) :: tau_H1_th_lo, tau_He1_th_lo, tau_He2_th_lo
    real(real64), intent(out) :: tau_H1_th_hi, tau_He1_th_hi, tau_He2_th_hi
    integer(int32), intent(in) :: Nl
    
    integer(int32) :: ii, ff


    if ( monotonic_decreasing( Edges, Nl+1 ) ) then

       ! in this case Edges[0] = outer most radius and 
       ! Edges[Nl] = inner most radius
       !--------------------------------------------------------------

       if ( inl == 0 ) then

          ! this is the outer most layer
          !-------------------------------------------------------
          
          tau_H1_th_hi = zero
          tau_He1_th_hi = zero
          tau_He2_th_hi = zero
          
          ii = 1
          ff = Nl-1
          tau_H1_th_lo = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_lo = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_lo = sum( dtau_He2_th(ii:ff) )
          
       else if ( inl == Nl-1 ) then

          ! this is the inner most layer
          !-------------------------------------------------------

          ii = 0
          ff = Nl-2
          tau_H1_th_hi = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_hi = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_hi = sum( dtau_He2_th(ii:ff) )          
          
          tau_H1_th_lo = zero
          tau_He1_th_lo = zero
          tau_He2_th_lo = zero          

       else

          ! this is any interior layer
          !-------------------------------------------------------
          ii = 0
          ff = inl-1
          tau_H1_th_hi = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_hi = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_hi = sum( dtau_He2_th(ii:ff) )          
          
          ii = inl+1
          ff = Nl-1
          tau_H1_th_lo = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_lo = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_lo = sum( dtau_He2_th(ii:ff) )          
          
       end if


    else if ( monotonic_increasing( Edges, Nl+1 ) ) then    

       ! in this case Edges[0] = inner most radius and 
       ! Edges[Nl] = outer most radius
       !--------------------------------------------------------------

       if ( inl == 0 ) then

          ! this is the inner most layer
          !-------------------------------------------------------
          
          tau_H1_th_lo = zero
          tau_He1_th_lo = zero
          tau_He2_th_lo = zero
          
          ii = 1
          ff = Nl-1
          tau_H1_th_hi = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_hi = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_hi = sum( dtau_He2_th(ii:ff) )
          
       else if ( inl == Nl-1 ) then

          ! this is the outer most layer
          !-------------------------------------------------------

          ii = 0
          ff = Nl-2
          tau_H1_th_lo = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_lo = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_lo = sum( dtau_He2_th(ii:ff) )          
          
          tau_H1_th_hi = zero
          tau_He1_th_hi = zero
          tau_He2_th_hi = zero          

       else

          ! this is any interior layer
          !-------------------------------------------------------
          ii = 0
          ff = inl-1
          tau_H1_th_lo = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_lo = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_lo = sum( dtau_He2_th(ii:ff) )          
          
          ii = inl+1
          ff = Nl-1
          tau_H1_th_hi = sum( dtau_H1_th(ii:ff) )
          tau_He1_th_hi = sum( dtau_He1_th(ii:ff) )
          tau_He2_th_hi = sum( dtau_He2_th(ii:ff) )          
          
       end if


    else
       write(*,*) 'Edges not monotonic!'
       stop

    end if
       


  end subroutine set_optical_depths



  !======================================================================
  ! reverse_arrays
  !
  ! Reverses the order of important arrays
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  !------------------------------------------------------
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   Tmprtr: temperature in each shell [K]
  !   nH: hydrogen number density in each shell [cm^-3]
  !   nHe: helium number density in each shell [cm^-3]
  !------------------------------------------------------
  !   H1i_src: H1 source photoionization rate 
  !   He1i_src: He1 source photoionization rate 
  !   He2i_src: He2 source photoionization rate 
  !   H1i_rec: H1 recombination photoionization rate 
  !   He1i_rec: He1 recombination photoionization rate 
  !   He2i_rec: He2 recombination photoionization rate 
  !   H1h_src: H1 source photoheating rate
  !   He1h_src: He1 source photoheating rate
  !   He2h_src: He2 source photoheating rate
  !   H1h_rec: H1 recombination photoheating rate
  !   He1h_rec: He1 recombination photoheating rate
  !   He2h_rec: He2 recombination photoheating rate
  !------------------------------------------------------
  !   Nl: number of shells
  !
  !======================================================================
  subroutine reverse_arrays( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, &
       H1i_src, He1i_src, He2i_src, &
       H1i_rec, He1i_rec, He2i_rec, &
       H1h_src, He1h_src, He2h_src, &
       H1h_rec, He1h_rec, He2h_rec, Nl  )

    real(real64), dimension(1:Nl), intent(inout) :: xH1, xH2
    real(real64), dimension(1:Nl), intent(inout) :: xHe1, xHe2, xHe3
    real(real64), dimension(1:Nl+1), intent(inout) :: Edges
    real(real64), dimension(1:Nl), intent(inout) :: nH, nHe
    real(real64), dimension(1:Nl), intent(inout) :: Tmprtr
    real(real64), dimension(1:Nl), intent(inout) :: H1i_src, H1i_rec
    real(real64), dimension(1:Nl), intent(inout) :: He1i_src, He1i_rec
    real(real64), dimension(1:Nl), intent(inout) :: He2i_src, He2i_rec
    real(real64), dimension(1:Nl), intent(inout) :: H1h_src, H1h_rec
    real(real64), dimension(1:Nl), intent(inout) :: He1h_src, He1h_rec
    real(real64), dimension(1:Nl), intent(inout) :: He2h_src, He2h_rec
    integer(int32), intent(in) :: Nl

    xH1 = xH1(Nl:1:-1)
    xH2 = xH2(Nl:1:-1)
    xHe1 = xHe1(Nl:1:-1)
    xHe2 = xHe2(Nl:1:-1)
    xHe3 = xHe3(Nl:1:-1)

    Edges = Edges(Nl+1:1:-1)
    nH = nH(Nl:1:-1)
    nHe = nHe(Nl:1:-1)
    Tmprtr = Tmprtr(Nl:1:-1)

    H1i_src = H1i_src(Nl:1:-1)
    He1i_src = He1i_src(Nl:1:-1)
    He2i_src = He2i_src(Nl:1:-1)

    H1h_src = H1h_src(Nl:1:-1)
    He1h_src = He1h_src(Nl:1:-1)
    He2h_src = He2h_src(Nl:1:-1)

    H1i_rec = H1i_rec(Nl:1:-1)
    He1i_rec = He1i_rec(Nl:1:-1)
    He2i_rec = He2i_rec(Nl:1:-1)

    H1h_rec = H1h_rec(Nl:1:-1)
    He1h_rec = He1h_rec(Nl:1:-1)
    He2h_rec = He2h_rec(Nl:1:-1)

  end subroutine reverse_arrays




  !======================================================================
  ! monotonic_decreasing
  !
  ! Checks if an array is monotonically decreasing. 
  !
  ! Input: 
  !   arr: an array
  !   N: size of array
  !
  !======================================================================
  function monotonic_decreasing(  arr, N  ) result( bool )

    real(real64), dimension(0:N-1), intent(in) :: arr
    integer(int32), intent(in) :: N
    logical :: bool
    integer(int32) :: i
    
    bool = .true.
    do i = 1, N-1
       if ( arr(i) > arr(i-1) ) then
          bool = .false.
          exit
       end if
    end do

  end function monotonic_decreasing

  !======================================================================
  ! monotonic_increasing
  !
  ! Checks if an array is monotonically increasing. 
  !
  ! Input: 
  !   arr: an array
  !   N: size of array
  !
  !======================================================================
  function monotonic_increasing(  arr, N  ) result( bool )

    real(real64), dimension(0:N-1), intent(in) :: arr
    integer(int32), intent(in) :: N
    logical :: bool
    integer(int32) :: i
    
    bool = .true.
    do i = 1, N-1
       if ( arr(i) < arr(i-1) ) then
          bool = .false.
          exit
       end if
    end do

  end function monotonic_increasing


  !======================================================================
  ! get_dr_and_r_c
  !
  ! Calculates radial distance through each layer and central radius of 
  ! each layer. 
  !
  ! Input: 
  !   Edges: radius of each shell
  !   Nl: number of layers
  ! 
  ! Output:
  !   dr: radial distance through each layer 
  !   r_c: central radius of each layer
  !
  !======================================================================
  subroutine get_dr_and_r_c(  Edges, dr, r_c, Nl ) 

    real(real64), dimension(0:Nl), intent(in) :: Edges
    real(real64), dimension(0:Nl-1), intent(out) :: dr
    real(real64), dimension(0:Nl-1), intent(out) :: r_c
    integer(int32), intent(in) :: Nl

    dr = abs( Edges(0:Nl-1) - Edges(1:Nl)  )
    r_c = ( Edges(0:Nl-1) + Edges(1:Nl) ) * half

  end subroutine get_dr_and_r_c



end module sphere_base
