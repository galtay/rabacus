!
! Routines useful to all planar gas distributions. 
!
!------------------------------------------------------------------------

module slab_base
  use types
  use physical_constants, only: ev_2_erg, pi

  use chem_cool_rates, only: get_kchem

  use photo_xsections, only: return_sigma_H1
  use photo_xsections, only: return_sigma_He1
  use photo_xsections, only: return_sigma_He2

  use source_background, only: bgnd_src_return_ionrate_shld_ei
  use source_background, only: bgnd_src_return_heatrate_shld_ei

  use source_plane, only: pln_src_return_ionrate_shld
  use source_plane, only: pln_src_return_ionrate_thin
  use source_plane, only: pln_src_return_heatrate_shld

!  use special_functions, only: E1
  use zhang_jin_f, only: e1xa

  use utils, only: utils_integrate_trap

  implicit none


  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: half = 0.5d0
  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: two = 2.0d0
  real(real64), parameter :: four = 4.0d0

  logical, parameter :: i_ray_H2 = .true.
  logical, parameter :: i_ray_He2 = .true.
  logical, parameter :: i_ray_He3 = .true.

contains


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
  subroutine set_optical_depths( &
       inl, dtau_H1_th, dtau_He1_th, dtau_He2_th, &
       tau_H1_th_lo, tau_He1_th_lo, tau_He2_th_lo, &
       tau_H1_th_hi, tau_He1_th_hi, tau_He2_th_hi, Nl )

    integer(int32), intent(in) :: inl
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_H1_th
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_He1_th
    real(real64), dimension(0:Nl-1), intent(in) :: dtau_He2_th
    real(real64), intent(out) :: tau_H1_th_lo, tau_He1_th_lo, tau_He2_th_lo
    real(real64), intent(out) :: tau_H1_th_hi, tau_He1_th_hi, tau_He2_th_hi
    integer(int32), intent(in) :: Nl
    
    integer(int32) :: ii, ff

    if ( inl == 0 ) then
       tau_H1_th_lo = zero
       tau_He1_th_lo = zero
       tau_He2_th_lo = zero          

       ii = 1
       ff = Nl-1
       tau_H1_th_hi = sum( dtau_H1_th(ii:ff) )
       tau_He1_th_hi = sum( dtau_He1_th(ii:ff) )
       tau_He2_th_hi = sum( dtau_He2_th(ii:ff) )

    else if ( inl == Nl-1 ) then
       ii = 0
       ff = Nl-2
       tau_H1_th_lo = sum( dtau_H1_th(ii:ff) )
       tau_He1_th_lo = sum( dtau_He1_th(ii:ff) )
       tau_He2_th_lo = sum( dtau_He2_th(ii:ff) )          

       tau_H1_th_hi = zero
       tau_He1_th_hi = zero
       tau_He2_th_hi = zero          

    else
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

  end subroutine set_optical_depths





  !======================================================================
  ! set_recomb_photons
  !
  ! Calculates the photoionization and photoheating rates in each layer due 
  ! to the recombination photons emitted from every other layer.  Recombination
  ! radiation is always emitted isotropically
  ! 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  ! ---------------------------------------------------------------
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH: hydrogen number density in each shell [cm^-3]
  !   nHe: helium number density in each shell [cm^-3]
  !   Tmprtr: temperature in each shell [K]
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
  !   Nl: number of shells
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

  subroutine set_recomb_photons( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_photo_fit, i_rate_fit, &
       H1i, He1i, He2i, H1h, He1h, He2h, Nl )

    ! Arguments
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1), intent(in) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(in) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl), intent(in) :: Edges
    real(real64), dimension(0:Nl-1), intent(in) ::  nH, nHe
    real(real64), dimension(0:Nl-1), intent(in) :: Tmprtr
    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th
    integer(int32), intent(in) :: i_photo_fit
    integer(int32), intent(in) :: i_rate_fit
    real(real64), dimension(0:Nl-1), intent(out) :: H1i, He1i, He2i
    real(real64), dimension(0:Nl-1), intent(out) :: H1h, He1h, He2h
    integer(int32), intent(in) :: Nl

    ! Local
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1) :: dl, z_c
    real(real64), dimension(0:Nl-1) :: dNH, dNHe

    real(real64), dimension(0:Nl-1) :: dtau_H1_th
    real(real64), dimension(0:Nl-1) :: dtau_He1_th
    real(real64), dimension(0:Nl-1) :: dtau_He2_th

    real(real64) :: dtau_H1_th_k
    real(real64) :: dtau_He1_th_k
    real(real64) :: dtau_He2_th_k

    integer(int32), parameter :: Nnu = 1
    real(real64), dimension(0:Nnu-1) :: E_eV
    real(real64) :: sigma_H1_at_H1_th
    real(real64) :: sigma_H1_at_He1_th
    real(real64) :: sigma_H1_at_He2_th
    real(real64) :: sigma_He1_at_He1_th
    real(real64) :: sigma_He1_at_He2_th
    real(real64) :: sigma_He2_at_He2_th

    real(real64) :: sigma_H1_ra_H1
    real(real64) :: sigma_H1_ra_He1
    real(real64) :: sigma_H1_ra_He2
    real(real64) :: sigma_He1_ra_He1
    real(real64) :: sigma_He1_ra_He2
    real(real64) :: sigma_He2_ra_He2


    real(real64), dimension(0:Nl-1) :: fcA_H2, fcA_He2, fcA_He3
    real(real64), dimension(0:Nl-1) :: reH2_A, reHe2_A, reHe3_A
    real(real64), dimension(0:Nl-1) :: reH2_B, reHe2_B, reHe3_B
    real(real64), dimension(0:Nl-1) :: reH2_1, reHe2_1, reHe3_1
    real(real64), dimension(0:Nl-1) :: ciH1, ciHe1, ciHe2

    real(real64), dimension(0:Nl-1) :: ne
    real(real64), dimension(0:Nl-1) :: em_H2, em_He2, em_He3
    real(real64), dimension(0:Nl-1) :: J_H2, J_He2, J_He3
    real(real64), dimension(0:Nl-1) :: J_H2_dum, J_He2_dum, J_He3_dum


    ! for integrating 
    real(real64), dimension(0:Nl) :: em_H2_e, em_He2_e, em_He3_e
    real(real64), dimension(0:Nl) :: tau_H2_lo, tau_He2_lo, tau_He3_lo
    real(real64), dimension(0:Nl) :: tau_H2_hi, tau_He2_hi, tau_He3_hi
    real(real64), dimension(0:Nl) :: J_H2_lo, J_He2_lo, J_He3_lo
    real(real64), dimension(0:Nl) :: J_H2_hi, J_He2_hi, J_He3_hi
    real(real64), dimension(0:Nl) :: J_H2_e, J_He2_e, J_He3_e
    real(real64), dimension(0:Nl) :: E1_H2_lo, E1_He2_lo, E1_He3_lo
    real(real64), dimension(0:Nl) :: E1_H2_hi, E1_He2_hi, E1_He3_hi

    integer(int32) :: k, i, nn
    real(real64) :: tau

    real(real64) :: sigma_H1, sigma_H1_ra
    real(real64) :: sigma_He1, sigma_He1_ra
    real(real64) :: sigma_He2, sigma_He2_ra



    ! geometry
    !---------------------------------------------
    dl = Edges(1:Nl) - Edges(0:Nl-1)
    z_c = Edges(0:Nl-1) + dl * half
    dNH = dl * nH
    dNHe = dl * nHe

    dtau_H1_th = dNH * xH1 * sigma_H1_th 
    dtau_He1_th = dNHe * xHe1 * sigma_He1_th
    dtau_He2_th = dNHe * xHe2 * sigma_He2_th

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

    ! calculate cross-section ratios
    !------------------------------------------------------------
    sigma_H1_ra_H1 = one
    sigma_H1_ra_He1 = sigma_H1_at_He1_th / sigma_H1_at_H1_th    
    sigma_H1_ra_He2 = sigma_H1_at_He2_th / sigma_H1_at_H1_th

    sigma_He1_ra_He1 = one
    sigma_He1_ra_He2 = sigma_He1_at_He2_th / sigma_He1_at_He1_th

    sigma_He2_ra_He2 = one



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

    ! calculate emission coefficients in layers
    ! photons / ( s cm^3 sr )
    !------------------------------------------------------------
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    em_H2 = ( reH2_1 * nH * xH2 * ne ) / (four * pi)
    em_He2 = ( reHe2_1 * nHe * xHe2 * ne ) / (four * pi)
    em_He3 = ( reHe3_1 * nHe * xHe3 * ne ) / (four * pi)


    ! calculate emission coefficients on edges
    ! photons / ( s cm^3 sr )
    !------------------------------------------------------------
    do k = 0, Nl
       if (k==0) then
          em_H2_e(k) = em_H2(k) * half 
          em_He2_e(k) = em_He2(k) * half 
          em_He3_e(k) = em_He3(k) * half 

       else if ( k==Nl ) then
          em_H2_e(k) = em_H2(k-1) * half 
          em_He2_e(k) = em_He2(k-1) * half 
          em_He3_e(k) = em_He3(k-1) * half 

       else
          em_H2_e(k) = ( em_H2(k-1) + em_H2(k) ) * half
          em_He2_e(k) = ( em_He2(k-1) + em_He2(k) ) * half
          em_He3_e(k) = ( em_He3(k-1) + em_He3(k) ) * half

       end if
    end do

    
    ! calculate optical depths below each edge 
    !------------------------------------------------
    do k = 0, Nl

       if ( k == 0 ) then
          
          tau_H2_lo(k) = zero
          tau_He2_lo(k) = zero
          tau_He3_lo(k) = zero

       else

          tau_H2_lo(k) = tau_H2_lo(k-1) + &
               dtau_H1_th(k-1) * sigma_H1_ra_H1
          
          tau_He2_lo(k) = tau_He2_lo(k-1) + &
               dtau_H1_th(k-1) * sigma_H1_ra_He1 + &
               dtau_He1_th(k-1) * sigma_He1_ra_He1

          tau_He3_lo(k) = tau_He3_lo(k-1) + &
               dtau_H1_th(k-1) * sigma_H1_ra_He2 + &
               dtau_He1_th(k-1) * sigma_He1_ra_He2 + &
               dtau_He2_th(k-1) * sigma_He2_ra_He2
          
       end if

    end do


    ! calculate right moving "ray", i.e. moving from
    ! low to high z
    !------------------------------------------------
    J_H2_lo = zero    
    J_He2_lo = zero
    J_He3_lo = zero

    do k = 1, Nl
       J_H2_lo(k) = J_H2_lo(k-1) + &
            em_H2(k-1) * dl(k-1) * e1xa( tau_H2_lo(k) - tau_H2_lo(k-1) )
       J_He2_lo(k) = J_He2_lo(k-1) + &
            em_He2(k-1) * dl(k-1) * e1xa( tau_He2_lo(k) - tau_He2_lo(k-1) )
       J_He3_lo(k) = J_He3_lo(k-1) + &
            em_He3(k-1) * dl(k-1) * e1xa( tau_He3_lo(k) - tau_He3_lo(k-1) )
    end do

    ! calculate left moving "ray", i.e. moving from
    ! high to low z
    !------------------------------------------------
    J_H2_hi = zero    
    J_He2_hi = zero
    J_He3_hi = zero

    do k = Nl-1, 0, -1
       J_H2_hi(k) = J_H2_hi(k+1) + &
            em_H2(k) * dl(k) * e1xa( tau_H2_lo(k+1) - tau_H2_lo(k) )
       J_He2_hi(k) = J_He2_hi(k+1) + &
            em_He2(k) * dl(k) * e1xa( tau_He2_lo(k+1) - tau_He2_lo(k) )
       J_He3_hi(k) = J_He3_hi(k+1) + &
            em_He3(k) * dl(k) * e1xa( tau_He3_lo(k+1) - tau_He3_lo(k) )
    end do


    ! add lo and hi to get J at each edge
    !------------------------------------------------
    J_H2_e = J_H2_lo + J_H2_hi
    J_He2_e = J_He2_lo + J_He2_hi
    J_He3_e = J_He3_lo + J_He3_hi
    
    ! average edge values to get layer values
    !------------------------------------------------
    J_H2 = ( J_H2_e(0:Nl-1) + J_H2_e(1:Nl) ) * half 
    J_He2 = ( J_He2_e(0:Nl-1) + J_He2_e(1:Nl) ) * half 
    J_He3 = ( J_He3_e(0:Nl-1) + J_He3_e(1:Nl) ) * half 

    ! four pi for calculating photo rates
    ! half for Eq. 36 in rabacus paper
    J_H2 = J_H2 * four * pi * half
    J_He2 = J_He2 * four * pi * half
    J_He3 = J_He3 * four * pi * half



    ! photo-ionization rates
    !----------------------------------------------
    H1i = J_He3 * sigma_H1_at_He2_th + &
          J_He2 * sigma_H1_at_He1_th + &
          J_H2  * sigma_H1_at_H1_th

    He1i = J_He3 * sigma_He1_at_He2_th + &
           J_He2 * sigma_He1_at_He1_th

    He2i = J_He3 * sigma_He2_at_He2_th


    ! photo-heating rates
    !----------------------------------------------
    H1h = J_He3 * ( E_He2_th - E_H1_th ) * sigma_H1_at_He2_th * eV_2_erg + &
          J_He2 * ( E_He1_th - E_H1_th ) * sigma_H1_at_He1_th * eV_2_erg

    He1h = J_He3 * ( E_He2_th - E_H1_th ) * sigma_He1_at_He2_th * eV_2_erg

    He2h = zero
    


  end subroutine set_recomb_photons







end module slab_base
