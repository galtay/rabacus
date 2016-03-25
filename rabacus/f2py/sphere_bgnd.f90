!
! Routines to calculate ionization and / or temperature equilibrium for 
! spherically symmetric gas distributions with uniform background sources
!
!------------------------------------------------------------------------


module sphere_bgnd

  use types
  use physical_constants, only: pi

  use ion_solver, only: solve_pce, solve_pcte
  use chem_cool_rates, only: get_kchem

  use photo_xsections, only: return_E_H1_th, return_sigma_H1
  use photo_xsections, only: return_E_He1_th, return_sigma_He1
  use photo_xsections, only: return_E_He2_th, return_sigma_He2

  use source_background, only: bgnd_src_return_ionrate_thin
  use source_background, only: bgnd_src_return_ionrate_shld
  use source_background, only: bgnd_src_return_heatrate_thin
  use source_background, only: bgnd_src_return_heatrate_shld

  use sphere_base, only: reverse_arrays
  use sphere_base, only: monotonic_decreasing
  use sphere_base, only: monotonic_increasing
  use sphere_base, only: set_recomb_photons_outward
  use sphere_base, only: set_recomb_photons_radial
  use sphere_base, only: set_recomb_photons_isotropic

  use geometry, only: m_for_ray, ds_segment
  use legendre_polynomial, only: p_quadrature_rule

  implicit none


  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: half = 0.5d0
  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: two = 2.0d0
  real(real64), parameter :: three = 3.0d0
  real(real64), parameter :: four = 4.0d0

  integer(int32), parameter :: MAX_ITER = 500
  real(real64), parameter :: TOL_FRAC = 1.0d-2

  
contains
  
  !======================================================================
  ! sphere_bgnd_solve 
  !
  ! Main driver for solution.  Shell refers to the Nl spherical shells. 
  !
  ! Input: 
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH: hydrogen number density in each shell [cm^-3]
  !   nHe: helium number density in each shell [cm^-3]
  !   Tmprtr: temperature in each shell [K]
  !   Nmu: number of polar angle samples
  ! -----------------------------------------------------------------
  !   bgnd_src_eV: energy samples for spectrum [eV]
  !   bgnd_src_shape: shape of spectrum polychromatic=[erg/s/Hz], 
  !                                     monochromatic=[erg/s]
  ! -----------------------------------------------------------------
  !   i_rec_meth: how to treat recombinations? {1=fixed, 2=outward, 3=radial}
  !   fixed_fcA: If i_rec_meth=1, constant caseA fraction
  ! -----------------------------------------------------------------
  !   i_photo_fit: photo_xsection fits {1=verner96}
  !   i_rate_fit: atomic rate fits {1=hg97}
  !   i_find_Teq: solve for eqlbrm T {0=no,1=yes} if no, Tmprtr is constant
  !   i_thin: return optically thin values
  ! ---------------------------------------------------------------
  !   em_H1_fac: multiplicative factor for H1 recomb emission
  !   em_He1_fac: multiplicative factor for He1 recomb emission
  !   em_He2_fac: multiplicative factor for He2 recomb emission
  ! -----------------------------------------------------------------
  !   z: redshift, only used if i_find_Teq=1
  !   Hz: hubble parameter [1/s] @ z (0.0 for no hubble cooling)
  !   tol: tolerance for all convergence tests
  !   Nl: number of shells
  !   Nnu: number of energy samples in spectrum
  !
  ! Output: 
  !   Tmprtr: temperature in each shell [K]
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
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
  !
  !======================================================================
  
  subroutine sphere_bgnd_solve( &
       Edges, nH, nHe, Tmprtr, Nmu, &
       bgnd_src_eV, bgnd_src_shape, &
       i_rec_meth, fixed_fcA, &
       i_photo_fit, i_rate_fit, &
       i_find_Teq, i_thin, &
       em_H1_fac, em_He1_fac, em_He2_fac, &
       z, Hz, tol, Nl, Nnu, &
       xH1, xH2, xHe1, xHe2, xHe3, &
       H1i_src, He1i_src, He2i_src, &
       H1i_rec, He1i_rec, He2i_rec, &
       H1h_src, He1h_src, He2h_src, &
       H1h_rec, He1h_rec, He2h_rec )


    ! arguments
    !--------------------------------------------------------    
    real(real64), dimension(0:Nl), intent(inout) :: Edges
    real(real64), dimension(0:Nl-1), intent(inout) :: nH, nHe
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr  
    integer(int32), intent(in) :: Nmu

    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_eV  
    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_shape 

    integer(int32), intent(in) :: i_rec_meth
    real(real64), intent(in) :: fixed_fcA

    integer(int32), intent(in) :: i_photo_fit   
    integer(int32), intent(in) :: i_rate_fit    
    integer(int32), intent(in) :: i_find_Teq
    integer(int32), intent(in) :: i_thin

    real(real64), intent(in) :: em_H1_fac, em_He1_fac, em_He2_fac

    real(real64), intent(in) :: z          
    real(real64), intent(in) :: Hz
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: Nl       
    integer(int32), intent(in) :: Nnu       

    real(real64), dimension(0:Nl-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(out) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl-1), intent(out) :: H1i_src, H1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1i_src, He1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2i_src, He2i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: H1h_src, H1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1h_src, He1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2h_src, He2h_rec


    ! local
    !--------------------------------------------------------    
    logical :: reversed
    real(real64) :: E_H1_th, E_He1_th, E_He2_th
    real(real64) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), dimension(0:Nnu-1) :: sigma_H1, sigma_He1, sigma_He2

    integer(int32) :: iter
    logical :: not_converged
    real(real64), dimension(0:Nl-1) :: ne
    real(real64), dimension(0:Nl-1) :: conv_old, conv_new, conv_change






    !=======================================
    ! Initialize 
    !=======================================

    ! we don't make any assumptions except monotonicity in the 
    ! Edges array.  here we make sure it is decreasing. 

    reversed = .false.

    if ( monotonic_increasing(  Edges, Nl+1 ) ) then

       reversed = .true.
       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i_src, He1i_src, He2i_src, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_src, He1h_src, He2h_src, &
            H1h_rec, He1h_rec, He2h_rec, Nl  )

    else if ( monotonic_decreasing(  Edges, Nl+1 ) ) then
       ! do nothing

    else
       write(*,*) 'Edges not monotonic!'
       stop

    end if

    ! ionization thresholds 
    !---------------------------------------------
    E_H1_th = return_E_H1_th( i_rate_fit )
    E_He1_th = return_E_He1_th( i_rate_fit )
    E_He2_th = return_E_He2_th( i_rate_fit )

    ! photoionization xsections at thresholds
    ! (scalar quantities)
    !---------------------------------------------
    sigma_H1_th  = sum( return_sigma_H1(  (/E_H1_th/), i_photo_fit, 1 ) )
    sigma_He1_th = sum( return_sigma_He1( (/E_He1_th/), i_photo_fit, 1 ) )
    sigma_He2_th = sum( return_sigma_He2( (/E_He2_th/), i_photo_fit, 1 ) )

    ! photoionization xsections
    ! (frequency dependent vectors)
    !---------------------------------------------
    sigma_H1 = return_sigma_H1( bgnd_src_eV, i_photo_fit, Nnu )
    sigma_He1 = return_sigma_He1( bgnd_src_eV, i_photo_fit, Nnu )
    sigma_He2 = return_sigma_He2( bgnd_src_eV, i_photo_fit, Nnu )


    !=======================================
    ! Set optically thin values 
    !=======================================
    call sphere_bgnd_optically_thin( &
         Edges, nH, nHe, Tmprtr, &
         bgnd_src_eV, bgnd_src_shape, &
         sigma_H1, sigma_He1, sigma_He2, &
         sigma_H1_th, sigma_He1_th, sigma_He2_th, &
         E_H1_th, E_He1_th, E_He2_th, &
         i_rec_meth, fixed_fcA, &
         i_find_Teq, z, Hz, i_rate_fit, tol, Nl, Nnu, &
         xH1, xH2, xHe1, xHe2, xHe3, &
         H1i_src, He1i_src, He2i_src, &
         H1i_rec, He1i_rec, He2i_rec, &
         H1h_src, He1h_src, He2h_src, &
         H1h_rec, He1h_rec, He2h_rec  ) 

    if ( i_thin == 1 ) then
       if ( reversed ) then
          call reverse_arrays( &
               xH1, xH2, xHe1, xHe2, xHe3, &
               Edges, nH, nHe, Tmprtr, &
               H1i_src, He1i_src, He2i_src, &
               H1i_rec, He1i_rec, He2i_rec, &
               H1h_src, He1h_src, He2h_src, &
               H1h_rec, He1h_rec, He2h_rec, Nl  )
       end if       
       return
    end if


    !=======================================
    ! Perform sweeps 
    !=======================================
    iter = 0
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    conv_old = ne
    not_converged = .true. 
    do while( not_converged )
          
       call sphere_bgnd_sweep( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, Nmu, &
            bgnd_src_eV, bgnd_src_shape, &
            sigma_H1, sigma_He1, sigma_He2, &
            sigma_H1_th, sigma_He1_th, sigma_He2_th, &
            E_H1_th, E_He1_th, E_He2_th, &
            i_rec_meth, fixed_fcA, i_find_Teq, z, Hz, &
            em_H1_fac, em_He1_fac, em_He2_fac, &
            i_photo_fit, i_rate_fit, tol, Nl, Nnu, &
            H1i_src, He1i_src, He2i_src, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_src, He1h_src, He2h_src, &
            H1h_rec, He1h_rec, He2h_rec )
       

       ! check convergence 
       !------------------------------------------------
       if ( iter > MAX_ITER ) then
          write(*,*) "MAX ITER in sphere_bgnd_solve"
          stop
       end if
       iter = iter + 1
       ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe      
    
       conv_new = ne
       conv_change = conv_new / conv_old - one
       if ( all( abs(conv_change) < tol ) ) then
          not_converged = .false.
       end if
       conv_old = conv_new


    end do

    ! return arrays in the same order they were input
    !-------------------------------------------------------
    if ( reversed ) then
       call reverse_arrays( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            H1i_src, He1i_src, He2i_src, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_src, He1h_src, He2h_src, &
            H1h_rec, He1h_rec, He2h_rec, Nl  )
    end if


  end subroutine sphere_bgnd_solve



  !======================================================================
  ! sphere_bgnd_optically_thin
  !
  ! Sets optically thin ionization fractions 
  !
  ! Input: 
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH: hydrogen number density in each shell [cm^-3]
  !   nHe: helium number density in each shell [cm^-3]
  !   Tmprtr: temperature in each shell [K]
  ! ---------------------------------------------------------------
  !   bgnd_src_eV: energy samples for spectrum [eV]
  !   bgnd_src_shape: shape of spectrum polychromatic=[erg/s/Hz], 
  !                                     monochromatic=[erg/s]
  ! ---------------------------------------------------------------
  !   sigma_H1: H1 photo_xsection at each energy in bgnd_src_eV [cm^2]
  !   sigma_He1: He1 photo_xsection at each energy in bgnd_src_eV [cm^2]
  !   sigma_He2: He2 photo_xsection at each energy in bgnd_src_eV [cm^2]
  ! ---------------------------------------------------------------
  !   sigma_H1_th: H1 photoion xsection at H1 threshold
  !   sigma_He1_th: He1 photoion xsection at He1 threshold
  !   sigma_He2_th: He2 photoion xsection at He2 threshold
  ! ---------------------------------------------------------------
  !   E_H1_th: H1 ionization threshold [eV]
  !   E_He1_th: He1 ionization threshold [eV]
  !   E_He2_th: He2 ionization threshold [eV]
  ! ---------------------------------------------------------------
  !   i_rec_meth: how to treat recombinations? {1=fixed, 2=outward, 3=radial}
  !   fixed_fcA: If i_rec_meth=1, constant caseA fraction
  ! ---------------------------------------------------------------
  !   i_find_Teq: solve for eqlibrm T {0=no,1=yes} if no, Tmprtr is constant
  !   z: redshift, only used if i_find_Teq=1
  !   Hz: hubble parameter [1/s] @ z (0 for no hubble cooling)
  !   i_rate_fit: atomic rate fits {1=hg97}
  ! ---------------------------------------------------------------
  !   tol: tolerance for all convergence tests
  !   Nl: number of shells
  !   Nnu: number of energy samples in spectrum
  !
  ! Output: 
  !   Tmprtr: temperature in each layer [K]
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
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
  !
  !======================================================================

  subroutine sphere_bgnd_optically_thin( &
       Edges, nH, nHe, Tmprtr, &
       bgnd_src_eV, bgnd_src_shape, &
       sigma_H1, sigma_He1, sigma_He2, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_rec_meth, fixed_fcA, i_find_Teq, z, Hz, &
       i_rate_fit, tol, Nl, Nnu, &
       xH1, xH2, xHe1, xHe2, xHe3, &
       H1i_src, He1i_src, He2i_src, &
       H1i_rec, He1i_rec, He2i_rec, &
       H1h_src, He1h_src, He2h_src, &
       H1h_rec, He1h_rec, He2h_rec )

    ! Arguments
    !-------------------------------------------------------------
    real(real64), dimension(0:Nl), intent(in) :: Edges
    real(real64), dimension(0:Nl-1), intent(in) ::  nH, nHe
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr    
    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_eV
    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_shape
    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_H1
    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_He1
    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_He2
    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th

    integer(int32), intent(in) :: i_rec_meth
    real(real64), intent(in) :: fixed_fcA
    integer(int32), intent(in) :: i_find_Teq
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    integer(int32), intent(in) :: i_rate_fit
    real(real64), intent(in) :: tol

    integer(int32), intent(in) :: Nl
    integer(int32), intent(in) :: Nnu

    real(real64), dimension(0:Nl-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(out) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:Nl-1), intent(out) :: H1i_src, H1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1i_src, He1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2i_src, He2i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: H1h_src, H1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1h_src, He1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2h_src, He2h_rec

    ! Local
    !---------------------------------------------------------------
    real(real64), dimension(0:Nl-1) :: fcA_H2, fcA_He2, fcA_He3
    real(real64), dimension(0:Nl-1) :: reH2, reHe2, reHe3
    real(real64), dimension(0:Nl-1) :: ciH1, ciHe1, ciHe2    


    ! get chemistry rates
    !---------------------------------------------
    if ( i_rec_meth == 1 ) then
       fcA_H2 = fixed_fcA
       fcA_He2 = fixed_fcA
       fcA_He3 = fixed_fcA
    else 
       fcA_H2 = one
       fcA_He2 = one
       fcA_He3 = one
    end if

    call get_kchem( Tmprtr, fcA_H2, fcA_He2, fcA_He3, &
         i_rate_fit, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, Nl )


    ! get source optically thin photo ionization/heating rates 
    !---------------------------------------------
    H1i_src(:) = bgnd_src_return_ionrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_H1, Nnu ) 
    He1i_src(:) = bgnd_src_return_ionrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_He1, Nnu ) 
    He2i_src(:) = bgnd_src_return_ionrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_He2, Nnu ) 
    
    H1h_src(:) = bgnd_src_return_heatrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_H1, E_H1_th, Nnu ) 
    He1h_src(:) = bgnd_src_return_heatrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_He1, E_He1_th, Nnu ) 
    He2h_src(:) = bgnd_src_return_heatrate_thin( &
         bgnd_src_eV, bgnd_src_shape, sigma_He2, E_He2_th, Nnu ) 




    ! solve for equilibrium with only source photons
    !------------------------------------------------

    if ( i_find_Teq == 1 ) then
          
       call solve_pcte( nH, nHe, &
            H1i_src, He1i_src, He2i_src, &
            H1h_src, He1h_src, He2h_src, &
            z, Hz, fcA_H2, fcA_He2, fcA_He3, i_rate_fit, &
            xH1, xH2, xHe1, xHe2, xHe3, Tmprtr, tol*TOL_FRAC, Nl )

    else

       call solve_pce( nH, nHe, reH2, reHe2, reHe3, &
            ciH1, ciHe1, ciHe2, H1i_src, He1i_src, He2i_src, &
            xH1, xH2, xHe1, xHe2, xHe3, tol*TOL_FRAC, Nl )
       
    end if


    ! if not solving for recombination photons were done
    !----------------------------------------------------------------
    H1i_rec = zero
    He1i_rec = zero
    He2i_rec = zero

    H1h_rec = zero
    He1h_rec = zero
    He2h_rec = zero


  end subroutine sphere_bgnd_optically_thin





  !======================================================================
  ! sphere_bgnd_sweep
  !
  ! Performs one sweep through the shells from small to large radii. 
  !
  ! Input: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  !------------------------------------------------------
  !   Edges: radius of the edges of the shells [cm] (Nl+1 entries)
  !   nH: hydrogen number density in each shell [cm^-3]
  !   nHe: helium number density in each shell [cm^-3]
  !   Tmprtr: temperature in each shell [K]
  !   Nmu: number of polar angle samples
  !------------------------------------------------------
  !   bgnd_src_eV: energy samples for spectrum [eV]
  !   bgnd_src_shape: shape of spectrum polychromatic=[erg/s/Hz], 
  !                                     monochromatic=[erg/s]
  !------------------------------------------------------
  !   sigma_H1: H1 photo_xsection at each energy in bgnd_src_eV [cm^2]
  !   sigma_He1: He1 photo_xsection at each energy in bgnd_src_eV [cm^2]
  !   sigma_He2: He2 photo_xsection at each energy in bgnd_src_eV [cm^2]
  !------------------------------------------------------
  !   sigma_H1_th: H1 photo_xsection at H1 ionization threshold [cm^2]
  !   sigma_He1_th: He1 photo_xsection at He1 ionization threshold [cm^2]
  !   sigma_He2_th: He2 photo_xsection at He2 ionization threshold [cm^2]
  !------------------------------------------------------
  !   E_H1_th: H1 ionization threshold [eV]
  !   E_He1_th: He1 ionization threshold [eV]
  !   E_He2_th: He2 ionization threshold [eV]
  !------------------------------------------------------
  !   i_rec_meth: how to treat recombinations? {1=fixed, 2=outward, 3=radial}
  !   fixed_fcA: If i_rec_meth=1, constant caseA fraction
  !------------------------------------------------------ 
  !   i_find_Teq: solve for eqlbrm T {0=no,1=yes} if no, Tmprtr is constant
  !   z: redshift, only used if i_find_Teq=1
  !   Hz: hubble parameter [1/s] @ z (0 for no hubble cooling)
  ! ---------------------------------------------------------------
  !   em_H1_fac: multiplicative factor for H1 recomb emission
  !   em_He1_fac: multiplicative factor for He1 recomb emission
  !   em_He2_fac: multiplicative factor for He2 recomb emission
  ! ---------------------------------------------------------------
  !   i_photo_fit: photo_xsection fits {1=verner96}
  !   i_rate_fit: atomic rate fits {1=hg97}
  !------------------------------------------------------ 
  !   tol: tolerance for all convergence tests
  !   Nl: number of shells
  !   Nnu: number of energy samples in spectrum
  !
  ! Output: 
  !   xH1: H1 ionization fraction = nH1 / nH
  !   xH2: H2 ionization fraction = nH2 / nH
  !   xHe1: He1 ionization fraction = nHe1 / nHe
  !   xHe2: He2 ionization fraction = nHe2 / nHe
  !   xHe3: He3 ionization fraction = nHe3 / nHe
  !   Tmprtr: temperature in each shell [K]
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
  !
  !======================================================================
  subroutine sphere_bgnd_sweep( &
       xH1, xH2, xHe1, xHe2, xHe3, &
       Edges, nH, nHe, Tmprtr, Nmu, &
       bgnd_src_eV, bgnd_src_shape, &
       sigma_H1, sigma_He1, sigma_He2, &
       sigma_H1_th, sigma_He1_th, sigma_He2_th, &
       E_H1_th, E_He1_th, E_He2_th, &
       i_rec_meth, fixed_fcA, i_find_Teq, z, Hz, &
       em_H1_fac, em_He1_fac, em_He2_fac, &
       i_photo_fit, i_rate_fit, tol, Nl, Nnu, &
       H1i_src, He1i_src, He2i_src, &
       H1i_rec, He1i_rec, He2i_rec, &
       H1h_src, He1h_src, He2h_src, &
       H1h_rec, He1h_rec, He2h_rec  )

    ! arguments
    !--------------------------------------------------------------
    real(real64), dimension(0:Nl-1), intent(inout) :: xH1, xH2
    real(real64), dimension(0:Nl-1), intent(inout) :: xHe1, xHe2, xHe3

    real(real64), dimension(0:Nl), intent(inout) :: Edges
    real(real64), dimension(0:Nl-1), intent(inout) :: nH, nHe
    real(real64), dimension(0:Nl-1), intent(inout) :: Tmprtr
    integer(int32), intent(in) :: Nmu

    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_eV
    real(real64), dimension(0:Nnu-1), intent(in) :: bgnd_src_shape

    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_H1
    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_He1
    real(real64), dimension(0:Nnu-1), intent(in) :: sigma_He2

    real(real64), intent(in) :: sigma_H1_th, sigma_He1_th, sigma_He2_th
    real(real64), intent(in) :: E_H1_th, E_He1_th, E_He2_th

    integer(int32), intent(in) :: i_rec_meth
    real(real64), intent(in) :: fixed_fcA
    integer(int32), intent(in) :: i_find_Teq
    real(real64), intent(in) :: z          
    real(real64), intent(in) :: Hz

    real(real64), intent(in) :: em_H1_fac, em_He1_fac, em_He2_fac

    integer(int32), intent(in) :: i_photo_fit
    integer(int32), intent(in) :: i_rate_fit
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: Nl
    integer(int32), intent(in) :: Nnu

    real(real64), dimension(0:Nl-1), intent(out) :: H1i_src, H1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1i_src, He1i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2i_src, He2i_rec
    real(real64), dimension(0:Nl-1), intent(out) :: H1h_src, H1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He1h_src, He1h_rec
    real(real64), dimension(0:Nl-1), intent(out) :: He2h_src, He2h_rec


    ! locals
    !--------------------------------------------------------------
    real(real64), dimension(0:Nl-1) :: fcA_H2, fcA_He2, fcA_He3

    real(real64), dimension(0:Nnu-1) :: sigma_H1_ra
    real(real64), dimension(0:Nnu-1) :: sigma_He1_ra
    real(real64), dimension(0:Nnu-1) :: sigma_He2_ra

    real(real64), dimension(0:Nnu-1) :: tau_H1, tau_He1, tau_He2

    real(real64), dimension(0:Nmu-1) :: H1i_mu, H1h_mu
    real(real64), dimension(0:Nmu-1) :: He1i_mu, He1h_mu
    real(real64), dimension(0:Nmu-1) :: He2i_mu, He2h_mu
    real(real64), dimension(0:Nmu-1) :: xi, wi

    real(real64) :: ds, dmu, mu
    real(real64) :: NH1, NHe1, NHe2
    real(real64) :: tau_H1_th, tau_He1_th, tau_He2_th

    real(real64), dimension(0:Nl-1) :: H1i, He1i, He2i
    real(real64), dimension(0:Nl-1) :: H1h, He1h, He2h

    real(real64) :: reH2, reHe2, reHe3
    real(real64) :: ciH1, ciHe1, ciHe2    

    integer(int32) :: il, jnl, m, j, imu, j0



    ! set static variables
    !===========================================================

    ! sigma ratios
    !---------------------------------------------
    sigma_H1_ra = sigma_H1 / sigma_H1_th
    sigma_He1_ra = sigma_He1 / sigma_He1_th
    sigma_He2_ra = sigma_He2 / sigma_He2_th


    ! case A fractions
    !---------------------------------------------
    if ( i_rec_meth == 1 ) then
       fcA_H2 = fixed_fcA
       fcA_He2 = fixed_fcA
       fcA_He3 = fixed_fcA
    else 
       fcA_H2 = one
       fcA_He2 = one
       fcA_He3 = one
    end if


    ! calculate effect of recombinationo photons at current state
    ! i_rec_meth = 1 is fixed recombination rates 
    ! i_rec_meth = 2 is outward only
    ! i_rec_meth = 3 is radial only 
    ! i_rec_meth = 4: isotropic
    !===========================================================

    if ( i_rec_meth == 1 ) then
       H1i_rec = zero
       He1i_rec = zero
       He2i_rec = zero
       H1h_rec = zero
       He1h_rec = zero
       He2h_rec = zero

    else if ( i_rec_meth == 2 ) then

       call set_recomb_photons_outward( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            sigma_H1_th, sigma_He1_th, sigma_He2_th, &
            E_H1_th, E_He1_th, E_He2_th, &
            i_photo_fit, i_rate_fit, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_rec, He1h_rec, He2h_rec, Nl )

    else if ( i_rec_meth == 3 ) then

       call set_recomb_photons_radial( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, &
            sigma_H1_th, sigma_He1_th, sigma_He2_th, &
            E_H1_th, E_He1_th, E_He2_th, &
            i_photo_fit, i_rate_fit, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_rec, He1h_rec, He2h_rec, Nl )

    else if ( i_rec_meth == 4 ) then

       call set_recomb_photons_isotropic( &
            xH1, xH2, xHe1, xHe2, xHe3, &
            Edges, nH, nHe, Tmprtr, Nmu, &
            sigma_H1_th, sigma_He1_th, sigma_He2_th, &
            E_H1_th, E_He1_th, E_He2_th, &
            i_photo_fit, i_rate_fit, &
            em_H1_fac, em_He1_fac, em_He2_fac, &
            H1i_rec, He1i_rec, He2i_rec, &
            H1h_rec, He1h_rec, He2h_rec, Nl )

    end if



    ! calculate xi and wi for Gauss-Legendre Quadrature
    !============================================================

    call p_quadrature_rule( Nmu, xi, wi )


    ! loop over layers and calculate photoionization
    ! and photoheating rates
    !============================================================

    ! begin parallel region
    !============================================================

    !$omp parallel private( il, imu, dmu, mu, m, j, j0, jnl, ds, &
    !$omp&  NH1, NHe1, NHe2, tau_H1_th, tau_He1_th, tau_He2_th, & 
    !$omp&  tau_H1, tau_He1, tau_He2, H1i_mu, He1i_mu, He2i_mu, &
    !$omp&  H1h_mu, He1h_mu, He2h_mu )

    !$omp do
    over_layers_1: do il = 0, Nl-1


       ! calculate rates in different directions
       !============================================================
       over_mu_1: do imu = 0, Nmu-1

          mu = xi(imu)
          m = m_for_ray( Edges, il, mu, Nl )

          if ( mu <= zero ) then
             j0 = il
          else
             j0 = 2*m - il
          end if

          NH1 = zero
          NHe1 = zero
          NHe2 = zero

          over_segements: do j = j0, 2*m

             jnl = m-abs(m-j)
             ds = ds_segment( Edges, il, mu, j, Nl )

             if ( j == j0 ) then
                ds = ds * half
             end if

             NH1 = NH1 + ds * nH(jnl) * xH1(jnl)
             NHe1 = NHe1 + ds * nHe(jnl) * xHe1(jnl)
             NHe2 = NHe2 + ds * nHe(jnl) * xHe2(jnl)

          end do over_segements

          tau_H1_th = NH1 * sigma_H1_th
          tau_He1_th = NHe1 * sigma_He1_th
          tau_He2_th = NHe2 * sigma_He2_th

          tau_H1 = tau_H1_th * sigma_H1_ra
          tau_He1 = tau_He1_th * sigma_He1_ra
          tau_He2 = tau_He2_th * sigma_He2_ra


          H1i_mu(imu) = bgnd_src_return_ionrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_H1, &
               tau_H1, tau_He1, tau_He2, Nnu ) 


          He1i_mu(imu) = bgnd_src_return_ionrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_He1, &
               tau_H1, tau_He1, tau_He2, Nnu ) 


          He2i_mu(imu) = bgnd_src_return_ionrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_He2, &
               tau_H1, tau_He1, tau_He2, Nnu ) 



          H1h_mu(imu) = bgnd_src_return_heatrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_H1, E_H1_th, &
               tau_H1, tau_He1, tau_He2, Nnu ) 

          He1h_mu(imu) = bgnd_src_return_heatrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_He1, E_He1_th, &
               tau_H1, tau_He1, tau_He2, Nnu ) 

          He2h_mu(imu) = bgnd_src_return_heatrate_shld( &
               bgnd_src_eV, bgnd_src_shape, sigma_He2, E_He2_th, &
               tau_H1, tau_He1, tau_He2, Nnu ) 


       end do over_mu_1


       ! integrate over directions
       !============================================================
       over_mu_2: do imu = 0, Nmu-1
          dmu = wi(imu)

          H1i_src(il) = H1i_src(il) + H1i_mu(imu) * dmu
          He1i_src(il) = He1i_src(il) + He1i_mu(imu) * dmu
          He2i_src(il) = He2i_src(il) + He2i_mu(imu) * dmu

          H1h_src(il) = H1h_src(il) + H1h_mu(imu) * dmu
          He1h_src(il) = He1h_src(il) + He1h_mu(imu) * dmu
          He2h_src(il) = He2h_src(il) + He2h_mu(imu) * dmu
       end do over_mu_2


       H1i_src(il) = H1i_src(il) * half 
       He1i_src(il) = He1i_src(il) * half 
       He2i_src(il) = He2i_src(il) * half 
       
       H1h_src(il) = H1h_src(il) * half 
       He1h_src(il) = He1h_src(il) * half 
       He2h_src(il) = He2h_src(il) * half 


       ! set total photoionization rates
       !---------------------------------------------------
       if ( i_rec_meth == 1 ) then
          H1i(il) = H1i_src(il) 
          He1i(il) = He1i_src(il) 
          He2i(il) = He2i_src(il) 
       else 
          H1i(il) = H1i_src(il) + H1i_rec(il) 
          He1i(il) = He1i_src(il) + He1i_rec(il) 
          He2i(il) = He2i_src(il) + He2i_rec(il) 
       end if

       ! set total photoheating rates
       !---------------------------------------------------
       if ( i_rec_meth == 1 ) then
          H1h(il) = H1h_src(il) 
          He1h(il) = He1h_src(il) 
          He2h(il) = He2h_src(il) 
       else 
          H1h(il) = H1h_src(il) + H1h_rec(il) 
          He1h(il) = He1h_src(il) + He1h_rec(il) 
          He2h(il) = He2h_src(il) + He2h_rec(il) 
       end if


       ! solve for equilibrium
       !------------------------------------------------
       if ( i_find_Teq == 1 ) then
          
          call solve_pcte( &
               nH(il), nHe(il), H1i(il), He1i(il), He2i(il), &
               H1h(il), He1h(il), He2h(il), z, Hz, &
               fcA_H2(il), fcA_He2(il), fcA_He3(il), i_rate_fit, &
               xH1(il), xH2(il), xHe1(il), xHe2(il), xHe3(il), &
               Tmprtr(il), tol*TOL_FRAC, 1 )
          
       else
          
          call get_kchem( &
               Tmprtr(il), fcA_H2(il), fcA_He2(il), fcA_He3(il), &
               i_rate_fit, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, 1 )
          
          call solve_pce( &
               nH(il), nHe(il), reH2, reHe2, reHe3, &
               ciH1, ciHe1, ciHe2, H1i(il), He1i(il), He2i(il), &
               xH1(il), xH2(il), xHe1(il), xHe2(il), xHe3(il), &
               tol*TOL_FRAC, 1 )
          
       end if


    end do over_layers_1
    !$omp end do

    !$omp end parallel 


                       

  end subroutine sphere_bgnd_sweep




end module sphere_bgnd
