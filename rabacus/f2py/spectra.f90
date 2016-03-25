! This module sets energy and shape arrays for various spectral types.
! 
! A spectrum is defined by two components E_eV and shape.  The first is 
! simply an array of the energies sampled for the spectrum.  The second is 
! the shape of the spectrum, but the units are determined by what kind of 
! source it is. 
!
!    E_eV: energies sampled [eV]
!
!    Monochromatic:
!       point source: 
!          shape: Lu [erg/s] 
!       plane source:
!          shape: Fu [erg/(s cm^2)]
!       background source:
!          shape: Inu [erg/(s cm^2 sr)]
!        
!    Polychromatic (thermal, powerlaw, HM12 ... )
!       point source: 
!          shape: dLu/dnu [erg/(s Hz)] 
!       plane source:
!          shape: dFu/dnu [erg/(s cm^2 Hz)]
!       background source:
!          shape: Inu [erg/(s cm^2 sr Hz)]
! 
! Types: 
!     Monochromatic: Single energy spectrum
!     Thermal: Blackbody spectrum
!     Powerlaw: Powerlaw spectrum 
!     HM12: Haardt and Madau 2012 
!
!--------------------------------------------------------------------------

module spectra
  use types
  use physical_constants, only: h_erg_s, h_eV_s, kb_eV_K, c_cm_s
  use photo_xsections
  use hm12, only: hm12_spectrum_at_z
  implicit none


contains


  ! Constructs an array of energies with the first entry equal to E_min and 
  ! the last entry equal to E_max. Units should be in eV. 
  !------------------------------------------------------------------------
  subroutine set_E_array( E_min, E_max, Npts, E_eV )
    real(real64), intent(in) :: E_min 
    real(real64), intent(in) :: E_max
    integer(int32), intent(in) :: Npts
    real(real64), dimension(0:Npts-1), intent(out) :: E_eV

    real(real64) :: log_E_min
    real(real64) :: log_E_max
    real(real64) :: d_logE
    integer(int32) :: i
    real(real64), parameter :: eps = 1.0e-12

    if ( Npts == 1 ) then

       if ( E_min /= E_max ) then
          write(*,*) 'E_min and E_max must be equal for monochromatic spectra'
          stop
       end if
       E_eV = E_min
       
    else

       log_E_min = log10( E_min ) + eps
       log_E_max = log10( E_max ) - eps
       
       d_logE = ( log_E_max - log_E_min ) / ( Npts-1 )
       do i = 0,Npts-1
          E_eV(i) = 10.0d0**(log_E_min + d_logE * i)
       end do
       
    end if

  end subroutine set_E_array


  ! For each element in a supplied list of energies, E_snap (the snap array), 
  ! find the element in E_eV_in that is closest to it and change the value in 
  ! E_eV_in to the value in E_snap.  This can be used (for example) to ensure 
  ! that the energy array includes values that are exactly at the ionization 
  ! energy thresholds.  Units should be in eV. 
  !------------------------------------------------------------------------
  function snap_E_array( E_snap, E_eV_in, n_snap, n_E ) result( E_eV_out )
    real(real64), dimension(0:n_snap-1), intent(in) :: E_snap
    real(real64), dimension(0:n_E-1), intent(in) :: E_eV_in
    integer(int32), intent(in) :: n_snap
    integer(int32), intent(in) :: n_E
    real(real64), dimension(0:n_E-1) :: E_eV_out

    real(real64), dimension(0:n_E-1) :: diff

    integer(int32) :: i_snap, indx

    E_eV_out = E_eV_in

    do i_snap = 0, n_snap-1

       diff = abs( E_eV_in - E_snap(i_snap) )
       indx = sum( minloc( diff ) ) - 1
       E_eV_out(indx) = E_snap(i_snap)

    end do

  end function snap_E_array



  ! Returns an un-normlized intensity array for monochromatic spectrum.
  ! 
  ! Input: 
  !    E_eV: Input array of energy samples [eV] (must be size=1)
  !
  ! Output: 
  !    shape: Unormalized array that gives spectral shape. (1.0)
  !------------------------------------------------------------------------
  subroutine set_shape_monochromatic( E_eV, nn, shape )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1), intent(out) :: shape
    
    if ( nn /= 1 ) then
       write(*,*) 'E_eV must have size=1 for monochromatic spectrum. ' 
       stop
    end if
    shape = 1.0d0
       
  end subroutine set_shape_monochromatic



  ! Returns an un-normlized intensity array in the shape of a black body 
  ! spectrum. 
  ! 
  ! Input: 
  !    E_eV: Input array of energy samples [eV]
  !    T_eff: Effective temperature for blackbody [K]. 
  !
  ! Output: 
  !    shape: Unormalized array that gives spectral shape.
  !------------------------------------------------------------------------
  subroutine set_shape_thermal( E_eV, T_eff, nn, shape )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), intent(in) :: T_eff
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1), intent(out) :: shape

    real(real64), dimension(0:nn-1) :: nu_Hz
    real(real64), dimension(0:nn-1) :: t1, t2, e1

    nu_Hz = E_eV / h_eV_s

    t1 = 2.0d0 * h_erg_s * nu_Hz**3 / c_cm_s**2
    e1 = E_eV / ( kb_eV_K * T_eff )
    t2 = exp(e1)
    
    shape = t1 / (t2-1.0d0)
       
  end subroutine set_shape_thermal



  ! Returns an un-normlized intensity array in the shape of a powerlaw
  ! spectrum. 
  ! 
  ! Input: 
  !    E_eV: Input array of energy samples [eV]
  !    alpha: Slope of spectrum, shape = (E/E[0])^alpha
  !
  ! Output: 
  !    shape: Unormalized array that gives spectral shape.
  !------------------------------------------------------------------------
  subroutine set_shape_powerlaw( E_eV, alpha, nn, shape )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), intent(in) :: alpha
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1), intent(out) :: shape
    
    shape = ( E_eV / E_eV(0) )**alpha
       
  end subroutine set_shape_powerlaw


  ! Returns an un-normlized intensity array in the shape of the HM12 
  ! model for the UV background
  ! 
  ! Input: 
  !    E_eV: Input array of energy samples [eV]
  !    z: desired redshift
  !
  ! Output: 
  !    shape: Unormalized array that gives spectral shape.
  !------------------------------------------------------------------------
  subroutine set_shape_hm12( E_eV, z, nn, shape )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    real(real64), intent(in) :: z
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1), intent(out) :: shape
    
    call hm12_spectrum_at_z( z, E_eV, shape, nn ) 
       
  end subroutine set_shape_hm12







end module spectra
