module utils
  use types
  implicit none


contains


  ! numerical integration using the trapezoid rule
  !---------------------------------------------------------
  function utils_integrate_trap( xx, yy, nn ) result( ss )
    real(real64), dimension(0:nn-1), intent(in) :: xx
    real(real64), dimension(0:nn-1), intent(in) :: yy
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-2) :: dx
    real(real64) :: ss

    dx = xx(1:nn-1) - xx(0:nn-2)
    ss = 0.5d0 * sum( ( yy(1:nn-1) + yy(0:nn-2) ) * dx )

  end function utils_integrate_trap


  ! linear interpolation 
  !
  ! xx and yy define the function to be interpolated
  ! xr are the requested x values
  ! yr are the interpolated y values
  !---------------------------------------------------------
  function utils_interp_lin_1d( xx, yy, xr, nn, nr ) result( yr )
    real(real64), dimension(0:nn-1), intent(in) :: xx
    real(real64), dimension(0:nn-1), intent(in) :: yy
    real(real64), dimension(0:nr-1), intent(in) :: xr
    integer(int32), intent(in) :: nn
    integer(int32), intent(in) :: nr
    real(real64), dimension(0:nr-1) :: yr

    real(real64) :: frac, dx, dy
    integer(int32) :: ii, ir, i_lo, i_hi


    ! check input
    !------------------------------
    if ( minval(xr) < xx(0) ) then
       write(*,*) 'intrp values out of range'
       stop
    end if

    if ( maxval(xr) > xx(nn-1) ) then
       write(*,*) 'intrp values out of range'
       stop
    end if


    ! loop over interpolation requests
    !----------------------------------
    do ir = 0, nr-1

       ! find bracketing indices
       !----------------------------------
       do ii = 0, nn-1
          if ( xx(ii) > xr(ir) ) then
             i_lo = ii - 1
             i_hi = i_lo + 1
             exit
          end if
       end do

       ! interpolate
       !----------------------------------
       dx = xx(i_hi) - xx(i_lo) 
       dy = yy(i_hi) - yy(i_lo) 

       frac = ( xr(ir) - xx(ii) ) / dx
       yr(ir) = yy(ii) + frac * dy

    end do
    

  end function utils_interp_lin_1d


end module utils



