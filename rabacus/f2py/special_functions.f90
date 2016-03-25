module special_functions
use slatec, only: DE1
use zhang_jin, only: e1xb
use zhang_jin_f, only: e1xa
use types
implicit none


integer(int32) :: table_E1_N
real(real64) :: table_E1_x_lo, table_E1_log_x_lo
real(real64) :: table_E1_x_hi, table_E1_log_x_hi
real(real64) :: table_E1_dlog_x, table_E1_dx
real(real64), allocatable :: table_E1_x(:), table_E1_log_x(:)
real(real64), allocatable :: table_E1_v(:), table_E1_log_v(:)

integer(int32) :: table_E2_N
real(real64) :: table_E2_x_lo, table_E2_log_x_lo
real(real64) :: table_E2_x_hi, table_E2_log_x_hi
real(real64) :: table_E2_dlog_x, table_E2_dx
real(real64), allocatable :: table_E2_x(:), table_E2_log_x(:)
real(real64), allocatable :: table_E2_v(:), table_E2_log_v(:)


logical :: use_lookup = .false.
logical :: intrp_lin = .true.



contains





! Exponential integral E1
!--------------------------------------------
function E1( x )
  real(real64), intent(in) :: x
  real(real64) :: E1

  if ( use_lookup ) then
     E1 = table_E1_interpolate( x )
  else
     E1 = e1xa(x)
     !  call e1xb( x, E1 )
     !  E1 = DE1(x)
  endif
     
end function E1



! Exponential integral E2
!--------------------------------------------
function E2( x ) 
  real(real64), intent(in) :: x
  real(real64) :: E2

  if ( use_lookup ) then
     E2 = table_E2_interpolate( x )
  else
     E2 = exp(-x) - x * E1(x)
  endif

end function E2


function E2xa( x ) 
  real(real64) :: x
  real(real64) :: E2xa
  real(real64) :: E1tmp
  E1tmp = e1xa( x )
  E2xa = exp(-x) - x * E1tmp
end function E2xa


function E2xb( x ) 
  real(real64) :: x
  real(real64) :: E2xb
  real(real64) :: E1tmp
  call e1xb( x, E1tmp )
  E2xb = exp(-x) - x * E1tmp
end function E2xb





! creates a lookup table for E1 
!--------------------------------------------
subroutine table_E1_create( N, x_lo, x_hi ) 
  integer(int32), intent(in) :: N
  real(real64), intent(in) :: x_lo, x_hi
  real(real64) :: log_x, x
  integer(int32) :: i


  if ( allocated(table_E1_x) ) then
     if ( N == size(table_E1_x) ) then
        return
     else
        deallocate(table_E1_x, table_E1_log_x)  
        deallocate(table_E1_v, table_E1_log_v)  
     end if
  end if

  table_E1_N = N

  table_E1_x_lo = x_lo
  table_E1_x_hi = x_hi
  table_E1_dx = ( table_E1_x_hi - table_E1_x_lo ) / ( N - 1 )

  table_E1_log_x_lo = log10( x_lo )
  table_E1_log_x_hi = log10( x_hi )
  table_E1_dlog_x = ( table_E1_log_x_hi - table_E1_log_x_lo ) / ( N - 1 )

  allocate( table_E1_x(0:N-1), table_E1_log_x(0:N-1) )
  allocate( table_E1_v(0:N-1), table_E1_log_v(0:N-1) )

  do i = 0, N-1
     if (intrp_lin) then
        x = table_E1_x_lo + table_E1_dx * i
        log_x = log10( x )
     else
        log_x = table_E1_log_x_lo + table_E1_dlog_x * i
        x = 10.0d0**log_x
     end if
     table_E1_x(i) = x
     table_E1_log_x(i) = log_x
     !
     call e1xb( x, table_E1_v(i) )
     !
  end do

  table_E1_log_v = log10( table_E1_v )
  
end subroutine table_E1_create


function table_E1_interpolate( x ) result( E1_intrp )
  real(real64) :: x
  real(real64) :: log_E1, E1_intrp
  integer(int32) :: ix
  real(real64) :: log_x, frac


  ! create lookup table if it doesn't exist
  !----------------------------------------------
  if ( .not. allocated( table_E1_x ) ) then
     call table_E1_create( 1024, 1.0d-6, 2.0d1 )
  end if


  if ( x >= table_E1_x_hi ) then
!     call e1xb( x, E1_intrp )
     E1_intrp = table_E1_v( table_E1_N-1 )

  else if ( x <= table_E1_x_lo ) then
     call e1xb( x, E1_intrp )

  else
     if ( intrp_lin ) then
        frac = (x - table_E1_x_lo) / table_E1_dx
        ix = int( frac )
        frac = frac - ix
        E1_intrp = table_E1_v(ix) + &
             ( table_E1_v(ix+1) - table_E1_v(ix) ) * frac

     else
        log_x = log10( x )
        frac = (log_x - table_E1_log_x_lo) / table_E1_dlog_x
        ix = int( frac )
        frac = frac - ix
        log_E1 = table_E1_log_v(ix) + &
             ( table_E1_log_v(ix+1) - table_E1_log_v(ix) ) * frac
        E1_intrp = 10.0d0**log_E1

     end if
     
  end if

end function table_E1_interpolate




! creates a lookup table for E2 
!--------------------------------------------
subroutine table_E2_create( N, x_lo, x_hi ) 
  integer(int32), intent(in) :: N
  real(real64), intent(in) :: x_lo, x_hi
  real(real64) :: log_x, x
  integer(int32) :: i


  if ( allocated(table_E2_x) ) then
     if ( N == size(table_E2_x) ) then
        return
     else
        deallocate(table_E2_x, table_E2_log_x)  
        deallocate(table_E2_v, table_E2_log_v)  
     end if
  end if

  table_E2_N = N

  table_E2_x_lo = x_lo
  table_E2_x_hi = x_hi
  table_E2_dx = ( table_E2_x_hi - table_E2_x_lo ) / ( N - 1 )

  table_E2_log_x_lo = log10( x_lo )
  table_E2_log_x_hi = log10( x_hi )
  table_E2_dlog_x = ( table_E2_log_x_hi - table_E2_log_x_lo ) / ( N - 1 )

  allocate( table_E2_x(0:N-1), table_E2_log_x(0:N-1) )
  allocate( table_E2_v(0:N-1), table_E2_log_v(0:N-1) )

  do i = 0, N-1
     if (intrp_lin) then
        x = table_E2_x_lo + table_E2_dx * i
        log_x = log10( x )
     else
        log_x = table_E2_log_x_lo + table_E2_dlog_x * i
        x = 10.0d0**log_x
     end if
     table_E2_x(i) = x
     table_E2_log_x(i) = log_x
     !
     table_E2_v(i) = E2xb(x)
     !
  end do

  table_E2_log_v = log10( table_E2_v )
  
end subroutine table_E2_create


function table_E2_interpolate( x ) result( E2_intrp )
  real(real64) :: x
  real(real64) :: log_E2, E2_intrp
  integer(int32) :: ix
  real(real64) :: log_x, frac


  ! create lookup table if it doesn't exist
  !----------------------------------------------
  if ( .not. allocated( table_E2_x ) ) then
     call table_E2_create( 1024, 1.0d-6, 2.0d1 )
  end if


  if ( x >= table_E2_x_hi ) then
!     call e1xb( x, E2_intrp )
     E2_intrp = table_E2_v( table_E2_N-1 )

  else if ( x <= table_E2_x_lo ) then
     E2_intrp = E2xb( x )

  else
     if ( intrp_lin ) then
        frac = (x - table_E2_x_lo) / table_E2_dx
        ix = int( frac )
        frac = frac - ix
        E2_intrp = table_E2_v(ix) + &
             ( table_E2_v(ix+1) - table_E2_v(ix) ) * frac

     else
        log_x = log10( x )
        frac = (log_x - table_E2_log_x_lo) / table_E2_dlog_x
        ix = int( frac )
        frac = frac - ix
        log_E2 = table_E2_log_v(ix) + &
             ( table_E2_log_v(ix+1) - table_E2_log_v(ix) ) * frac
        E2_intrp = 10.0d0**log_E2

     end if

  end if

end function table_E2_interpolate




end module special_functions
