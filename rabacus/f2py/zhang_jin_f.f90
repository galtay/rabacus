!
! more useful routines from 
!  http://people.sc.fsu.edu/~jburkardt/f_src/special_functions
!

module zhang_jin_f
implicit none



contains



function e1xa ( x ) result( e1 )

!*****************************************************************************80
!
!! E1XA computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
  implicit none

  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ) :: e1

  real ( kind = 8 ) es1
  real ( kind = 8 ) es2


  if ( x == 0.0D+00 ) then

    e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

    e1 = - log ( x ) + (((( &
        1.07857D-03 * x &
      - 9.76004D-03 ) * x &
      + 5.519968D-02 ) * x &
      - 0.24991055D+00 ) * x &
      + 0.99999193D+00 ) * x &
      - 0.57721566D+00

  else

    es1 = ((( x &
      +  8.5733287401D+00 ) * x &
      + 18.059016973D+00  ) * x &
      +  8.6347608925D+00 ) * x &
      +  0.2677737343D+00

    es2 = ((( x &
      +  9.5733223454D+00 ) * x &
      + 25.6329561486D+00 ) * x &
      + 21.0996530827D+00 ) * x &
      +  3.9584969228D+00

    e1 = exp ( - x ) / x * es1 / es2

  end if

  return
end function e1xa


function e2xa ( x ) result( e2 )

!*****************************************************************************80
!
!! E2XA computes the exponential integral E2(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
  implicit none

  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ) :: e1
  real ( kind = 8 ) :: e2

  real ( kind = 8 ) es1
  real ( kind = 8 ) es2


  if ( x == 0.0D+00 ) then

    e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

    e1 = - log ( x ) + (((( &
        1.07857D-03 * x &
      - 9.76004D-03 ) * x &
      + 5.519968D-02 ) * x &
      - 0.24991055D+00 ) * x &
      + 0.99999193D+00 ) * x &
      - 0.57721566D+00

  else

    es1 = ((( x &
      +  8.5733287401D+00 ) * x &
      + 18.059016973D+00  ) * x &
      +  8.6347608925D+00 ) * x &
      +  0.2677737343D+00

    es2 = ((( x &
      +  9.5733223454D+00 ) * x &
      + 25.6329561486D+00 ) * x &
      + 21.0996530827D+00 ) * x &
      +  3.9584969228D+00

    e1 = exp ( - x ) / x * es1 / es2

  end if

  e2 = exp(-x) - x * e1

  return
end function e2xa




subroutine e1xb ( x, e1 )

!*****************************************************************************80
!
!! E1XB computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
  implicit none

  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ), intent(out) :: e1

  real ( kind = 8 ) ga
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) t0


  if ( x == 0.0D+00 ) then

    e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

    e1 = 1.0D+00
    r = 1.0D+00

    do k = 1, 25
      r = -r * k * x / ( k + 1.0D+00 )**2
      e1 = e1 + r
      if ( abs ( r ) <= abs ( e1 ) * 1.0D-15 ) then
        exit
      end if
    end do

    ga = 0.5772156649015328D+00
    e1 = - ga - log ( x ) + x * e1

  else

    m = 20 + int ( 80.0D+00 / x )
    t0 = 0.0D+00
    do k = m, 1, -1
      t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
    end do
    t = 1.0D+00 / ( x + t0 )
    e1 = exp ( -x ) * t

  end if

  return
end subroutine e1xb








subroutine enxa ( n, x, en )

!*****************************************************************************80
!
!! ENXA computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EN(0:N), the function values.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n
  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ), intent(out) :: en(0:n)

  real ( kind = 8 ) e1
  real ( kind = 8 ) ek
  integer ( kind = 4 ) k


  en(0) = exp ( - x ) / x 
  call e1xb ( x, e1 )

  en(1) = e1
  do k = 2, n
    ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
    en(k) = ek
    e1 = ek
  end do

  return
end subroutine enxa



subroutine enxb ( n, x, en )

!*****************************************************************************80
!
!! ENXB computes the exponential integral En(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    10 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) EN(0:N), the function values.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: n
  real ( kind = 8 ), intent(in) :: x
  real ( kind = 8 ), intent(out) :: en(0:n)

  real ( kind = 8 ) ens
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) ps
  real ( kind = 8 ) r
  real ( kind = 8 ) rp
  real ( kind = 8 ) s
  real ( kind = 8 ) s0
  real ( kind = 8 ) t
  real ( kind = 8 ) t0

  if ( x == 0.0D+00 ) then

    en(0) = 1.0D+300
    en(1) = 1.0D+300
    do k = 2, n
      en(k) = 1.0D+00 / ( k - 1.0D+00 )
    end do
    return

  else if ( x <= 1.0D+00 ) then

    en(0) = exp ( - x ) / x
    do l = 1, n
      rp = 1.0D+00
      do j = 1, l - 1
        rp = - rp * x / j
      end do
      ps = -0.5772156649015328D+00
      do m = 1, l - 1
        ps = ps + 1.0D+00 / m
      end do
      ens = rp * ( - log ( x ) + ps )
      s = 0.0D+00
      do m = 0, 20
        if ( m /= l - 1 ) then
          r = 1.0D+00
          do j = 1, m
            r = - r * x / j
          end do
          s = s + r / ( m - l + 1.0D+00 )
          if ( abs ( s - s0 ) < abs ( s ) * 1.0D-15 ) then
            exit
          end if
          s0 = s
        end if
      end do

      en(l) = ens - s

    end do

  else

    en(0) = exp ( - x ) / x
    m = 15 + int ( 100.0D+00 / x )
    do l = 1, n
      t0 = 0.0D+00
      do k = m, 1, -1
        t0 = ( l + k - 1.0D+00 ) / ( 1.0D+00 + k / ( x + t0 ) )
      end do
      t = 1.0D+00 / ( x + t0 )
      en(l) = exp ( - x ) * t
    end do

  end if

  return
end subroutine enxb




end module zhang_jin_f
