module legendre_polynomial
implicit none

contains

subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end subroutine imtqlx


subroutine p_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! P_POLYNOMIAL evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ), intent(in) ::  m
  integer ( kind = 4 ), intent(in) ::  n

  integer ( kind = 4 ) i
  real ( kind = 8 ), intent(out) ::  v(m,0:n)
  real ( kind = 8 ), intent(in) ::  x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
  end do
 
  return
end subroutine p_polynomial

subroutine p_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
!
!  First terms:
!
!     1
!     0     1
!    -1/2   0      3/2
!     0    -3/2    0     5/2
!     3/8   0    -30/8   0     35/8
!     0    15/8    0   -70/8    0     63/8
!    -5/16  0    105/16  0   -315/16   0    231/16
!     0   -35/16   0   315/16   0   -693/16   0    429/16
!
!     1.00000
!     0.00000  1.00000
!    -0.50000  0.00000  1.50000
!     0.00000 -1.50000  0.00000  2.5000
!     0.37500  0.00000 -3.75000  0.00000  4.37500
!     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
!    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
!     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the 
!    Legendre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ), intent(in) ::  n

  real ( kind = 8 ), intent(out) ::  c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n <= 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0:i-2) =          real (   - i + 1, kind = 8 ) * c(i-2,0:i-2) &
                        / real (     i,     kind = 8 )
    c(i,1:i) = c(i,1:i) + real ( i + i - 1, kind = 8 ) * c(i-1,0:i-1) &
                        / real (     i,     kind = 8 )
  end do
 
  return
end subroutine p_polynomial_coefficients


subroutine p_polynomial_prime ( m, n, x, vp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VP(M,0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = 4 ), intent(in) ::  m
  integer ( kind = 4 ), intent(in) ::  n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ), intent(out) ::  vp(m,0:n)
  real ( kind = 8 ), intent(in) ::  x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00
  vp(1:m,0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0D+00
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
    vp(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * ( v(1:m,i-1) &
                                                   + x(1:m) * vp(1:m,i-1) ) &
                - real (     i - 1, kind = 8 ) *   vp(1:m,i-2)               ) &
                / real (     i,     kind = 8 )
 
  end do
 
  return
end subroutine p_polynomial_prime


subroutine p_polynomial_prime2 ( m, n, x, vpp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!    P"(0,X) = 0
!    P"(1,X) = 0
!    P"(N,X) = ( (2*N-1)*(2*P(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VPP(M,0:N), the second derivative of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = 4 ), intent(in) ::  m
  integer ( kind = 4 ), intent(in) ::  n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) vp(m,0:n)
  real ( kind = 8 ), intent(out) ::  vpp(m,0:n)
  real ( kind = 8 ), intent(in) ::  x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00
  vp(1:m,0) = 0.0D+00
  vpp(1:m,0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0D+00
  vpp(1:m,1) = 0.0D+00
 
  do i = 2, n
 
    v(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
      - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
      / real (     i,     kind = 8 )
 
    vp(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * ( v(1:m,i-1) + x(1:m) * vp(1:m,i-1) ) &
      - real (     i - 1, kind = 8 ) *   vp(1:m,i-2)               ) &
      / real (     i,     kind = 8 )

    vpp(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * ( 2.0D+00 * vp(1:m,i-1) &
                                         + x(1:m) * vpp(1:m,i-1) ) &
      - real (     i - 1, kind = 8 ) *   vpp(1:m,i-2)               ) &
      / real (     i,     kind = 8 )

  end do
 
  return
end subroutine p_polynomial_prime2


subroutine p_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUES returns values of the Legendre polynomials P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.4062500000000000D+00, &
    -0.3359375000000000D+00, &
     0.1577148437500000D+00, &
     0.3397216796875000D+00, &
     0.2427673339843750D-01, &
    -0.2799186706542969D+00, &
    -0.1524540185928345D+00, &
     0.1768244206905365D+00, &
     0.2212002165615559D+00, &
     0.0000000000000000D+00, &
    -0.1475000000000000D+00, &
    -0.2800000000000000D+00, &
    -0.3825000000000000D+00, &
    -0.4400000000000000D+00, &
    -0.4375000000000000D+00, &
    -0.3600000000000000D+00, &
    -0.1925000000000000D+00, &
     0.8000000000000000D-01, &
     0.4725000000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end subroutine p_polynomial_values


subroutine p_polynomial_zeros ( nt, t )

!*****************************************************************************80
!
!! P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), the zeros.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), intent(out) :: t(nt)
  real ( kind = 8 ) wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( 2.0D+00 )

  call imtqlx ( nt, t, bj, wts )

  return
end subroutine p_polynomial_zeros



subroutine p_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights
!    of the rule.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), intent(out) :: t(nt)
  real ( kind = 8 ), intent(out) :: wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1) = sqrt ( 2.0D+00 )
  wts(2:nt) = 0.0D+00

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2

  return
end subroutine p_quadrature_rule



end module legendre_polynomial
