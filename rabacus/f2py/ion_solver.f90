module ion_solver
  use types
  use chem_cool_rates, only: get_kchem, get_kcool, get_heat, get_cool
  implicit none

  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: half = 0.5d0
  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: two = 2.0d0
  real(real64), parameter :: four = 4.0d0


  integer(int32), parameter :: MAX_ITER = 5000
  real(real64), parameter :: TFLOOR = 7.5d3
  real(real64), parameter :: TMAX = 1.0d5


! overload scalar and vector versions of functions
!----------------------------------------------------------
  interface fpc
     module procedure fpc_s, fpc_v
  end interface fpc

  interface check_nan
     module procedure check_nan_s, check_nan_v
  end interface check_nan

  interface analytic_H1
     module procedure analytic_H1_s, analytic_H1_v
  end interface analytic_H1

  interface implicit_ionfracs
     module procedure implicit_ionfracs_s, implicit_ionfracs_v
  end interface implicit_ionfracs

  interface solve_ce
     module procedure solve_ce_s, solve_ce_v
  end interface solve_ce

  interface solve_pce
     module procedure solve_pce_s, solve_pce_v
  end interface solve_pce

  interface solve_pcte
     module procedure solve_pcte_s, solve_pcte_v
  end interface solve_pcte

  interface get_dudT
     module procedure get_dudT_s, get_dudT_v
  end interface get_dudT


contains



!====================================================================
!====================================================================
!
! Floating point correction 
!
!====================================================================
!====================================================================


  ! Enforces ionization fractions sum to one
  !----------------------------------------------------------------------
  subroutine fpc_s( xH1, xH2, xHe1, xHe2, xHe3, nn )
    real(real64), intent(inout) :: xH1, xH2
    real(real64), intent(inout) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn
    real(real64) :: total
    total = xH1 + xH2
    xH1 = xH1 / total
    xH2 = xH2 / total    
    total = xHe1 + xHe2 + xHe3
    xHe1 = xHe1 / total
    xHe2 = xHe2 / total
    xHe3 = xHe3 / total
  end subroutine fpc_s

  subroutine fpc_v( xH1, xH2, xHe1, xHe2, xHe3, nn )
    real(real64), dimension(0:nn-1), intent(inout) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(inout) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: total
    total = xH1 + xH2
    xH1 = xH1 / total
    xH2 = xH2 / total    
    total = xHe1 + xHe2 + xHe3
    xHe1 = xHe1 / total
    xHe2 = xHe2 / total
    xHe3 = xHe3 / total
  end subroutine fpc_v


!====================================================================
!====================================================================
!
! Check for NaN
!
!====================================================================
!====================================================================


  ! check nan
  !----------------------------------------------------------------------
  function check_nan_s( xH1, xH2, xHe1, xHe2, xHe3, tag, nn ) result( has_nan )
    real(real64), intent(in) :: xH1, xH2
    real(real64), intent(in) :: xHe1, xHe2, xHe3
    character(1024), intent(in) :: tag
    integer(int32), intent(in) :: nn
    logical :: has_nan

    has_nan = .false. 
    if ( isnan(xH1) ) has_nan = .true.
    if ( isnan(xH2) ) has_nan = .true.
    if ( isnan(xHe1) ) has_nan = .true.
    if ( isnan(xHe2) ) has_nan = .true.
    if ( isnan(xHe3) ) has_nan = .true.

  end function check_nan_s

  function check_nan_v( xH1, xH2, xHe1, xHe2, xHe3, tag, nn ) result( has_nan )
    real(real64), dimension(0:nn-1), intent(in) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(in) :: xHe1, xHe2, xHe3
    character(1024), intent(in) :: tag
    integer(int32), intent(in) :: nn
    logical :: has_nan

    has_nan = .false. 
    if ( any(isnan(xH1)) ) has_nan = .true.
    if ( any(isnan(xH2)) ) has_nan = .true.
    if ( any(isnan(xHe1)) ) has_nan = .true.
    if ( any(isnan(xHe2)) ) has_nan = .true.
    if ( any(isnan(xHe3)) ) has_nan = .true.

  end function check_nan_v




  ! Calculates analytic neutral hydrogen fraction.
  !
  ! Input: 
  !    nH:    hydrogen number density
  !    y:     ne = ( xH2 + y ) * nH
  !    reH2:  recombination rates
  !    ciH1:  collisional ionization rates
  !    H1i:   photoionization rates
  !
  ! Output: 
  !    xH1:   hydrogen ionization fractions
  !
  ! Units: 
  !    number densities [1/cm^3]
  !    chemistry rates [cm^3/s]
  !    photoionization rates [1/s]
  !
  !----------------------------------------------------------------------
  subroutine analytic_H1_s( nH, y, reH2, ciH1, H1i, xH1, nn )
    real(real64), intent(in) :: nH 
    real(real64), intent(in) :: y
    real(real64), intent(in) :: reH2
    real(real64), intent(in) :: ciH1
    real(real64), intent(in) :: H1i
    real(real64), intent(out) :: xH1
    integer(int32), intent(in) :: nn

    real(real64) :: RR, QQ, PP, q, dd

    RR = (ciH1 + reH2) * nH
    QQ = -( H1i + reH2 * nH + RR * (one + y) )
    PP = reH2 * nH * (one + y) 
    dd = QQ*QQ - four*RR*PP
    if ( dd < zero ) then
       dd = zero
    end if
    
    ! QQ is always negative in this case
    ! q = -0.5 * (QQ + np.sign(QQ) * np.sqrt(QQ*QQ - 4.0*RR*PP))
    q = -half * (QQ - sqrt(dd))
    xH1 = PP/q
        
  end subroutine analytic_H1_s


  subroutine analytic_H1_v( nH, y, reH2, ciH1, H1i, xH1, nn )
    real(real64), dimension(0:nn-1), intent(in) :: nH 
    real(real64), dimension(0:nn-1), intent(in) :: y
    real(real64), dimension(0:nn-1), intent(in) :: reH2
    real(real64), dimension(0:nn-1), intent(in) :: ciH1
    real(real64), dimension(0:nn-1), intent(in) :: H1i
    real(real64), dimension(0:nn-1), intent(out) :: xH1
    integer(int32), intent(in) :: nn

    real(real64), dimension(0:nn-1) :: RR, QQ, PP, q, dd

    RR = (ciH1 + reH2) * nH
    QQ = -( H1i + reH2 * nH + RR * (one + y) )
    PP = reH2 * nH * (one + y) 
    dd = QQ*QQ - four*RR*PP
    where ( dd < zero )
       dd = zero
    end where

    ! QQ is always negative in this case
    ! q = -0.5 * (QQ + np.sign(QQ) * np.sqrt(QQ*QQ - 4.0*RR*PP))
    q = -half * (QQ - sqrt(dd))
    xH1 = PP/q
        
  end subroutine analytic_H1_v





  ! Calculates ionization fractions given electron number density,  
  ! chemistry rates, and photoionization rates.  
  !
  ! Input: 
  !    ne:                 electron number density
  !    reH2, reHe2, reHe2: recombination rates
  !    ciH1, ciHe1, ciHe2: collisional ionization rates
  !    H1i, He1i, He2i:    photoionization rates
  !
  ! Output: 
  !    xH1, xH2:           hydrogen ionization fractions
  !    xHe1, xHe2, xHe3:   helium ionization fractions
  !
  !
  ! Units: 
  !    number densities [1/cm^3]
  !    chemistry rates [cm^3/s]
  !    photoionization rates [1/s]
  !
  !----------------------------------------------------------------------
  subroutine implicit_ionfracs_s( ne, reH2, reHe2, reHe3, &
       ciH1, ciHe1, ciHe2, H1i, He1i, He2i, &
       xH1, xH2, xHe1, xHe2, xHe3, nn ) 
    real(real64), intent(in) :: ne 
    real(real64), intent(in) :: reH2, reHe2, reHe3 
    real(real64), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), intent(in) :: H1i, He1i, He2i
    real(real64), intent(out) :: xH1, xH2
    real(real64), intent(out) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn
    real(real64) :: R2, CG1, DH
    real(real64) :: R3, CG2, DHe

    ! solve hydrogen
    !----------------------------------------------
    R2 = reH2 * ne
    CG1 = ciH1 * ne + H1i

    DH = R2 + CG1
    xH1 = R2 / DH
    xH2 = CG1 / DH

    ! solve helium
    !----------------------------------------------
    R2 = reHe2 * ne
    R3 = reHe3 * ne
    CG1 = ciHe1 * ne + He1i
    CG2 = ciHe2 * ne + He2i

    DHe = (R2 * R3) + (R3 * CG1) + (CG1 * CG2)
    xHe1 = R2  * R3  / DHe
    xHe2 = R3  * CG1 / DHe
    xHe3 = CG1 * CG2 / DHe

    ! floating point correction
    !----------------------------------------------
    call fpc( xH1, xH2, xHe1, xHe2, xHe3, nn )

  end subroutine implicit_ionfracs_s


  subroutine implicit_ionfracs_v( ne, reH2, reHe2, reHe3, &
       ciH1, ciHe1, ciHe2, H1i, He1i, He2i, &
       xH1, xH2, xHe1, xHe2, xHe3, nn ) 
    real(real64), dimension(0:nn-1), intent(in) :: ne 
    real(real64), dimension(0:nn-1), intent(in) :: reH2, reHe2, reHe3 
    real(real64), dimension(0:nn-1), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), dimension(0:nn-1), intent(in) :: H1i, He1i, He2i
    real(real64), dimension(0:nn-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(out) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn

    real(real64), dimension(0:nn-1) :: R2, CG1, DH
    real(real64), dimension(0:nn-1) :: R3, CG2, DHe


    ! solve hydrogen
    !----------------------------------------------
    R2 = reH2 * ne
    CG1 = ciH1 * ne + H1i

    DH = R2 + CG1
    xH1 = R2 / DH
    xH2 = CG1 / DH

    ! solve helium
    !----------------------------------------------
    R2 = reHe2 * ne
    R3 = reHe3 * ne
    CG1 = ciHe1 * ne + He1i
    CG2 = ciHe2 * ne + He2i

    DHe = (R2 * R3) + (R3 * CG1) + (CG1 * CG2)
    xHe1 = R2  * R3  / DHe
    xHe2 = R3  * CG1 / DHe
    xHe3 = CG1 * CG2 / DHe

    ! floating point correction
    !----------------------------------------------
    call fpc( xH1, xH2, xHe1, xHe2, xHe3, nn )

  end subroutine implicit_ionfracs_v




  ! Calculates collisional ionization equilibrium given chemistry rates.
  !
  ! Input: 
  !    reH2, reHe2, reHe2: recombination rates
  !    ciH1, ciHe1, ciHe2: collisional ionization rates
  !
  ! Output: 
  !    xH1, xH2:           hydrogen ionization fractions
  !    xHe1, xHe2, xHe3:   helium ionization fractions
  !
  !
  ! Units: 
  !    chemistry rates [cm^3/s]
  !
  !----------------------------------------------------------------------
  subroutine solve_ce_s( reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
       xH1, xH2, xHe1, xHe2, xHe3, nn )
    real(real64), intent(in) :: reH2, reHe2, reHe3
    real(real64), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), intent(out) :: xH1, xH2
    real(real64), intent(out) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn
    real(real64) :: DD
    character(1024) :: tag

    ! solve hydrogen
    !-------------------------------------
    DD = reH2 + ciH1
    xH1 = reH2 / DD
    xH2 = ciH1 / DD

    ! solve helium
    !-------------------------------------
    DD = reHe2 * reHe3 + &
         reHe3 * ciHe1 + &
         ciHe1 * ciHe2
    xHe1 = reHe2 * reHe3 / DD
    xHe2 = reHe3 * ciHe1 / DD
    xHe3 = ciHe1 * ciHe2 / DD

    ! floating point correction
    !----------------------------------------------
    call fpc( xH1, xH2, xHe1, xHe2, xHe3, nn )

    tag = 'solve_ce_s'
    if ( check_nan( xH1, xH2, xHe1, xHe2, xHe3, tag, nn ) ) then
       write(*,*) trim(tag)
       stop
    end if

  end subroutine solve_ce_s


  subroutine solve_ce_v( reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
       xH1, xH2, xHe1, xHe2, xHe3, nn )
    real(real64), dimension(0:nn-1), intent(in) :: reH2, reHe2, reHe3
    real(real64), dimension(0:nn-1), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), dimension(0:nn-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(out) :: xHe1, xHe2, xHe3
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: DD
    character(1024) :: tag

    ! solve hydrogen
    !-------------------------------------
    DD = reH2 + ciH1
    xH1 = reH2 / DD
    xH2 = ciH1 / DD

    ! solve helium
    !-------------------------------------
    DD = reHe2 * reHe3 + &
         reHe3 * ciHe1 + &
         ciHe1 * ciHe2
    xHe1 = reHe2 * reHe3 / DD
    xHe2 = reHe3 * ciHe1 / DD
    xHe3 = ciHe1 * ciHe2 / DD

    ! floating point correction
    !----------------------------------------------
    call fpc( xH1, xH2, xHe1, xHe2, xHe3, nn )

    tag = 'solve_ce_v'
    if ( check_nan( xH1, xH2, xHe1, xHe2, xHe3, tag, nn ) ) then
       write(*,*) trim(tag)
       stop
    end if

  end subroutine solve_ce_v




  ! Calculates photo-collisional ionization equilibrium given the number 
  ! density of hydrogen and helium, chemistry rates, and photoionization 
  ! rates.  
  !
  ! Input: 
  !    nH, nHe:            hydrogen/helium number density
  !    reH2, reHe2, reHe2: recombination rates
  !    ciH1, ciHe1, ciHe2: collisional ionization rates
  !    H1i, He1i, He2i:    photoionization rates
  !    tol:                abs(1 - ne_new/ne_old) < tol
  !    nn:                 always 1 for scalar solver
  !
  ! Output: 
  !    xH1, xH2:           hydrogen ionization fractions
  !    xHe1, xHe2, xHe3:   helium ionization fractions
  !
  !
  ! Units: 
  !    number densities [1/cm^3]
  !    chemistry rates [cm^3/s]
  !    photoionization rates [1/s]
  !
  !----------------------------------------------------------------------
  subroutine solve_pce_s( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
       H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn )

    real(real64), intent(in) :: nH, nHe
    real(real64), intent(in) :: reH2, reHe2, reHe3
    real(real64), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), intent(in) :: H1i, He1i, He2i
    real(real64), intent(out) :: xH1, xH2
    real(real64), intent(out) :: xHe1, xHe2, xHe3
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn

    real(real64) :: xH1_old, xHe1_old, xHe3_old
    real(real64) :: ne_min, ne_max
    real(real64) :: ne_left, ne_right
    real(real64) :: ne_old, ne, y
    real(real64) :: err, ratio_ne, ratio_x
    integer(int32) :: iter
    character(1024) :: tag


    ! calculate collisional ionization eq. 
    !-------------------------------------
    call solve_ce( reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         xH1, xH2, xHe1, xHe2, xHe3, nn )

    ! calculate minimum and maximum possible ne 
    !-------------------------------------
    ne_min = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe
    ne_max = nH + two * nHe

    ! initialize brackets 
    !-------------------------------------
    ne_left = ne_min
    ne_right = ne_max

    ! calculate analytic estimate with photoionizations
    !--------------------------------------
    y = ( xHe2 + two * xHe3 ) * nHe / nH
    call analytic_H1( nH, y, reH2, ciH1, H1i, xH1, nn )    
    xH2 = one-xH1
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe

    ! initialize old variables
    !--------------------------------------
    ne_old = ne
    xH1_old = xH1
    xHe1_old = xHe1
    xHe3_old = xHe3

    ! initialize error array and iteration counter
    !--------------------------------------
    iter = 0
    err = huge(one) 

    do while( err > tol )

       ! calculate ionization fractions from ne
       !--------------------------------------
       call implicit_ionfracs( ne, reH2, reHe2, reHe3, &
            ciH1, ciHe1, ciHe2, H1i, He1i, He2i,  &
            xH1, xH2, xHe1, xHe2, xHe3, nn )
       
       ! calculate new ne from ionization fractions 
       !--------------------------------------
       ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe       

       ! calculate ratio and error
       !--------------------------------------       
       ratio_ne = ne / ne_old
       ratio_x = max( xH1/xH1_old, xHe1/xHe1_old, xHe3/xHe3_old ) 
       err = abs( ratio_x - one )

       ! make new guess for ne
       !--------------------------------------
       if ( ratio_ne > one ) then
          ne = ( ne_old + ne_right ) * half
          ne_left = ne_old
       else
          ne = ( ne_left + ne_old ) * half
          ne_right = ne_old
       end if

       ! update old variables
       !--------------------------------------
       ne_old = ne
       xH1_old = xH1
       xHe1_old = xHe1
       xHe3_old = xHe3

       ! increment
       !--------------------------------------
       iter = iter + 1

       ! check max iterations 
       !--------------------------------------
       if ( iter > MAX_ITER ) then
          write(*,*) 'solve_pce'
          write(*,*) 'iter > MAX_ITER'
          write(*,*) 'err: ', err
          write(*,*) 'nH, nHe, ne: ', nH, nHe, ne
          write(*,*) 'i:  ', H1i, He1i, He2i
          write(*,*) 'ci: ', ciH1, ciHe1, ciHe2
          write(*,*) 're: ', reH2, reHe2, reHe3
          write(*,*) 
          stop
       end if

    end do

    if ( nH <= zero ) then
       xH1 = zero
       xH2 = zero
    end if

    if ( nHe <= zero ) then
       xHe1 = zero
       xHe2 = zero
       xHe3 = zero
    end if


  end subroutine solve_pce_s


  subroutine solve_pce_v( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
       H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn  )

    real(real64), dimension(0:nn-1), intent(in) :: nH, nHe
    real(real64), dimension(0:nn-1), intent(in) :: reH2, reHe2, reHe3
    real(real64), dimension(0:nn-1), intent(in) :: ciH1, ciHe1, ciHe2
    real(real64), dimension(0:nn-1), intent(in) :: H1i, He1i, He2i
    real(real64), dimension(0:nn-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(out) :: xHe1, xHe2, xHe3
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn

    real(real64), dimension(0:nn-1) :: xH1_old, xHe1_old, xHe3_old
    real(real64), dimension(0:nn-1) :: ne_min, ne_max
    real(real64), dimension(0:nn-1) :: ne_left, ne_right
    real(real64), dimension(0:nn-1) :: ne_old, ne, y
    real(real64), dimension(0:nn-1) :: err, ratio_ne, ratio_x
    integer(int32) :: iter



    ! calculate collisional ionization eq. 
    !-------------------------------------
    call solve_ce( reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         xH1, xH2, xHe1, xHe2, xHe3, nn )

    ! calculate minimum and maximum possible ne
    !--------------------------------------------
    ne_min = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe
    ne_max = nH + two * nHe

    ! initialize brackets
    !--------------------------------------------
    ne_left = ne_min
    ne_right = ne_max

    ! calculate analytic estimate with photoionizations 
    !--------------------------------------
    y = ( xHe2 + two * xHe3 ) * nHe / nH
    call analytic_H1( nH, y, reH2, ciH1, H1i, xH1, nn )    
    xH2 = one-xH1
    ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe
 
    ! initialize "old" variables
    !--------------------------------------
    ne_old = ne
    xH1_old = xH1
    xHe1_old = xHe1
    xHe3_old = xHe3
    
    ! initialize error array and iteration counter
    !--------------------------------------
    iter = 0
    err = huge(one)

    do while( any( err > tol ) )

       ! calculate ionization fractions from ne
       !--------------------------------------
       call implicit_ionfracs( ne, reH2, reHe2, reHe3, &
            ciH1, ciHe1, ciHe2, H1i, He1i, He2i,  &
            xH1, xH2, xHe1, xHe2, xHe3, nn )
       
       ! calculate new ne from ionization fractions 
       !--------------------------------------
       ne = xH2 * nH + ( xHe2 + two * xHe3 ) * nHe       

       ! calculate ratio and error
       !--------------------------------------
       ratio_ne = ne / ne_old
       ratio_x = max( xH1/xH1_old, xHe1/xHe1_old, xHe3/xHe3_old ) 
       err = abs( ratio_x - one )

       ! make new guesses for ne
       !--------------------------------------
       if ( any( ratio_ne > one ) ) then
          where( ratio_ne > one )
             ne = ( ne_old + ne_right ) * half
             ne_left = ne_old
          end where
       end if

       if ( any( ratio_ne < one ) ) then
          where( ratio_ne < one )
             ne = ( ne_left + ne_old ) * half
             ne_right = ne_old
          end where
       end if

       ! update old variables
       !--------------------------------------
       ne_old = ne
       xH1_old = xH1
       xHe1_old = xHe1
       xHe3_old = xHe3


       ! increment
       !--------------------------------------
       iter = iter + 1

       ! check max iterations 
       !--------------------------------------
       if ( iter > MAX_ITER ) then
          write(*,*) 'solve_pce'
          write(*,*) 'iter > MAX_ITER'
          write(*,*) 'err: ', minval(err), maxval(err)
          write(*,*) 'tol: ', tol
          write(*,*) 'any(err>tol): ', any(err > tol)
          stop
       endif

    end do


    if ( any(nH <= zero) ) then
       where ( nH <= zero )
          xH1 = zero
          xH2 = zero
       end where
    end if

    if ( any(nHe <= zero) ) then
       where ( nHe <= zero )
          xHe1 = zero
          xHe2 = zero
          xHe3 = zero
       end where
    end if


  end subroutine solve_pce_v



  ! Calculates photo-collisional-thermal ionization/temperature equilibrium 
  ! given the number density of hydrogen and helium, photoionization rates 
  ! and photoheating rates.  
  !
  ! Input: 
  !    nH, nHe:            hydrogen/helium number density
  !    H1i, He1i, He2i:    photoionization rates
  !    H1h, He1h, He2h:    photoheating rates
  !    z:                  redshift (for compton cooling)
  !    Hz:                 hubble parameter [1/s] @ z (0 for no hubble cooling)
  !    fcA_H2:             case A fraction for H2 (1.0=caseA, 0.0=caseB)
  !    fcA_He2:            case A fraction for He2 (1.0=caseA, 0.0=caseB)
  !    fcA_He3:            case A fraction for He3 (1.0=caseA, 0.0=caseB)
  !    irate:              rate fit [1 = hui gnedin 97]
  !    tol:                abs(1 - ne_new/ne_old) < tol
  !    nn:                 array size
  !
  ! Output: 
  !    xH1, xH2:           hydrogen ionization fractions
  !    xHe1, xHe2, xHe3:   helium ionization fractions
  !    T:                  equilibrium temperature
  !
  !
  ! Units: 
  !    number densities [1/cm^3]
  !    photoionization rates [1/s]
  !    photoheating rates [erg/s]
  !
  !----------------------------------------------------------------------
  subroutine solve_pcte_s( nH, nHe, H1i, He1i, He2i, H1h, He1h, He2h, z, Hz, &
       fcA_H2, fcA_He2, fcA_He3, irate, xH1, xH2, xHe1, xHe2, xHe3, T, &
       tol, nn )

    real(real64), intent(in) :: nH, nHe 
    real(real64), intent(in) :: H1i, He1i, He2i 
    real(real64), intent(in) :: H1h, He1h, He2h 
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), intent(out) :: xH1, xH2
    real(real64), intent(out) :: xHe1, xHe2, xHe3
    real(real64), intent(out) :: T 
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn


    real(real64) :: Tleft, Tright, Told
    real(real64) :: reH2, reHe2, reHe3
    real(real64) :: ciH1, ciHe1, ciHe2

    real(real64) :: dudT, dudTmin, dudTmax
    real(real64) :: err
    integer(int32) :: iter

    
    ! verify that we can bracket the equilibrium temperature 
    !-------------------------------------
    Tleft = TFLOOR
    Tright = TMAX

    ! bracket the equilibrium temperature.  Teq is bracketed if 
    ! dudT at Tmin is positive and dudT and Tmax is negative
    !
    ! Note that for some density / rates combinations the equilibrium
    ! temperature will be below the temperature floor (and dudT at Tmin
    ! will be negative).  For these entries we simply return the 
    ! equilibrium ionization state at the temperature floor.  This will 
    ! only happen in the case of low temperatures and low photoionization
    ! and/or heating rates. 
    !================================================================

        
    ! dudT for Tmin (=Tleft at the moment)
    !-----------------------------------
    call get_dudT( nH, nHe, Tleft, H1i, He1i, He2i, H1h, He1h, He2h, &
         z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudTmin, tol, nn )


    ! if dudT at Tmin is negative, Teq < Tmin and we simply calculate 
    ! the ionization state at Tmin
    !-----------------------------------
    if ( dudTmin < zero ) then

       call get_kchem( Tleft, fcA_H2, fcA_He2, fcA_He3, irate, &
            reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, nn )

       call solve_pce( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
            H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn )

       T = Tleft

       return

    end if

    ! dudT for Tmax (=Tright at the moment)
    !-----------------------------------
    call get_dudT( nH, nHe, Tright, H1i, He1i, He2i, H1h, He1h, He2h, &
         z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudTmax, tol, nn )


    ! check we have the T bracketed for all values 
    !-----------------------------------
    if ( dudTmin < zero .or. dudTmax > zero ) then
       write(*,*) ' !!!! Teq not bracketed !!!! '
       write(*,*) 'Tleft, dudTmin: ', Tleft, dudTmin
       write(*,*) 'Tright, dudTmax: ', Tright, dudTmax
       stop
    endif


    ! choose an initial guess for T
    !-----------------------------------
    T = ( log10( Tleft ) + log10( Tright ) ) * half
    T = 10.0d0**T

    ! initialize error array and iteration counter
    !--------------------------------------
    err = huge(one)
    iter = 0

    do while( err > tol )

       ! calculate dudT with new T
       !--------------------------------------
       call get_dudT( nH, nHe, T, H1i, He1i, He2i, H1h, He1h, He2h, &
            z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudT, tol, nn )

       ! if dudT is negative, Teq is between Tleft and T
       ! if dudT is positive, Teq is between T and Tright
       !--------------------------------------
       Told = T

       if ( dudT < zero ) then
          Tright = T
          T = (Tleft + T) * half
       else
          Tleft = T
          T = (T + Tright) * half
       end if 

       ! check err against tol
       !--------------------------------------
       err = abs( one - T / Told )

       ! increment
       !--------------------------------------
       iter = iter + 1

       ! check max iterations 
       !--------------------------------------
       if ( iter > MAX_ITER ) then
          write(*,*) 'solve_pcte'
          write(*,*) 'iter > MAX_ITER'
          stop
       endif

    end do

    ! solve for ionization fractions at final temperature
    !-----------------------------------
    call get_kchem( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, nn )

    call solve_pce( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn )

    if ( nH <= zero ) then
       xH1 = zero
       xH2 = zero
    end if

    if ( nHe <= zero ) then
       xHe1 = zero
       xHe2 = zero
       xHe3 = zero
    end if

 
  end subroutine solve_pcte_s



  subroutine solve_pcte_v( nH, nHe, H1i, He1i, He2i, H1h, He1h, He2h, z, Hz, &
       fcA_H2, fcA_He2, fcA_He3, irate, xH1, xH2, xHe1, xHe2, xHe3, T, &
       tol, nn )

    real(real64), dimension(0:nn-1), intent(in) :: nH, nHe 
    real(real64), dimension(0:nn-1), intent(in) :: H1i, He1i, He2i 
    real(real64), dimension(0:nn-1), intent(in) :: H1h, He1h, He2h 
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), dimension(0:nn-1), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), dimension(0:nn-1), intent(out) :: xH1, xH2
    real(real64), dimension(0:nn-1), intent(out) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:nn-1), intent(out) :: T 
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn

    logical, dimension(0:nn-1) :: mask
    logical, dimension(0:nn-1) :: b_floor
    integer(int32) :: n_floor

    real(real64), dimension(0:nn-1) :: Tleft, Tright, Told

    real(real64), dimension(0:nn-1) :: reH2, reHe2, reHe3
    real(real64), dimension(0:nn-1) :: ciH1, ciHe1, ciHe2

    real(real64), dimension(0:nn-1) :: dudT, dudTmin, dudTmax
    real(real64), dimension(0:nn-1) :: err
    integer(int32) :: iter, i
    


    ! verify that we can bracket the equilibrium temperature 
    !-------------------------------------
    Tleft = TFLOOR
    Tright = TMAX

    ! bracket the equilibrium temperature.  Teq is bracketed if 
    ! dudT at Tmin is positive and dudT and Tmax is negative
    !
    ! Note that for some density / rates combinations the equilibrium
    ! temperature will be below the temperature floor (and dudT at Tmin
    ! will be negative).  For these entries we simply return the 
    ! equilibrium ionization state at the temperature floor.  This will 
    ! only happen in the case of low temperatures and low photoionization
    ! and/or heating rates. 
    !================================================================
        

    ! dudT for Tmin (=Tleft at the moment)
    !-----------------------------------
    call get_dudT( nH, nHe, Tleft, H1i, He1i, He2i, H1h, He1h, He2h, &
         z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudTmin, tol, nn )


    ! mark entries for which dudT at Tmin is negative
    !---------------------------------------------------
    b_floor = .false.
    do i = 0, nn-1
       if ( dudTmin(i) < zero ) then
          b_floor(i) = .true.
       end if
    end do
    n_floor = count( b_floor )


    ! dudT for Tmax (=Tright at the moment)
    !-----------------------------------
    call get_dudT( nH, nHe, Tright, H1i, He1i, He2i, H1h, He1h, He2h, &
         z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudTmax, tol, nn )


    ! check we have the Teq bracketed for all non Tfloor entries 
    !-----------------------------------
    mask = (dudTmin < zero) .and. (.not. b_floor)
    if ( any(mask) ) then
       write(*,*) ' !!!! solve_pcte_v: dudTmin < 0 !!!! '
       stop
    end if

    mask = (dudTmax > zero) .and. (.not. b_floor)
    if ( any(mask) ) then
       write(*,*) 'dudTmax: ', dudTmax
       write(*,*) ' !!!! solve_pcte_v: dudTmax > 0 !!!! '
       stop
    end if


    ! choose an initial guess for T
    !-----------------------------------
    where( b_floor )
       T = Tleft
    elsewhere
       T = ( log10( Tleft ) + log10( Tright ) ) * half
       T = 10.0d0**T
    end where


    ! initialize error array and iteration counter
    !--------------------------------------
    err = huge(one)
    iter = 0


    do while( any( err > tol ) )

       ! calculate dudT with new T
       !--------------------------------------
       call get_dudT( nH, nHe, T, H1i, He1i, He2i, H1h, He1h, He2h, &
            z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudT, tol, nn )

       ! if dudT is negative, Teq is between Tleft and T
       ! if dudT is positive, Teq is between T and Tright
       ! leave entries with T=T_floor alone
       !--------------------------------------
       Told = T

       mask = (dudT < zero) .and. (.not. b_floor)
       where( mask )
          Tright = T
          T = (Tleft + T) * half
       end where

       mask = (dudT > zero) .and. (.not. b_floor)
       where ( mask )
          Tleft = T
          T = (T + Tright) * half
       end where

       ! check err against tol
       !--------------------------------------
       err = abs( one - T / Told )

       ! increment
       !--------------------------------------
       iter = iter + 1

       ! check max iterations 
       !--------------------------------------
       if ( iter > MAX_ITER ) then
          write(*,*) 'solve_pcte'
          write(*,*) 'iter > MAX_ITER'
          stop
       endif

    end do

    ! solve for ionization fractions at final temperature
    !-----------------------------------
    call get_kchem( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, nn )

    call solve_pce( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn )


    if ( any(nH <= zero) ) then
       where ( nH <= zero )
          xH1 = zero
          xH2 = zero
       end where
    end if

    if ( any(nHe <= zero) ) then
       where ( nHe <= zero )
          xHe1 = zero
          xHe2 = zero
          xHe3 = zero
       end where
    end if

 
  end subroutine solve_pcte_v




  ! Calculates dudT given the number density of hydrogen and helium, 
  ! temperature, photoionization rates, and photoheating rates. 
  !
  ! Input: 
  !    nH, nHe:            hydrogen/helium number density
  !    T:                  temperature
  !    H1i, He1i, He2i:    photoionization rates
  !    H1h, He1h, He2h:    photoheating rates
  !    z:                  redshift (for compton cooling)
  !    Hz:                 hubble parameter [1/s] @ z (0 for no hubble cooling)
  !    fcA_H2:             case A fraction for H2 (1.0=caseA, 0.0=caseB)
  !    fcA_He2:            case A fraction for He2 (1.0=caseA, 0.0=caseB)
  !    fcA_He3:            case A fraction for He3 (1.0=caseA, 0.0=caseB)
  !    irate:              rate fit [1 = hui gnedin 97]
  !    tol:                abs(1 - ne_new/ne_old) < tol
  !    nn:                 array size
  !
  ! Output: 
  !    dudT:               du/dT
  !
  !
  ! Units: 
  !    number densities [1/cm^3]
  !    temperature [K]
  !    photoionization rates [1/s]
  !    photoheating rates [erg/s]
  !
  !----------------------------------------------------------------------
  subroutine get_dudT_s( nH, nHe, T, H1i, He1i, He2i, H1h, He1h, He2h, &
       z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudT, tol, nn )
    real(real64), intent(in) :: nH, nHe, T
    real(real64), intent(in) :: H1i, He1i, He2i
    real(real64), intent(in) :: H1h, He1h, He2h
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), intent(out) :: dudT
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn


    real(real64) :: reH2, reHe2, reHe3
    real(real64) :: ciH1, ciHe1, ciHe2
    real(real64) :: recH2, recHe2, recHe3
    real(real64) :: cicH1, cicHe1, cicHe2
    real(real64) :: cecH1, cecHe2
    real(real64) :: bremss, compton
    real(real64) :: xH1, xH2
    real(real64) :: xHe1, xHe2, xHe3
    real(real64) :: heat, cool


    call get_kchem( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, nn )

    call get_kcool( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         recH2, recHe2, recHe3, cicH1, cicHe1, cicHe2, &
         cecH1, cecHe2, bremss, compton, nn )


    call solve_pce( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn  )
    
    call get_heat( nH, nHe, H1h, He1h, He2h, xH1, xHe1, xHe2, heat, nn )
    
    call get_cool( nH, nHe, T, recH2, recHe2, recHe3, &
         cicH1, cicHe1, cicHe2, cecH1, cecHe2, bremss, compton, &
         xH1, xH2, xHe1, xHe2, xHe3, z, Hz, cool, nn )
    
    dudT = heat - cool 


  end subroutine get_dudT_s



  subroutine get_dudT_v( nH, nHe, T, H1i, He1i, He2i, H1h, He1h, He2h, &
       z, Hz, fcA_H2, fcA_He2, fcA_He3, irate, dudT, tol, nn )
    real(real64), dimension(0:nn-1), intent(in) :: nH, nHe, T
    real(real64), dimension(0:nn-1), intent(in) :: H1i, He1i, He2i
    real(real64), dimension(0:nn-1), intent(in) :: H1h, He1h, He2h
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), dimension(0:nn-1), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), dimension(0:nn-1), intent(out) :: dudT
    real(real64), intent(in) :: tol
    integer(int32), intent(in) :: nn


    real(real64), dimension(0:nn-1) :: reH2, reHe2, reHe3
    real(real64), dimension(0:nn-1) :: ciH1, ciHe1, ciHe2
    real(real64), dimension(0:nn-1) :: recH2, recHe2, recHe3
    real(real64), dimension(0:nn-1) :: cicH1, cicHe1, cicHe2
    real(real64), dimension(0:nn-1) :: cecH1, cecHe2
    real(real64), dimension(0:nn-1) :: bremss, compton
    real(real64), dimension(0:nn-1) :: xH1, xH2
    real(real64), dimension(0:nn-1) :: xHe1, xHe2, xHe3
    real(real64), dimension(0:nn-1) :: heat, cool


    call get_kchem( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, nn )

    call get_kcool( T, fcA_H2, fcA_He2, fcA_He3, irate, &
         recH2, recHe2, recHe3, cicH1, cicHe1, cicHe2, &
         cecH1, cecHe2, bremss, compton, nn )


    call solve_pce( nH, nHe, reH2, reHe2, reHe3, ciH1, ciHe1, ciHe2, &
         H1i, He1i, He2i, xH1, xH2, xHe1, xHe2, xHe3, tol, nn )
    
    call get_heat( nH, nHe, H1h, He1h, He2h, xH1, xHe1, xHe2, heat, nn )
    
    call get_cool( nH, nHe, T, recH2, recHe2, recHe3, &
         cicH1, cicHe1, cicHe2, cecH1, cecHe2, bremss, compton, &
         xH1, xH2, xHe1, xHe2, xHe3, z, Hz, cool, nn )
    
    dudT = heat - cool 


  end subroutine get_dudT_v





end module ion_solver





