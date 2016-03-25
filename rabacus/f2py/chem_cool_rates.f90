  !-----------------------------------------------------------
  ! RETURNS SETS OF RATES
  !
  !  This is the main file for all chemistry and cooling rates.
  !  It wraps calls to the various fitting formula presented by
  !  different groups.
  ! 
  !    irate = 1 :  Hui Gnedin 97
  !
  !-----------------------------------------------------------

module chem_cool_rates
  use types
  use hui_gnedin_97
  use physical_constants, only: kb_erg_K
  implicit none

  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: ten = 1.0d1

  interface get_kchem
     module procedure get_kchem_s, get_kchem_v
  end interface get_kchem

  interface get_kcool
     module procedure get_kcool_s, get_kcool_v
  end interface get_kcool

  interface get_heat
     module procedure get_heat_s, get_heat_v
  end interface get_heat

  interface get_cool
     module procedure get_cool_s, get_cool_v
  end interface get_cool

contains



  ! Return all chemistry rates.
  ! fcA is used to linearly interpolate between caseA and caseB
  ! (fcA = 1.0 yields case A rates, fcA = 0.0 yields case B rates) 
  !-------------------------------------------------------------------
  subroutine get_kchem_s( &
       Tin, fcA_H2, fcA_He2, fcA_He3, &
       irate, &
       reH2, reHe2, reHe3, &
       ciH1, ciHe1, ciHe2, &
       nn )

    real(real64), intent(in) :: Tin
    real(real64), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), intent(out) :: reH2, reHe2, reHe3
    real(real64), intent(out) :: ciH1, ciHe1, ciHe2
    integer(int32), intent(in) :: nn
    real(real64) :: T(1)
    real(real64) :: reH2a, reH2b, reHe2a, reHe2b, reHe3a, reHe3b

    T = Tin

    if ( irate == 1 ) then

       reH2a = sum( hg97_reH2a(T,nn) )
       reH2b = sum( hg97_reH2b(T,nn) )
       reH2 = ten**( log10(reH2a) * fcA_H2 + log10(reH2b) * (one-fcA_H2) )

       reHe2a = sum( hg97_reHe2a(T,nn) )
       reHe2b = sum( hg97_reHe2b(T,nn) )
       reHe2 = ten**( log10(reHe2a) * fcA_He2 + log10(reHe2b) * (one-fcA_He2) )
       reHe2 = reHe2 + sum( hg97_reHe2di(T,nn) )       
       
       reHe3a = sum( hg97_reHe3a(T,nn) )
       reHe3b = sum( hg97_reHe3b(T,nn) )
       reHe3 = ten**( log10(reHe3a) * fcA_He3 + log10(reHe3b) * (one-fcA_He3) )

       ciH1 = sum( hg97_ciH1(T,nn) )
       ciHe1 = sum( hg97_ciHe1(T,nn) )
       ciHe2 = sum( hg97_ciHe2(T,nn) )

    else

       write(*,*) ' get_kchem' 
       write(*,*) ' irate not recognized ' 
       write(*,*) ' irate = ', irate

    endif


  end subroutine get_kchem_s



  subroutine get_kchem_v( &
       T, fcA_H2, fcA_He2, fcA_He3, &
       irate, &
       reH2, reHe2, reHe3, &
       ciH1, ciHe1, ciHe2, &
       nn )

    real(real64), dimension(0:nn-1), intent(in) :: T
    real(real64), dimension(0:nn-1), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), dimension(0:nn-1), intent(out) :: reH2, reHe2, reHe3
    real(real64), dimension(0:nn-1), intent(out) :: ciH1, ciHe1, ciHe2
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: reH2a, reH2b
    real(real64), dimension(0:nn-1) :: reHe2a, reHe2b
    real(real64), dimension(0:nn-1) :: reHe3a, reHe3b

    if ( irate == 1 ) then

       reH2a = hg97_reH2a(T,nn)
       reH2b = hg97_reH2b(T,nn)
       reH2 = ten**( log10(reH2a) * fcA_H2 + log10(reH2b) * (one-fcA_H2) )

       reHe2a = hg97_reHe2a(T,nn)
       reHe2b = hg97_reHe2b(T,nn)
       reHe2 = ten**( log10(reHe2a) * fcA_He2 + log10(reHe2b) * (one-fcA_He2) )
       reHe2 = reHe2 + hg97_reHe2di(T,nn) 

       reHe3a = hg97_reHe3a(T,nn)
       reHe3b = hg97_reHe3b(T,nn)
       reHe3 = ten**( log10(reHe3a) * fcA_He3 + log10(reHe3b) * (one-fcA_He3) )

       ciH1 = hg97_ciH1(T,nn)
       ciHe1 = hg97_ciHe1(T,nn)
       ciHe2 = hg97_ciHe2(T,nn)

    else

       write(*,*) ' get_kchem' 
       write(*,*) ' irate not recognized ' 
       write(*,*) ' irate = ', irate

    endif


  end subroutine get_kchem_v


  ! Return all cooling rates.
  ! fcA is used to linearly interpolate between caseA and caseB
  ! (fcA = 1.0 yields case A rates, fcA = 0.0 yields case B rates) 
  !-----------------------------------------------------------
  subroutine get_kcool_s( Tin, &
       fcA_H2, fcA_He2, fcA_He3, irate, &
       recH2, recHe2, recHe3, &
       cicH1, cicHe1, cicHe2, &
       cecH1, cecHe2, &
       bremss, compton, nn )

    real(real64), intent(in) :: Tin
    real(real64), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), intent(out) :: recH2, recHe2, recHe3
    real(real64), intent(out) :: cicH1, cicHe1, cicHe2
    real(real64), intent(out) :: cecH1, cecHe2
    real(real64), intent(out) :: bremss, compton
    integer(int32), intent(in) :: nn
    real(real64) :: T(1)
    real(real64) :: recH2a, recH2b, recHe2a, recHe2b, recHe3a, recHe3b

    T = Tin

    if ( irate == 1 ) then

       recH2a = sum( hg97_recH2a(T,nn) )
       recH2b = sum( hg97_recH2b(T,nn) )
       recH2 = ten**( log10(recH2a) * fcA_H2 + &
            log10(recH2b) * (one-fcA_H2) )

       recHe2a = sum( hg97_recHe2a(T,nn) )
       recHe2b = sum( hg97_recHe2b(T,nn) )
       recHe2 = ten**( log10(recHe2a) * fcA_He2 + &
            log10(recHe2b) * (one-fcA_He2) )
       recHe2 = recHe2 + sum( hg97_recHe2di(T,nn) )       

       recHe3a = sum( hg97_recHe3a(T,nn) )
       recHe3b = sum( hg97_recHe3b(T,nn) )
       recHe3 = ten**( log10(recHe3a) * fcA_He3 + &
            log10(recHe3b) * (one-fcA_He3) )
       
       cicH1 = sum( hg97_cicH1(T,nn) )
       cicHe1 = sum( hg97_cicHe1(T,nn) )
       cicHe2 = sum( hg97_cicHe2(T,nn) )
       
       cecH1 = sum( hg97_cecH1(T,nn) )
       cecHe2 = sum( hg97_cecHe2(T,nn) )
       
       bremss = sum( hg97_bremss(T,nn) )
       compton = sum( hg97_compton(T,nn) )

    else

       write(*,*) ' get_kcool' 
       write(*,*) ' irate nor recognized ' 
       write(*,*) ' irate = ', irate

    endif
    

  end subroutine get_kcool_s



  subroutine get_kcool_v( T, fcA_H2, fcA_He2, fcA_He3, irate, &
                        recH2, recHe2, recHe3, &
                        cicH1, cicHe1, cicHe2, &
                        cecH1, cecHe2, &
                        bremss, compton, nn )

    real(real64), dimension(0:nn-1), intent(in) :: T
    real(real64), dimension(0:nn-1), intent(in) :: fcA_H2, fcA_He2, fcA_He3
    integer(int32), intent(in) :: irate
    real(real64), dimension(0:nn-1), intent(out) :: recH2, recHe2, recHe3
    real(real64), dimension(0:nn-1), intent(out) :: cicH1, cicHe1, cicHe2
    real(real64), dimension(0:nn-1), intent(out) :: cecH1, cecHe2
    real(real64), dimension(0:nn-1), intent(out) :: bremss, compton
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: recH2a, recH2b
    real(real64), dimension(0:nn-1) :: recHe2a, recHe2b
    real(real64), dimension(0:nn-1) :: recHe3a, recHe3b

    if ( irate == 1 ) then

       recH2a = hg97_recH2a(T,nn)
       recH2b = hg97_recH2b(T,nn)
       recH2 = ten**( log10(recH2a) * fcA_H2 + &
            log10(recH2b) * (one-fcA_H2) )

       recHe2a = hg97_recHe2a(T,nn)
       recHe2b = hg97_recHe2b(T,nn)
       recHe2 = ten**( log10(recHe2a) * fcA_He2 + &
            log10(recHe2b) * (one-fcA_He2) )
       recHe2 = recHe2 + hg97_recHe2di(T,nn) 

       recHe3a = hg97_recHe3a(T,nn)
       recHe3b = hg97_recHe3b(T,nn)
       recHe3 = ten**( log10(recHe3a) * fcA_He3 + &
            log10(recHe3b) * (one-fcA_He3) )
       
       cicH1 = hg97_cicH1(T,nn)
       cicHe1 = hg97_cicHe1(T,nn)
       cicHe2 = hg97_cicHe2(T,nn)
       
       cecH1 = hg97_cecH1(T,nn)
       cecHe2 = hg97_cecHe2(T,nn)
       
       bremss = hg97_bremss(T,nn)
       compton = hg97_compton(T,nn)

    else

       write(*,*) ' get_kcool' 
       write(*,*) ' irate nor recognized ' 
       write(*,*) ' irate = ', irate

    endif
    

  end subroutine get_kcool_v


  ! Return heating rate
  !-----------------------------------------------------------
  subroutine get_heat_s( nH, nHe, H1h, He1h, He2h, xH1, xHe1, xHe2, heat, nn )
    real(real64), intent(in) :: nH, nHe
    real(real64), intent(in) :: H1h, He1h, He2h
    real(real64), intent(in) :: xH1, xHe1, xHe2
    real(real64), intent(out) :: heat
    integer(int32), intent(in) :: nn
    heat = xH1 * nH * H1h + xHe1 * nHe * He1h + xHe2 * nHe * He2h
  end subroutine get_heat_s

  subroutine get_heat_v( nH, nHe, H1h, He1h, He2h, xH1, xHe1, xHe2, heat, nn )
    real(real64), dimension(0:nn-1), intent(in) :: nH, nHe
    real(real64), dimension(0:nn-1), intent(in) :: H1h, He1h, He2h
    real(real64), dimension(0:nn-1), intent(in) :: xH1, xHe1, xHe2
    real(real64), dimension(0:nn-1), intent(out) :: heat
    integer(int32), intent(in) :: nn
    heat = xH1 * nH * H1h + xHe1 * nHe * He1h + xHe2 * nHe * He2h
  end subroutine get_heat_v
    

  ! Return cooling rate
  !-----------------------------------------------------------
  subroutine get_cool_s( nH, nHe, T, recH2, recHe2, recHe3, &
       cicH1, cicHe1, cicHe2, cecH1, cecHe2, bremss, compton, &
       xH1, xH2, xHe1, xHe2, xHe3, z, Hz, cool, nn )
    real(real64), intent(in) :: nH, nHe, T
    real(real64), intent(in) :: recH2, recHe2, recHe3
    real(real64), intent(in) :: cicH1, cicHe1, cicHe2
    real(real64), intent(in) :: cecH1, cecHe2
    real(real64), intent(in) :: bremss, compton
    real(real64), intent(in) :: xH1, xH2, xHe1, xHe2, xHe3
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), intent(out) :: cool
    integer(int32), intent(in) :: nn
    real(real64) :: ne
    real(real64) :: Tcmb
    real(real64) :: cc

    ! two-particle processes
    !--------------------------------------------------
    ne = xH2 * nH + (xHe2 + 2.0d0 * xHe3) * nHe

    cool = recH2  * ne * xH2  * nH  + &
           recHe2 * ne * xHe2 * nHe + &
           recHe3 * ne * xHe3 * nHe + &
           cicH1  * ne * xH1  * nH  + &
           cicHe1 * ne * xHe1 * nHe + &
           cicHe2 * ne * xHe2 * nHe + &
           cecH1  * ne * xH1  * nH  + &
           cecHe2 * ne * xHe2 * nHe

    ! compton cooling
    !--------------------------------------------------
    Tcmb = 2.725d0 * (1.0d0 + z)
    cool = cool + compton * (1.0d0+z)**4 * (T-Tcmb) * ne

    ! bremsstrahlung cooling
    !--------------------------------------------------
    cool = cool + bremss * ne * ( nH * xH2 + nHe * xHe2 + 4.0d0 * nHe * xHe3 )

    ! hubble cooling
    !--------------------------------------------------
    cool = cool + 3.0d0 * Hz * kb_erg_K * T * (nH + nHe + ne)

  end subroutine get_cool_s



  subroutine get_cool_v( nH, nHe, T, recH2, recHe2, recHe3, &
       cicH1, cicHe1, cicHe2, cecH1, cecHe2, bremss, compton, &
       xH1, xH2, xHe1, xHe2, xHe3, z, Hz, cool, nn )
    real(real64), dimension(0:nn-1), intent(in) :: nH, nHe, T
    real(real64), dimension(0:nn-1), intent(in) :: recH2, recHe2, recHe3
    real(real64), dimension(0:nn-1), intent(in) :: cicH1, cicHe1, cicHe2
    real(real64), dimension(0:nn-1), intent(in) :: cecH1, cecHe2
    real(real64), dimension(0:nn-1), intent(in) :: bremss, compton
    real(real64), dimension(0:nn-1), intent(in) :: xH1, xH2, xHe1, xHe2, xHe3
    real(real64), intent(in) :: z
    real(real64), intent(in) :: Hz
    real(real64), dimension(0:nn-1), intent(out) :: cool
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: ne
    real(real64) :: Tcmb

    ! two-particle processes
    !--------------------------------------------------
    ne = xH2 * nH + (xHe2 + 2.0d0 * xHe3) * nHe

    cool = recH2  * ne * xH2  * nH  + &
           recHe2 * ne * xHe2 * nHe + &
           recHe3 * ne * xHe3 * nHe + &
           cicH1  * ne * xH1  * nH  + &
           cicHe1 * ne * xHe1 * nHe + &
           cicHe2 * ne * xHe2 * nHe + &
           cecH1  * ne * xH1  * nH  + &
           cecHe2 * ne * xHe2 * nHe

    ! compton cooling
    !--------------------------------------------------
    Tcmb = 2.725d0 * (1.0d0 + z)
    cool = cool + compton * (1.0d0+z)**4 * (T-Tcmb) * ne

    ! bremsstrahlung cooling
    !--------------------------------------------------
    cool = cool + bremss * ne * ( nH * xH2 + nHe * xHe2 + 4.0d0 * nHe * xHe3 )

    ! hubble cooling
    !--------------------------------------------------
    cool = cool + 3.0d0 * Hz * kb_erg_K * T * (nH + nHe + ne)

  end subroutine get_cool_v



 

end module chem_cool_rates
