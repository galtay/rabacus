module hui_gnedin_97
  use types
  use physical_constants
  implicit none

  real(real64), parameter :: T_H1 = 1.57807d5  !< H1  ion thresh [K]
  real(real64), parameter :: T_He1 = 2.85335d5 !< He1 ion thresh [K]
  real(real64), parameter :: T_He2 = 6.31515d5 !< He2 ion thresh [K]


contains


  !-----------------------------------------------------------
  ! CHEMISTRY
  !-----------------------------------------------------------


!> HII recomb rate case A/B [cm^3/s] - Ferland3_1e9
!-----------------------------------------------------------
  function hg97_reH2a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 1.269d-13
    real(real64), parameter :: t2 = 1.503d0
    real(real64), parameter :: t3 = 0.522d0
    real(real64), parameter :: t4 = 0.470d0
    real(real64), parameter :: t5 = 1.923d0
    lambda = 2.d0 * T_H1 / T
    rate = t1 * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_reH2a

  function hg97_reH2b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 2.753d-14
    real(real64), parameter :: t2 = 1.500d0
    real(real64), parameter :: t3 = 2.740d0
    real(real64), parameter :: t4 = 0.407d0
    real(real64), parameter :: t5 = 2.242d0
    lambda = 2.d0 * T_H1 / T
    rate = t1 * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_reH2b

!> HeII recomb rate case A/B [cm^3/s] - BurgessSeaton5e3_5e5
!-----------------------------------------------------------
  function hg97_reHe2a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    lambda = 2.d0 * T_He1 / T
    rate = 3.0d-14 * lambda**0.654d0 
  end function hg97_reHe2a

  function hg97_reHe2b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    lambda = 2.d0 * T_He1 / T
    rate = 1.26d-14 * lambda**0.750d0 
  end function hg97_reHe2b

!> HeII dielectronic recombination [cm^3 s^-1] - AP3e4_1e6 
!-----------------------------------------------------------
  function hg97_reHe2di(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    lambda = 2.0d0 * T_He2  / T
    rate = 1.90d-3 * T**(-1.5d0) * exp(-0.75d0 * lambda * 0.5d0) * &
         ( 1.0d0 + 0.3d0 * exp(-0.15d0 * lambda * 0.5d0) )
  end function hg97_reHe2di

!> HeIII recomb rate case A/B [cm^3/s] - Ferland3_1e9
!-----------------------------------------------------------
  function hg97_reHe3a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 2.0d0 * 1.269d-13
    real(real64), parameter :: t2 = 1.503d0
    real(real64), parameter :: t3 = 0.522d0
    real(real64), parameter :: t4 = 0.470d0
    real(real64), parameter :: t5 = 1.923d0
    lambda = 2.d0 * T_He2 / T
    rate = t1 * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_reHe3a

  function hg97_reHe3b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 2.0d0 * 2.753d-14
    real(real64), parameter :: t2 = 1.500d0
    real(real64), parameter :: t3 = 2.740d0
    real(real64), parameter :: t4 = 0.407d0
    real(real64), parameter :: t5 = 2.242d0
    lambda = 2.d0 * T_He2 / T
    rate = t1 * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_reHe3b

!> HI collisional ionization rate [cm^3/s] - Lotz1e4_1e9
!-----------------------------------------------------------
  function hg97_ciH1(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 21.11d0
    real(real64), parameter :: t2 = -1.089d0
    real(real64), parameter :: t3 = 0.354d0
    real(real64), parameter :: t4 = 0.874d0
    real(real64), parameter :: t5 = 1.101d0
    lambda = 2.d0 * T_H1 / T
    rate = t1 * T**(-1.5d0) * exp(-lambda*0.5d0) * lambda**(t2) / &
         ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_ciH1

!> HeI collisional ionization rate [cm^3/s] - Lotz1e4_1e9
!-----------------------------------------------------------
  function hg97_ciHe1(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 32.38d0
    real(real64), parameter :: t2 = -1.146d0
    real(real64), parameter :: t3 = 0.416d0
    real(real64), parameter :: t4 = 0.987d0
    real(real64), parameter :: t5 = 1.056d0
    lambda = 2.d0 * T_He1 / T
    rate = t1 * T**(-1.5d0) * exp(-lambda*0.5d0) * lambda**(t2) / &
         ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_ciHe1

!> HeII collisional ionization rate [cm^3/s] - Lotz1e4_1e9
!-----------------------------------------------------------
  function hg97_ciHe2(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 19.95d0
    real(real64), parameter :: t2 = -1.089d0
    real(real64), parameter :: t3 = 0.553d0
    real(real64), parameter :: t4 = 0.735d0
    real(real64), parameter :: t5 = 1.275d0
    lambda = 2.d0 * T_He2 / T
    rate = t1 * T**(-1.5d0) * exp(-lambda*0.5d0) * lambda**(t2) / &
         ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_ciHe2



  !-----------------------------------------------------------
  ! COOLING
  !-----------------------------------------------------------



!> HII recomb cooling rate case A/B [erg cm^3/s] - Ferland3_1e9
!-------------------------------------------------------------
  function hg97_recH2a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 1.778d-29
    real(real64), parameter :: t2 = 1.965d0
    real(real64), parameter :: t3 = 0.541d0
    real(real64), parameter :: t4 = 0.502d0
    real(real64), parameter :: t5 = 2.697d0
    lambda = 2.d0 * T_H1 / T
    rate = t1 * T * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_recH2a


  function hg97_recH2b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 3.435d-30
    real(real64), parameter :: t2 = 1.970d0
    real(real64), parameter :: t3 = 2.250d0
    real(real64), parameter :: t4 = 0.376d0
    real(real64), parameter :: t5 = 3.720d0
    lambda = 2.d0 * T_H1 / T
    rate = t1 * T * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_recH2b



!> HeII recomb cooling rate case A/B [erg cm^3/s] - BurgessSeaton5e3_5e5
!----------------------------------------------------------------------
  function hg97_recHe2a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = kb_erg_K * T * hg97_reHe2a(T,nn)
  end function hg97_recHe2a

  function hg97_recHe2b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = kb_erg_K * T * hg97_reHe2b(T,nn)
  end function hg97_recHe2b


!>  He dielectronic recombination cooling [erg cm^3 s^-1]
!---------------------------------------------------------
  function hg97_recHe2di(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = 0.75d0 * kb_erg_K * T_He2 * hg97_reHe2di(T,nn)
  end function hg97_recHe2di


!> HeIII recomb cooling rate case A/B [erg cm^3/s] - Ferland3_1e9
!---------------------------------------------------------------
  function hg97_recHe3a(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 8.0d0 * 1.778d-29 
    real(real64), parameter :: t2 = 1.965d0
    real(real64), parameter :: t3 = 0.541d0 
    real(real64), parameter :: t4 = 0.502d0
    real(real64), parameter :: t5 = 2.697d0
    lambda = 2.d0 * T_He2 / T
    rate = t1 * T * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_recHe3a


  function hg97_recHe3b(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    real(real64), parameter :: t1 = 8.0d0 * 3.435d-30 
    real(real64), parameter :: t2 = 1.970d0
    real(real64), parameter :: t3 = 2.250d0 
    real(real64), parameter :: t4 = 0.376d0
    real(real64), parameter :: t5 = 3.720d0
    lambda = 2.d0 * T_He2 / T
    rate = t1 * T * lambda**(t2) / ( 1.d0 + (lambda/t3)**t4 )**t5
  end function hg97_recHe3b




!> HI collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
!---------------------------------------------------------------
  function hg97_cicH1(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = kb_erg_K * T_H1 * hg97_ciH1(T,nn)
  end function hg97_cicH1


!> HeI collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
!---------------------------------------------------------------
  function hg97_cicHe1(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = kb_erg_K * T_He1 * hg97_ciHe1(T,nn)
  end function hg97_cicHe1


!> HeII collisional ionization cooloing rate [erg cm^3/s] - Lotz1e4_1e9
!---------------------------------------------------------------
  function hg97_cicHe2(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate 
    rate = kb_erg_K * T_He2 * hg97_ciHe2(T,nn)
  end function hg97_cicHe2



!> HI collisional excitation cooling [erg cm^3 s^-1] - Black5e3_5e5  
!----------------------------------------------------------------------------
  function hg97_cecH1(T,nn)  result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    lambda = 2.0d0 * T_H1  / T
    rate = 7.5d-19  * exp(-0.75d0*lambda*0.5d0) / (1.0d0 + sqrt(T*1.0d-5))
  end function hg97_cecH1


!> HeII collisional excitation cooling [erg cm^3 s^-1] - Black5e3_5e5  
!----------------------------------------------------------------------------
  function hg97_cecHe2(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: lambda, rate
    lambda = 2.0d0 * T_He2  / T
    rate = 5.54d-17 * (1.0d0/T)**(0.397d0) * &
         exp(-0.75d0*lambda*0.5d0) / (1.0d0 + sqrt(T*1.0d-5))
  end function hg97_cecHe2



!> Compton cooling rate
!-----------------------------------------------------------
  function hg97_compton(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: rate
    rate = 5.65d-36
  end function hg97_compton

!> Bremsstrahlung cooling rate
!-----------------------------------------------------------
  function hg97_bremss(T,nn) result(rate)
    real(real64), dimension(0:nn-1) :: T
    integer(int32) :: nn
    real(real64), dimension(0:nn-1) :: gff, rate
    gff = 1.1d0 + 0.34d0 * exp( -(5.5d0 - log10(T))**2.0d0/3.0d0 )
    rate = 1.43d-27 * sqrt(T) * gff
  end function hg97_bremss





end module hui_gnedin_97

