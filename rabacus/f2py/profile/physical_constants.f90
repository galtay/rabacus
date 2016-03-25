module physical_constants
  use types
  implicit none
 
  ! mathematical constants 
  !------------------------------------------------------- 
  real(real64), parameter :: pi = 3.141592653589793d+00
 
  ! physical constants 
  !------------------------------------------------------- 
  real(real64), parameter :: kb_erg_K = 1.380648800000000d-16
  real(real64), parameter :: kb_eV_K = 8.617332573209022d-05
  real(real64), parameter :: c_cm_s = 2.997924580000000d+10
  real(real64), parameter :: c_A_s = 2.997924580000000d+18
  real(real64), parameter :: h_eV_s = 4.135667603369524d-15
  real(real64), parameter :: h_erg_s = 6.626069570000000d-27
 
  ! conversion factors 
  !------------------------------------------------------- 
  real(real64), parameter :: eV_2_erg = 1.602176530000000d-12
  real(real64), parameter :: erg_2_eV = 6.241509479607719d+11
 
end module physical_constants
