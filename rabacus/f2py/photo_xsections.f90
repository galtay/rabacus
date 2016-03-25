  !-----------------------------------------------------------
  ! RETURNS PHOTOIONIZATION CROSS-SECTIONS
  !
  !  This is the main file for all photoionization x-sections. 
  !  It wraps calls to the various fitting formula presented by
  !  different groups.
  ! 
  !    ifit = 1 :  Verner 96
  !
  !-----------------------------------------------------------

module photo_xsections
  use types
  use verner_96
  implicit none

  contains


  ! Return ionization energy for H1
  !
  ! Input: 
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    E_H1_th: H1 ionization energy [eV]
  !
  !------------------------------------------------------------------------
  function return_E_H1_th( ifit ) result( E_H1_th )
    integer(int32), intent(in) :: ifit
    real(real64) :: E_H1_th

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       E_H1_th = sum( v96_return_Eth( (/1/), (/1/), 1 ) )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_E_H1_th


  ! Return ionization energy for He1
  !
  ! Input: 
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    E_He1_th: He1 ionization energy [eV]
  !
  !------------------------------------------------------------------------
  function return_E_He1_th( ifit ) result( E_He1_th )
    integer(int32), intent(in) :: ifit
    real(real64) :: E_He1_th

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       E_He1_th = sum( v96_return_Eth( (/2/), (/2/), 1 ) )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_E_He1_th


  ! Return ionization energy for He2
  !
  ! Input: 
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    E_He2_th: He2 ionization energy [eV]
  !
  !------------------------------------------------------------------------
  function return_E_He2_th( ifit ) result( E_He2_th )
    integer(int32), intent(in) :: ifit
    real(real64) :: E_He2_th

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       E_He2_th = sum( v96_return_Eth( (/2/), (/1/), 1 ) )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_E_He2_th







  ! Return photoionization cross-sections for H1 given an array of
  ! energies and a fit type
  !
  ! Input: 
  !    E_arr: energies where we want the cross-section [eV]
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    sigma_H1: H1 X-section at requested energies [cm^2]
  !
  !------------------------------------------------------------------------
  function return_sigma_H1( E_eV, ifit, nn ) result( sigma_H1 )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    integer(int32), intent(in) :: ifit
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: sigma_H1

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       sigma_H1 = v96_return_sigma( 1, 1, E_eV, nn )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_sigma_H1



  ! Return photoionization cross-sections for He1 given an array of
  ! energies and a fit type
  !
  ! Input: 
  !    E_eV: energies where we want the cross-section [eV]
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    sigma_He1: He1 X-section at requested energies [cm^2]
  !
  !------------------------------------------------------------------------
  function return_sigma_He1( E_eV, ifit, nn ) result( sigma_He1 )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    integer(int32), intent(in) :: ifit
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: sigma_He1

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       sigma_He1 = v96_return_sigma( 2, 2, E_eV, nn )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_sigma_He1



  ! Return photoionization cross-sections for He2 given an array of
  ! energies and a fit type
  !
  ! Input: 
  !    E_eV: energies where we want the cross-section [eV]
  !    ifit: fit type [1=verner96]
  !
  ! Output: 
  !    sigma_He2: He2 X-section at requested energies [cm^2]
  !
  !------------------------------------------------------------------------
  function return_sigma_He2( E_eV, ifit, nn ) result( sigma_He2 )
    real(real64), dimension(0:nn-1), intent(in) :: E_eV
    integer(int32), intent(in) :: ifit
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: sigma_He2

    ! For verner 96 fit
    !-----------------------------------------
    if ( ifit == 1 ) then       
       sigma_He2 = v96_return_sigma( 2, 1, E_eV, nn )
    else
       write(*,*) 'fit not recognized' 
    end if

  end function return_sigma_He2





 

end module photo_xsections
