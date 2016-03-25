module verner_96
  use types
  implicit none
  
  integer(int32), parameter :: N_entries = 185
  integer(int32), dimension(0:N_entries-1) :: Z_table    ! Atomic Number
  integer(int32), dimension(0:N_entries-1) :: N_table    ! Electron Number
  real(real64), dimension(0:N_entries-1) :: Eth_table    ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: Emax_table   ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: E0_table     ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: sigma0_table ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: ya_table     ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: P_table      ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: yw_table     ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: y0_table     ! Fit parameters
  real(real64), dimension(0:N_entries-1) :: y1_table     ! Fit parameters
  

  logical :: table_read = .false. 

  
contains


  ! reads in fit parameters to module variables
  !--------------------------------------------------------------------
  subroutine v96_read_table( ) 
    character(1024) :: fname
    integer(int32) :: lun
    integer(int32) :: i
    character(1024) :: envar

    call get_environment_variable( 'RABACUS_DIR', envar )

    lun = 21
    fname = trim(envar) // '/atomic/verner/photox/photo.dat'

    open( unit=lun, file=fname, action='read' ) 
    do i = 0, N_entries-1
       read(lun,*) Z_table(i), N_table(i), Eth_table(i), Emax_table(i), &
            E0_table(i), sigma0_table(i), ya_table(i), P_table(i), &
            yw_table(i), y0_table(i), y1_table(i)
    end do
    close( unit=lun )

    sigma0_table = sigma0_table * 1.0d-18

  end subroutine v96_read_table



  ! returns ionization energy 
  !
  ! Input: 
  !     Z: Atomic number ( H1 = 1, He1 = 2, He2 = 2 ... )
  !     N: Electron number ( H1 = 1, He1 = 2, He2 = 1 ... )
  !
  ! Output: 
  !     Eth: Ionization energy [eV]
  !--------------------------------------------------------------------
  function v96_return_Eth( Z, N, nn ) result( Eth )
    integer(int32), dimension(0:nn-1), intent(in) :: Z
    integer(int32), dimension(0:nn-1), intent(in) :: N
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: Eth

    integer(int32) :: i, j

    ! read table if we haven't already
    !--------------------------------
    if ( .not. table_read ) call v96_read_table()

    if ( any(Z < 1) .or. any(Z > 26) ) then
       write(*,*) 'all Z must be between 1 and 26'
       write(*,*) 'Z = ', Z
       stop
    end if
    
    if ( any(N < 1) .or. any(N > Z) ) then
       write(*,*) 'all N must be between 1 and Z'
       write(*,*) 'Z = ', Z
       write(*,*) 'N = ', N
       stop
    end if

    do i = 0,nn-1
       do j = 0,N_entries-1         
          if ( Z(i) == Z_table(j) .and. N(i) == N_table(j) ) then
             Eth(i) = Eth_table(j)
          end if          
       end do       
    end do


  end function v96_return_Eth
  

  ! returns a fit to photoionization x-section
  !
  ! Input: 
  !     Z: Atomic number ( H1 = 1, He1 = 2, He2 = 2 ... )
  !     N: Electron number ( H1 = 1, He1 = 2, He2 = 1 ... )
  !     E: Energies where we want the X-section [eV]
  !
  ! Output: 
  !     sigma: Photoionization X-section [cm^2]
  !--------------------------------------------------------------------
  function v96_return_sigma( Z, N, E, nn ) result( sigma )
    integer(int32), intent(in) :: Z
    integer(int32), intent(in) :: N    
    real(real64), dimension(0:nn-1), intent(in) :: E
    integer(int32), intent(in) :: nn
    real(real64), dimension(0:nn-1) :: sigma

    real(real64), dimension(0:nn-1) :: x, y, t1, t2, t3, F
    real(real64) :: Eth

    integer(int32) :: i


    ! read table if we haven't already
    !--------------------------------
    if ( .not. table_read ) call v96_read_table()

    ! find requested ion and construct fit
    !--------------------------------
    do i = 0,N_entries-1

       ! if found correct ion
       !--------------------------------
       if ( Z == Z_table(i) .and. N == N_table(i) ) then

          Eth = Eth_table(i)
          x = E / E0_table(i) - y0_table(i)
          y = sqrt( x*x + y1_table(i)*y1_table(i) )
          t1 = (x-1.0d0)*(x-1.0d0) + yw_table(i)*yw_table(i)
          t2 = y**(0.5d0 * P_table(i) - 5.5d0)
          t3 = (1.0d0 + sqrt(y/ya_table(i)))**(-P_table(i))

          F = t1 * t2 * t3

          sigma = sigma0_table(i) * F

       end if
    end do

    ! zero out sigma below threshold
    !--------------------------------
    where( E < Eth )
       sigma = 0.0d0
    end where

    
  end function v96_return_sigma


  
end module verner_96
