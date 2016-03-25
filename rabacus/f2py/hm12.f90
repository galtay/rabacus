module hm12
  use types
  use physical_constants, only: c_A_s, h_eV_s
  use utils, only: utils_interp_lin_1d
  implicit none

  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: ten = 1.0d1

  ! photorates table
  !-------------------------------------------------------
  integer(int32), parameter :: pr_Nz = 59
  real(real64), dimension(0:pr_Nz-1) :: pr_z             ! redshift
  real(real64), dimension(0:pr_Nz-1) :: pr_l1pz          ! log10(1+z)
  real(real64), dimension(0:pr_Nz-1) :: pr_H1i           ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_H1h           ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_He1i          ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_He1h          ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_He2i          ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_He2h          ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_Comptonh      ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_H1i       ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_H1h       ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_He1i      ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_He1h      ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_He2i      ! [1/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_He2h      ! [eV/s]
  real(real64), dimension(0:pr_Nz-1) :: pr_log_Comptonh  ! [eV/s]
  logical :: pr_table_read = .false. 

  ! uvb table
  !-------------------------------------------------------
  integer(int32), parameter :: uvb_Nz = 60             ! number of redshifts
  integer(int32), parameter :: uvb_Nnu = 575            ! number of wavelengths
  real(real64), dimension(0:uvb_Nz-1) :: uvb_z         ! redshift 
  real(real64), dimension(0:uvb_Nz-1) :: uvb_l1pz      ! log10(1+z)
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_lam_A     ! wavelength [A]
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_nu_Hz     ! frequency [Hz]
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_E_eV      ! energy [eV]
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_log_lam_A     ! wavelength [A]
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_log_nu_Hz     ! frequency [Hz]
  real(real64), dimension(0:uvb_Nnu-1) :: uvb_log_E_eV      ! energy [eV]
  real(real64), dimension(0:uvb_Nnu-1, 0:uvb_Nz-1) :: uvb_Inu ! Spec. Intensity
  real(real64), dimension(0:uvb_Nnu-1, 0:uvb_Nz-1) :: uvb_log_Inu 
  logical :: uvb_table_read = .false.


contains


  ! reads in UVB.out table
  !--------------------------------------------------------------------
  subroutine hm12_uvb_read() 
    character(1024) :: fname
    integer(int32) :: lun
    integer(int32) :: i
    character(1024) :: envar

    integer(int32), parameter :: Nskip = 20
    logical, dimension(0:uvb_Nnu-1, 0:uvb_Nz-1) :: mask
    real(real64) :: Inu_min

    call get_environment_variable( 'RABACUS_DIR', envar )

    lun = 21
    fname = trim(envar) // '/uv_bgnd/hm12_dat/UVB.out'

    open( unit=lun, file=fname, action='read' ) 

    ! read header
    do i = 0, Nskip-1
       read(lun,*) 
    end do

    ! read redshifts
    read(lun,*) uvb_z

    ! read wavelengths and specific intensities
    do i = 0, uvb_Nnu-1
       read(lun,*) uvb_lam_A(i), uvb_Inu(i,:)
    end do

    close( unit=lun )

    ! reverse ordering of wavelengths so that energy goes from low to high
    uvb_lam_A = uvb_lam_A( uvb_Nnu-1:0:-1 )
    
    ! reverse ordering of specific intensity
    do i = 0, uvb_Nz
       uvb_Inu(:,i) = uvb_Inu( uvb_Nnu-1:0:-1, i )
    end do

    uvb_nu_Hz = c_A_s / uvb_lam_A
    uvb_E_eV = uvb_nu_Hz * h_eV_s

    mask = uvb_Inu > 0.0
    Inu_min = minval( uvb_Inu, mask )
    where ( mask )
       uvb_log_Inu = log10( uvb_Inu )
    elsewhere
       uvb_log_Inu = log10( Inu_min )
    end where

    uvb_l1pz = log10( uvb_z + one )

    uvb_log_lam_A = log10( uvb_lam_A )
    uvb_log_nu_Hz = log10( uvb_nu_Hz )
    uvb_log_E_eV = log10( uvb_E_eV )

  end subroutine hm12_uvb_read


  ! reads in photorates.out table
  !--------------------------------------------------------------------
  subroutine hm12_photorates_read() 
    character(1024) :: fname
    integer(int32) :: lun
    integer(int32) :: i
    character(1024) :: envar

    integer(int32), parameter :: Nskip = 24

    call get_environment_variable( 'RABACUS_DIR', envar )

    lun = 21
    fname = trim(envar) // '/uv_bgnd/hm12_dat/photorates.out'

    open( unit=lun, file=fname, action='read' ) 
    do i = 0, Nskip-1
       read(lun,*) 
    end do
    do i = 0, pr_Nz-1
       read(lun,*) pr_z(i), pr_H1i(i), pr_H1h(i), pr_He1i(i), pr_He1h(i), &
            pr_He2i(i), pr_He2h(i), pr_Comptonh(i)
    end do
    close( unit=lun )

    pr_l1pz = log10( pr_z + one )

    pr_log_H1i = log10( pr_H1i )
    pr_log_H1h = log10( pr_H1h )

    pr_log_He1i = log10( pr_He1i )
    pr_log_He1h = log10( pr_He1h )

    pr_log_He2i = log10( pr_He2i )
    pr_log_He2h = log10( pr_He2h )

    pr_log_Comptonh = log10( pr_Comptonh )

  end subroutine hm12_photorates_read



  ! interpolates values from photorates.out table
  ! we interpolate log(1+z) and log(Q)
  !--------------------------------------------------------------------
  subroutine hm12_photorates_at_z( &
       z, H1i, H1h, He1i, He1h, He2i, He2h, Comptonh )
    real(real64), intent(in) :: z
    real(real64), intent(out) :: H1i, H1h
    real(real64), intent(out) :: He1i, He1h
    real(real64), intent(out) :: He2i, He2h
    real(real64), intent(out) :: Comptonh

    integer(int32) :: iz, iz_lo, iz_hi
    real(real64) :: l1pz, frac, dz, dQ

    ! read table if we haven't already
    !--------------------------------
    if ( .not. pr_table_read ) call hm12_photorates_read()

    ! check input redshift is in bounds
    !--------------------------------
    if ( z < pr_z(0) ) then
       write(*,*) ' z out of bounds low '
       stop
    end if

    if ( z > pr_z(pr_Nz-1) ) then
       write(*,*) ' z out of bounds high '
       stop
    end if

    ! find bracketing indices and fractions
    !--------------------------------
    do iz = 0, pr_Nz-1
       if ( z < pr_z(iz) ) then
          iz_hi = iz
          iz_lo = iz_hi - 1
          exit
       end if
       iz_hi = pr_Nz-1
       iz_lo = iz_hi - 1
    end do

    l1pz = log10( one + z )
    dz = pr_l1pz(iz_hi) - pr_l1pz(iz_lo)
    frac = ( l1pz - pr_l1pz(iz_lo) ) / dz

    dQ = pr_log_H1i(iz_hi) - pr_log_H1i(iz_lo)
    H1i = ten**( pr_log_H1i(iz_lo) + dQ * frac )
    
    dQ = pr_log_H1h(iz_hi) - pr_log_H1h(iz_lo)
    H1h = ten**( pr_log_H1h(iz_lo) + dQ * frac )

    dQ = pr_log_He1i(iz_hi) - pr_log_He1i(iz_lo)
    He1i = ten**( pr_log_He1i(iz_lo) + dQ * frac )

    dQ = pr_log_He1h(iz_hi) - pr_log_He1h(iz_lo)
    He1h = ten**( pr_log_He1h(iz_lo) + dQ * frac )

    dQ = pr_log_He2i(iz_hi) - pr_log_He2i(iz_lo)
    He2i = ten**( pr_log_He2i(iz_lo) + dQ * frac )

    dQ = pr_log_He2h(iz_hi) - pr_log_He2h(iz_lo)
    He2h = ten**( pr_log_He2h(iz_lo) + dQ * frac )

    dQ = pr_log_Comptonh(iz_hi) - pr_log_Comptonh(iz_lo)
    Comptonh = ten**( pr_log_Comptonh(iz_lo) + dQ * frac )

  end subroutine hm12_photorates_at_z



  ! interpolates values from UVB.out table
  !--------------------------------------------------------------------
  subroutine hm12_spectrum_at_z( z, E_eV, Inu, nn )
    real(real64), intent(in) :: z     ! requested redshift
    real(real64), dimension(0:nn-1), intent(in) :: E_eV  ! requested samples
    real(real64), dimension(0:nn-1), intent(out) :: Inu  ! returned Inu 
    integer(int32), intent(in) :: nn

    integer(int32) :: iz, iz_lo, iz_hi, iE, iE_lo, iE_hi, Nclip, i
    real(real64) :: l1pz, z_frac, dz
    real(real64) :: E_min, E_max

    real(real64), dimension(0:nn-1) :: Inu_lo, Inu_hi, log_E_eV, log_Inu
    real(real64), dimension(0:uvb_Nnu-1) :: dQ

    real(real64), dimension(:), allocatable :: dlog_Inu_clip, log_Inu_clip
    real(real64), dimension(:), allocatable :: log_E_eV_clip

    ! read table if we haven't already
    !--------------------------------
    if ( .not. uvb_table_read ) call hm12_uvb_read()

    ! check input redshift is in bounds
    !--------------------------------
    if ( z < uvb_z(0) ) then
       write(*,*) ' z out of bounds low '
       stop
    end if

    if ( z > uvb_z(pr_Nz-1) ) then
       write(*,*) ' z out of bounds high '
       stop
    end if


    ! find bracketing energy indices and fractions
    !--------------------------------
    E_min = minval( E_eV )
    do i = 0, uvb_Nnu-1
       if ( uvb_E_eV(i) > E_min ) then
          iE_lo = i - 1
          exit
       end if
    end do

    E_max = maxval( E_eV )
    do i = uvb_Nnu-1, 1, -1
       if ( uvb_E_eV(i) < E_max ) then
          iE_hi = i + 1
          exit
       end if
    end do


    ! create clipped energy array
    !--------------------------------
    Nclip = iE_hi - iE_lo + 1
    allocate( log_E_eV_clip(0:Nclip-1) ) 
    log_E_eV_clip = uvb_log_E_eV(iE_lo:iE_hi)


    ! find bracketing redshift indices and fractions
    !--------------------------------
    do iz = 0, uvb_Nz-1
       if ( z < uvb_z(iz) ) then
          iz_hi = iz
          iz_lo = iz_hi - 1
          exit
       end if
       iz_hi = uvb_Nz-1
       iz_lo = iz_hi - 1
    end do

    l1pz = log10( one + z )
    dz = uvb_l1pz(iz_hi) - uvb_l1pz(iz_lo)
    z_frac = ( l1pz - uvb_l1pz(iz_lo) ) / dz

    ! create clipped Inu array
    !--------------------------------
    allocate( log_Inu_clip(0:Nclip-1), dlog_Inu_clip(0:Nclip-1) ) 

    dlog_Inu_clip = uvb_log_Inu(iE_lo:iE_hi, iz_hi) - &
                    uvb_log_Inu(iE_lo:iE_hi, iz_lo) 

    log_Inu_clip = uvb_log_Inu(iE_lo:iE_hi, iz_lo) + &
         dlog_Inu_clip * z_frac


    ! interpolate clipped array onto requested E_eV
    !--------------------------------
    log_E_ev = log10( E_eV )
    log_Inu = utils_interp_lin_1d( log_E_eV_clip, log_Inu_clip, log_E_eV, &
         Nclip, nn)
    Inu = 10.0d0**log_Inu


  end subroutine hm12_spectrum_at_z


  
  
end module hm12
