program test_slab_bgnd
use iso_fortran_env
#ifdef UseOMP
use omp_lib
#endif
use slab_bgnd, only: slab_bgnd_solve
use special_functions, only: table_E1_create
use special_functions, only: table_E2_create
implicit none


integer(int32) :: Nl
real(real64) :: Lslab 
real(real64) :: dl
real(real64), allocatable :: Edges(:), nH(:), nHe(:), T(:)

integer(int32) :: Nnu
real(real64), allocatable :: E_eV(:), Inu(:)

integer(int32) :: i_rec_meth
real(real64) :: fixed_fcA
integer(int32) :: i_photo_fit
integer(int32) :: i_rate_fit
integer(int32) :: i_find_Teq
integer(int32) :: i_thin
real(real64) :: z
real(real64) :: tol

real(real64), allocatable :: xH1(:), xH2(:)
real(real64), allocatable :: xHe1(:), xHe2(:), xHe3(:)
real(real64), allocatable :: H1i_src(:), He1i_src(:), He2i_src(:)
real(real64), allocatable :: H1i_rec(:), He1i_rec(:), He2i_rec(:)
real(real64), allocatable :: H1h_src(:), He1h_src(:), He2h_src(:)
real(real64), allocatable :: H1h_rec(:), He1h_rec(:), He2h_rec(:)

integer(int32) :: i
integer(int32) :: lun
character(1024) :: fname

integer(int32) :: start, finish
integer(int32) :: count_rate
integer(int32) :: nthreads, ithread

! initialize slab
!----------------------------------------------
Nl = 1024
allocate( xH1(0:Nl-1), xH2(0:Nl-1) )
allocate( xHe1(0:Nl-1), xHe2(0:Nl-1), xHe3(0:Nl-1) )
allocate( Edges(0:Nl), nH(0:Nl-1), nHe(0:Nl-1), T(0:Nl-1) )
allocate( H1i_src(0:Nl-1), He1i_src(0:Nl-1), He2i_src(0:Nl-1) )
allocate( H1i_rec(0:Nl-1), He1i_rec(0:Nl-1), He2i_rec(0:Nl-1) )
allocate( H1h_src(0:Nl-1), He1h_src(0:Nl-1), He2h_src(0:Nl-1) )
allocate( H1h_rec(0:Nl-1), He1h_rec(0:Nl-1), He2h_rec(0:Nl-1) )

Lslab = 1.3177850552122683d+22
dl = Lslab / Nl
do i = 0,Nl
   Edges(i) = i * dl
end do
nH = 8.218831761053482d-3
nHe = 6.765635287054012d-4
T = 1.0d4



! read in spectrum
!----------------------------------------------
Nnu = 128
allocate( E_eV(0:Nnu-1), Inu(0:Nnu-1) )

lun = 20
fname = 'spectrum_E_eV.txt'
open( unit=lun, file=fname, action='read' )
do i = 0, Nnu-1
   read(lun,*) E_eV(i) 
end do
fname = 'spectrum_Inu.txt'
open( unit=lun, file=fname, action='read' )
do i = 0, Nnu-1
   read(lun,*) Inu(i) 
end do


! set other input variables
!----------------------------------------------
i_rec_meth = 2
fixed_fcA = 1.0d0
i_photo_fit = 1
i_rate_fit = 1
i_find_Teq = 1
i_thin = 0
z = 3.0d0
tol = 1.0d-5

#ifdef UseOMP
!$omp parallel 
nthreads = omp_get_num_threads()
ithread = omp_get_thread_num()
if (ithread==0) then
   write(*,*) 'working with nthreads = ', nthreads
end if
!$omp end parallel
#endif

call system_clock(start, count_rate)

call slab_bgnd_solve( &
     Edges, nH, nHe, T, &
     E_eV, Inu, &
     i_rec_meth, fixed_fcA, &   
     i_photo_fit, i_rate_fit, i_find_Teq, &
     i_thin, z, tol, Nl, Nnu, &
     xH1, xH2, xHe1, xHe2, xHe3, &
     H1i_src, He1i_src, He2i_src, &
     H1i_rec, He1i_rec, He2i_rec, &
     H1h_src, He1h_src, He2h_src, &
     H1h_rec, He1h_rec, He2h_rec  )

call system_clock(finish)
print '("Time = ",f6.3," seconds.")', real(finish-start)/real(count_rate)

write(*,*) 'xH1(Nl-1): ', xH1(Nl-1)
write(*,*) 'xH1(Nl/2): ', xH1(Nl/2)



end program test_slab_bgnd
