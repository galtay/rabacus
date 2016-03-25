!
! Routines to perform 2-D computational geometry
!
!------------------------------------------------------------------------

module geometry
  use types
  use m_mrgrnk
  use physical_constants, only: pi
  implicit none


  real(real64), parameter :: zero = 0.0d0
  real(real64), parameter :: half = 0.5d0
  real(real64), parameter :: one = 1.0d0
  real(real64), parameter :: two = 2.0d0
  real(real64), parameter :: four = 4.0d0

contains



   
  ! shells are indexed from 0 - Nl starting from the outer most.
  ! layers are indexed from 0 - Nl-1 starting from the outer most. 
  ! segments are labeled from 0 to 2m starting from the left.
  ! rays are characterized by the layer of closest approach.
  ! this assumes that radii are measured from the center and 
  ! that Edges is monotonically decreasing (i.e. outer most first)
  !
  ! Edges: radii of all shells (Nl+1 shells and Nl layers)
  ! i_s_lyr: index of layer containg ray origin. 
  ! mu: cos(theta) for ray where theta is angle relative to radial
  ! j: index of segment to calculate ds for
  ! Nl: number of layers
  !----------------------------------------------------------------
  function m_for_ray( Edges, i_s_lyr, mu, Nl ) result( m )
    real(real64), dimension(0:Nl), intent(in) :: Edges
    integer(int32), intent(in) :: i_s_lyr
    real(real64), intent(in) :: mu
    integer(int32), intent(in) :: Nl
    integer(int32) :: m

    integer(int32) :: i_sh_1
    integer(int32) :: i_sh_2
    real(real64) :: r_ray
    real(real64) :: p_ray
    integer(int32) :: i


    ! calculate radius of ray start point
    !-----------------------------------------------------
    i_sh_1 = i_s_lyr
    i_sh_2 = i_s_lyr + 1
    r_ray = ( Edges(i_sh_1) + Edges(i_sh_2) ) * half

    ! calculate impact parameter of ray
    ! note this doesn't care if mu is negative
    !-----------------------------------------------------
    p_ray = r_ray * sqrt( one - mu*mu )

    ! calculate layer of closest approach
    ! this is the layer that contains the point on the
    ! ray closest to the center of the circle. 
    !-----------------------------------------------------
    do i = 1,Nl
       if ( Edges(i) <= p_ray ) then
          m = i-1
          exit
       end if
    end do

  end function m_for_ray



  ! shells are indexed from 0 - Nl starting from the outer most.
  ! layers are indexed from 0 - Nl-1 starting from the outer most. 
  ! segments are labeled from 0 to 2m starting from the left.
  ! rays are characterized by the layer of closest approach.
  ! this assumes that radii are measured from the center and 
  ! that Edges is monotonically decreasing (i.e. outer most first)
  !
  ! Edges: radii of all shells (Nl+1 shells and Nl layers)
  ! i_s_lyr: index of layer containg ray origin. 
  ! mu: cos(theta) for ray where theta is angle relative to radial
  ! j: index of segment to calculate ds for
  ! Nl: number of layers
  !----------------------------------------------------------------
  function ds_segment( Edges, i_s_lyr, mu, j, Nl ) result( ds )
    real(real64), dimension(0:Nl), intent(in) :: Edges
    integer(int32), intent(in) :: i_s_lyr
    real(real64), intent(in) :: mu
    integer(int32), intent(in) :: j
    integer(int32), intent(in) :: Nl
    real(real64) :: ds

    integer(int32) :: i_sh_1
    integer(int32) :: i_sh_2
    real(real64) :: r_ray     ! radius of ray start
    real(real64) :: p_ray     ! impact parameter of ray

    integer(int32) :: m, i
    integer(int32) :: i1, i2
    real(real64) :: s1
    real(real64) :: s2
    integer(int32) :: n


    ! calculate radius of ray start point
    !-----------------------------------------------------
    i_sh_1 = i_s_lyr
    i_sh_2 = i_s_lyr + 1
    r_ray = ( Edges(i_sh_1) + Edges(i_sh_2) ) * half

    ! calculate impact parameter of ray
    ! note this doesn't care if mu is negative
    !-----------------------------------------------------
    p_ray = r_ray * sqrt( one - mu*mu )

    ! calculate layer of closest approach
    ! this is the layer that contains the point on the
    ! ray closest to the center of the circle. 
    !-----------------------------------------------------
    do i = 1,Nl
       if ( Edges(i) <= p_ray ) then
          m = i-1
          exit
       end if
    end do

    ! calculate the number of segments on the ray (n)
    !-----------------------------------------------------
    n = 2 * m + 1

    ! check input
    !-----------------------------------------------------
    if ( j < 0 .or. j > n-1 ) then
       write(*,*) 'j out of range in ds_segment'
       write(*,*) 'j,n: ', j,n
       stop
    end if

    ! calculate path length
    !-----------------------------------------------------
    if ( j < m ) then
       i1 = j
       i2 = j+1
       s1 = sqrt( Edges(i1)*Edges(i1) - p_ray*p_ray )
       s2 = sqrt( Edges(i2)*Edges(i2) - p_ray*p_ray )
       ds = abs( s2-s1 )

    else if ( j == m ) then
       i1 = j
       s1 = sqrt( Edges(i1)*Edges(i1) - p_ray*p_ray )
       ds = 2*s1

    else if ( j > m ) then
       i1 = m - abs(m-j)
       i2 = i1 + 1
       s1 = sqrt( Edges(i1)*Edges(i1) - p_ray*p_ray )
       s2 = sqrt( Edges(i2)*Edges(i2) - p_ray*p_ray )
       ds = abs( s2-s1 )

    end if


  end function ds_segment




end module geometry
