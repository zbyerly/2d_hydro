subroutine getrtheta(nx,ny,x,y,r,theta)
  implicit none
  integer :: nx,ny,i,j
  double precision :: pi
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)
  double precision :: origin_x,origin_y,adjusted_x,adjusted_y


  !%%%%%%%%%%%%%%%
  ! Changing the origin in this subroutine does NOT work!
  !%%%%%%%%%%%%%%

  origin_x = 0d0!(x(nx/2,1)/2d0)
  origin_y = 0d0!(y(ny/2,1)/2d0)

!  origin_x = 0.6d0
!  origin_y = 0d0


  pi = 4d0*atan(1d0)

  !$OMP PARALLEL DO PRIVATE(i,j,adjusted_x,adjusted_y)
  do j=1,ny
     do i=1,nx
        adjusted_x = x(i,j)-origin_x
        adjusted_y = y(i,j)-origin_y
        r(i,j) = dsqrt(adjusted_x**2+adjusted_y**2)
 !       print*,adjusted_y,adjusted_x,pi
        !        theta(i,j) = -atan2(adjusted_y,-adjusted_x)+pi
        theta(i,j) = datan2(adjusted_y,adjusted_x)
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine getrtheta
