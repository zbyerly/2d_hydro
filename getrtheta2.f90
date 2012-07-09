subroutine getrtheta2(nx,ny,x,y,r,theta)
  implicit none
  integer :: nx,ny,i,j
  double precision, parameter :: pi = 3.14195d0
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)
  double precision :: origin_x,origin_y,adjusted_x,adjusted_y

  origin_x = 0d0
  origin_y = 0d0

  !$OMP PARALLEL DO PRIVATE(i,j,adjusted_x,adjusted_y)
  do j=1,ny
     do i=1,nx
        adjusted_x = x(i,j)-origin_x
        adjusted_y = y(i,j)-origin_y
        r(i,j) = sqrt(adjusted_x**2+adjusted_y**2)
 !       print*,adjusted_y,adjusted_x,pi
        theta(i,j) = -atan2(adjusted_y,-adjusted_x)+pi
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine getrtheta2
