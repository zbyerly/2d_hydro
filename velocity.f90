subroutine velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y,omega_grid,x,y)
  implicit none
  integer :: i,j,nx,ny
  double precision :: rho(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: v_x(nx,ny),v_y(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)


  double precision :: omega_grid, v_grid_x(nx,ny), v_grid_y(nx,ny)

!  do j=1,ny
!     do i=1,nx
!        v_x(i,j) = y(i,j)
!     end do
!  end do
  call getrtheta(nx,ny,x,y,r,theta)


  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx

!        v_x(i,j) = mom_x(i,j)/rho(i,j)
!        v_y(i,j) = mom_y(i,j)/rho(i,j)

        v_grid_x(i,j) = omega_grid*y(i,j)
        v_grid_y(i,j) = -omega_grid*x(i,j)

        v_x(i,j) = mom_x(i,j)/rho(i,j)-v_grid_x(i,j)
        v_y(i,j) = mom_y(i,j)/rho(i,j)-v_grid_y(i,j)


     end do
  end do
  !$OMP END PARALLEL DO

end subroutine velocity
