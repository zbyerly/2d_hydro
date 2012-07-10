subroutine velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y)
  implicit none
  integer :: i,j,nx,ny
  double precision :: rho(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: v_x(nx,ny),v_y(nx,ny)


  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        v_x(i,j) = mom_x(i,j)/rho(i,j)
        v_y(i,j) = mom_y(i,j)/rho(i,j)
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine velocity
