subroutine reconstruct(nx,ny,u,u_left_x,u_right_x,u_left_y,u_right_y)
  implicit none

  double precision :: u(nx,ny)
  double precision :: u_left_x(nx,ny),u_right_x(nx,ny)
  double precision :: u_left_y(nx,ny),u_right_y(nx,ny)
  double precision :: s1,s2,u_slope,minmod
  integer :: i,j,nx,ny

  u_left_x = 1d-9
  u_left_y = 1d-9
  u_right_x = 1d-9
  u_right_y = 1d-9

!$OMP PARALLEL DO PRIVATE(i,j,u_slope,s1,s2)
  do j=2,ny-1
     do i=2,nx-1        

        s1 = u(i+1,j) - u(i,j)
        s2 = u(i,j) - u(i-1,j)
        u_slope = minmod(s1,s2)
        
        u_right_x(i,j) = u(i,j) - 0.5d0*u_slope
        u_left_x(i+1,j) = u(i,j) + 0.5d0*u_slope
        
        s1 = u(i,j+1) - u(i,j)
        s2 = u(i,j) - u(i,j-1)
        u_slope = minmod(s1,s2)
        
        u_right_y(i,j) = u(i,j) - 0.5d0*u_slope
        u_left_y(i,j+1) = u(i,j) + 0.5d0*u_slope
     end do
  end do
!$OMP END PARALLEL DO

end subroutine reconstruct
