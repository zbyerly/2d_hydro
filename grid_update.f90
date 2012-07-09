subroutine grid_update(nx,ny,u,delta_u,u_0,beta,dt)
  implicit none
  integer :: nx,ny,i,j
  double precision :: u(nx,ny),delta_u(nx,ny),u_0(nx,ny),beta,dt

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx

        u(i,j) = beta*(u(i,j) + delta_u(i,j)*dt)
        u(i,j) = u(i,j) + (1d0-beta)*u_0(i,j)

     end do
  end do

end subroutine grid_update
