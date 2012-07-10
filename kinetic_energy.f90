subroutine kinetic_energy(nx,ny,rho,mom_x,mom_y,e_kinetic,omega_grid)
  implicit none
  include 'variables.h'
  double precision :: v_x(nx,ny),v_y(nx,ny)
  double precision :: rho(nx,ny),mom_x(nx,ny)
  double precision :: mom_y(nx,ny),e_kinetic(nx,ny)
  double precision :: omega_grid, x(nx,ny),y(nx,ny)

  call velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y,omega_grid,x,y)

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        e_kinetic(i,j) = 0.5d0*rho(i,j)*(v_x(i,j)*v_x(i,j)+v_y(i,j)*v_y(i,j))
!        if (e_kinetic(i,j) .lt. 1d-9) then
!           e_kinetic(i,j) = 1d-9
!        end if
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine kinetic_energy
