subroutine internal_energy(nx,ny,rho,mom_x,mom_y,etot,e_internal)
  implicit none
  integer :: nx,ny,j,i
  double precision :: rho(nx,ny),mom_x(nx,ny)
  double precision :: mom_y(nx,ny),etot(nx,ny),e_internal(nx,ny),e_kinetic(nx,ny)

  call kinetic_energy(nx,ny,rho,mom_x,mom_y,e_kinetic)
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        e_internal(i,j) = etot(i,j) - e_kinetic(i,j)
        if (e_internal(i,j) .lt. 0d0) then
           e_internal(i,j) = 0d0
        end if
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine internal_energy
