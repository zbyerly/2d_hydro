subroutine flux(nx,ny,mom_x,mom_y,rho,u,u_flux,dir)
  implicit none
  double precision :: v_x(nx,ny),v_y(nx,ny),u_flux(nx,ny),u(nx,ny)
  double precision :: rho(nx,ny),mom_y(nx,ny),mom_x(nx,ny)
  integer :: i,j,nx,ny,dir
  
  call velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y)

  !2/22/2011 2:14pm changed 'ny-1' to 'ny' and 'nx-1' to 'nx'
  
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx

        if (dir .eq. 0) then
           u_flux(i,j)   = v_x(i,j)*u(i,j)
        else 
           u_flux(i,j)   = v_y(i,j)*u(i,j)
        end if
        
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine flux

