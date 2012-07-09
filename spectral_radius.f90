subroutine spectral_radius(nx,ny,rho,mom_x,mom_y,etot,gamma,tau,sr_x,sr_y)
  implicit none
  include 'variables.h'
  double precision :: r(nx,ny),theta(nx,ny),rho(nx,ny),mom_x(nx,ny)
  double precision :: mom_y(nx,ny),etot(nx,ny)
  double precision :: sr_x(nx,ny),sr_y(nx,ny),x(nx,ny),y(nx,ny)
  double precision :: v_x(nx,ny),v_y(nx,ny),c_s(nx,ny),tau(nx,ny)

  call soundspeed(nx,ny,rho,mom_x,mom_y,etot,gamma,tau,c_s)  
  call velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y)

!  sr_x = c_s!abs(v_x) + c_s
!  sr_y = c_s!abs(v_y) + c_s

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        sr_x(i,j) = abs(v_x(i,j)) + c_s(i,j)
        sr_y(i,j) = abs(v_y(i,j)) + c_s(i,j)
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine spectral_radius
