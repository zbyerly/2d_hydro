subroutine get_pressure(nx,ny,rho,mom_x,mom_y,etot,gamma,tau,pressure)

  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),mom_x(nx,ny)
  double precision :: mom_y(nx,ny),etot(nx,ny),e_internal(nx,ny)
  double precision :: pressure(nx,ny),tau(nx,ny)

  double precision :: x_min,x_max,y_min,y_max
  include 'params.h'

!  print*,'getting internal energy'
  call internal_energy(nx,ny,rho,mom_x,mom_y,etot,e_internal)

!  print*,'calculating pressure'
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx

        pressure(i,j) = kappa*rho(i,j)**gamma
        
        if ( e_internal(i,j)  .gt. 1e-3*etot(i,j)) then
!           pressure(i,j) = (gamma-1d0)*e_internal(i,j)
        else           
!           print*,i,j,tau(i,j)
!           pressure(i,j) = (gamma-1d0)*(tau(i,j)**gamma)
!           print*,i,j,tau(i,j)
        end if
!        if (pressure(i,j) .lt. 1d-9) then
!           print*,i,j,e_internal(i,j)
!        end if        
        
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine get_pressure
