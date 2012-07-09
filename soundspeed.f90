subroutine soundspeed(nx,ny,rho,mom_x,mom_y,etot,gamma,tau,c_s)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),mom_x(nx,ny)
  double precision :: mom_y(nx,ny),etot(nx,ny)
  double precision :: c_s(nx,ny),pressure(nx,ny),tau(nx,ny)


  call get_pressure(nx,ny,rho,mom_x,mom_y,etot,gamma,tau,pressure)
  
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
!        print*,i,j,rho(i,j)
        c_s(i,j) = sqrt(gamma*pressure(i,j)/rho(i,j))
        if (c_s(i,j) .lt. 0) then
           print*,i,j,c_s(i,j),pressure(i,j),rho(i,j)
           call exit(0)
        end if
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine soundspeed
