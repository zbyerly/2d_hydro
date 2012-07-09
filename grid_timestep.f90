subroutine grid_timestep(nx,ny,dx,dy,kappa,gamma,cfl_factor,&
     rho,tau,etot,mom_x,mom_y,dt_cfl)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny),phi(nx,ny),x(nx,ny),y(nx,ny)
  double precision :: sr_x(nx,ny),sr_y(nx,ny)
  double precision :: dt_cfl,dt_x(nx,ny),dt_y(nx,ny)

  dt_cfl = 100d0

!  print*,'calling spectral_radius'
  call spectral_radius(nx,ny,rho,mom_x,mom_y,etot&
       ,gamma,tau,sr_x,sr_y)


  ! what is the # of ghost zones? j=2,ny-1 means 1 ghost zone here?
  do j=2,ny-1
     do i=2,nx-1
        dt_x(i,j) = cfl_factor*dx/sr_x(i,j)
        dt_y(i,j) = cfl_factor*dy/sr_y(i,j)
        if (dt_x(i,j) .lt. dt_cfl) then
           dt_cfl = dt_x(i,j)
        end if
        if (dt_y(i,j) .lt. dt_cfl) then
           dt_cfl = dt_y(i,j)
        end if
!        write(42,*) i,j,dt_x(i,j),dt_y(i,j)
     end do
  end do
!  call exit(0)
  
  

end subroutine grid_timestep
