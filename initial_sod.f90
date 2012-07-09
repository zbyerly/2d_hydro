subroutine initial_sod(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
       rho,tau,etot,mom_A,mom_B,x,y,phi)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny),mom_A(nx,ny),mom_B(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny),phi(nx,ny)
  double precision :: e_kinetic(nx,ny),z(nx,ny)
!  integer, parameter :: n=10000
!  double precision :: rho_1D(n),phi_1D(n),dr
  double precision :: x_min,x_max,y_min,y_max,e_internal(nx,ny)
!  integer :: radial_index

  !read in parameters
  !gamma,kappa,x_min,x_max,y_min,y_max,rho_floor,mom_geom,tau_floor,endtime
  include 'params.h'
  
  dx = (x_max - x_min)/dble(Nx)
  dy = (y_max - y_min)/dble(Ny)  
    
  gamma = 1.4d0

  do j=1,Ny
     do i=1,Nx

        x(i,j) = (i-0.5d0)*dx+x_min
        y(i,j) = (j-0.5d0)*dy+y_min        
        
        z(i,j) = -x(i,j)-y(i,j)

        if (y(i,j) .lt. -1d-1) then
!        if (x(i,j) .lt. 6.9d0) then    !off axis
           rho(i,j) = 1d0
           e_internal(i,j) = 2.5d0
        else
           rho(i,j) = 0.125d0
           e_internal(i,j) = 0.25d0
        end if

        if ( rho(i,j) .lt. rho_floor) then
           rho(i,j) = rho_floor 
        end if

!        tau(i,j) = rho(i,j)*(kappa/(gamma-1))**(1/gamma)      
        
        
        mom_A(i,j) = 0d0
        mom_B(i,j) = 0d0  

     end do
  end do

  
  call mom_cyl2cart(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y)

  call kinetic_energy(nx,ny,rho,mom_x,mom_y,e_kinetic)

  do j=1,Ny
     do i=1,Nx
!        e_internal = tau(i,j)**gamma
        etot(i,j) = e_internal(i,j)+e_kinetic(i,j)
        tau(i,j) = (e_internal(i,j))**(1d0/gamma)         
     end do
  end do

end subroutine initial_sod
