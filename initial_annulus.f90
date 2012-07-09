subroutine initial_annulus(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
       rho,tau,etot,mom_A,mom_B,x,y,phi)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny),mom_A(nx,ny),mom_B(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny),phi(nx,ny)
  double precision :: e_kinetic(nx,ny),v_r(nx,ny),v_theta(nx,ny)
  double precision :: mom_r(nx,ny),mom_theta(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  integer, parameter :: n=100000
  double precision :: rho_1D(n),phi_1D(n),omega_1D(n),einternal_1D(n),dr
  double precision :: x_min,x_max,y_min,y_max,e_internal
  double precision :: slope
  integer :: radial_index

  !read in parameters
  !gamma,kappa,x_min,x_max,y_min,y_max,rho_floor,tau_floor,endtime
  include 'params.h'

  dx = (x_max - x_min)/dble(Nx)
  dy = (y_max - y_min)/dble(Ny)

  do j=1,Ny
     do i=1,Nx

        x(i,j) = (i-0.5d0)*dx+x_min
        y(i,j) = (j-0.5d0)*dy+y_min

     end do
  end do

  call getrtheta(nx,ny,x,y,r,theta)
        
  call radial_initial_annulus(n,dr,rho_1D,phi_1D,omega_1D)

  open(99,file='initial.dat') 

  do j=1,ny
     do i=1,nx

        radial_index = floor(r(i,j)/dr)+1

        phi(i,j) = phi_1D(radial_index)
        rho(i,j) = rho_1D(radial_index)
        
        v_r(i,j) = 0d0
        v_theta(i,j) = omega_1D(radial_index)*r(i,j) 
  

        if ( rho(i,j) .lt. rho_floor) then
           rho(i,j) = rho_floor 
!           v_r(i,j) = 0d0
!           v_theta(i,j) = 0d0 
        end if

        tau(i,j) = rho(i,j)*(kappa/(gamma-1))**(1/gamma)      
      ! v_theta is the linear velocity in the theta direction
        
     end do
  end do

  open(99,file='tau_initial.dat')
  do j=1,ny
     do i=1,nx
        write(99,*) x(i,j),y(i,j),tau(i,j)
!,phi(i,j)-phi(i,ny-j+1)
!x(i,j)-x(i,ny-j+1)
     end do
  end do
  close(99)
!  call exit(0)

  call getrtheta(nx,ny,x,y,r,theta)


  mom_r = v_r*rho
  mom_theta = (v_theta*rho)*r  

  mom_A = mom_r
  mom_B = mom_theta


! defining mom_x,mom_y from mom_A,mom_B 
  call mom_cyl2cart(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y)



  call kinetic_energy(nx,ny,rho,mom_x,mom_y,e_kinetic)

  do j=1,Ny
     do i=1,Nx
        e_internal = tau(i,j)**gamma
        etot(i,j) = e_internal+e_kinetic(i,j)
     end do
  end do

end subroutine initial_annulus
