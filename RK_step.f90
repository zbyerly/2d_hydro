subroutine RK_step(nx,ny,dx,dy,kappa,gamma,dt,rho_floor,&
          rho,tau,etot,mom_A,mom_B,mom_x,mom_y,phi,x,y,flux_in,alpha,recons)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny)
  double precision :: mom_A(nx,ny),mom_B(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: phi(nx,ny),x(nx,ny),y(nx,ny)
  double precision :: rho_0(nx,ny),tau_0(nx,ny),etot_0(nx,ny)
  double precision :: mom_A_0(nx,ny),mom_B_0(nx,ny)
  double precision :: mom_x_0(nx,ny),mom_y_0(nx,ny)
  double precision :: beta,r(nx,ny),theta(nx,ny),e_internal(nx,ny)
  double precision :: Floor_plus_tot,flux_out,flux_in
  double precision :: alpha
  integer :: recons
!  print*,nx,ny
 
!  print*,'alpha=',alpha


!  print*,'storing grid'
  call grid_store(nx,ny,&
       rho,tau,etot,mom_A,mom_B,mom_x,mom_y,&
       rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0)
!  print*,'done storing grid'
  


  beta = 1d0
  call RK_substep(nx,ny,dx,dy,kappa,gamma,dt,beta,rho_floor,&
       rho,tau,etot,mom_A,mom_B,mom_x,mom_y,x,y,phi,&
       rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0,flux_out,recons)
  call momentum_sync(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y,alpha)
  
  beta = 0.25d0
  call RK_substep(nx,ny,dx,dy,kappa,gamma,dt,beta,rho_floor,&
       rho,tau,etot,mom_A,mom_B,mom_x,mom_y,x,y,phi,&
       rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0,flux_out,recons)
  call momentum_sync(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y,alpha)

  
  beta = 2d0/3d0
  call RK_substep(nx,ny,dx,dy,kappa,gamma,dt,beta,rho_floor,&
       rho,tau,etot,mom_A,mom_B,mom_x,mom_y,x,y,phi,&
       rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0,flux_out,recons)
  call momentum_sync(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y,alpha)


  Floor_plus_tot = 0d0
  ! density floor
  !  print*,'density floor'
!!!!!  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        if (rho(i,j) .lt. rho_floor) then
           Floor_plus_tot = Floor_plus_tot + (rho_floor-rho(i,j))*dx*dy
           rho(i,j) = rho_floor           
        end if
     end do
  end do
!!!!!  !$OMP END PARALLEL DO
!  print*, 'floor_plus =',Floor_plus_tot

  
  !update tau from etot under certain conditions
  
  call internal_energy(nx,ny,rho,mom_x,mom_y,etot,e_internal)
!!!  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        if (e_internal(i,j) .gt. 1d-1*etot(i,j)) then
           tau(i,j) = (e_internal(i,j))**(1d0/gamma)
        end if
     end do
  end do
!!!  !$OMP END PARALLEL DO

end subroutine RK_step
  
