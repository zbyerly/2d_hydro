subroutine RK_substep(nx,ny,dx,dy,kappa,gamma,dt,beta,rho_floor,&
       rho,tau,etot,mom_A,mom_B,mom_x,mom_y,x,y,phi,&
       rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0,flux_out_tot,recons,omega_grid)
  implicit none
  include 'variables.h'
  include 'grid.h'
  double precision :: beta, r_here, omega_grid
  integer :: recons
  double precision :: Flux_out_tot,floor_plus_tot
!  print*,'succesfully entered RK_substep subroutine, beta =',beta
!  print*,'reconstructing grid...'
  call grid_reconstruct(nx,ny,rho_floor,&
     rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y,&
     tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y,&
     etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y,&
     mom_A,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y,&
     mom_B,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y,&
     mom_x,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y,&
     mom_y,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y,&
     x,x_left_x,x_right_x,x_left_y,x_right_y,&
     y,y_left_x,y_right_x,y_left_y,y_right_y,&
     phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y,recons)
  !print*,'finished reconstructing grid'

!  print*,'grid_fluxes...'
  call grid_fluxes(nx,ny,gamma,dx,dy,Flux_out_tot,&
     rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y,&
     tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y,&
     etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y,&
     mom_A,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y,&
     mom_B,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y,&
     mom_x,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y,&
     mom_y,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y,&
     x,x_left_x,x_right_x,x_left_y,x_right_y,&
     y,y_left_x,y_right_x,y_left_y,y_right_y,&
     phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y,&
     rho_Fx,tau_Fx,etot_Fx,mom_A_Fx,mom_B_Fx,mom_x_Fx,mom_y_Fx,&
     rho_Fy,tau_Fy,etot_Fy,mom_A_Fy,mom_B_Fy,mom_x_Fy,mom_y_Fy,omega_grid)
  !print*,'finished grid fluxes'
  
!  print*,'grid_timedifs...'
  call grid_timedifs(nx,ny,gamma,dx,dy,&
     rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y,&
     tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y,&
     etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y,&
     mom_A,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y,&
     mom_B,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y,&
     mom_x,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y,&
     mom_y,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y,&
     x,x_left_x,x_right_x,x_left_y,x_right_y,&
     y,y_left_x,y_right_x,y_left_y,y_right_y,&
     phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y,&
     rho_Fx,tau_Fx,etot_Fx,mom_A_Fx,mom_B_Fx,mom_x_Fx,mom_y_Fx,&
     rho_Fy,tau_Fy,etot_Fy,mom_A_Fy,mom_B_Fy,mom_x_Fy,mom_y_Fy,&
     delta_rho,delta_tau,delta_etot,delta_mom_A,delta_mom_B,delta_mom_x,delta_mom_y,omega_grid)
  !print*,'finished grid_timedifs'

!  print*,'beta =',beta,'dt =',dt

!  print*,'grid_update...'
  call grid_update(nx,ny,rho,delta_rho,rho_0,beta,dt)
  call grid_update(nx,ny,tau,delta_tau,tau_0,beta,dt)
  call grid_update(nx,ny,etot,delta_etot,etot_0,beta,dt)
  call grid_update(nx,ny,mom_A,delta_mom_A,mom_A_0,beta,dt)
  call grid_update(nx,ny,mom_B,delta_mom_B,mom_B_0,beta,dt)
  call grid_update(nx,ny,mom_x,delta_mom_x,mom_x_0,beta,dt)
  call grid_update(nx,ny,mom_y,delta_mom_y,mom_y_0,beta,dt)
    
!  print*,'boundaries...'
  call boundaries(nx,ny,rho)
  call boundaries(nx,ny,tau)
  call boundaries(nx,ny,etot)
  call boundaries(nx,ny,mom_A)
  call boundaries(nx,ny,mom_B)
  call boundaries(nx,ny,mom_x)
  call boundaries(nx,ny,mom_y)

  

  Flux_out_tot = Flux_out_tot*dt
 
!  if (beta .eq. 1d0) then
!     print*, 'flux_out =',Flux_out_tot
!  end if

  ! floor everything within a certain radius
  do j=1,ny
     do i=1,nx
        r_here = dsqrt(x(i,j)**2d0+y(i,j)**2d0)
        
        if (r_here .lt. 1d-5) then
           rho(i,j) = rho_floor
           mom_A(i,j) = 0d0
           mom_B(i,j) = 0d0
           mom_x(i,j) = 0d0
           mom_y(i,j) = 0d0
           etot(i,j) = rho_floor           
        end if
     end do
  end do

end subroutine RK_substep
