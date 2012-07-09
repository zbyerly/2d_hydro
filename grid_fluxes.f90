subroutine grid_fluxes(nx,ny,gamma,dx,dy,Flux_out_tot,&
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
     rho_Fy,tau_Fy,etot_Fy,mom_A_Fy,mom_B_Fy,mom_x_Fy,mom_y_Fy)
  implicit none
  include 'variables.h'
  include 'grid.h'
  integer :: dir
  double precision :: pressure(nx,ny),temp(nx,ny),sr_xL(nx,ny),sr_yL(nx,ny)
  double precision :: sr_x(nx,ny),sr_y(nx,ny),a(nx,ny),sr_xR(nx,ny),sr_yR(nx,ny)
  double precision :: rho_fluxl(nx,ny),tau_fluxl(nx,ny),etot_fluxl(nx,ny)
  double precision :: mom_a_fluxl(nx,ny),mom_b_fluxr(nx,ny)
  double precision :: mom_a_fluxr(nx,ny), mom_b_fluxl(nx,ny)
  double precision :: mom_x_fluxl(nx,ny),mom_y_fluxr(nx,ny)
  double precision :: mom_x_fluxr(nx,ny), mom_y_fluxl(nx,ny)
  double precision :: rho_fluxr(nx,ny),tau_fluxr(nx,ny),etot_fluxr(nx,ny)

  double precision :: Flux_out_tot

  !xdir
  dir = 0
  !left_side


  call get_pressure(nx,ny,rho_left_x,mom_x_left_x,&
       mom_y_left_x,etot_left_x,gamma,tau_left_x,pressure)
!  print*,'done getting pressure'

!  print*,'computing fluxes...'
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,rho_left_x  ,rho_fluxL  ,dir)
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,tau_left_x  ,tau_fluxL  ,dir)
  temp = etot_left_x+pressure
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,temp        ,etot_fluxL ,dir)
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,mom_A_left_x,mom_A_fluxL,dir)
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,mom_B_left_x,mom_B_fluxL,dir)

  !cartesian momentum
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,mom_x_left_x,mom_x_fluxL,dir)
  call flux(nx,ny,mom_x_left_x,mom_y_left_x,rho_left_x,mom_y_left_x,mom_y_fluxL,dir)



!  print*,'right_side'
!  print*,'getting pressure...'
  call get_pressure(nx,ny,rho_right_x,mom_x_right_x,&
       mom_y_right_x,etot_right_x,gamma,tau_right_x,pressure)
!  print*,'done getting pressure'


  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,rho_right_x  ,rho_fluxR  ,dir)
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,tau_right_x  ,tau_fluxR  ,dir)
  temp = etot_right_x+pressure
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,temp         ,etot_fluxR ,dir)
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,mom_A_right_x,mom_A_fluxR,dir)
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,mom_B_right_x,mom_B_fluxR,dir)

  !cartesian momentum
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,mom_x_right_x,mom_x_fluxR,dir)
  call flux(nx,ny,mom_x_right_x,mom_y_right_x,rho_right_x,mom_y_right_x,mom_y_fluxR,dir)

  
  !print*,'getting spectral radius...'
  call spectral_radius(nx,ny,rho_left_x,&
       mom_x_left_x,mom_y_left_x,etot_left_x,gamma,tau_left_x,sr_xL,sr_y)
  call spectral_radius(nx,ny,rho_right_x,&
       mom_x_right_x,mom_y_right_x,etot_right_x,gamma,tau_right_x,sr_xR,sr_y)
  !print*,'done getting spectral radius'


  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        a(i,j) = max(sr_xL(i,j),sr_xR(i,j))
  
        rho_Fx(i,j)   = -0.5d0*a(i,j)*(rho_right_x(i,j)-rho_left_x(i,j))
        tau_Fx(i,j)   = -0.5d0*a(i,j)*(tau_right_x(i,j)-tau_left_x(i,j))
        etot_Fx(i,j)  = -0.5d0*a(i,j)*(etot_right_x(i,j)-etot_left_x(i,j))
        mom_A_Fx(i,j) = -0.5d0*a(i,j)*(mom_A_right_x(i,j)-mom_A_left_x(i,j))
        mom_B_Fx(i,j) = -0.5d0*a(i,j)*(mom_B_right_x(i,j)-mom_B_left_x(i,j))
        mom_x_Fx(i,j) = -0.5d0*a(i,j)*(mom_x_right_x(i,j)-mom_x_left_x(i,j))
        mom_y_Fx(i,j) = -0.5d0*a(i,j)*(mom_y_right_x(i,j)-mom_y_left_x(i,j))
        
        rho_Fx(i,j)   = rho_Fx(i,j)   + 0.5d0*(rho_fluxR(i,j)+rho_fluxL(i,j))
        tau_Fx(i,j)   = tau_Fx(i,j)   + 0.5d0*(tau_fluxR(i,j)+tau_fluxL(i,j))
        etot_Fx(i,j)  = etot_Fx(i,j)  + 0.5d0*(etot_fluxR(i,j)+etot_fluxL(i,j))
        mom_A_Fx(i,j) = mom_A_Fx(i,j) + 0.5d0*(mom_A_fluxR(i,j)+mom_A_fluxL(i,j))
        mom_B_Fx(i,j) = mom_B_Fx(i,j) + 0.5d0*(mom_B_fluxR(i,j)+mom_B_fluxL(i,j))
        mom_x_Fx(i,j) = mom_x_Fx(i,j) + 0.5d0*(mom_x_fluxR(i,j)+mom_x_fluxL(i,j))
        mom_y_Fx(i,j) = mom_y_Fx(i,j) + 0.5d0*(mom_y_fluxR(i,j)+mom_y_fluxL(i,j))

     end do
  end do
  !$OMP END PARALLEL DO

!  print*,'done with xdirection'


  !ydir--------------
  dir = 1
  !left_side

  call get_pressure(nx,ny,rho_left_y,mom_x_left_y,&
       mom_y_left_y,etot_left_y,gamma,tau_left_y,pressure)
!  print*,'done getting pressure'

!  print*,'computing fluxes...'
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,rho_left_y  ,rho_fluxL  ,dir)
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,tau_left_y  ,tau_fluxL  ,dir)
  temp = etot_left_y+pressure
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,temp        ,etot_fluxL ,dir)
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,mom_A_left_y,mom_A_fluxL,dir)
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,mom_B_left_y,mom_B_fluxL,dir)

  !cartesian momentum
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,mom_x_left_y,mom_x_fluxL,dir)
  call flux(nx,ny,mom_x_left_y,mom_y_left_y,rho_left_y,mom_y_left_y,mom_y_fluxL,dir)



!  print*,'getting pressure...'
  call get_pressure(nx,ny,rho_right_y,mom_x_right_y,&
       mom_y_right_y,etot_right_y,gamma,tau_right_y,pressure)
!  print*,'done getting pressure'


  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,rho_right_y  ,rho_fluxR  ,dir)
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,tau_right_y  ,tau_fluxR  ,dir)
  temp = etot_right_y+pressure
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,temp         ,etot_fluxR ,dir)
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,mom_A_right_y,mom_A_fluxR,dir)
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,mom_B_right_y,mom_B_fluxR,dir)

  !cartesian momentum
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,mom_x_right_y,mom_x_fluxR,dir)
  call flux(nx,ny,mom_x_right_y,mom_y_right_y,rho_right_y,mom_y_right_y,mom_y_fluxR,dir)

  
!  print*,'done with right side'
  call spectral_radius(nx,ny,rho_left_y,&
       mom_x_left_y,mom_y_left_y,etot_left_y,gamma,tau_left_y,sr_xL,sr_yL)
  call spectral_radius(nx,ny,rho_right_y,&
       mom_x_right_y,mom_y_right_y,etot_right_y,gamma,tau_right_y,sr_xR,sr_yR)
  
!  print*,'done with sr'

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        a(i,j) = max(sr_yL(i,j),sr_yR(i,j))
        
        rho_Fy(i,j)   = -0.5d0*a(i,j)*(rho_right_y(i,j)-rho_left_y(i,j))
        tau_Fy(i,j)   = -0.5d0*a(i,j)*(tau_right_y(i,j)-tau_left_y(i,j))
        etot_Fy(i,j)  = -0.5d0*a(i,j)*(etot_right_y(i,j)-etot_left_y(i,j))
        mom_A_Fy(i,j) = -0.5d0*a(i,j)*(mom_A_right_y(i,j)-mom_A_left_y(i,j))
        mom_B_Fy(i,j) = -0.5d0*a(i,j)*(mom_B_right_y(i,j)-mom_B_left_y(i,j))
        mom_x_Fy(i,j) = -0.5d0*a(i,j)*(mom_x_right_y(i,j)-mom_x_left_y(i,j))
        mom_y_Fy(i,j) = -0.5d0*a(i,j)*(mom_y_right_y(i,j)-mom_y_left_y(i,j))

        rho_Fy(i,j)   = rho_Fy(i,j)   + 0.5d0*(rho_fluxR(i,j)+rho_fluxL(i,j))
        tau_Fy(i,j)   = tau_Fy(i,j)   + 0.5d0*(tau_fluxR(i,j)+tau_fluxL(i,j))
        etot_Fy(i,j)  = etot_Fy(i,j)  + 0.5d0*(etot_fluxR(i,j)+etot_fluxL(i,j))
        mom_A_Fy(i,j) = mom_A_Fy(i,j) + 0.5d0*(mom_A_fluxR(i,j)+mom_A_fluxL(i,j))
        mom_B_Fy(i,j) = mom_B_Fy(i,j) + 0.5d0*(mom_B_fluxR(i,j)+mom_B_fluxL(i,j))
        mom_x_Fy(i,j) = mom_x_Fy(i,j) + 0.5d0*(mom_x_fluxR(i,j)+mom_x_fluxL(i,j))
        mom_y_Fy(i,j) = mom_y_Fy(i,j) + 0.5d0*(mom_y_fluxR(i,j)+mom_y_fluxL(i,j))

     end do
  end do
  !$OMP END PARALLEL DO

 ! print*,'done with gridfluxes'

  !diagnostics for the fluxes out of the edge of the grid
!!$

  Flux_out_tot = 0d0
  !left/right
  do j=3,ny-2
     !left
     i=3
     Flux_out_tot = Flux_out_tot - rho_Fx(i,j)*dy

     rho_Fx(i,j) = 0d0
     tau_Fx(i,j) = 0d0
     etot_Fx(i,j) = 0d0
     mom_A_Fx(i,j) = 0d0
     mom_B_Fx(i,j) = 0d0
     mom_x_Fx(i,j) = 0d0
     mom_y_Fx(i,j) = 0d0

     !right
     i=nx-1
     Flux_out_tot = Flux_out_tot + rho_Fx(i,j)*dy

     rho_Fx(i,j) = 0d0
     tau_Fx(i,j) = 0d0
     etot_Fx(i,j) = 0d0
     mom_A_Fx(i,j) = 0d0
     mom_B_Fx(i,j) = 0d0
     mom_x_Fx(i,j) = 0d0
     mom_y_Fx(i,j) = 0d0

  end do
  
  !top/bottom
  do i=3,nx-2
     !bottom
     j=3
     Flux_out_tot = Flux_out_tot - rho_Fy(i,j)*dx

     rho_Fy(i,j) = 0d0
     tau_Fy(i,j) = 0d0
     etot_Fy(i,j) = 0d0
     mom_A_Fy(i,j) = 0d0
     mom_B_Fy(i,j) = 0d0
     mom_x_Fy(i,j) = 0d0
     mom_y_Fy(i,j) = 0d0

     !top
     j=ny-1
     Flux_out_tot = Flux_out_tot + rho_Fy(i,j)*dx

     rho_Fy(i,j) = 0d0
     tau_Fy(i,j) = 0d0
     etot_Fy(i,j) = 0d0
     mom_A_Fy(i,j) = 0d0
     mom_B_Fy(i,j) = 0d0
     mom_x_Fy(i,j) = 0d0
     mom_y_Fy(i,j) = 0d0
      end do

!  Flux_out_tot = Flux_out_tot
!  print*,Flux_out_tot
  
!  print*,'done with gridfluxes'
end subroutine grid_fluxes
