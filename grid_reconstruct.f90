subroutine grid_reconstruct(nx,ny,rho_floor,&
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
  implicit none
  integer :: recons
  include 'variables.h'
  include 'grid.h'

  double precision :: mom_A_2(nx,ny),mom_B_2(nx,ny),mom_x_2(nx,ny),mom_y_2(nx,ny)

  call getrtheta(nx,ny,x,y,r,theta)

  !divide conserved quantities by rho
  mom_A_2=mom_A/rho
  mom_B_2=mom_B/rho/r
  mom_x_2=mom_x/rho
  mom_y_2=mom_y/rho

  if (recons .eq. 1) then
!     print*, 'using ppm!'
     call reconstruct_ppm(nx,ny,rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y)
     call reconstruct_ppm(nx,ny,tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y)
     call reconstruct_ppm(nx,ny,etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y)
     call reconstruct_ppm(nx,ny,mom_A_2,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y)
     call reconstruct_ppm(nx,ny,mom_B_2,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y)
     call reconstruct_ppm(nx,ny,mom_x_2,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y)
     call reconstruct_ppm(nx,ny,mom_y_2,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y)
  else if (recons .eq. 0) then
!     print*, 'using minmod!'
     call reconstruct(nx,ny,rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y)
     call reconstruct(nx,ny,tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y)
     call reconstruct(nx,ny,etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y)
     call reconstruct(nx,ny,mom_A_2,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y)
     call reconstruct(nx,ny,mom_B_2,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y)
     call reconstruct(nx,ny,mom_x_2,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y)
     call reconstruct(nx,ny,mom_y_2,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y)
  else
     print*, 'error invalid reconstruction!'
     call exit(1)
  end if


  
  do j=2,ny
     do i=2,nx
        x_left_x(i,j) = ( x(i,j) + x(i-1,j) )/2d0
        x_right_x(i,j) = ( x(i,j) + x(i-1,j) )/2d0
        x_left_y(i,j) = ( x(i,j) + x(i,j-1) )/2d0
        x_right_y(i,j) =  ( x(i,j) + x(i,j-1) )/2d0

        y_left_x(i,j) = ( y(i,j) + y(i-1,j) )/2d0
        y_right_x(i,j) = ( y(i,j) + y(i-1,j) )/2d0
        y_left_y(i,j) = ( y(i,j) + y(i,j-1) )/2d0
        y_right_y(i,j) = ( y(i,j) + y(i,j-1) )/2d0
     end do
  end do

  do i=1,nx
     j=1
     x_left_x(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_right_x(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_left_y(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_right_y(i,j) =  ( x(i,j) + x(i,j) )/2d0
     
     y_left_x(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_right_x(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_left_y(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_right_y(i,j) = ( y(i,j) + y(i,j) )/2d0
  end do

  do j=1,ny
     i=1
     x_left_x(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_right_x(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_left_y(i,j) = ( x(i,j) + x(i,j) )/2d0
     x_right_y(i,j) =  ( x(i,j) + x(i,j) )/2d0
     
     y_left_x(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_right_x(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_left_y(i,j) = ( y(i,j) + y(i,j) )/2d0
     y_right_y(i,j) = ( y(i,j) + y(i,j) )/2d0
  end do

  !multiplying reconstructed values by rho to return to conserved variables

  mom_A_left_x = mom_A_left_x*rho_left_x
  mom_A_right_x = mom_A_right_x*rho_right_x
  mom_A_left_y = mom_A_left_y*rho_left_y
  mom_A_right_y = mom_A_right_y*rho_right_y


  call getrtheta(nx,ny,x_left_x,y_left_x,r_left_x,theta)
  call getrtheta(nx,ny,x_right_x,y_right_x,r_right_x,theta)
  call getrtheta(nx,ny,x_left_y,y_left_y,r_left_y,theta)
  call getrtheta(nx,ny,x_right_y,y_right_y,r_right_y,theta)

  mom_B_left_x  = mom_B_left_x  * rho_left_x  * r_left_x
  mom_B_right_x = mom_B_right_x * rho_right_x * r_right_x
  mom_B_left_y  = mom_B_left_y  * rho_left_y  * r_left_y
  mom_B_right_y = mom_B_right_y * rho_right_y * r_right_y

  mom_x_left_x = mom_x_left_x*rho_left_x
  mom_x_right_x = mom_x_right_x*rho_right_x
  mom_x_left_y = mom_x_left_y*rho_left_y
  mom_x_right_y = mom_x_right_y*rho_right_y

  mom_y_left_x = mom_y_left_x*rho_left_x
  mom_y_right_x = mom_y_right_x*rho_right_x
  mom_y_left_y = mom_y_left_y*rho_left_y
  mom_y_right_y = mom_y_right_y*rho_right_y


!  call reconstruct(nx,ny,x,x_left_x,x_right_x,x_left_y,x_right_y)
!  call reconstruct(nx,ny,y,y_left_x,y_right_x,y_left_y,y_right_y)
!  call reconstruct(nx,ny,phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y)

  
  !density floor
 !!!! !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx
        rho_left_x(i,j)  = max(rho_left_x(i,j),rho_floor)
        rho_left_y(i,j)  = max(rho_left_y(i,j),rho_floor)
        rho_right_x(i,j) = max(rho_right_x(i,j),rho_floor)
        rho_right_y(i,j) = max(rho_right_y(i,j),rho_floor)        
     end do
  end do
!!!!!  !$OMP END PARALLEL DO


end subroutine grid_reconstruct
