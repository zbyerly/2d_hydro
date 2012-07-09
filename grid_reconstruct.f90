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
     phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y)
  implicit none
  include 'variables.h'
  include 'grid.h'

  !divide conserved quantities by rho
!  tau = tau/rho
!  etot = etot/rho
!  mom_A = mom_A/rho
!  mom_B = mom_B/rho

  call reconstruct(nx,ny,rho,rho_left_x,rho_right_x,rho_left_y,rho_right_y)
  call reconstruct(nx,ny,tau,tau_left_x,tau_right_x,tau_left_y,tau_right_y)
  call reconstruct(nx,ny,etot,etot_left_x,etot_right_x,etot_left_y,etot_right_y)
  call reconstruct(nx,ny,mom_A,mom_A_left_x,mom_A_right_x,mom_A_left_y,mom_A_right_y)
  call reconstruct(nx,ny,mom_B,mom_B_left_x,mom_B_right_x,mom_B_left_y,mom_B_right_y)
  call reconstruct(nx,ny,mom_x,mom_x_left_x,mom_x_right_x,mom_x_left_y,mom_x_right_y)
  call reconstruct(nx,ny,mom_y,mom_y_left_x,mom_y_right_x,mom_y_left_y,mom_y_right_y)

  !multiplying reconstructed values by rho to return to conserved variables
!!$  tau_left_x = tau_left_x*rho_left_x
!!$  tau_right_x = tau_right_x*rho_right_x
!!$  tau_left_y = tau_left_y*rho_left_y
!!$  tau_right_y = tau_right_y*rho_right_y
!!$
!!$  etot_left_x = etot_left_x*rho_left_x
!!$  etot_right_x = etot_right_x*rho_right_x
!!$  etot_left_y = etot_left_y*rho_left_y
!!$  etot_right_y = etot_right_y*rho_right_y
!!$
!!$  mom_A_left_x = mom_A_left_x*rho_left_x
!!$  mom_A_right_x = mom_A_right_x*rho_right_x
!!$  mom_A_left_y = mom_A_left_y*rho_left_y
!!$  mom_A_right_y = mom_A_right_y*rho_right_y
!!$
!!$  mom_B_left_x = mom_B_left_x*rho_left_x
!!$  mom_B_right_x = mom_B_right_x*rho_right_x
!!$  mom_B_left_y = mom_B_left_y*rho_left_y
!!$  mom_B_right_y = mom_B_right_y*rho_right_y
!!$  


!  call reconstruct(nx,ny,x,x_left_x,x_right_x,x_left_y,x_right_y)
!  call reconstruct(nx,ny,y,y_left_x,y_right_x,y_left_y,y_right_y)
!  call reconstruct(nx,ny,phi,phi_left_x,phi_right_x,phi_left_y,phi_right_y)

  
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
