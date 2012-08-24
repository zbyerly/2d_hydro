subroutine grid_timedifs(nx,ny,gamma,dx,dy,&
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

  implicit none
  include 'variables.h'
  include 'grid.h'
  
  double precision :: v_x(nx,ny),v_y(nx,ny)
  double precision :: pressure_left(nx,ny),pressure_right(nx,ny)

  double precision :: delta_mom_X1(nx,ny),delta_mom_Y1(nx,ny)
  double precision :: omega_grid

  double precision :: temp1(nx,ny),temp2(nx,ny),temp3(nx,ny),temp4(nx,ny)

  delta_rho = 0d0
  delta_tau = 0d0
  delta_etot = 0d0
  delta_mom_A = 0d0
  delta_mom_B = 0d0
  delta_mom_X = 0d0
  delta_mom_Y = 0d0
  delta_mom_X1 = 0d0
  delta_mom_Y1 = 0d0


  !Advection
  ! 2/22/2011 2:30pm changed 'ny-1' to 'ny' and for x
  ! for all the loops in here. should be 'j=3,ny-2' perhaps?
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=2,ny-1
     do i=2,nx-1
        !xdir
        delta_rho(i,j) = delta_rho(i,j)-(rho_Fx(i+1,j)-rho_Fx(i,j))/dx
        delta_tau(i,j) = delta_tau(i,j)-(tau_Fx(i+1,j)-tau_Fx(i,j))/dx
        delta_etot(i,j) = delta_etot(i,j)-(etot_Fx(i+1,j)-etot_Fx(i,j))/dx
        delta_mom_A(i,j) = delta_mom_A(i,j)-(mom_A_Fx(i+1,j)-mom_A_Fx(i,j))/dx
        delta_mom_B(i,j) = delta_mom_B(i,j)-(mom_B_Fx(i+1,j)-mom_B_Fx(i,j))/dx
        delta_mom_x(i,j) = delta_mom_x(i,j)-(mom_x_Fx(i+1,j)-mom_x_Fx(i,j))/dx
        delta_mom_y(i,j) = delta_mom_y(i,j)-(mom_y_Fx(i+1,j)-mom_y_Fx(i,j))/dx
        !ydir
        delta_rho(i,j) = delta_rho(i,j)-(rho_Fy(i,j+1)-rho_Fy(i,j))/dy
        delta_tau(i,j) = delta_tau(i,j)-(tau_Fy(i,j+1)-tau_Fy(i,j))/dy
        delta_etot(i,j) = delta_etot(i,j)-(etot_Fy(i,j+1)-etot_Fy(i,j))/dy
        delta_mom_A(i,j) = delta_mom_A(i,j)-(mom_A_Fy(i,j+1)-mom_A_Fy(i,j))/dy
        delta_mom_B(i,j) = delta_mom_B(i,j)-(mom_B_Fy(i,j+1)-mom_B_Fy(i,j))/dy     
        delta_mom_x(i,j) = delta_mom_x(i,j)-(mom_x_Fy(i,j+1)-mom_x_Fy(i,j))/dy
        delta_mom_y(i,j) = delta_mom_y(i,j)-(mom_y_Fy(i,j+1)-mom_y_Fy(i,j))/dy     
     end do
  end do
  !$OMP END PARALLEL DO


  call velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y,omega_grid,x,y)

!!  !$OMP PARALLEL DO PRIVATE(i,j)
!  do j=1,ny
!     do i=1,nx        
!        mom_x(i,j) = rho(i,j)*v_x(i,j)
!        mom_y(i,j) = rho(i,j)*v_y(i,j)
!     end do
!  end do
!!  !$OMP END PARALLEL DO

  ! pressure on the left will be at x_right_x(i,j)
  ! pressure on the right will be at x_left_x(i+1,j)
  call get_pressure(nx,ny,rho_right_x,&
       mom_x_right_x,mom_y_right_x,etot_right_x,gamma,tau_right_x,pressure_right)
  call get_pressure(nx,ny,rho_left_x,&
       mom_x_left_x,mom_y_left_x,etot_left_x,gamma,tau_left_x,pressure_left)

  !setting pressure to zero at boundaries
  !left/right
  do j=3,ny-2
     !left
     i=3

     pressure_left(i,j) = 0d0
     pressure_right(i,j) = 0d0

     !right
     i=nx-1

     pressure_left(i,j) = 0d0
     pressure_right(i,j) = 0d0

  end do


  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=2,ny-1
     do i=2,nx-1
        !source
        delta_mom_X1(i,j) = delta_mom_X1(i,j) + (0.5d0/dx)*( pressure_left(i,j)+pressure_right(i,j)-pressure_right(i+1,j)-pressure_left(i+1,j) )
        delta_mom_X1(i,j) = delta_mom_X1(i,j) - (0.5d0/dx)*(rho(i,j)*( phi(i+1,j) - phi(i-1,j)))

!+(0.5d0/dx)*(&
!             ( pressure_left(i,j)+pressure_right(i,j)-pressure_right(i+1,j)-pressure_left(i+1,j) )&
!             -rho(i,j)*( phi(i+1,j) - phi(i-1,j) ) &
!             )                
!!$        if ( (i .eq. 60) .and. (j .eq. 60) ) then
!!$           print*,' '
!!$           print*, (0.5d0/dx)*( pressure_left(i,j)+pressure_right(i,j)-pressure_right(i+1,j)-pressure_left(i+1,j) )
!!$           print*, (0.5d0/dx)*(rho(i,j)*( phi(i+1,j) - phi(i-1,j)))
!!$           print*, delta_mom_X(i,j)
!!$        end if
           delta_etot(i,j) = delta_etot(i,j) - (mom_x(i,j)/dx)*&
             0.5d0*(phi(i+1,j)-phi(i-1,j))
     end do
  end do
  !$OMP END PARALLEL DO

  ! pressure on the left will be at x_right_x(i,j)

  ! pressure on the right will be at x_left_x(i+1,j)


  call get_pressure(nx,ny,rho_right_y,&
       mom_x_right_y,mom_y_right_y,etot_right_y,gamma,tau_right_y,pressure_right)
  call get_pressure(nx,ny,rho_left_y,&
       mom_x_left_y,mom_y_left_y,etot_left_y,gamma,tau_left_y,pressure_left)

  !setting pressure to zero at boundaries
  do i=3,nx-2
     !bottom
     j=3

     pressure_left(i,j) = 0d0
     pressure_right(i,j) = 0d0

     !top
     j=ny-1

     pressure_left(i,j) = 0d0
     pressure_right(i,j) = 0d0

  end do


  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=2,ny-1
     do i=2,nx-1
        !source
        delta_mom_Y1(i,j) = delta_mom_Y1(i,j) + (0.5d0/dy)*( pressure_left(i,j)+pressure_right(i,j)-pressure_right(i,j+1)-pressure_left(i,j+1) )
        delta_mom_Y1(i,j) = delta_mom_Y1(i,j) - (0.5d0/dy)*(rho(i,j)*( phi(i,j+1) - phi(i,j-1)))


!        delta_mom_Y(i,j) = delta_mom_Y(i,j)+(0.5d0/dy)*(&
!             ( pressure_left(i,j)+pressure_right(i,j)-pressure_right(i,j+1)-pressure_left(i,j+1) ) -&
!             rho(i,j)*(phi(i,j+1) - phi(i,j-1)) &
!             )

        delta_etot(i,j) = delta_etot(i,j) - (mom_y(i,j)/dy)*&
             0.5d0*(phi(i,j+1)-phi(i,j-1))
     end do
  end do
  !$OMP END PARALLEL DO

  call getrtheta(nx,ny,x,y,r,theta)

!  temp1 = -0.5d0*r*omega_grid**2d0
!  temp2 = 0d0

!  call mom_cyl2cart(nx,ny,temp1,temp2,temp3,temp4,x,y)

!  do j=1,ny
!     do i=1,nx
!        write(98,*) i, j, temp1(i,j),temp4(i,j)
!     end do
!  end do
!  call exit(0)

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,ny
     do i=1,nx

           delta_mom_X(i,j) = delta_mom_X(i,j) + delta_mom_X1(i,j)
           delta_mom_Y(i,j) = delta_mom_Y(i,j) + delta_mom_Y1(i,j)
           delta_mom_A(i,j) = delta_mom_A(i,j) + delta_mom_X1(i,j)*&
                (x(i,j)/r(i,j)) + delta_mom_Y1(i,j)*(y(i,j)/r(i,j))
           delta_mom_B(i,j) = delta_mom_B(i,j) + r(i,j)*(delta_mom_Y1(i,j)*&
                (x(i,j)/r(i,j)) - delta_mom_X1(i,j)*(y(i,j)/r(i,j)))
           
           delta_mom_A(i,j) = delta_mom_A(i,j) + rho(i,j)*( (mom_B(i,j)/rho(i,j))**2d0 )/( r(i,j)**3d0 )

!           delta_mom_B(i,j) = delta_mom_B(i,j) - rho(i,j)*( (-y(i,j)*v_x(i,j)+x(i,j)*v_y(i,j))*(x(i,j)*v_x(i,j)+y(i,j)*v_y(i,j)) )/( r(i,j)**3d0)

           if (omega_grid .gt. 1d-9) then
              delta_mom_X(i,j) = delta_mom_X(i,j) - omega_grid*mom_Y(i,j)
              delta_mom_Y(i,j) = delta_mom_Y(i,j) + omega_grid*mom_X(i,j)
           end if

     end do
  end do
  !$OMP END PARALLEL DO



end subroutine grid_timedifs
