subroutine driver(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
       rho,tau,etot,mom_A,mom_B,x,y,phi,alpha,recons,omega_grid)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny)
  double precision :: mom_A(nx,ny),mom_B(nx,ny),phi(nx,ny),x(nx,ny),y(nx,ny)  
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: beta,dt_cfl,J_tot_last,M_tot_last,flux_in,flux_in_tot
  double precision :: out_time_last,output_freq

  double precision :: alpha, omega_grid

  integer :: outcount,output_yes,timeint,quit,recons

  timestep = 0
  time = 0d0
  
  flux_in_tot = 0d0

  out_time_last = 0d0
  !output frequency for rotating polytrope:
!  output_freq = 1d-1

  !output freq for annulus
  output_freq = 6.28318530718d0/omega_grid/1d2
!  output_freq = 6.28318530718d0/omega_grid/1d6
  outcount = 0
  print*,'output_freq=',output_freq

! defining mom_x,mom_y from mom_A,mom_B 
  call mom_cyl2cart(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y)


  !output initial data
  open(30,file='J_tot.dat') 
  open(31,file='diag.dat')
  call output(nx,ny,rho,x,y,timestep,time,gamma,&
       phi,mom_A,mom_B,mom_x,mom_y,etot,tau,output_freq)
 ! call exit(0)
  print*,'endtime = ',endtime
  do while (time .lt. endtime)
!  do while (timestep .lt. 100)
     
     !     print*,'calling timestep'
     call grid_timestep(nx,ny,dx,dy,kappa,gamma,cfl_factor,&
          rho,tau,etot,mom_x,mom_y,dt_cfl,omega_grid,x,y)
     
          print*,'starting timestep',timestep
     if (timestep .lt. 2) then
        dt_last = 1d-16
     end if
     
     if ( (dt_cfl .gt. dt_last*1.5d0) .and. (time .lt. output_freq) )then
        dt = dt_last*1.5d0
     else         
        dt = dt_cfl
     end if

     if (dt .gt. (outcount*output_freq-time)) then
        dt = outcount*output_freq - time
        output_yes = 1
        outcount = outcount + 1
     end if

     !     print*,'calling diagnostic'

     !     if (mod(timestep,10) .eq. 0) then
     call diagnostic(nx,ny,dx,dy,kappa,gamma,&
          rho,tau,etot,mom_A,mom_B,x,y,time,dt,J_tot_last,M_tot_last,timestep)
     !     end if

     call RK_step(nx,ny,dx,dy,kappa,gamma,dt,rho_floor,&
          rho,tau,etot,mom_A,mom_B,mom_x,mom_y,phi,x,y,flux_in,alpha,recons,omega_grid)
     
     timestep = timestep + 1
     time = time + dt     

     !       !         !            !             !          ! 
     !     print*,timestep,dt,dt_cfl,time,dt_last
     dt_last = dt
       
     call checknan(nx,ny,rho,mom_A,mom_B,mom_x,mom_y,etot,tau,quit)

     if (quit .eq. 1) then
        output_yes = 1
     end if

!        call output(nx,ny,rho,x,y,t
!     timeint = floor(time/output_freq)
     timeint = outcount
     if ( output_yes .eq. 1 ) then
        call output(nx,ny,rho,x,y,timestep,timeint,gamma,&
             phi,mom_A,mom_B,mom_x,mom_y,etot,tau,omega_grid)
        out_time_last = time
     end if
     output_yes = 0

     if (quit .eq. 1) then
        print*,'exiting...'
        call exit(0)
     end if

  end do

!     timeint = floor(time/output_freq)
!     timeint = timestep
!     if ( output_yes .eq. 1 ) then
!        call output(nx,ny,rho,x,y,timestep,timeint,gamma,&
!             phi,mom_A,mom_B,mom_x,mom_y,etot,tau)
!        out_time_last = time
!     end if


  close(31)
  close(30)

end subroutine driver
