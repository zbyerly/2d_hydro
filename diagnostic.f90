subroutine diagnostic(nx,ny,dx,dy,kappa,gamma,&
     rho,tau,etot,mom_A,mom_B,x,y,time,dt,J_tot_last,M_tot_last,timestep)
  implicit none
  include 'variables.h'
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny),mom_A(nx,ny),&
       mom_B(nx,ny),x(nx,ny),y(nx,ny)
  double precision :: mom_r(nx,ny),mom_theta(nx,ny),J_tot,r(nx,ny),theta(nx,ny)
  double precision :: J_tot_last,J_tot_dot,M_tot,E_tot,E_kin_tot,e_kinetic(nx,ny)
  double precision :: M_tot_last,M_tot_dot
  double precision :: mom_r_tot,max_ang_mom
  double precision :: com_x,com_y,r_here
  
!  call getrtheta(nx,ny,x,y,r,theta)
!  call kinetic_energy(nx,ny,r,theta,rho,mom_x,mom_y,e_kinetic,x,y)

           mom_r = mom_A
           mom_theta = mom_B

  J_tot = 0d0
  M_tot = 0d0
  E_tot = 0d0
!  E_kin_tot = 0d0
  mom_r_tot = 0d0
  com_x = 0d0
  com_y = 0d0

  max_ang_mom = 0d0

  do j=3,ny-2
     do i=3,nx-2
        J_tot = J_tot + mom_theta(i,j)*dx*dy
        M_tot = M_tot + rho(i,j)*dx*dy
        com_x = x(i,j)*rho(i,j)*dx*dy
        com_y = y(i,j)*rho(i,j)*dx*dy        
        E_tot = E_tot + etot(i,j)*dx*dy
!        E_kin_tot = E_kin_tot + e_kinetic(i,j)*dx*dy
        mom_r_tot = mom_r_tot + mom_r(i,j)

        if (mom_theta(i,j) .gt. max_ang_mom) then
           max_ang_mom = mom_theta(i,j)
        end if

     end do
  end do

  com_x = com_x/M_tot
  com_y = com_y/M_tot


!  J_tot_dot = (J_tot-J_tot_last)/dt
!  M_tot_dot = (M_tot-M_tot_last)

  write(30,*) time, J_tot!, J_tot_dot
  write(31,*) timestep, time, dt
  write(32,*) time, E_tot
  write(33,*) time, M_tot
!  write(34,*) time, E_kin_tot
  write(35,*) time, mom_r_tot
  write(36,*) time, mom_r_tot + J_tot
  write(37,*) time, com_x, com_y
  write(38,*) time, timestep, max_ang_mom


  J_tot_last = J_tot
  M_tot_last = M_tot

!  print*,'timestep  dt__           E_kin_tot         time            J_tot'
!  write(*,"(I7,D15.4,D15.4,D15.4,D15.4)"),timestep,dt,E_kin_tot,time,J_tot

!  print*,'max_ang_mom =',max_ang_mom

end subroutine diagnostic
