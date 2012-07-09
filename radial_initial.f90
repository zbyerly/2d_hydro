subroutine radial_initial(n,dr,rho,phi_1,omega)
  implicit none
  integer :: i,j,integrate,n,nstar
  
  double precision :: pi 
  double precision :: xi(n),r(n),nmax,dxi,dr,dj,ximax
  double precision :: rho(n),M_inc(n),pressure(n),phi(n),balance(n),phi_1(n)
  
  double precision :: kappa,G,r_star,theta_here,condition,constant,c_1,c_2
  integer :: n_star
  
  double precision :: inertia_moment
  double precision :: J_0(n),J_1_xi_1,J_0_xi_1,M_inc_1(n),balance_1(n)
  double precision :: dp_dr, dphi_dr, dphi1_dr,rho_here,r_here,xi_here,xi_1
  double precision :: theta(n)

  double precision :: dphi_dxi_1

  double precision :: omega,beta

  pi = 4d0*atan(1d0)

!  omega = 6.2831d-1
  beta =  (omega*omega)/(2d0*pi)

  print*,'omega=',omega

  !manually setting beta = omega^2/2*pi
!  beta = 0.05
  
  print*,'beta=',beta

  kappa=1d0

  ximax = 10d0
  dxi = ximax/dble(n)
   
  dr = dxi*sqrt((2d0*kappa)/(4d0*pi))
!  print*,'dr = ',dr,'dxi =',dxi
  
  integrate = 1

!  print*,'starting radial integrals'
!  print*,'computing bessel functions'

  i=0
  condition = 100
  do while (condition .gt. 0d0)
     i=i+1
     xi(i) = (i-1)*dxi
     r(i) = (i-1)*dr
     ! integrate
     J_0(i) = 0d0
     dj = pi/1000d0
     do j=1,1000
        theta_here=j*dj
        J_0(i) = J_0(i) + cos(xi(i)*sin(theta_here))*dj
     end do
     J_0(i) = (1d0/pi)*J_0(i)
     theta(i) = (1d0-beta)*J_0(i)+beta
     condition = theta(i)
  end do

  nstar = i-1
  
  !get a better zero for J_0 (xi_1)
!  xi_1 = -(J_0(nstar)*dxi)/(J_0(nstar+1)-J_0(nstar))+xi(nstar)
  xi_1 = -(theta(nstar)*dxi)/(theta(nstar+1)-theta(nstar))+xi(nstar)

  print*,'xi_1=',xi_1

  ! find J_1
  J_1_xi_1=0d0
  J_0_xi_1=0d0
  do j=1,1000
     theta_here=j*dj
     J_1_xi_1 = J_1_xi_1 + cos(theta_here - xi_1*sin(theta_here))*dj
     J_0_xi_1 = J_0_xi_1 + cos(xi_1*sin(theta_here))*dj
  end do
  J_1_xi_1 = J_1_xi_1/pi
  J_0_xi_1 = J_0_xi_1/pi

  
  do i=nstar+1,n
     xi(i) = (i-0.5d0)*dxi
     r(i) = (i-0.5d0)*dr
     theta(i) = 0d0
  end do

!  print*,'nstar = ',nstar,'xistar = ',xi_1

!  print*,'computing profiles from bessel fucntions'
  do i=1,n
!     rho(i) = (1-beta)J_0(i)+beta
!     rho(i) = J_0(i)
     rho(i) = theta(i)
     pressure(i) = kappa*(rho(i))**2d0
  end do

  rho_here = 0.5d0*(rho(1)+rho(2))
  M_inc(1) = 2d0*pi*rho_here*r(1)*dr
  do i=2,n-1
     rho_here = 0.5d0*(rho(i+1)+rho(i))
     M_inc(i) = M_inc(i-1)+2d0*pi*rho_here*r(i)*dr
  end do

  inertia_moment=0d0
  do i=1,nstar
!     M_inc_1(i) = kappa*xi(i)*J_1(i)
!     phi_1(i) = -2d0*kappa*J_0(i)
!     phi_1(i) = -2d0*(theta(i) - 1d0 - beta*xi(i)*xi(i)/4d0)
     phi_1(i) = 2d0*(1d0-beta)*J_0(i)+2d0*beta-2d0-beta*xi(i)*xi(i)/2d0
     phi_1(i) = -phi_1(i)
     inertia_moment = inertia_moment + rho(i)*(r(i)**3)*dr
  end do
!  c_1 = 2d0*kappa*xi_1*J_1_xi_1
!  c_1 = -2d0*xi_1*J_1_xi_1-beta*xi_1*xi_1
  c_1 = -2d0*(1d0-beta)*xi_1*J_1_xi_1-beta*xi_1*xi_1
  c_2 = 2d0*(1d0-beta)*J_0_xi_1+2d0*beta-2d0-beta*xi_1*xi_1/2d0-c_1*log(xi_1)

  dphi_dxi_1 = -1d0*(1d0-beta)*J_1_xi_1
  print*,'dphi_dxi_1=',dphi_dxi_1

  do i=nstar+1,n
!     M_inc_1(i) = M_inc(nstar)
!     phi_1(i) = c_1*log(xi(i)/xi_1) - 2d0*(theta(i) -1d0  - beta*xi_1*xi_1/4d0)
     phi_1(i) = c_1*log(xi(i))+c_2
     phi_1(i) = -phi_1(i)
  end do

  inertia_moment = 2*pi*inertia_moment

  print*,'intertiamoment=',inertia_moment

!  phi(1)=0d0
!  phi_1(1)=0d0
!  do i=1,n-1
!     xi_here = 0.5d0*(xi(i+1)+xi(i))
!     phi_1(i+1) = phi_1(i)+ (2d0/xi_here)*M_inc_1(i)*dxi
!     phi(i+1) = phi(i) + (2d0/xi_here)*M_inc(i)*dxi
!  end do
  
  ! look at radial pressure / potential balance
!  balance(1) = 0d0
!  balance_1(1) = 0d0
!  do i=2,n
!     balance(i) = ((pressure(i)-pressure(i-1)) + 0.5d0*(rho(i)+rho(i-1))*&
!          (phi(i)   - phi(i-1))     )/(pressure(i) - pressure(i-1))     
!     balance_1(i) = ((pressure(i)-pressure(i-1)) + 0.5d0*(rho(i)+rho(i-1))*&
!          (phi_1(i) - phi_1(i-1)) )/(pressure(i) - pressure(i-1))     
!  end do

  open(15,file='initial_radial_profile.dat')
  do i=1,n-1
     dp_dr = (pressure(i+1)-pressure(i))/dr
     dphi_dr = ((phi(i+1)-phi(i))/dr)
     dphi1_dr = ((phi_1(i+1)-phi_1(i))/dr)
     rho_here = 0.5d0*(rho(i+1)+rho(i))
     write(15,*) xi(i),pressure(i),phi_1(i)
  end do
  close(15)
  
!  call exit(0)
    
end subroutine radial_initial
