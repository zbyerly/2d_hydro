subroutine radial_initial_annulus(n,dr,rho,phi,omega,omega_k)
  implicit none
  integer :: i,j,integrate,n,nstar
  
  double precision :: pi 
  double precision :: r(n),dr,rmax
  double precision :: rho(n),phi(n),omega(n)
  double precision :: kappa,G,eps,eta(n),rho_over_tau,little_g

  double precision :: M_c,a,q,omega_k,r_0

  pi = 4d0*atan(1d0)

  kappa = 1d0
  G     = 1d0

  rmax = 2.2d-4

  dr = rmax/dble(n)

  eps = 0.133d0
  q = 2d0
  r_0 = 9.486d-5
  a = r_0*eps

  little_g = a*sqrt(2d0*pi*G/kappa)

  

!  M_c = r_0/(eps*eps)

!  omega_k = sqrt(G*M_c/(r_0**3d0))

  omega_k = 2d0/(eps*r_0)
  M_c = (omega_k**2d0)*(r_0**3d0)

  print*,'eps=',eps
  print*,'a=',a
  print*,'M_c=',M_c
  print*,'omega_k=',omega_k


  do i=1,n
     r(i) = (i-0.5d0)*dr
     eta(i) = (r(i) - r_0)/a
     phi(i) = -G*M_c/r(i)
     omega(i) = omega_k*(r_0/r(i))**q
     
     rho_over_tau=(omega_k*omega_k)/(4d0*pi*G)

     if ( (r(i) .gt. r_0-a) .and. (r(i) .lt. r_0+a)) then        
!        rho(i) = (2d0*q-3d0)*(omega_k**2d0)*(a*a-(r(i)-r_0)**2d0)/(4d0*kappa)
        rho(i) = rho_over_tau*( cos(little_g*eta(i))/cos(little_g) - 1d0 )
     else
        rho(i) = 0d0
     end if

     if ( r(i) .lt. r_0-1.2d0*a) then
        omega(i) = omega_k*(r_0/(r_0-1.2d0*a))**q
     end if

  end do

  open(15,file='initial_radial_profile.dat')
  do i=1,n-1
     write(15,*) r(i),phi(i),omega(i)
  end do
  close(15)
    
end subroutine radial_initial_annulus
