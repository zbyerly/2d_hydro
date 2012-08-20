subroutine radial_initial_annulus(n,dr,rho,phi,omega,omega_k)
  implicit none
  integer :: i,j,integrate,n,nstar
  
  double precision :: pi 
  double precision :: r(n),dr,rmax
  double precision :: rho(n),phi(n),omega(n)
  double precision :: kappa,G,eps,eta(n),rho_over_tau,little_g

  double precision :: M_c,a,q,omega_k,r_outer,r_inner,c_1,c_2,temp

  double precision :: r_rho_max,rho_max

  pi = 4d0*atan(1d0)

  kappa = 1d0
  G     = 1d0

  rmax = 2.2d-4

  dr = rmax/dble(n)

  M_c = 2d-2

  !assuming q=2

  eps = 0.133d0
!  eps = 0.9d0
  r_outer = 1.0747d-4
  r_inner = r_outer*(1d0-eps)/(1d0+eps)
  print*,'r_outer=',r_outer
  print*,'r_inner=',r_inner


  if (r_inner .lt. 0d0) then
     print*,'r_inner less than zero!'
  end if
  
  C_2 = sqrt(2d0*G*M_c*R_inner*R_outer/(R_inner+R_outer))
  C_1 = 0.5d0*(C_2/R_inner)**2d0-G*M_c/R_inner
  
!  omega_k = 2d0/(eps*r_0)
!  M_c = (omega_k**2d0)*(r_0**3d0)


!  M_c = r_0/(eps*eps)

!  omega_k = sqrt(G*M_c/(r_0**3d0))

  print*,'eps=',eps
  print*,'c1=',c_1
  print*,'c2=',c_2
  print*,'M_c=',M_c

  
  rho_max = 0d0

  do i=1,n
     r(i) = (i-0.5d0)*dr
!     eta(i) = (r(i) - r_0)/a
     phi(i) = -G*M_c/r(i)
     omega(i) = c_2/(r(i)**2d0)
     
!     rho_over_tau=(omega_k*omega_k)/(4d0*pi*G)

     if ( (r(i) .gt. r_inner) .and. (r(i) .lt. r_outer)) then        
!        rho(i) = (2d0*q-3d0)*(omega_k**2d0)*(a*a-(r(i)-r_0)**2d0)/(4d0*kappa)
!        rho(i) = rho_over_tau*( cos(little_g*eta(i))/cos(little_g) - 1d0 )
        rho(i) = (0.5d0/kappa)*(C_1+G*M_c/r(i)-0.5d0*(C_2/r(i))**2d0)
     else
        rho(i) = 0d0
     end if

     if ( r(i) .lt. r_inner) then
        omega(i) = C_2/(r_inner**2)
     end if

     if (rho(i) .gt. rho_max) then
        rho_max = rho(i)
        r_rho_max = r(i)
        omega_k = omega(i)
     end if

  end do

  

  print*,'r_rho_max =',r_rho_max
  print*,'omega_k=',omega_k

  open(15,file='initial_radial_profile.dat')
  do i=1,n-1
     write(15,*) r(i),rho(i),omega(i)
  end do
  close(15)
    

end subroutine radial_initial_annulus
