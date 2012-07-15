program main
  implicit none

  integer, parameter :: nx=200,ny=200  
  double precision :: gamma,kappa,cfl_factor,endtime,time,alpha
  double precision :: dx,dy
  double precision :: rho_floor
  integer :: timestep, recons, rotation
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny)
  double precision :: mom_A(nx,ny),mom_B(nx,ny),phi(nx,ny),x(nx,ny),y(nx,ny)

  double precision :: omega_k, omega_grid
  double precision :: e_internal, e_kinetic(nx,ny)

  !openmp variables
  integer :: numthreads,OMP_GET_NUM_THREADS

  character(len=32) :: mom_geom
  character(len=32) :: reconstruction
  character(len=32) :: rotate

  call get_command_argument(1 , mom_geom)
  call get_command_argument(2 , reconstruction)
!  reconstruction = "minmod"
  call get_command_argument(3 , rotate)
 
 

  if (mom_geom .eq. 'cyl') then
     alpha = 0d0
  else if (mom_geom .eq. 'cart') then
     alpha = 1d0
  else if (mom_geom .eq. 'half') then
     alpha = 0.5d0
  else if (mom_geom .eq. 'quarter_cyl') then
     alpha = 0.25d0
  else if (mom_geom .eq. 'quarter_cart') then
     alpha = 0.75d0
  else
     print*,'valid cmdline arguments are cyl or cart'
  end if

  print*,'alpha =',alpha

  if (reconstruction .eq. 'minmod') then
     print*,'using minmod'
     recons = 0
  else if (reconstruction .eq. 'ppm') then
     print*,'using ppm'
     recons = 1
  else
!     print*,'valid cmdline arguments are minmod or ppm'
  end if

  if (rotate .eq. 'norotate') then
     print*,'rotation off'
     rotation = 0
  else if (rotate .eq. 'rotate') then
     print*,'rotation on'
     rotation = 1
  else
!     print*,'valid cmdline arguments are minmod or ppm'
     print*,'rotation off by default...'
     rotation = 0
  end if




  call initial_annulus(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
       rho,tau,etot,mom_A,mom_B,x,y,phi,omega_k)

!  call initial_sod(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
!       rho,tau,etot,mom_A,mom_B,x,y,phi)

  if (rotation .eq. 1) then
     omega_grid = omega_k
  else
     omega_grid = 0d0
  end if


  call driver(nx,ny,dx,dy,kappa,gamma,cfl_factor,endtime,rho_floor,&
       rho,tau,etot,mom_A,mom_B,x,y,phi,alpha,recons,omega_grid)

end program main
