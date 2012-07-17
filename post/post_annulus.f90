program post
  implicit none
  
  integer, parameter :: nx = 200
  integer, parameter :: ny = 200
  integer, parameter :: rbins = 100
  integer, parameter :: n_theta = 128
  integer :: i,j,k
  double precision :: pi,d_theta
  double precision :: gamma,kappa,cfl_factor,endtime
  double precision :: dx,dy
  double precision :: rho_floor
  double precision :: time,dt,dt_last, dt_here
  integer :: mom_geom,timestep
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny),mom_A(nx,ny),mom_B(nx,ny)
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny),phi(nx,ny)
  double precision :: e_kinetic(nx,ny),e_internal(nx,ny),pressure(nx,ny)
  double precision :: v_x(nx,ny),v_y(nx,ny)
  double precision :: mom_theta(nx,ny),mom_r(nx,ny),J_enc(rbins),M_enc(rbins)
  double precision :: E_tot(nx,ny), sr_x(nx,ny), sr_y(nx,ny), dt_x(nx,ny)
  double precision :: dt_y(nx,ny),c_s(nx,ny)
  double precision :: x_min,x_max,y_min,y_max
  integer :: timeint
  character filename*40  
  
  double precision :: M_tot,r_here
  double precision :: rmax,dr
  integer :: r_index,r_count(rbins)
  double precision :: rho_theta(n_theta)
  double precision :: rho_1D(rbins),tau_over_rho_1D(ny),pressure_1D(ny),vel_1D(ny)
  double precision :: Arho_1D(ny),Atau_over_rho_1D(ny),Apressure_1D(ny),Avel_1D(ny)
  double precision :: spec_ang_mom_1d(rbins)
  double precision :: omega_r(rbins)

  double precision :: i_double,j_double
  double precision :: rho_LL,rho_UL,rho_LR,rho_UR
  double precision :: LL_dist,UL_dist,LR_dist,UR_dist,total_dist
  double precision :: theta_here, x_here, y_here,r_0
  integer :: i_minus,i_plus,j_minus,j_plus

  double precision :: total_mom_B,rho_here,transform(n_theta)


!  double precision :: M_1,M_2,M_3,J_1,J_2,J_3

  double precision :: x_head,x_0,c_1,w_3,c_3,x_tail,t,w_4,x_contact,c_5,p_4,p_5,x_shock,W,w_2,exponent

  cfl_factor = 0.4d0

  gamma = 2d0
  kappa = 1d0
  
  x_min = -1.5d-4
  x_max = 1.5d-4
  y_min = -1.5d-4
  y_max = 1.5d-4
  
!  x_min = -1.5d0
!  x_max = 1.5d0
!  y_min = -1.5d0
!  y_max = 1.5d0

  rho_floor = 1d-9
  mom_geom = 0

  dx = (x_max - x_min)/dble(Nx)
  dy = (y_max - y_min)/dble(Ny)

!&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&  LOOP OVER TIMES  &&&
!&&&&&&&&&&&&&&&&&&&&&&&&&
  
  open(11,file='mode1.dat')
  open(12,file='mode2.dat')
  open(13,file='mode3.dat')
  open(14,file='mode4.dat')

  
  do timeint=0,5000,1
     
     print*,timeint
     
     write(filename,'(A10,i7.7,A4)') './data/rho',timeint,'.dat'
     open(10,FILE=filename,status="old")

     rmax = 1.5d0
     dr = 1.5d0/dble(rbins)
     
     do i=1,rbins
        rho_1d(i) = 0d0
        spec_ang_mom_1d(i) = 0d0
        omega_r(i) = 0d0
     end do
          
     !actually read in the data
     do j=1,ny
        do i=1,nx
           read(10,*) x(i,j), y(i,j), rho(i,j)

           r(i,j) = sqrt(x(i,j)*x(i,j) + y(i,j)*y(i,j))
!           theta_here = datan2(y(i,j),x(i,j))           

        end do
     end do

     close(10)

     pi = 4d0*atan(1d0)
     d_theta = (2d0*pi)/n_theta
     r_0 = 9.486d-5
!     print*,d_theta,r_0,dx,dy
     do k=1,n_theta
        theta_here = (dble(k)-0.5d0)*d_theta
        x_here = r_0*cos(theta_here)
        y_here = r_0*sin(theta_here)
        i_double = (x_here-x_min)/dx
        j_double = (y_here-y_min)/dy
        
        i_minus = floor(i_double)
        i_plus = i_minus+1
        j_minus = floor(j_double)
        j_plus = j_minus+1

        rho_LL = rho(i_minus,j_minus)
        LL_dist = dsqrt((dble(i_minus)-i_double)**2d0+(dble(j_minus)-j_double)**2d0)
        rho_UL = rho(i_minus,j_plus)
        UL_dist = dsqrt((dble(i_minus)-i_double)**2d0+(dble(j_plus)-j_double)**2d0)
        rho_LR = rho(i_plus,j_minus)
        LR_dist = dsqrt((dble(i_plus)-i_double)**2d0+(dble(j_minus)-j_double)**2d0)
        rho_UR = rho(i_plus,j_plus)
        UR_dist = dsqrt((dble(i_plus)-i_double)**2d0+(dble(j_plus)-j_double)**2d0)
        
        total_dist = LL_dist + UL_dist + LR_dist + UR_dist

        rho_here = rho_LL*(LL_dist/total_dist)+rho_UL*(UL_dist/total_dist)+rho_LR*(LR_dist/total_dist)+rho_UR*(UR_dist/total_dist)

        rho_theta(k) = rho_here

!        write(30,*) theta_here, rho_here

     end do

     do i=1,n_theta
        transform(i) = 0d0
        do j=1,n_theta 
           transform(i) = transform(i) + rho_theta(j)*(cos(2d0*pi*(i-1)*(j-1)/dble(n_theta))) 
        end do
        transform(i) = transform(i)/dble(n_theta)
     end do
        
     write(11,*) timeint, transform(1)/transform(0)
     write(11,*) timeint, transform(2)/transform(0)
     write(11,*) timeint, transform(3)/transform(0)
     write(11,*) timeint, transform(4)/transform(0)


  end do

  !close the data files
  close(11)           
  close(12)           
  close(13)           
  close(14)
  
end program post
