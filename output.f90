subroutine output(nx,ny,rho,x,y,timestep,timeint,gamma,&
     phi,mom_A,mom_B,mom_x,mom_y,etot,tau,omega_grid)
  implicit none
  include 'variables.h'
  include 'grid.h'
!  double precision :: rho(nx,ny),x(nx,ny),y(nx,ny),phi(nx,ny)
!  double precision :: mom_A(nx,ny),mom_B(nx,ny),etot(nx,ny)
  double precision :: e_kinetic(nx,ny),e_internal(nx,ny),pressure(nx,ny)
! double precision :: solution_rho(nx)
  double precision :: solution_rho(ny),pressure_left(nx,ny),pressure_right(nx,ny)
  double precision :: e_internal2(nx,ny),v_x(nx,ny),v_y(nx,ny)
  double precision :: mom_theta(nx,ny),mom_r(nx,ny),out_time_last
  double precision :: omega_grid
  
  double precision :: r_here,rho_slice(150),omega_slice(150),output_freq
  integer :: slice_count(150),index_here

  integer :: timeint

  character filename*40  


  ! this is a comment
  
  call velocity(nx,ny,mom_x,mom_y,rho,v_x,v_y,omega_grid,x,y)


!  e_internal2 = tau**gamma
  
  !     call analytic(time,nx,x,dx,solution_rho)
!  call analytic(time,ny,y,dy,solution_rho)

  write(*,*) "outputting data... timeint =",timeint

!  if ( mod(timeint,100) .eq. 0 ) then
     write(filename,'(A10,i7.7,A4)') './data/rho',timeint,'.dat'
     open(10,FILE=filename)

!     write(filename,'(A16,i7.7,A4)') './data/mom_theta',timeint,'.dat'
!     open(11,FILE=filename)

!     write(filename,'(A12,i7.7,A4)') './data/mom_A',timeint,'.dat'
!     open(12,FILE=filename)

!     write(filename,'(A12,i7.7,A4)') './data/mom_B',timeint,'.dat'
!     open(13,FILE=filename)

!     write(filename,'(A11,i7.7,A4)') './data/etot',timeint,'.dat'
!     open(14,FILE=filename)

!     write(filename,'(A12,i7.7,A4)') './data/mom_x',timeint,'.dat'
!     open(15,FILE=filename)

!     write(filename,'(A12,i7.7,A4)') './data/mom_y',timeint,'.dat'
!     open(16,FILE=filename)

!     write(filename,'(A11,i7.7,A4)') './data/pres',timeint,'.dat'
!     open(17,FILE=filename)

!     write(filename,'(A10,i7.7,A4)') './data/tau',timeint,'.dat'
!     open(18,FILE=filename)

!     write(filename,'(A12,i7.7,A4)') './data/mom_R',timeint,'.dat'
!     open(19,FILE=filename)

     write(filename,'(A10,i7.7,A4)') './data/v_x',timeint,'.dat'
     open(20,FILE=filename)

     write(filename,'(A10,i7.7,A4)') './data/v_y',timeint,'.dat'
     open(21,FILE=filename)

     do j=1,Ny
!     j=Ny/2
        do i=1,Nx
!        i = nx/2
           write(10,*) x(i,j), y(i,j), rho(i,j)
!           write(10,*) x(i,j), rho(i,j),solution_rho(i)
!           write(11,*) x(i,j), y(i,j), rho(i,j) - solution_rho(j)
!           write(11,*) x(i,j), y(i,j), etot(i,j)
!           write(12,*) x(i,j), y(i,j), mom_theta(i,j)/(rho(i,j)*(x(i,j)*x(i,j)+y(i,j)*y(i,j)))
!           write(12,*) x(i,j), y(i,j), mom_A(i,j)
!           write(13,*) x(i,j), y(i,j), mom_B(i,j)
!           write(14,*) x(i,j), y(i,j), etot(i,j)
!           write(15,*) x(i,j), y(i,j), mom_x(i,j)
!           write(16,*) x(i,j), y(i,j), mom_y(i,j)
!           write(17,*) x(i,j), y(i,j), pressure(i,j)
!           write(18,*) x(i,j), y(i,j), tau(i,j)
!           write(19,*) x(i,j), y(i,j), mom_r(i,j)
           write(20,*) x(i,j), y(i,j), v_x(i,j)
           write(21,*) x(i,j), y(i,j), v_y(i,j)
        end do
     end do

     close(10)           
!     close(11)           
!     close(12)           
!     close(13)           
!     close(14)           
!     close(15)           
!     close(16)           
!     close(17)           
!     close(18)           
!     close(19)           
     close(20)           
     close(21)           

 

  
end subroutine output
