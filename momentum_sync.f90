subroutine momentum_sync(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y,alpha)
  implicit none
  include 'variables.h'
  
  double precision :: mom_A(nx,ny),mom_B(nx,ny),mom_x(nx,ny),mom_y(nx,ny)
  double precision :: mom_x1(nx,ny),mom_y1(nx,ny)
  double precision :: mom_x2(nx,ny),mom_y2(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)

  double precision :: alpha


  call getrtheta(nx,ny,x,y,r,theta)

  !mom_x1,mom_y1 is cartesian momentum from the cyl_conserved

  call mom_cyl2cart(nx,ny,mom_A,mom_B,mom_x1,mom_y1,x,y)


!  print*,'alpha=',alpha

!  do i=1,nx
!     do j=1,ny
        
        ! setting alpha = 1 uses on cart mom

!        if (r(i,j) .lt. 0.015d0) then
!        if (r(i,j) .lt. 0.03d0) then
!        if (r(i,j) .lt. 0.06d0) then
!           alpha(i,j) = 1d0
!        else
!           alpha(i,j) = 0d0
!        end if

!     end do
!  end do

  !mom_x2,mom_y2 is the cartesian momentum after being synced

  mom_x2 = (1d0-alpha)*mom_x1+alpha*mom_x
  mom_y2 = (1d0-alpha)*mom_y1+alpha*mom_y
  
  mom_x = mom_x2
  mom_y = mom_y2

  call mom_cart2cyl(nx,ny,mom_A,mom_B,mom_x,mom_y,x,y)

  
!!$  open(71,FILE='mom_A_sync')
!!$  open(72,FILE='mom_B_sync')
!!$  open(73,FILE='mom_x_sync')
!!$  open(74,FILE='mom_y_sync')
!  open(75,FILE='alpha')

!  do i=1,nx
!     do j=1,ny
!!$        write(71,*) x(i,j),y(i,j),mom_A(i,j)
!!$        write(72,*) x(i,j),y(i,j),mom_B(i,j)
!!$        write(73,*) x(i,j),y(i,j),mom_x(i,j)
!!$        write(74,*) x(i,j),y(i,j),mom_y(i,j)
!        write(75,*) x(i,j),y(i,j),alpha(i,j)
!     end do
!  end do

!!$  close(71)
!!$  close(72)
!!$  close(73)
!!$  close(74)
!  close(75)
!  call exit(0)


end subroutine momentum_sync
