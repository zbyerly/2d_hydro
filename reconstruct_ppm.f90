
subroutine reconstruct_ppm(nx,ny,u,u_left_x,u_right_x,u_left_y,u_right_y)
  implicit none

  double precision :: u(nx,ny)
  double precision :: u_left_x(nx,ny),u_right_x(nx,ny)
  double precision :: u_left_y(nx,ny),u_right_y(nx,ny)
  double precision :: s1,s2,u_slope(nx,ny),minmod,tmp1,tmp2,condition
  integer :: i,j,nx,ny

  u_left_x = 1d-9
  u_left_y = 1d-9
  u_right_x = 1d-9
  u_right_y = 1d-9

  !!left/right       
  do j=1,ny
     do i=2,nx-1       
        s1 = u(i+1,j) - u(i,j)
        s2 = u(i,j) - u(i-1,j)
        u_slope(i,j) = minmod(s1,s2)
     end do
  end do

  do j=1,ny
     do i=2,nx-2
        u_right_x(i,j)=( u(i,j) + u(i+1,j) ) * 0.5d0 + (u_slope(i,j)-u_slope(i+1,j))*(1d0/6d0)
        u_left_x(i,j)=u_right_x(i-1,j)
     end do
  end do

  do j=1,ny
     do i=3,nx-2
        tmp1 = u_right_x(i,j) - u_left_x(i,j)
        tmp2 = u_right_x(i,j) + u_left_x(i,j)
        condition = ( u_right_x(i,j)-u(i,j) )*( u(i,j)-u_left_x(i,j) )
        if (condition .le. 0d0) then           
           u_right_x(i,j) = u(i,j)
           u_left_x(i,j) = u(i,j)
        else if (tmp1*(u(i,j)-0.5d0*tmp2) .gt. (1d0/6d0)*tmp1*tmp1) then           
           u_left_x(i,j) = 3d0*u(i,j)-2d0*u_right_x(i,j)           
        else if (-(1d0/6d0)*tmp1*tmp1 .gt. tmp1*(u(i,j)-0.5d0*tmp2)) then           
           u_right_x(i,j) = 3d0*u(i,j)-2d0*u_left_x(i,j)          
        end if
       
     end do
  end do

  do j=1,ny
     do i= nx-2,4,-1
        u_right_x(i,j) = u_right_x(i-1,j)
     end do
  end do

  !!up/down *************************   
  do j=2,ny-1
     do i=1,nx       
        s1 = u(i,j+1) - u(i,j)
        s2 = u(i,j) - u(i,j-1)
        u_slope(i,j) = minmod(s1,s2)
     end do
  end do

  do j=2,ny-2
     do i=1,nx
        u_right_y(i,j)=( u(i,j) + u(i,j+1) ) * 0.5d0 + (u_slope(i,j)-u_slope(i,j+1))*(1d0/6d0)
        u_left_y(i,j)=u_right_y(i,j-1)
     end do
  end do

  do j=3,ny-2
     do i=1,nx
        tmp1 = u_right_y(i,j) - u_left_y(i,j)
        tmp2 = u_right_y(i,j) + u_left_y(i,j)
        condition = ( u_right_y(i,j)-u(i,j) )*( u(i,j)-u_left_y(i,j) )
        if (condition .le. 0d0) then           
           u_right_y(i,j) = u(i,j)
           u_left_y(i,j) = u(i,j)
        else if (tmp1*(u(i,j)-0.5d0*tmp2) .gt. (1d0/6d0)*tmp1*tmp1) then           
           u_left_y(i,j) = 3d0*u(i,j)-2d0*u_right_y(i,j)           
        else if (-(1d0/6d0)*tmp1*tmp1 .gt. tmp1*(u(i,j)-0.5d0*tmp2)) then           
           u_right_y(i,j) = 3d0*u(i,j)-2d0*u_left_y(i,j)          
        end if
       
     end do
  end do

  do j= ny-2,4,-1
     do i=1,nx
        u_right_y(i,j) = u_right_y(i,j-1)
     end do
  end do



end subroutine reconstruct_ppm
