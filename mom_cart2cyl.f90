subroutine mom_cart2cyl(nx,ny,mom_r,mom_theta,mom_x,mom_y,x,y)
  implicit none
  
  integer :: nx,ny,i,j
  double precision :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: mom_r(nx,ny),mom_theta(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)

  call getrtheta(nx,ny,x,y,r,theta)

  mom_r = mom_x*(x/r) + mom_y*(y/r)
  mom_theta = r*(mom_y*(x/r) - mom_x*(y/r))

end subroutine mom_cart2cyl
