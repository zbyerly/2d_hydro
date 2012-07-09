subroutine mom_cyl2cart(nx,ny,mom_r,mom_theta,mom_x,mom_y,x,y)
  implicit none
  
  integer :: nx,ny,i,j
  double precision, intent(in)   :: mom_r(nx,ny),mom_theta(nx,ny)
  double precision  :: mom_x(nx,ny),mom_y(nx,ny)
  double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny)

  call getrtheta(nx,ny,x,y,r,theta)
  
  mom_x = mom_r*(x/r) - mom_theta*(y/r)/r
  mom_y = mom_r*(y/r) + mom_theta*(x/r)/r

end subroutine mom_cyl2cart
