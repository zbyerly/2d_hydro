subroutine grid_store(nx,ny,&
     rho,tau,etot,mom_A,mom_B,mom_x,mom_y,&
     rho_0,tau_0,etot_0,mom_A_0,mom_B_0,mom_x_0,mom_y_0)
  implicit none
  integer :: nx,ny
  double precision :: rho(nx,ny),tau(nx,ny),etot(nx,ny)
  double precision :: mom_A(nx,ny),mom_B(nx,ny)
  double precision :: rho_0(nx,ny),tau_0(nx,ny),etot_0(nx,ny)
  double precision :: mom_A_0(nx,ny),mom_B_0(nx,ny)

  double precision :: mom_x(nx,ny),mom_y(nx,ny),mom_x_0(nx,ny),mom_y_0(nx,ny)

  rho_0 = rho
  tau_0 = tau
  etot_0 = etot
  mom_A_0 = mom_A
  mom_B_0 = mom_B
  mom_x_0 = mom_x
  mom_y_0 = mom_y


end subroutine grid_store
