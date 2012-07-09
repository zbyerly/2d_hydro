!state / change in state / last state
  double precision :: rho(nx,ny),delta_rho(nx,ny),rho_0(nx,ny)
  double precision :: tau(nx,ny),delta_tau(nx,ny),tau_0(nx,ny)
  double precision :: etot(nx,ny),delta_etot(nx,ny),etot_0(nx,ny)
  double precision :: mom_A(nx,ny),delta_mom_A(nx,ny),mom_A_0(nx,ny)
  double precision :: mom_B(nx,ny),delta_mom_B(nx,ny),mom_B_0(nx,ny)
  double precision :: mom_x(nx,ny),delta_mom_x(nx,ny),mom_x_0(nx,ny)
  double precision :: mom_y(nx,ny),delta_mom_y(nx,ny),mom_y_0(nx,ny)
  
  ! L/R sides and F of left face (x-direction)
    double precision :: rho_left_x(nx,ny),rho_right_x(nx,ny),rho_Fx(nx,ny)
    double precision :: tau_left_x(nx,ny),tau_right_x(nx,ny),tau_Fx(nx,ny)
    double precision :: etot_left_x(nx,ny),etot_right_x(nx,ny),etot_Fx(nx,ny)
    double precision :: mom_A_left_x(nx,ny),mom_A_right_x(nx,ny),mom_A_Fx(nx,ny)
    double precision :: mom_B_left_x(nx,ny),mom_B_right_x(nx,ny),mom_B_Fx(nx,ny)
    double precision :: mom_x_left_x(nx,ny),mom_x_right_x(nx,ny),mom_x_Fx(nx,ny)
    double precision :: mom_y_left_x(nx,ny),mom_y_right_x(nx,ny),mom_y_Fx(nx,ny)
    
    ! L/R sides and F of left face (y-direction)
    double precision :: rho_left_y(nx,ny),rho_right_y(nx,ny),rho_Fy(nx,ny)
    double precision :: tau_left_y(nx,ny),tau_right_y(nx,ny),tau_Fy(nx,ny)
    double precision :: etot_left_y(nx,ny),etot_right_y(nx,ny),etot_Fy(nx,ny)
    double precision :: mom_A_left_y(nx,ny),mom_A_right_y(nx,ny),mom_A_Fy(nx,ny)
    double precision :: mom_B_left_y(nx,ny),mom_B_right_y(nx,ny),mom_B_Fy(nx,ny)
    double precision :: mom_x_left_y(nx,ny),mom_x_right_y(nx,ny),mom_x_Fy(nx,ny)
    double precision :: mom_y_left_y(nx,ny),mom_y_right_y(nx,ny),mom_y_Fy(nx,ny)
    
    
    !static variables
    double precision :: x(nx,ny),y(nx,ny),r(nx,ny),theta(nx,ny),phi(nx,ny)
  
  double precision :: x_left_x(nx,ny),y_left_x(nx,ny),r_left_x(nx,ny)
    double precision :: x_right_x(nx,ny),y_right_x(nx,ny),r_right_x(nx,ny)
    double precision :: theta_left_x(nx,ny),phi_left_x(nx,ny)
    double precision :: theta_right_x(nx,ny),phi_right_x(nx,ny)
  
  double precision :: x_left_y(nx,ny),y_left_y(nx,ny),r_left_y(nx,ny)
    double precision :: x_right_y(nx,ny),y_right_y(nx,ny),r_right_y(nx,ny)
    double precision :: theta_left_y(nx,ny),phi_left_y(nx,ny)
    double precision :: theta_right_y(nx,ny),phi_right_y(nx,ny)
