
subroutine analytic(t,Ncells,x,dx,solution_rho)
  implicit none
  
  ! This analytic solution only valid for wave travelling to the right
  ! and with x_min = -0.5d0  x_max = 0.5d0
  ! and with x_0 = -0.1d0, gamma = 1.4
  !
  ! but, with any resolution

  integer :: Ncells,i
  double precision :: t,dx
  double precision :: x(Ncells), rho(Ncells)
  double precision :: c_s(Ncells), pressure(Ncells)
  double precision :: solution_rho(Ncells)

  
  double precision :: gamma, exponent
  double precision :: x_0,c_1,w_3,c_3,w_4,W,c_5,p_4,p_5,w_2
  double precision :: x_min,x_max
  double precision :: x_1,x_3,x_4,x_5
  double precision :: c_1_byerly,w_3_byerly,c_3_byerly,w_4_byerly,c_5_byerly,p_5_byerly,p_4_byerly
  double precision :: x_head,x_tail,x_contact,x_shock
  integer :: N_1, N_4, N_5, N_3
  
  x_0 = -0.1d0
  x_min = -0.5d0
  x_max = 0.5d0
  gamma = 1.4d0

  ! the idea here is that we will use values from Motl's paper for 
  ! w's, c's, p's, then afterwards check to make sure the values are consist.

  !sound speed of undisturbed gas on the left

  !  head of rarefaction wave speed 
  c_1 = 1.183d0 !Motl's value
  x_head = x_0 - c_1*t
  
  w_3 = 0.9274d0 !Motl's value
  c_3 = 0.9978d0 !Motl's value

  !  tail of rarefaction wave speed
  x_tail = x_0 + (w_3 - c_3)*t

  w_4 = 0.9274d0 !Motl's value

  x_contact = x_0 + w_4*t
  
  c_5 = 1.058d0 !Motl's value
  p_4 = 0.3031d0 !Motl's value
  p_5 = 0.1d0 !Motl's value
  
  W = c_5*sqrt(1 + (gamma+1)*(p_4 - p_5)/(2d0*gamma*p_5))

  x_shock = x_0 + W*t

  !Now, compare values:
  
!!$  x_1 = (x_head+(-0.5d0))/2d0 
!!$  N_1 = floor((x_1-x_min)/dx)
!!$  
!!$  x_3 = (x_tail+x_contact)/2d0
!!$  N_3 = floor((x_3-x_min)/dx)
!!$		
!!$  x_4 = (x_contact+x_shock)/2d0
!!$  N_4 = floor((x_4-x_min)/dx)
!!$  
!!$  x_5 = (x_shock+x_max)/2d0
!!$  N_5 = floor((x_5-x_min)/dx)

!!$  c_1_byerly = c_s(N_1)
!!$
!!$  c_3_byerly = c_s(N_3)
!!$  w_3_byerly = momentum_cell(N_3)/rho(N_3)
!!$
!!$  w_4_byerly = momentum_cell(N_4)/rho(N_4)
!!$  p_4_byerly = pressure(N_4)
!!$
!!$  c_5_byerly = c_s(N_5)
!!$  p_5_byerly = pressure(N_5)

!!$  print*,'Motl,byerly'
!!$  print*,c_1,c_1_byerly
!!$  print*,c_3,c_3_byerly
!!$  print*,w_3,w_3_byerly
!!$  print*,w_4,w_4_byerly
!!$  print*,p_4,p_4_byerly
!!$  print*,c_5,c_5_byerly
!!$  print*,p_5,p_5_byerly
!!$  print*,'x_head,x_tail,x_contact,x_shock'
!!$  print*,x_head,x_tail,x_contact,x_shock
!!$  print*,'x_1,x_3,x_4,x_5'
!!$  print*,x_1,x_3,x_4,x_5

do i=1,Ncells
   if (x(i) .lt. x_head) then
      solution_rho(i)= 1.0d0 
   else if (x(i) .lt. x_tail) then 
      w_2 = 2d0*(c_1+(x(i)-x_0)/t)/(gamma+1d0)
      exponent = 2d0/(gamma-1d0)
      solution_rho(i)= (1d0-(gamma-1d0)*w_2/(2d0*c_1))**exponent
   else if (x(i) .lt. x_contact) then
      solution_rho(i)= 0.4263d0
   else if (x(i) .lt. x_shock) then
      solution_rho(i)= 0.2656d0
   else 
      solution_rho(i)= 0.125d0
   end if
   
end do

end subroutine analytic
