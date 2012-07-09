subroutine checknan(nx,ny,rho,mom_A,mom_B,mom_x,mom_y,etot,tau,quit)
  implicit none
  include 'variables.h'
  
  double precision :: mom_A(nx,ny),mom_B(nx,ny),mom_x(nx,ny),mom_y(nx,ny)
  double precision :: rho(nx,ny),etot(nx,ny),tau(nx,ny)

  integer :: quit

  do i=1,nx
     do j=1,ny
        if ( isnan(rho(i,j)) ) then
           print*,'ERROR, rho = NaN',i,j
           quit = 1
        end if

        if ( isnan(mom_A(i,j)) ) then
           print*,'ERROR, mom_A = NaN',i,j
           quit = 1
        end if

        if ( isnan(mom_B(i,j)) ) then
           print*,'ERROR, mom_B = NaN',i,j
           quit = 1
        end if

        if ( isnan(mom_x(i,j)) ) then
           print*,'ERROR, mom_x = NaN',i,j
           quit = 1
        end if

        if ( isnan(mom_y(i,j)) ) then
           print*,'ERROR, mom_y = NaN',i,j
           quit = 1
        end if

        if ( isnan(etot(i,j)) ) then
           print*,'ERROR, etot = NaN',i,j
           quit = 1
        end if

        if ( isnan(tau(i,j)) ) then
           print*,'ERROR, tau = NaN',i,j
           quit = 1
        end if

     end do
  end do
  
end subroutine checknan
