subroutine boundaries(nx,ny,u)
  implicit none
  
  integer ::nx,ny,i,j
  double precision :: u(nx,ny)
  double precision :: avg
 
  !$OMP PARALLEL DO PRIVATE(j)
  do j=1,ny
     u(1,j) = u(3,j)
     u(2,j) = u(3,j)
     u(nx,j)   = u(nx-2,j)
     u(nx-1,j) = u(nx-2,j)
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE(i)
  do i=1,nx
     u(i,1) = u(i,3)
     u(i,2) = u(i,3)
     u(i,ny)   = u(i,ny-2)
     u(i,ny-1) = u(i,ny-2)     
  end do
  !$OMP END PARALLEL DO

!  avg = ( u(nx/2+1,ny/2+1) + u(nx/2+1,ny/2-1) + u(nx/2-1,ny/2+1) + u(nx/2-1,ny/2-1) )/4d0

!  u(nx/2+1,ny/2+1) = avg
!  u(nx/2+1,ny/2-1) = avg
!  u(nx/2-1,ny/2+1) = avg
!  u(nx/2-1,ny/2-1) = avg


end subroutine boundaries
