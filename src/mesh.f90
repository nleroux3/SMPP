! Reads the inputs to create the mesh
subroutine mesh
use Declarations 
real(dp) :: delta_x, &    ! Horizontal space step [m]
        &   delta_y       ! Vertical space ste [m]
integer :: choice_nb_layers ! Variable to decide if vertical space step is equal to layer thicknesses (=1) or not (=0)
real(dp), allocatable :: y(:), &        ! Variable used to calculate yn
                     &   thickness(:)   ! Variable used to calculate yn

choice_nb_layers = 0

open(7,file='inputs.txt') 
read(7,*) L1
read(7,*) L2
read(7,*) N
read(7,*) M
read(7,*) 
read(7,*) 
read(7,*) beta
read(7,*) nb_layers
allocate (thickness(nb_layers))
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
read(7,*)
do j = 1,nb_layers
   read(7,*) thickness(j)
enddo

close (7) 

if (M .eq. 999) then   ! use layer thicknesses for vertical mesh resolution
   M = nb_layers + 1
   choice_nb_layers = 1 
endif

!================== Allocate mesh variables =================================
allocate (y(M),xn(N,M),yn(N,M),xv(N-1,M-1),yv(N-1,M-1), Vol(N-1,M-1))
allocate (Sx_west(N-1,M-1),Sy_west(N-1,M-1),Sx_east(N-1,M-1),Sy_east(N-1,M-1))
allocate (Sx_north(N-1,M-1),Sy_north(N-1,M-1),Sx_south(N-1,M-1),Sy_south(N-1,M-1)) 


y(1) = 0._dp

if (choice_nb_layers .eq. 0) then
   delta_y = L2 / real(M-1,dp)
   do j = 1,M-1
      y(j+1) = y(j) + delta_y
   enddo
else 
   do j=1,M-1
      y(j+1) = y(j)+thickness(j)
   enddo
endif

do j = 1,M
   xn(1,j) = dsin(beta)*y(j)
   yn(1,j) = dcos(beta)*y(j)
enddo

delta_x = L1 / real(N-1,dp) 

do j = 1,M
   do i = 2,N
      xn(i,j) = xn(1,j) + real((i-1),dp) * delta_x * dcos(beta)
      yn(i,j) = yn(1,j) - dsin(beta)*real(i-1,dp)*delta_x
   enddo
enddo
do j = 1,M-1
   do i = 1,N-1
      xv(i,j) = (xn(i,j)+xn(i+1,j)+xn(i,j+1)+xn(i+1,j+1))*0.25_dp
      yv(i,j) = (yn(i,j)+yn(i+1,j)+yn(i,j+1)+yn(i+1,j+1))*0.25_dp
      Vol(i,j) = (yn(i,j+1)-yn(i,j))*delta_x/dcos(beta)
      Sx_west(i,j) = dabs(yn(i,j+1)-yn(i,j))
      Sy_west(i,j) = dabs(xn(i,j+1)-xn(i,j))
      Sx_east(i,j) = dabs(yn(i+1,j+1)-yn(i+1,j))
      Sy_east(i,j) = dabs(xn(i+1,j+1)-xn(i+1,j))
      Sx_north(i,j) = dabs(yn(i+1,j+1)-yn(i,j+1))
      Sy_north(i,j) = dabs(xn(i+1,j+1)-xn(i,j+1))
      Sx_south(i,j) = dabs(yn(i+1,j)-yn(i,j))
      Sy_south(i,j) = dabs(xn(i+1,j)-xn(i,j))
    enddo
enddo


! Check if sum of all small volumes is equal to rectangle volume
write(*,*) 'Sum of all small volumes [m2] = ', sum(Vol)
write(*,*) 'Total volume [m2] = ', L1 * L2

deallocate(y, thickness)

return
end

