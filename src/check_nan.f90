! Checks if model diverged and writes some outputs to try to locate error
subroutine check_nan(time, mass_conservation)
use Declarations
implicit none
real(dp), intent(in) :: time, &             ! Time since beginning of simulation [s]
                    &   mass_conservation   ! Mass conservation index
integer :: flag = 0


!=============== Check if saturation is greater than 1 or theta less than 0 =======================
do j = 1,M-1
   do i = 1,N-1
      if (S(i,j) .gt. 1._dp .or. water_content(i,j,1) .lt. 0.) then
         flag = 1
         exit
      endif
   enddo
enddo
  
!============ If water content = NaN or  flag == 1 ---> stop program and plot some outputs ==============
if (isnan(sum(water_content(:,1:M-1,1))) .or. flag .eq. 1 ) then
   open(99,file='./outputs/water_nan.dat')
   write(99,*) water_content(:,1:M-1,1)
   close(99)

   open(99,file='./outputs/P_nan.dat')
   write(99,*) P(:,1:M-1)
   close(99)

   open(99,file='./outputs/Fgrav_nan.dat')
   write(99,*) Fgrav(:,1:M-1)
   close(99)

   open(99,file='./outputs/Fgravn_nan.dat')
   write(99,*) Fgravn(:,1:M-1)
   close(99)

   open(99,file='./outputs/Fgravs_nan.dat')
   write(99,*) Fgravs(:,1:M-1)
   close(99)

   open(99,file='./outputs/K_south_nan.dat')
   write(99,*) K_south(:,1:M-1)
   close(99)

   open(99,file='./outputs/K_north_nan.dat')
   write(99,*) K_north(:,1:M-1)
   close(99)

   open(99,file='./outputs/Q_in_nan.dat')
   write(99,*) Qin
   close(99)

   open(99,file='./outputs/Sy_north.dat')
   write(99,*) Sy_north(:,1:M-1)
   close(99)

   open(99,file='./outputs/Dt.dat')
   write(99,*) Dt
   close(99)

   open(99,file='./outputs/grain_classic.dat')
   write(99,*) grain_classic(:,1:M-1)
   close(99)
   
   write(*,*) 'NaN'

   Vol = 0._dp

   call write_outputs(time,mass_conservation)
   stop 

endif

end subroutine
