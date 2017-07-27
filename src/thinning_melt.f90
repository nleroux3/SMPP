! Estimates snowpack thinning and melt at the surface
subroutine thinning_melt(Tss_ini,time,mass_conservation)
use declarations
implicit none
real(dp), intent(in) :: time, &                 ! Time since beginning of simulation [s]
                    &   mass_conservation, &    ! Mass conservation parameter (=1 if mass is conserved)
                    &   Tss_ini(N-1)            ! Initial snow surface temperature at the beginning of the time step
real(dp) :: Vn(N-1), &                          ! Velocity of melting surface [m/s]
        &   delta_x                             ! Spatial step in x-direction
integer :: flag                                 ! Dummy variable used in this script identifying if top layers merge

delta_x = L1/real(N-1,dp)



!============ Computation of velocity of moving boundary =================
where (dabs(Hn) .gt. 0.0_dp) 
   Tss= dabs(yn(1:N-1,M)-yn(1:N-1,M-1))/(2.0_dp*k_eff_north(:,M-1))*Hn+T(:,M-1,1)
   where (Tss .ge. 0.0_dp .and. Tss_ini .lt. 0.0_dp)  ! not all Hn is used to melt, some is used to increase the surface temperature up to 0C
      Tss = 0.0_dp
      Hn = Hn - (0._dp - T(:,M-1,1)) * 2._dp * k_eff_north(:,M-1) / dabs(yn(1:N-1,M)-yn(1:N-1,M-1))
      Vn = Hn/(Lf*dry_density(:,M-1))
   end where
   where (Tss .ge. 0.0_dp .and. Tss_ini .ge. 0.0_dp) 
      Tss = 0.0_dp
      Vn = Hn/(Lf*dry_density(:,M-1))
   end where
   where (Tss .lt. 0._dp)
      Vn = 0._dp
   endwhere
elsewhere  ! Hn .eq. 0
   Vn = 0.0_dp 
end where

!============= Update mesh =========================
if (sum(Vn) .gt. 0.0_dp) then

   ! Infiltration rate 
   Qin = Vn*(dry_density(:,M-1)/rho_w+water_content(:,M-1,1))  

   ! Change of snowdepth
   L2 = L2-Vn(1)*Dt 
   
   ! Moving mesh
   xn(1,M) = xn(1,M)-Vn(1)*Dt*dsin(beta)
   yn(1,M) = yn(1,M)-Vn(1)*Dt*dcos(beta)  
   
   xn(N,M) = xn(1,M)+real(N-1,dp)*delta_x*dcos(beta)
   yn(N,M) = yn(N,M)-Vn(N-1)*Dt*dcos(beta)
   
   xn(2:N-1,M) = xn(1:N-2,M)+delta_x*dcos(beta)
   yn(2:N-1,M) = yn(2:N-1,M)-(Vn(1:N-2)+Vn(2:N-1))*0.5_dp*Dt*dcos(beta)
  
   j = M-1
   xv(1:N-1,j) = (xn(1:N-1,j)+xn(2:N,j)+xn(1:N-1,j+1)+xn(2:N,j+1))*0.25_dp
   yv(1:N-1,j) = (yn(1:N-1,j)+yn(2:N,j)+yn(1:N-1,j+1)+yn(2:N,j+1))*0.25_dp
   Sx_west(1:N-1,j) = dabs(yn(1:N-1,j+1)-yn(1:N-1,j))
   Sy_west(1:N-1,j) = dabs(xn(1:N-1,j+1)-xn(1:N-1,j))
   Sx_east(1:N-1,j) = dabs(yn(2:N,j+1)-yn(2:N,j))
   Sy_east(1:N-1,j) = dabs(xn(2:N,j+1)-xn(2:N,j))
   Sx_north(1:N-1,j) = dabs(yn(2:N,j+1)-yn(1:N-1,j+1))
   Sy_north(1:N-1,j) = dabs(xn(2:N,j+1)-xn(1:N-1,j+1))
   Vol(1:N-1,j) = dabs(yn(1:N-1,j+1)-yn(1:N-1,j))*delta_x/dcos(beta)

else

   Qin = 0.0_dp

endif

!========= Melt upper layers if yn(i,j+1) <= yn(i,j) + 1mm ===============
flag=0
do i=1,N
   if (yn(i,M) .le. yn(i,M-1)+0.001_dp .and. M .gt. 2) then
      flag = 1
   endif
enddo


if (flag .eq. 1) then

   yn(:,M-1) = yn(:,M)
   xn(:,M-1) = xn(:,M)

   Vn = Vn + 1e-3_dp / Dt

   Qin = Vn*(dry_density(:,M-1)/rho_w+water_content(:,M-1,1))  
   M = M-1
endif


!======== Stop model if snowdepth <= 2cm =================
if (yn(1,M)-yn(1,1).le. 0.02_dp) then

   outflow = outflow + density(:,M-1) / rho_w * Vol(:,M-1) + Qin
   yn(:,M) = yn(:,1)
   xn(:,M) = xn(:,1)

   xv(:,M-1) = 0.0_dp
   yv(:,M-1) = 0.0_dp
   Sx_west(:,M-1) = 0.0_dp
   Sy_west(:,M-1) = 0.0_dp
   Sx_east(:,M-1) = 0.0_dp
   Sy_east(:,M-1) = 0.0_dp
   Sx_north(:,M-1) = 0.0_dp
   Sy_north(:,M-1) = 0.0_dp
   Sx_south(:,M-1) = 0.0_dp
   Sy_south(:,M-1) = 0.0_dp
   Vol(:,M-1) = 0.0_dp

   call write_outputs(time,mass_conservation)

 endif   





return
end
