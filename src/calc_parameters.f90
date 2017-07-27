! Estimates of the thermal and hydraulic parameters
subroutine calc_parameters
USE Declarations
implicit none
real(dp) :: viscosity, &  ! Kinematic viscosity of water [m2/s]
        &   k_w           ! Thermal conductivity of water [W/(K m)]

viscosity = 1.787e-6_dp

do j = 1,M-1
   do i = 1,N-1

      S(i,j) = max((water_content(i,j,1)-irreducible)/(porosity(i,j)-irreducible),0._dp)  

      k_w = (1.32_dp + 0.00559_dp * T(i,j,1) - 0.0000263*T(i,j,1)**2)/0.2389_dp

      !================ Effective thermal conductivity, (Calonne et al., 2011) [W/m K] ===============================
      k_eff(i,j) = (2.5e-6_dp*(dry_density(i,j)*dry_density(i,j))-1.23e-4_dp*dry_density(i,j)+0.024_dp) * &
                 & (1._dp - water_content(i,j,1)) + k_w * water_content(i,j,1)


      !================== Hydraulic conductivity of saturated snow (Calonne et al., 2012)  =====================
      Ks(i,j) = 9.81_dp / viscosity * 3.0_dp * (grain_opt(i,j)*0.5_dp)**2 &
         & * dexp(-0.013_dp*dry_density(i,j))        

      !================== Van Genuchten parameters  =========================================
      alpha(i,j) = 4.4e6_dp * (grain_classic(i,j)/dry_density(i,j))**(0.98_dp)    ! (Yamaguchi et al., 2012)
      nn(i,j) = 1.0_dp + 2.7e-3_dp * (grain_classic(i,j)/dry_density(i,j))**(-0.61_dp)
      mm(i,j) = 1.0_dp - 1.0_dp/nn(i,j)

      alpha_wetting(i,j) = 2._dp * alpha(i,j)
        
      !================ Hydraulic conductivity for unsaturated snow ===============================
      !============================= Mualem-VG model ===========================================
      K(i,j) = Ks(i,j)*dsqrt(S(i,j))*(1.0_dp-(1.0_dp-S(i,j)** &
         & (1.0_dp/mm(i,j)))**mm(i,j))**2

   enddo
enddo



! ================= Estimation of the hydraulic conductivity by arithmetic mean ====================== 
! =================== and arithmetic mean for effective thermal conductivity at the cell interfaces ==
if (M .gt. 2) then ! more than one snow layer
   do j = 1,M-1
      do i = 1,N-1
         if (i .ne. 1 .and. i .ne. N-1 .and. j .ne. 1 .and. j .ne. M-1) then

            !====================== Middle cells ========================
            K_south(i,j) = (K(i,j)+K(i,j-1))*(0.5_dp)
            K_north(i,j) = (K(i,j)+K(i,j+1))*(0.5_dp)
            K_east(i,j) = (K(i,j)+K(i+1,j))*(0.5_dp)
            K_west(i,j) = (K(i-1,j)+K(i,j))*(0.5_dp)

            k_eff_west(i,j) = (0.5_dp)*(k_eff(i,j)+k_eff(i-1,j))
            k_eff_east(i,j) = (0.5_dp)*(k_eff(i,j)+k_eff(i+1,j))
            k_eff_north(i,j) = (0.5_dp)*(k_eff(i,j)+k_eff(i,j+1))
            k_eff_south(i,j) = (0.5_dp)*(k_eff(i,j)+k_eff(i,j-1))
         endif

         if (i .ne. 1 .and. i .ne. N-1 .and. j .eq. 1) then  
            !============ Lower boundary except west and east cells ============
            K_east(i,1) = (K(i,1)+K(i+1,1))*(0.5_dp)
            K_west(i,1) = (K(i-1,1)+K(i,1))*(0.5_dp)
            K_north(i,1) = (K(i,1)+K(i,2))*(0.5_dp)
            K_south(i,1) = K(i,1)

            k_eff_west(i,1) = (0.5_dp)*(k_eff(i,1)+k_eff(i-1,1))
            k_eff_east(i,1) = (0.5_dp)*(k_eff(i,1)+k_eff(i+1,1))
            k_eff_north(i,1) = (0.5_dp)*(k_eff(i,1)+k_eff(i,2))
            k_eff_south(i,1) = k_eff(i,1)
         endif

         if (i .ne. 1 .and. i .ne. N-1 .and. j .eq. M-1) then   
            !======= upper boundary except west and east cells ================
            K_north(i,M-1) = K(i,M-1)
            K_south(i,M-1) = (K(i,M-1)+K(i,M-2))*(0.5_dp)
            K_east(i,M-1) = (K(i,M-1)+K(i+1,M-1))*(0.5_dp)
            K_west(i,M-1) = (K(i-1,M-1)+K(i,M-1))*(0.5_dp)

            k_eff_west(i,M-1) = (0.5_dp)*(k_eff(i,M-1)+k_eff(i-1,M-1))
            k_eff_east(i,M-1) = (0.5_dp)*(k_eff(i,M-1)+k_eff(i+1,M-1))
            k_eff_north(i,M-1) = k_eff(i,M-1)
            k_eff_south(i,M-1) = (0.5_dp)*(k_eff(i,M-1)+k_eff(i,M-2))
         endif

         if (j .ne. 1 .and. j .ne. M-1 .and. i .eq. N-1) then 
            !============ Right boundary except north and south cells ============
            if (choice_lat_BC .eq. 0) then
               K_east(N-1,j) = (K(N-1,j)+K(1,j))*0.5_dp
            else
               K_east(N-1,j) = K(N-1,j)
            endif
            K_west(N-1,j) = (K(N-2,j)+K(N-1,j))*(0.5_dp)
            K_south(N-1,j) = (K(N-1,j)+K(N-1,j-1))*(0.5_dp)
            K_north(N-1,j) = (K(N-1,j)+K(N-1,j+1))*(0.5_dp)

            k_eff_west(N-1,j) = (0.5_dp)*(k_eff(N-1,j)+k_eff(N-2,j))
            k_eff_east(N-1,j) = k_eff(N-1,j)
            k_eff_north(N-1,j) = (0.5_dp)*(k_eff(N-1,j)+k_eff(N-1,j+1))
            k_eff_south(N-1,j) = (0.5_dp)*(k_eff(N-1,j)+k_eff(N-1,j-1))
         endif

         if (j .ne. 1 .and. j .ne. M-1 .and. i .eq. 1) then
            !========= Left boundary except north and south cells ==================
            K_east(1,j) = (K(1,j)+K(2,j))*(0.5_dp)
            if (choice_lat_BC .eq. 0) then
               K_west(1,j) = (K(N-1,j)+K(1,j))*0.5_dp
            else
               K_west(1,j) = K(1,j)
            endif
            K_south(1,j) = (K(1,j)+K(1,j-1))*(0.5_dp)
            K_north(1,j) = (K(1,j)+K(1,j+1))*(0.5_dp)

            k_eff_west(1,j) = k_eff(1,j)
            k_eff_east(1,j) = (0.5_dp)*(k_eff(1,j)+k_eff(2,j))
            k_eff_north(1,j) = (0.5_dp)*(k_eff(1,j)+k_eff(1,j+1))
            k_eff_south(1,j) = (0.5_dp)*(k_eff(1,j)+k_eff(1,j-1))
         endif

         if (j .eq. 1 .and. i .eq. 1) then
            !=========== Lower left cell =========================
            K_east(1,1) = (K(1,1)+K(2,1))*(0.5_dp)
            K_north(1,1) = (K(1,1)+K(1,2))*(0.5_dp)
            if (choice_lat_BC .eq. 0) then
               K_west(1,1) = (K(N-1,1)+K(1,1))*0.5_dp
            else
               K_west(1,1) = K(1,1)
            endif
            K_south(1,1) = K(1,1)

            k_eff_west(1,1) = k_eff(1,1)
            k_eff_east(1,1) = (0.5_dp)*(k_eff(1,1)+k_eff(2,1))
            k_eff_north(1,1) = (0.5_dp)*(k_eff(1,1)+k_eff(1,2))
            k_eff_south(1,1) = k_eff(1,1)
         endif

        if (j .eq. M-1 .and. i .eq. 1) then
            !========== Upper left cell ====================
            K_south(1,M-1) = (K(1,M-1)+K(1,M-2))*(0.5_dp)
            K_east(1,M-1) = (K(1,M-1)+K(2,M-1))*(0.5_dp)
            if (choice_lat_BC .eq. 0) then
               K_west(1,M-1) = (K(N-1,M-1)+K(1,M-1))*0.5_dp
            else
               K_west(1,M-1) = K(1,M-1)
            endif
            K_north(1,M-1) = K(1,M-1)

            k_eff_west(1,M-1) = k_eff(1,M-1)
            k_eff_east(1,M-1) = (0.5_dp)*(k_eff(1,M-1)+k_eff(2,M-1))
            k_eff_north(1,M-1) = k_eff(1,M-1)
            k_eff_south(1,M-1) = (0.5_dp)*(k_eff(1,M-1)+k_eff(1,M-2))
         endif

         if (j .eq. 1 .and. i .eq. N-1) then 
            !=========== Lower right cell=======================
            K_north(N-1,1) = (K(N-1,1)+K(N-1,2))*(0.5_dp)
            K_west(N-1,1) = (K(N-2,1)+K(N-1,1))*(0.5_dp)
            if (choice_lat_BC .eq. 0) then
               K_east(N-1,1) = (K(N-1,1)+K(1,1))*(0.5_dp)
            else
               K_east(N-1,1) = K(N-1,1)
            endif
            K_south(N-1,1) = K(N-1,1)

            k_eff_west(N-1,1) = (0.5_dp)*(k_eff(N-1,1)+k_eff(N-2,1))
            k_eff_east(N-1,1) = k_eff(N-1,1)
            k_eff_north(N-1,1) = (0.5_dp)*(k_eff(N-1,1)+k_eff(N-1,2))
            k_eff_south(N-1,1) = k_eff(N-1,1)
         endif

         if (j .eq. M-1 .and. i .eq. N-1) then
            !======= Upper right cell =================
            K_south(N-1,M-1) = (K(N-1,M-1)+K(N-1,M-2))*(0.5_dp)
            K_west(N-1,M-1) = (K(N-2,M-1)+K(N-1,M-1))*(0.5_dp)
            if (choice_lat_BC .eq. 0) then
               K_east(N-1,M-1) = (K(N-1,M-1)+K(1,M-1))*(0.5_dp)
            else
               K_east(N-1,M-1) = K(N-1,M-1)
            endif
            K_north(N-1,M-1) = K(N-1,M-1)

            k_eff_west(N-1,M-1) = (0.5_dp)*(k_eff(N-1,M-1)+k_eff(N-2,M-1))
            k_eff_east(N-1,M-1) = k_eff(N-1,M-1)
            k_eff_north(N-1,M-1) = k_eff(N-1,M-1)
            k_eff_south(N-1,M-1) = (0.5_dp)*(k_eff(N-1,M-1)+k_eff(N-1,M-2))
         endif
      enddo
   enddo

else   ! M = 2, i.e. one layer
   !======== right cell =====================
   K_north(N-1,1) = K(N-1,1)
   K_west(N-1,1) = (K(N-2,1)+K(N-1,1))*(0.5_dp)
   if (choice_lat_BC .eq. 0) then
      K_east(N-1,1) = (K(N-1,1)+K(1,1))*0.5_dp
   else
      K_east(N-1,1) = K(N-1,1)
   endif
   K_south(N-1,1) = K(N-1,1)

   k_eff_west(N-1,1) = (0.5_dp)*(k_eff(N-1,1)+k_eff(N-2,1))
   k_eff_east(N-1,1) = k_eff(N-1,1)
   k_eff_north(N-1,1) = k_eff(N-1,1)
   k_eff_south(N-1,1) = k_eff(N-1,1)

   ! ============= left cell ========================
   K_east(1,1) = (K(1,1)+K(2,1))*(0.5_dp)
   K_north(1,1)= K(1,1)
   if (choice_lat_BC .eq. 0) then
      K_west(1,1) = (K(N-1,1)+K(1,1))*0.5_dp
   else
      K_west(1,1) = K(1,1)
   endif
   K_south(1,1) = K(1,1)

   k_eff_west(1,1) = k_eff(1,1)
   k_eff_east(1,1) = (0.5_dp)*(k_eff(1,1)+k_eff(2,1))
   k_eff_north(1,1) = k_eff(1,1)
   k_eff_south(1,1) = k_eff(1,1)

   !================ Middle cells ===================== 
   K_east(2:N-2,1) = (K(2:N-2,1)+K(3:N-1,1))*(0.5_dp)
   K_west(2:N-2,1) = (K(1:N-3,1)+K(2:N-2,1))*(0.5_dp)
   K_north(2:N-2,1) = K(2:N-2,1)
   K_south(2:N-2,1) = K(2:N-2,1)

   k_eff_west(2:N-2,1) = (0.5_dp)*(k_eff(2:N-2,1)+k_eff(1:N-3,1))
   k_eff_east(2:N-2,1) = (0.5_dp)*(k_eff(2:N-2,1)+k_eff(3:N-1,1))
   k_eff_north(2:N-2,1) = k_eff(2:N-2,1)
   k_eff_south(2:N-2,1) = k_eff(2:N-2,1)

endif


RETURN
end
