!Calculates the time step to assure model stability
subroutine timestep
use Declarations
real(dp) :: dt1, &      ! Dummy time step [s]
        &   dt2, &      ! Dummy time step [s] 
        &   dt3, &      ! Dummy time step [s]
        &   D_south, &  ! Diffusivity coefficient for water flow in y-direction (down) [m2/s]
        &   D_east, &   ! Diffusivity coefficient for water flow in x-direction [m2/s]
        &   D_north, &  ! Diffusivity coefficient for water flow in y-direction (up) [m2/s]
        &   uu, &       ! Convective velocity in water flow equation [m/s]
        &   C, &        ! Parameter in Richards equation
        &   thermal_diffusivity_south, & ! Diffusivity coefficient for heat equation in y-direction (down)
        &   thermal_diffusivity_north, & ! Diffusivity coefficient for heat equation in y-direction (up)
        &   thermal_diffusivity_east     ! Diffusivity coefficien
real(dp) :: Co                           ! Courant number, should be lower than 1 

Dt2 = Dt_ini
dt3 = Dt_ini

if (beta .le. 0.0001) then !  flat snowpack
      do j = 1,M-1
         do i = 1,N-1
            Co = 0.01_dp

            !===================== CFL heat transfer ==================================
            thermal_diffusivity_south = k_eff_south(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))

            thermal_diffusivity_north = k_eff_north(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))

               
            thermal_diffusivity_east = k_eff_east(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))

            if (j .ne. 1 .and. i .ne. N-1  .and. j .ne. M-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) , &
                  & thermal_diffusivity_east * Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) ))
            elseif (j .eq. 1  .and. i .ne. N-1 ) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)  , &
                  & thermal_diffusivity_east * Sx_east(i,j)/ dabs(xv(i+1,j)-xv(i,j)), &
                  & thermal_diffusivity_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) ))
            elseif (j .eq. M-1  .and. i.ne. N-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_north * Sy_north(i,j)/ (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                  & thermal_diffusivity_east * Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)), &
                  & thermal_diffusivity_south * Sy_south(i,j)/ dabs(yv(i,j)-yv(i,j-1))))
            elseif (j .ne. 1 .and. i .eq. N-1  .and. j .ne. M-1) then
               dt3 = Co * Vol(i,j) *0.5_dp / (thermal_diffusivity_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)))
            elseif (j .eq. 1  .and. i.eq. N-1 ) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp), &
                 & thermal_diffusivity_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) ))
            elseif (j .eq. M-1  .and. i .eq. N-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_north * Sy_north(i,j)/ (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                 & thermal_diffusivity_south * Sy_south(i,j)/ dabs(yv(i,j)-yv(i,j-1))) )
            endif




            !======================= CFL water flow ========================================

           if (S(i,j) .gt. 0.) then
               uu = (2.0_dp * Ks(i,j) * ((1.0_dp-S(i,j)**(1.0_dp/mm(i,j)))**mm(i,j) - 1.0_dp) * (1.0_dp-S(i,j)** &
                  & (1.0_dp/mm(i,j)))**(mm(i,j)-1.0_dp) * S(i,j)**0.5_dp * S(i,j)**(1.0_dp/mm(i,j)-1.0_dp)) / &
                  & (irreducible-porosity(i,j)) - (Ks(i,j) * ((1.0_dp-S(i,j)**(1.0_dp/mm(i,j)))**mm(i,j) - 1.0_dp)**2.0_dp) / &
                  & (2.0_dp*(irreducible-porosity(i,j)) * S(i,j)**0.5_dp)
            else    
               uu = 0._dp
            endif
 
            
            if (wrc(i,j) .eq. 1 .or. wrc(i,j) .eq. 4 ) then
               C = (porosity(i,j)-irreducible) * alpha_wetting(i,j)**nn(i,j) * nn(i,j) * mm(i,j) * dabs(P(i,j))** &
                  &   (nn(i,j)-1.0_dp) * (1.0_dp + dabs(alpha_wetting(i,j)*P(i,j))**nn(i,j))**(-mm(i,j) - 1.0_dp)
            else
               C = (porosity(i,j)-irreducible) * alpha(i,j)**nn(i,j) * nn(i,j) * mm(i,j) * dabs(P(i,j))** & 
                  & (nn(i,j)-1.0_dp) * (1.0_dp + dabs(alpha(i,j)*P(i,j))**nn(i,j))**(-mm(i,j) - 1.0_dp)
            endif

            D_south = K_south(i,j) / C ! grid Fourier number El-Kadi and Ling, 1993
            D_east =  K_east(i,j) / C ! grid Fourier number El-Kadi and Ling, 1993
            D_north = K_north(i,j) / C ! grid Fourier number El-Kabi and Ling, 1993

            if (j .ne. 1 .and. j .ne. M-1) then
               if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) , &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)) , &
                     & uu * Sy_south(i,j) *0.5_dp )) 
               else
                  dt1 = Co * Vol(i,j) *0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) , &
                     & uu * Sy_south(i,j) * 0.5_dp))
               endif
            elseif (j .eq. 1) then
               if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)) , &
                     & uu * Sy_south(i,j) * 0.5_dp, D_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) ))  
               else 
                  dt1 = Co * Vol(i,j) *0.5_dp / (max(D_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                     & uu * Sy_south(i,j) * 0.5_dp , D_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)))) 
               endif
            elseif (j .eq. M-1) then
              if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) , &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)) , &
                     & uu * Sy_south(i,j) * 0.5_dp, D_north * Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)))  
              else 
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_north * Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                     & uu * Sy_south(i,j) * 0.5_dp, D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)))) 
              endif
            endif   
              Dt2 = min(+Dt2,+dt1,+dt3)

          enddo
      enddo
  
     

else ! ==================== BETA .gt. 0 =======================================

      do j = 1,M-1
         do i = 1,N-1
             
             Co = 0.9_dp

            !===================== CFL heat transfer ==================================
            thermal_diffusivity_south = k_eff_south(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))

            thermal_diffusivity_north = k_eff_north(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))
               
            thermal_diffusivity_east = k_eff_east(i,j) / (rho_i*Cpi*(1.0_dp-porosity(i,j)) + &
               & water_content(i,j,1)*rho_w*Cpw + &
               & rho_a*Cpa*(porosity(i,j)-water_content(i,j,1)))

            if (j .ne. 1 .and. i .ne. N-1  .and. j .ne. M-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) , &
                  & thermal_diffusivity_east * Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) ))
            elseif (j .eq. 1  .and. i .ne. N-1 ) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                  & thermal_diffusivity_east * Sx_east(i,j)/ dabs(xv(i+1,j)-xv(i,j)) , &
                  & thermal_diffusivity_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) ))
            elseif (j .eq. M-1  .and. i.ne. N-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_north * Sy_north(i,j)/ (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                  & thermal_diffusivity_east * Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)), &
                  &  thermal_diffusivity_south * Sy_south(i,j)/ dabs(yv(i,j)-yv(i,j-1))))
            elseif (j.ne. 1 .and. i.eq. N-1  .and. j.ne.M-1) then
               dt3 = Co * Vol(i,j) *0.5_dp / (thermal_diffusivity_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)))
            elseif (j .eq. 1  .and. i.eq. N-1 ) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp), &
                  & thermal_diffusivity_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j))))
            elseif (j .eq. M-1  .and. i .eq. N-1) then
               dt3 = Co * Vol(i,j) * 0.5_dp / (max(thermal_diffusivity_north * Sy_north(i,j)/ (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &  
                  & thermal_diffusivity_south * Sy_south(i,j)/ dabs(yv(i,j)-yv(i,j-1))))
            endif


            !======================= CFL water flow ========================================

            if (S(i,j) .gt. 0.) then
               uu = (2.0_dp * Ks(i,j) * ((1.0_dp-S(i,j)**(1.0_dp/mm(i,j)))**mm(i,j) - 1.0_dp) * (1.0_dp-S(i,j)** &
                  & (1.0_dp/mm(i,j)))**(mm(i,j)-1.0_dp) * S(i,j)**0.5_dp * S(i,j)**(1.0_dp/mm(i,j)-1.0_dp)) / &
                  & (irreducible-porosity(i,j)) - (Ks(i,j) * ((1.0_dp-S(i,j)**(1.0_dp/mm(i,j)))**mm(i,j) - 1.0_dp)**2.0_dp) / &
                  & (2.0_dp*(irreducible-porosity(i,j)) * S(i,j)**0.5_dp)
            else    !S(i,j)=1e-6 is used to avoid model divergence
               uu = 0._dp
            endif


            if (wrc(i,j) .eq. 1 .or. wrc(i,j) .eq. 4) then
               C = (porosity(i,j)-irreducible) * 2._dp * alpha(i,j)**nn(i,j) * nn(i,j) * mm(i,j) * dabs(P(i,j))** &
                  &   (nn(i,j)-1.0_dp) * (1.0_dp + dabs(2._dp * alpha(i,j)*P(i,j))**nn(i,j))**(-mm(i,j) - 1.0_dp)
            else
               C = (porosity(i,j)-irreducible) * alpha(i,j)**nn(i,j) * nn(i,j) * mm(i,j) * dabs(P(i,j))** & 
                  & (nn(i,j)-1.0_dp) * (1.0_dp + dabs(alpha(i,j)*P(i,j))**nn(i,j))**(-mm(i,j) - 1.0_dp)
            endif
            D_south = K_south(i,j) / C ! grid Fourier number El-Kadi and Ling, 1993
            D_east =  K_east(i,j) / C ! grid Fourier number El-Kadi and Ling, 1993
            D_north = K_north(i,j) / C ! grid Fourier number El-Kabi and Ling, 1993

            if (j .ne. 1 .and. j .ne. M-1) then
               if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)), &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)), &
                     & uu * Sy_south(i,j) *0.5_dp )) 
               else
                  dt1 = Co * Vol(i,j) *0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)), &
                     & uu * Sy_south(i,j) * 0.5_dp))
               endif
            elseif (j .eq. 1) then
              if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp), &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)), &
                     & uu * Sy_south(i,j) * 0.5_dp , &
                     & D_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j))))  
              else 
                  dt1 = Co * Vol(i,j) *0.5_dp / (max(D_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) , &
                     & uu * Sy_south(i,j) * 0.5_dp , &
                     & D_north * Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) )) 
              endif
            elseif (j .eq. M-1) then
              if (i .lt. N-1) then
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_south * Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)), &
                     & D_east * Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)), &
                     & uu * Sy_south(i,j) * 0.5_dp, & 
                     &  D_north * Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)))  
              else 
                  dt1 = Co * Vol(i,j) * 0.5_dp / (max(D_north * Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp), &
                     & uu * Sy_south(i,j) * 0.5_dp, & 
                     & D_south * Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)) ) 
              endif
            endif   

              Dt2 = min(+Dt2,+dt1,+dt3)

          enddo
      enddo

endif


Dt = Dt2
Dt = max(Dt, 1e-3_dp)


RETURN
end
