! Computes water content and temperature every time step
subroutine SMPP
use declarations
implicit none
real(dp) :: mass_conservation, &              ! Mass conservation coefficient (=1 if mass is conserved)
        &   volume_water, &                   ! Volume of liquid water within the whole snowpack [m3]
        &   time , &                          ! Time since beginning of simulation [s]
        &   Tss_ini(N-1), &                   ! Initial snow surface temperature at each time step
        &   volume_water_ini, &               ! Initial volume of liquid water at time = 0 [m3]
        &   volume_water_frozen               ! Volume of liquid water that refreezes during each time step [m3]

volume_water_frozen = 0._dp
mass_conservation = 0._dp
time = 0._dp
outflow  = 0.0_dp 
sum_fluxes = 0.0_dp  
volume_water_ini = sum(water_content(:,1:M-1,1) * Vol(:,1:M-1))

!======= Convert classical grain size into optical grain size ==============
call conversion_grain

!======== Create the initial output files ==============================
call write_outputs(time,mass_conservation)

!============================ Main iterations ==========================================
do while (time .lt. tf)

   Tss_ini = Tss

   !========= Convert classical grain size into optical grain size ================
   call conversion_grain

   !============ Compute of hydraulic and thermal parameters ====================
   call calc_parameters

   !============ Compute time step (Dt) with CFL condition ===========================
   call timestep

   time = time+Dt

   !========= Compute of snowpack thinning and infiltration rate for melt case ==================
   call thinning_melt(Tss_ini,time,mass_conservation)

   !========= Compute mass and energy fluxes at the boundaries ====================
   call boundaries 

   !========= Compute mass and energy fluxes in internal cells ====================
   if (M .gt. 2) then
      if (beta .gt. 0.0001) then  ! sloping snowpack
         do j = 2,M-2
            do i = 2,N-2
               !=========== Mass flux due to water pressure gradient ===================
               Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
                  & ( Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)) + &
                  & Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j))) 

               Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
                  & ( Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) + &
                  & Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j))) 

               Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
                  & ( Sx_south(i,j) / dabs(xv(i,j-1)-xv(i,j)) + &
                  & Sy_south(i,j) / dabs(yv(i,j-1)-yv(i,j))) 
   
               Fpresn(i,j) = K_north(i,j) * (P(i,j+1)-P(i,j)) * &
                  & ( Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)) + &
                  & Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j))) 


               !============ Mass flux due to gravity force =============================
               Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)                                     
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)                                       
               Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)              
               Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

               !=========================== Heat flux =================================== 
               Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
                  & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
                  & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))

               Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
                  & (Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) + &
                  & Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)))

               Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
                  & (Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)) + &
                  & Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)))
   
               Fcondw(i,j) = k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) * &
                  & (Sy_west(i,j) / dabs(yv(i-1,j)-yv(i,j)) + &
                  & Sx_west(i,j) / dabs(xv(i-1,j)-xv(i,j)))

            enddo    
         enddo
      else   ! beta=0, i.e. flat ground
         do j = 2,M-2
            do i = 2,N-2
               !================ Mass flux due to water pressure gradients =========================
               Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * &
                  & (P(i-1,j)-P(i,j)) / dabs(xv(i,j)-xv(i-1,j))

               Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * &
                  & (P(i+1,j)-P(i,j)) / dabs(xv(i+1,j)-xv(i,j)) 

               Fpress(i,j) = Sy_south(i,j) * K_south(i,j) * &
                  & (P(i,j-1)-P(i,j)) / dabs(yv(i,j)-yv(i,j-1))

               Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * &
                  & (P(i,j+1)-P(i,j)) / dabs(yv(i,j+1)-yv(i,j))

   
               !================== Mass flux due to gravity force ===========================
               Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)                                     
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)                                       
               Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)              
               Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

               !====================== Heat fluxes =========================================
               Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * &
                  & (T(i,j-1,1)-T(i,j,1)) / dabs(yv(i,j)-yv(i,j-1))

               Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * &
                  & (T(i,j+1,1)-T(i,j,1)) / dabs(yv(i,j)-yv(i,j+1))

               Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * &
                  & (T(i+1,j,1)-T(i,j,1)) / dabs(xv(i,j)-xv(i+1,j))

               Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * &
                  & (T(i-1,j,1)-T(i,j,1)) / dabs(xv(i,j)-xv(i-1,j))

            enddo
         enddo
      endif
   endif


   !========== Condition to avoid mass flow from dry to wet cells ==========
   call condition_Darcy

   !========== Compute mass and heat fluxes coming into each cell
   do j = 1,M-1
      do i = 1,N-1
            Fpres(i,j) = Fpresn(i,j) + Fpress(i,j) + Fpresw(i,j) + Fprese(i,j)
            Fgrav(i,j) = Fgravn(i,j) + Fgravs(i,j) + Fgravw(i,j) + Fgrave(i,j)
            Fcond(i,j) = Fcondw(i,j) + Fconde(i,j) + Fconds(i,j) + Fcondn(i,j) 
      enddo
   enddo

   !============= Solve mass and heat equations using FVM =========================
   call solver_mass_heat(volume_water_frozen)
   
   !============= Check if water content computed is NaN ===========================
   call check_nan(time, mass_conservation)

   !============= Compute total volume of liquid water and mass conservation ================
   sum_fluxes = sum_fluxes + (sum(Fgravs(:,1) + Fgravn(:,M-1)) + sum(Fgrave(N-1,1:M-1) + Fgravw(1,1:M-1))) * Dt
   volume_water = sum(water_content(:,1:M-1,1) * Vol(:,1:M-1)) - volume_water_ini + volume_water_frozen
   mass_conservation = volume_water / sum_fluxes

 
   !========================= Write ouputs in files ===================================
   call write_outputs(time,mass_conservation)


enddo

!================== write the final outputs ===========================
call write_outputs(time,mass_conservation)

return
end
