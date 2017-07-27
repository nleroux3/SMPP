! Estimates the refreezing of liquid water content based on Illangasekare et al., 1990
subroutine refreezing(volume_water_frozen, theta_1, theta_2, porosity, Vol, T_1, T_2, density, dry_density)
use declarations, only : rho_i, Cpi, rho_w, Cpw, rho_a, Cpa, dp, Lf
implicit none
real(dp), intent(inout) :: volume_water_frozen, &     ! Variable used to compute mass conservation [m3]
                       &   dry_density, &             ! Dry snow density [kg/m3]
                       &   density, &                 ! Bulk snow density [kg/m3]
                       &   porosity                   ! Porosity
real(dp), intent(in):: Vol, &             ! Volume of a cell [m3]
                   &   theta_1, &         ! Initial water content (before refreezing)
                   &   T_1                ! Initial temperature before refreezing [C]
real(dp), intent(out):: theta_2, &        ! Water content after refreezing
                    &   T_2               ! Temperature after refreezing [C]
real(dp) :: ice_content, & ! Ice content of each cell [m3/m3]
        &   rhoCp, &       ! rho*Cp for each snow cell [J/K/m3] 
        &   freeze_max     ! Maximum liquid water content that can refreeze, same units as water content



if (theta_1 .gt. 0._dp) then

   rhoCp = rho_i * Cpi * (1._dp - porosity) + rho_w * Cpw * theta_1 + &
      & rho_a * Cpa * (porosity - theta_1)

   freeze_max = (-rhoCp*T_1)/(Lf*rho_w)   ! same units as water content

   if (freeze_max .lt. theta_1) then  ! not all liquid water content refreezes

     T_2 = 0._dp
     theta_2 = theta_1-freeze_max
     volume_water_frozen = volume_water_frozen + freeze_max * Vol
     ice_content = (1.0_dp-porosity) + freeze_max*rho_w/rho_i   
     porosity = 1.0_dp - ice_content
     density = ice_content*rho_i + rho_w*theta_1 + rho_a * (porosity - theta_1)   
     dry_density = ice_content*rho_i  

   elseif (freeze_max .ge. theta_1) then ! all liquid water content refreezes,  i.e. freeze = water_content
      theta_2 = 0.0_dp
      volume_water_frozen = volume_water_frozen + theta_1 * Vol
      T_2 = theta_1* Lf * rho_w / rhoCp + T_1
      ice_content = (1.0_dp-porosity) + theta_1*rho_w/rho_i 
      porosity = 1.0_dp-ice_content
      density = ice_content*rho_i + rho_a * porosity
      dry_density = ice_content*rho_i

   endif
else 
   T_2 = T_1
   theta_2 = theta_1
endif
 

end subroutine
