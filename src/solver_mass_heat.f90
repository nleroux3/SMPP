! Solves Richards equation and heat equation using FVM
subroutine solver_mass_heat(volume_water_frozen)
USE Declarations
implicit none
real(dp) :: rhoCp, &      ! rho*Cp for each snow cell [J/K/m3] !
        &   P_new, &      ! New pressure computed at the end of the time step [m]
        &   theta, &      ! New water content after solving Richards equation
        &   theta_plus, & ! Extra water content if theta > porosity
        &   theta_plus_refreeze ! Extra water content if theta > porosity after refreezing
real(dp) :: volume_water_frozen ! Volume of water that refroze (for mass conservation) [m3]


theta_plus_refreeze = 0._dp
theta_plus = 0._dp

!============== Internal melt and update of snow properties ==============
do j = M-1,1,-1
   do i = 1,N-1

      rhoCp = rho_i * Cpi * (1._dp - porosity(i,j)) + rho_w * Cpw * water_content(i,j,1) + &
         & rho_a * Cpa * (porosity(i,j) - water_content(i,j,1))
 

     !============= Solving heat transfer equation ============================
      T(i,j,2) = T(i,j,1) + Dt/Vol(i,j) * Fcond(i,j) / rhoCp
      T(i,j,1) = T(i,j,2)

     !================ Solving mass flow equation  ==============================
     theta = water_content(i,j,1) + Dt/Vol(i,j) * (Fgrav(i,j)+Fpres(i,j))

     !============== UPDATE WATER CONTENT AND DENSITY ===========================
     density(i,j) = density(i,j) - rho_w*water_content(i,j,1)  + rho_w*theta

     !============= Compute refreezing of liquid water content ===========================
     call refreezing(volume_water_frozen, theta, water_content(i,j,2), &
             & porosity(i,j), Vol(i,j), T(i,j,1), T(i,j,2), density(i,j), dry_density(i,j))
      
      T(i,j,1) = T(i,j,2)
      theta_plus_refreeze = max(0._dp, water_content(i,j,2)-porosity(i,j))
      water_content(i,j,2) = min(water_content(i,j,2),porosity(i,j))

     !============================= Hysteresis ===============================
     call hysteresis(wrc(i,j), P_new, P(i,j), alpha(i,j), alpha_wetting(i,j),&
             &  porosity(i,j), water_content(i,j,1), water_content(i,j,2), nn(i,j), mm(i,j), Pdi(i,j), Pid(i,j), &
             &  theta_id(i,j,:), theta_di(i,j,:), grain_classic(i,j), theta_plus,theta_s(i,j,:), &
             &  theta_r(i,j,:), order_scanning(i,j))

     !========== If theta > porosity, excess theta goes to cell below through Fgravs(i,j) and Fgravn(i,j-1) ============
     if (theta_plus .gt. 0._dp)  then
        if (j .gt. 2) then
            Fgravs(i,j) = Fgravs(i,j) -  theta_plus * Vol(i,j)  / Dt
            Fgrav(i,j) = Fgravn(i,j) + Fgravs(i,j) + Fgravw(i,j) + Fgrave(i,j)

            Fgravn(i,j-1) = - Fgravs(i,j)
            Fgrav(i,j-1) = Fgravn(i,j-1) + Fgravs(i,j-1) + Fgravw(i,j-1) + Fgrave(i,j-1)
        else
            Fgravs(i,j) = Fgravs(i,j) -  theta_plus * Vol(i,j)/Dt 
            Fgrav(i,j) = Fgravn(i,j) + Fgravs(i,j) + Fgravw(i,j) + Fgrave(i,j)
        endif
     endif  

     !========== If theta > porosity after refreezing, excess theta goes to cell below through Fgravs(i,j) and Fgravn(i,j-1) ============
     if (theta_plus_refreeze .gt. 0._dp)  then
        if (j .gt. 2) then
            Fgravs(i,j) = Fgravs(i,j) -  theta_plus_refreeze * Vol(i,j)  / Dt
            Fgrav(i,j) = Fgravn(i,j) + Fgravs(i,j) + Fgravw(i,j) + Fgrave(i,j)

            Fgravn(i,j-1) = - Fgravs(i,j)
            Fgrav(i,j-1) = Fgravn(i,j-1) + Fgravs(i,j-1) + Fgravw(i,j-1) + Fgrave(i,j-1)
        else
            Fgravs(i,j) = Fgravs(i,j) -  theta_plus_refreeze * Vol(i,j)/Dt 
            Fgrav(i,j) = Fgravn(i,j) + Fgravs(i,j) + Fgravw(i,j) + Fgrave(i,j)
        endif
     endif  

     if (j .eq. 1) then
         outflow(i) = outflow(i) - Fgravs(i,1)*Dt
     endif

     !=========== Update variables ==============================
     water_content(i,j,1) = water_content(i,j,2)
     S(i,j) = max((water_content(i,j,1)-irreducible)/(porosity(i,j)-irreducible),0.0_dp) 
     P(i,j) = P_new

  end do
end do


end subroutine
