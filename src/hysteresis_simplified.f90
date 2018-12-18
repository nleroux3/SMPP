! Computes the water pressure knowing the water content (Huang et al., 2015)
subroutine hysteresis(wrc, P_new, P, alpha, alpha_wetting, porosity, theta_1, theta_2, nn, mm, Pdi, &
         & Pid, theta_id, theta_di, grain_classic, theta_plus,theta_s, theta_r, order)
use declarations, only : dp, &              ! Precision of float numbers
                     &   irreducible        ! Irreducible water content [-]
implicit none
real(dp), intent(in) :: alpha, &            ! Parameter in van Genuchten model (for drying boundary curve) [m-1]
                    &   alpha_wetting, &    ! Parameter in van Genuchten model (for wettng boundary curve) [m-1]
                    &   porosity, &         ! Porosity of a cell [-]
                    &   nn, &               ! Parameter in van Genuchten model
                    &   mm, &               ! Parameter in van Genuchten model
                    &   theta_1, &          ! Initial water content [-]
                    &   grain_classic       ! Grain diameter [m]
real(dp), intent(inout) :: P, &             ! Head pressure [m]
                       &   theta_2, &       ! Water content output [-]
                       &   Pdi, &           ! Head pressure at reversal point between drying and imbibition [m]
                       &   Pid, &           ! Head pressure at reversal point between imbibition and drying [m]
                       &   theta_id(12), & ! Water content at reversal point between imbibition and drying [-] 
                       &   theta_di(12), & ! Water content at reversal point between drying and imbibition [-]
                       &   theta_r(12), &  ! Variable used in hysteresis model [-]
                       &   theta_s(12)     ! Variable used in hysteresis model [-]
real(dp), intent(out) :: theta_plus, &      ! Extra water content to avoid model instability if theta > porosity [-]
                     &   P_new              ! New computed head pressure [m] 
real(dp) :: P_w, &                          ! Dummy variable used in this script to compute pressure on wetting or drying boundary curve [m]
        &   Sd_Pdi, &                       ! Saturation on drying boundary curve at reversal point between drying and imbibition
        &   Sd_Pid, &                       ! Saturation on drying boundary curve at reversal point between imbibition and drying
        &   Si_Pdi, &                       ! Saturation on wetting boundary curve at reversal point between drying and imbibition               
        &   Si_Pid, &                       ! Saturation on wetting boundary curve at reversal point between imbibition and drying
        &   theta_w, &                       ! Dummy variable used in this script to compute water content on wetting boundary curve [-]
        &   P_d
integer :: wrc_ini, &                       ! Variable to know on which water retention curve the cell is initally
       &   order                   ! Variable to know the order of the water retention curve
integer, intent(inout) :: wrc               ! Variable to know on which water retention curve the cell is at the end of the script

wrc_ini = wrc


if (wrc_ini .eq. 0) then  ! on water entry pressure 

 
   theta_w = irreducible + (porosity - irreducible) * &
          &  (1._dp + (-alpha_wetting * &
          & P)**nn)**(-mm) ! Theta on main wetting curve
   
   if (theta_2 .gt. theta_w) then ! Moves to main wetting curve
      
      order = 1
      wrc = 1

      theta_s(order) = porosity
      theta_r(order) = irreducible
      theta_di(order) = irreducible
      theta_id(order) = porosity

      P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
         &   **(1.0_dp/nn)) / alpha_wetting  

   else  ! Stays at water entry pressure
      P_new = P
   endif

endif


if (theta_2 .le. irreducible) then  ! If initially wrc > 0 (not at water entry pressure), moves to water entry pressure
    wrc_ini = 0 ! To avoid going to another loop below
    wrc = 0
    P_new = -(0.0437_dp / grain_classic * 1e-3_dp + 0.01074_dp) 
endif

if (wrc_ini .eq. 1) then  ! On main wetting WRC
   order = 1

   Pdi = -huge(1._dp)
   Pid = P
   theta_id(order) = porosity 
   theta_di(order) = irreducible

   if (theta_2 .lt. theta_1) then ! WRC moves to drying scanning curve (id : imbibition to dyring)

      order = order + 1

      wrc = 3

      Pdi = -huge(1._dp)
      Pid = P
      theta_id(order) = theta_1 
      theta_di(order) = irreducible

      Sd_Pdi = 0._dp
      Sd_Pid = (1._dp + (-alpha * Pid)**nn)**(-mm)

      theta_s(order) = (theta_di(order)*(1._dp - Sd_Pid) - theta_id(order) &
              & * (1._dp - Sd_Pdi))  /(Sd_Pdi - Sd_Pid)

      theta_r(order) = (Sd_Pdi * theta_s(order) - theta_di(order)) / (Sd_Pdi - 1._dp)

      P_new = - (((theta_2 - theta_r(order))/(theta_s(order) &
              & - theta_r(order)))**(-1._dp/mm)-1._dp)**(1._dp/nn) /alpha
       

   else ! stays on main wetting curve

      if (theta_2 .gt. (porosity-1e-6_dp)) then ! On main drying curve

         theta_plus = theta_2 - (porosity -1e-6_dp)
         theta_2 = porosity-1e-6_dp
         wrc = 2 ! Moves to main drying curve below
         order = 1

         Pdi = -huge(1._dp)
         Pid = 0._dp

         theta_di(order) = irreducible
         theta_id(order) = porosity

         theta_s(order) = porosity
         theta_r(order) = irreducible
         P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
            &   **(1.0_dp/nn)) / alpha  

      else

         theta_s(order) = porosity
         theta_r(order) = irreducible
         theta_di(order) = irreducible
         theta_id(order) = porosity

         P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
            &   **(1.0_dp/nn)) / alpha_wetting  

      endif
   endif
endif



if (wrc_ini .eq. 3) then  ! On drying scanning curve

      if (theta_2 .gt. (porosity-1e-6_dp)) then ! On main drying curve

         theta_plus = theta_2 - (porosity -1e-6_dp)
         theta_2 = porosity-1e-6_dp
         wrc = 2 ! Moves to main drying curve below
         order = 1

         Pdi = -huge(1._dp)
         Pid = 0._dp

         theta_di(order) = irreducible
         theta_id(order) = porosity

         theta_s(order) = porosity
         theta_r(order) = irreducible
         P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
            &   **(1.0_dp/nn)) / alpha  

      else

         P_new = - (((theta_2 - theta_r(order))/(theta_s(order) &
              & - theta_r(order)))**(-1._dp/mm)-1._dp)**(1._dp/nn)/alpha
       

         P_w = - (((theta_2 - irreducible)/(porosity - irreducible))**(-1._dp/mm)-1._dp)**(1._dp/nn)&
              & /alpha_wetting ! Pressure on the main drying curve


         if (dabs(P_new) .le. dabs(P_w) .or. theta_2 > theta_id(order)) then ! Moves to main drying curve
            wrc = 1
            order = 1
            P_new = P_w
         endif

      endif
endif

if (wrc_ini .eq. 4) then  ! On wetting scanning scanning


      if (theta_2 .gt. (porosity-1e-6_dp)) then ! On main drying curve

         theta_plus = theta_2 - (porosity -1e-6_dp)
         theta_2 = porosity-1e-6_dp
         wrc = 2 ! Moves to main drying curve below
         order = 1

         Pdi = -huge(1._dp)
         Pid = 0._dp

         theta_di(order) = irreducible
         theta_id(order) = porosity

         theta_s(order) = porosity
         theta_r(order) = irreducible
         P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
            &   **(1.0_dp/nn)) / alpha  

      else


         P_new = - (((theta_2 - theta_r(order))/(theta_s(order) &
                 & - theta_r(order)))**(-1._dp/mm)-1._dp)**(1._dp/nn) /alpha_wetting

         P_w = - (((theta_2 - irreducible)/(porosity - irreducible))**(-1._dp/mm)-1._dp)**(1._dp/nn)&
              & /alpha_wetting ! Pressure on the main wetting curve

         P_d = - (((theta_2 - irreducible)/(porosity - irreducible))**(-1._dp/mm)-1._dp)**(1._dp/nn)&
              & /alpha ! Pressure on the main drying curve

         if (dabs(P_new) .le. dabs(P_w)) then ! Moves to main wetting curve
            wrc = 1 
            order = 1
            P_new = P_w
          elseif (dabs(P_new) .ge. dabs(P_d) .or. theta_2 .le. theta_di(order)) then ! Moves to main wetting curve
            wrc = 2 
            order = 1
            P_new = P_d
         endif

      endif

endif

if (wrc_ini .eq. 2) then  ! On main drying curve

   order = 1

   Pdi = -huge(1._dp)
   Pid = 0._dp
   theta_id(order) = porosity 
   theta_di(order) = irreducible

   if (theta_2 .gt. theta_1) then ! Moves to wetting scanning (di: drainage to imbibition)

      if (theta_2 .gt. (porosity-1e-6_dp)) then ! On main drying curve

         theta_plus = theta_2 - (porosity -1e-6_dp)
         theta_2 = porosity-1e-6_dp
         wrc = 2 ! Moves to main drying curve below
         order = 1

         Pdi = -huge(1._dp)
         Pid = 0._dp

         theta_di(order) = irreducible
         theta_id(order) = porosity

         theta_s(order) = porosity
         theta_r(order) = irreducible
         P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
            &   **(1.0_dp/nn)) / alpha  

      else

         order = order + 1
         wrc = 4

         Pdi = P
         Pid = 0._dp
         theta_di(order) = theta_1 
         theta_id(order) = porosity

         Si_Pdi = min((1._dp + (-alpha_wetting * Pdi)**nn)**(-mm), 0.9999_dp)
         Si_Pid = (1._dp + (-alpha_wetting * Pid)**nn)**(-mm)


         theta_s(order) = (theta_di(order)*(1._dp - Si_Pid) - theta_id(order) &
              & * (1._dp - Si_Pdi))  / min(Si_Pdi - Si_Pid, -1e-30_dp)

         theta_r(order) = (Si_Pdi * theta_s(order) - theta_di(order)) / (Si_Pdi - 1._dp)

         P_new = - (((theta_2 - theta_r(order))/(theta_s(order) &
                 & - theta_r(order)))**(-1._dp/mm)-1._dp)**(1._dp/nn) /alpha_wetting

         P_w = - (((theta_2 - irreducible)/(porosity - irreducible))**(-1._dp/mm)-1._dp)**(1._dp/nn)&
               & /alpha_wetting ! Pressure on main wetting curve

         if (dabs(P_new) .le. dabs(P_w)) then ! Moves to main wetting curve
            wrc = 1
            order = 1 
            P_new = P_w
         endif  

      endif

   else  ! Stays on main drying curve

      theta_di(order) = irreducible
      theta_id(order) = porosity

      theta_s(order) = porosity
      theta_r(order) = irreducible
      P_new = -((((theta_2 - irreducible)/(porosity-irreducible))**(-1.0_dp/mm) - 1.0_dp) &
         &   **(1.0_dp/nn)) / alpha  

   endif
endif



end subroutine
