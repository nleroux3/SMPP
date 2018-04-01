! Calculates mass and energy fluxes at the boundaries of the domain
subroutine boundaries
USE declarations
implicit none



   if (M .gt. 2) then  ! more than one horizontal layer
      if (beta .gt. 0.0001)then  ! sloping snowpack
         do j = 2,M-2

            !========================================================================================================
            !=============================== West Boundary(except lower and north cells) ===============================
            !========================================================================================================
            i=1
            !=================== Mass flux due to water pressure gradients ===============================
            Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) + &
               & Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)))
            Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
               & ( Sx_south(i,j) / dabs(xv(i,j-1)-xv(i,j)) + &
               & Sy_south(i,j) / dabs(yv(i,j-1)-yv(i,j)))
            Fpresn(i,j) = K_north(i,j) * (P(i,j+1)-P(i,j)) * &
               & ( Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)) + &
               & Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)))

            if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
                Fpresw(i,j) = 0.0_dp
                Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
                Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                   & ( Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                   & (xv(i,j ) - xn(i,j)) ) + &
                   & Sy_west(i,j) / dabs( yv(N-1,j)-yv(i,j)))
                Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif

            !================ Mass flux due to gravity force ==============================
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

            !====================== Heat flux ==================================
            Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
               & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
               & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
            Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
               & (Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) + &
               & Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)))
            Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
               & (Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)) + &
               & Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)))
            Fcondw(i,j) = 0.0_dp



            !========================================================================================================
            !====================================  East Boundary (except upper and lower cells) ======================
            !========================================================================================================
            i = N-1

            !=========== Mass flux due to water pressure gradients ================================
            Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
               & (Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))) + &
               & Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))) )
            Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
               & (Sy_south(i,j) / (dabs(yv(i,j)-yv(i,j-1))) + &
               & Sx_south(i,j) / (dabs(xv(i,j)-xv(i,j-1))))
            Fpresn(i,j) = K_north(i,j) * (P(i,j+1)-P(i,j)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j+1)-yv(i,j))) + &
               & Sx_north(i,j) / (dabs(xv(i,j+1)-xv(i,j))))


            if (choice_lat_BC .eq. 2) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            elseif (choice_lat_BC .eq. 1) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
                  & ( Sx_east(i,j) / dabs((xn(N,j)-xv(i,j)) + &
                  & (xv(1,j) - xn(1,j)) ) + &
                  & Sy_east(i,j) / dabs( yv(i,j)-yv(1,j)) )
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            endif

            !=========== Mass flux due to gravity force =================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !====================== heat flux ==========================================
            Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
               & (Sy_south(i,j) / (dabs(yv(i,j)-yv(i,j-1))) + &
               & Sx_south(i,j) / (dabs(xv(i,j)-xv(i,j-1))))
            Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j)-yv(i,j+1))) + &
               & Sx_north(i,j) / (dabs(xv(i,j)-xv(i,j+1))))
            Fconde(i,j) = 0.0_dp
            Fcondw(i,j) = k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) * &
               & (Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))) + &
               & Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))))

         enddo
      else   ! flat snowpack (beta == 0)
         do j = 2,M-2

            !========================================================================================================
            !=============================== West Boundary(except lower and north cells) ===============================
            !========================================================================================================
            i=1
            !============ Mass flux due to WATER PRESSURE GRADIENTS =================
            Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * &
               & (P(i+1,j)-P(i,j)) / (dabs(xv(i+1,j)-xv(i,j)))
            Fpress(i,j) = Sy_south(i,j) * K_south(i,j) * &
               & (P(i,j-1)-P(i,j)) / (dabs(yv(i,j)-yv(i,j-1)))
            Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * &
               & (P(i,j+1)-P(i,j)) / (dabs(yv(i,j+1)-yv(i,j)))


             if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
               Fpresw(i,j) = 0.0_dp
               Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
                Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                   & (Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                   & (xv(i,j ) - xn(i,j)) ))
                Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif
            !=================== Mass flux due to gravity force =======================
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

            !====================== Heat flux ======================================
            Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * &
               & (T(i,j-1,1)-T(i,j,1)) / (dabs(yv(i,j)-yv(i,j-1)))
            Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * &
               & (T(i,j+1,1)-T(i,j,1)) / (dabs(yv(i,j)-yv(i,j+1)))
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * &
               & (T(i+1,j,1)-T(i,j,1)) / (dabs(xv(i,j)-xv(i+1,j)))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * &
                  & (T(N-1,j,1)-T(i,j,1)) / dabs(xv(i+1,j)-xv(i,j))

            !========================================================================================================
            !====================================  East Boundary (except upper and lower cells) ======================
            !========================================================================================================
            i = N-1
            !========= Mass flux due to WATER PRESSURE GRADIENTS --------------------------------
            Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) &
               & / (dabs(xv(i,j)-xv(i-1,j)))
            Fpress(i,j) = Sy_south(i,j) * K_south(i,j)*(P(i,j-1)-P(i,j)) &
               & / (dabs(yv(i,j)-yv(i,j-1)))
            Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * (P(i,j+1)-P(i,j)) &
               & / (dabs(yv(i,j)-yv(i,j+1)))

            if (choice_lat_BC .eq. 2) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            elseif (choice_lat_BC .eq. 1) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
                  & ( Sx_east(i,j) / dabs((xn(N,j)-xv(i,j)) + &
                  & (xv(1,j) - xn(1,j)) ))
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            endif

            !=========== Mass flux due to gravity force =================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !============= Heat flux ===========================================
            Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1))&
               & / (dabs(yv(i,j)-yv(i,j-1)))
            Fcondn(i,j) = Sy_south(i,j) * k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1))&
               & / (dabs(yv(i,j)-yv(i,j+1)))
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * &
                  & (T(1,j,1)-T(i,j,1)) / dabs(xv(i,j)-xv(i-1,j))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1))&
               & / (dabs(xv(i,j)-xv(i-1,j)))

         enddo

      endif
   endif



   !==================================== Lower left cell ========================================
   i=1
   if (M .gt. 2) then ! more than two horizontal snow layers
      j=1
      if (beta .gt. 0.0001)then   ! sloping snowpack

            !====================== Mass flux due to water pressure gradient ======================
            Fpress(i,j) = 0.0_dp
            Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) + &
               & Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)))
            Fpresn(i,j) = K_north(i,j) * (P(i,j+1)-P(i,j)) * &
               & ( Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)) + &
               & Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)))

            if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
               Fpresw(i,j) = 0.0_dp
               Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
                Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                   & ( Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                   & (xv(i,j) - xn(i,j)) ) + &
                   & Sy_west(i,j) / dabs( yv(N-1,j)-yv(i,j)))
                Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif

            !=============== Mass flux due to gravity force ==================
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !=================== Heat flux ===============================
            if (index_ground .eq. 0) then  ! use ground temperature
               Fconds(i,j) = k_eff_south(i,j)*(Tground-T(i,j,1)) * &
                  & (Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
                  &  Sx_south(i,j) / dabs(xv(i,j)-(xn(i+1,j)+xn(i,j))*0.5_dp))
            else ! use ground heat flux
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif
            Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
               & (Sy_north(i,j) / dabs(yv(i,j+1)-yv(i,j)) + &
               & Sx_north(i,j) / dabs(xv(i,j+1)-xv(i,j)))
            Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
               & (Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)) + &
               & Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)))
            Fcondw(i,j) = 0.0_dp

      else   ! flat snowpack (beta == 0)

            !============ Mass flux due to water pressure gradients =================================
            Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * (P(i+1,j)-P(i,j)) / (dabs(xv(i+1,j)-xv(i,j)))
            Fpress(i,j) = 0.0_dp
            Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * (P(i,j+1)-P(i,j)) / (dabs(yv(i,j+1)-yv(i,j)))

            if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
               Fpresw(i,j) = 0.0_dp
               Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                  & (Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                  & (xv(i,j) - xn(i,j)) ))
               Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif

            !============ Mass flux due to gravity force ===================================
            Fgrave(i,j) = -K_east(i,j)*Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j)*Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j)*Sy_south(i,j)

            !====================== Heat flux ==============================================
            if (index_ground .eq. 0) then ! temperature at soil-snow
               Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (Tground-T(i,j,1)) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
            else
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif
            Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) / (dabs(yv(i,j)-yv(i,j+1)))
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) / (dabs(xv(i,j)-xv(i+1,j)))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * &
                  & (T(N-1,j,1)-T(i,j,1)) / dabs(xv(i+1,j)-xv(i,j))

      endif
   endif




   !====================== Upper left cells =====================
   i=1
   if (M .gt. 2) then ! more than one layer
      j=M-1
      if (beta .gt. 0.0001) then

            !=========== Mass flux due to water pressure gradients =======================
            Fprese(i,j) = K_east(i,j)*(P(i+1,j)-P(i,j)) * &
               & (Sx_east(i,j) / (dabs(xv(i+1,j)-xv(i,j))) + &
               & Sy_east(i,j) / (dabs(yv(i+1,j)-yv(i,j))))
            Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
               & (Sy_south(i,j) / (dabs(yv(i,j)-yv(i,j-1))) + &
               & Sx_south(i,j) / (dabs(xv(i,j)-xv(i,j-1))))
            Fpresn(i,j) = 0.0_dp

            if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
               Fpresw(i,j) = 0.0_dp
               Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                  & ( Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                  & (xv(i,j) - xn(i,j))) + &
                  & Sy_west(i,j) / dabs( yv(N-1,j)-yv(i,j)) )
               Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif

            !============= Mass flux due to gravity force ======================
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !======================== heat flux =========================================
            Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
               & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
               & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
            Fcondn(i,j) = k_eff_north(i,j) * (Tss(i)-T(i,j,1)) * &
               & (Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
               & Sx_north(i,j) / dabs(xv(i,j)-(xn(i,j+1)+xn(i+1,j+1))*0.5_dp))
            Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
               & (Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)) + &
               & Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)))
            Fcondw(i,j) = 0.0_dp

      else   ! flat snowpack (beta == 0)

            !================= Mass flux due to water pressure gradients ================================
            Fprese(i,j) = Sx_east(i,j)*K_east(i,j)*(P(i+1,j)-P(i,j))/(dabs(xv(i+1,j)-xv(i,j)))
            Fpress(i,j) = Sy_south(i,j)*K_south(i,j)*(P(i,j-1)-P(i,j))/(dabs(yv(i,j)-yv(i,j-1)))
            Fpresn(i,j) = 0.0_dp


            if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
               Fpresw(i,j) = 0.0_dp
               Fgravw(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
                  & (Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
                  & (xv(i,j) - xn(i,j)) ))
               Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            endif

            !=========== Mass flux due to gravity force ================================
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !======================= heat flux ==========================================
            Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1))/(dabs(yv(i,j)-yv(i,j-1)))
            Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (Tss(i)-T(i,j,1))/(dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1))/(dabs(xv(i,j)-xv(i+1,j)))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * &
                  & (T(N-1,j,1)-T(i,j,1)) / dabs(xv(i+1,j)-xv(i,j))

      endif
   endif


   if (M .gt. 2) then ! more than two horizontal layers
      if (beta .gt. 0.0001) then !sloping snowpack
         do i = 2,N-2

            !==================================================================================================
            !================================ North Boundary (except west and east cells) =============================
            !==================================================================================================
            j=M-1

            !=========== Mass flux due to water pressure gradients =================================
            Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
               & (Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)) + &
               & Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)))
            Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
               & (Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) + &
               & Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)))
            Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
               & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
               & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
            Fpresn(i,j)=0.0_dp

            !=========== Mass flux due to gravity force ================================
            Fgravw(i,j) = K_west(i,j)*Sy_west(i,j)
            Fgrave(i,j) = -K_east(i,j)*Sy_east(i,j)
            Fgravn(i,j) = (Qin(i)+rain(i))*Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j)*Sy_south(i,j)

            !====================== heat flux ===========================================
            Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
               & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
               & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
            Fcondn(i,j) = k_eff_north(i,j) * (Tss(i)-T(i,j,1)) * &
               & (Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
               & Sx_north(i,j) / dabs(xv(i,j)-(xn(i,j+1)+xn(i+1,j+1))*0.5_dp))
            Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
               & (Sy_east(i,j) / dabs(yv(i,j)-yv(i+1,j)) + &
               & Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)))
            Fcondw(i,j) = k_eff_west(i,j)*(T(i-1,j,1)-T(i,j,1)) * &
               & (Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)) + &
               & Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)))


            !================================================================================================
            !========================== South boundary (except left and right cells) ===============================
            !================================================================================================
            j=1
            !=========== Mass flux due to water pressure gradients ===============================
             Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
               & (Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))) + &
               & Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))))
            Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
               & (Sx_east(i,j) / (dabs(xv(i+1,j)-xv(i,j))) + &
               & Sy_east(i,j) / (dabs(yv(i+1,j)-yv(i,j))))
            Fpress(i,j) = 0.0_dp
            Fpresn(i,j) = K_north(i,j) * (P(i,j+1)-P(i,j)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j+1)-yv(i,j))) + &
               & Sx_north(i,j) / (dabs(xv(i,j+1)-xv(i,j))))


            !=========== Mass flux due to gravity force =================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !====================== heat flux ======================================
            if (index_ground .eq. 0) then  ! use temperature at the interface snow-soil
               Fconds(i,j) = k_eff_south(i,j) * (Tground-T(i,j,1)) * &
                  & (Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
                  & Sx_south(i,j) / (dabs(xv(i,j)-(xn(i+1,j)+xn(i,j))*0.5_dp)) )
            else
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif

            Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j)-yv(i,j+1))) + &
               & Sx_north(i,j) / (dabs(xv(i,j)-xv(i,j+1))))
            Fconde(i,j) = k_eff_east(i,j)*(T(i+1,j,1)-T(i,j,1)) * &
               & (Sy_east(i,j) / (dabs(yv(i,j)-yv(i+1,j))) + &
               & Sx_east(i,j) / (dabs(xv(i,j)-xv(i+1,j))) )
            Fcondw(i,j) = k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) * &
               & (Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))) + &
               & Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))))
         enddo

      else   ! flat snowpack (beta == 0)

         do i = 2,N-2


            !==================================================================================================
            !================================ North Boundary (except west and east cells) =============================
            !==================================================================================================
            j=M-1
            !=========== Mass flux due to water pressure gradients ================================
            Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) / &
                 & dabs(xv(i,j)-xv(i-1,j))
            Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * (P(i+1,j)-P(i,j)) / &
                 & dabs(xv(i+1,j)-xv(i,j))
            Fpress(i,j) = Sy_south(i,j) * K_south(i,j) * (P(i,j-1)-P(i,j)) / &
                 & dabs(yv(i,j)-yv(i,j-1))
            Fpresn(i,j) = 0.0_dp

            !=========== Mass flux due to gravity force ================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !====================== heat flux ===========================================
            Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) / &
               & dabs(yv(i,j)-yv(i,j-1))
            Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (Tss(i)-T(i,j,1)) / &
               & (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) / &
               & dabs(xv(i,j)-xv(i+1,j))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) / &
               & dabs(xv(i,j)-xv(i-1,j))


            !================================================================================================
            !========================== South boundary (except left and right cells) ===============================
            !================================================================================================
            j=1
            !=========== Mass flux due to water pressure gradients ================================
            Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) / (dabs(xv(i,j)-xv(i-1,j)))
            Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * (P(i+1,j)-P(i,j)) / (dabs(xv(i+1,j)-xv(i,j)))
            Fpress(i,j) = 0.0_dp
            Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * (P(i,j+1)-P(i,j)) / (dabs(yv(i,j+1)-yv(i,j)))

            !=========== Mass flux due to gravity force =================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !===================== heat flux ============================
            if (index_ground .eq. 0) then ! use of temperature at the interface soil-snow
               Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (Tground-T(i,j,1)) &
                  & / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
            else
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif
            Fcondn(i,j) = Sy_south(i,j) * k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) &
               & / (dabs(yv(i,j)-yv(i,j+1)))
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) &
               & / (dabs(xv(i,j)-xv(i+1,j)))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) &
               & / (dabs(xv(i+1,j)-xv(i,j)))

         enddo
      endif
   endif




   !===================================== Lower right cell ===================================

   i=N-1
   if (M .gt. 2) then ! more than one horizontal snow layer
      j=1
      if (beta .gt. 0.0001) then ! sloping snowpack

            !=========== Mass flux due to water pressure gradients =========================
            Fpresw(i,j) = K_west(i,j)*(P(i-1,j)-P(i,j)) * &
               & (Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))) + &
               & Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))))
            Fpress(i,j) = 0.0_dp
            Fpresn(i,j) = K_north(i,j)*(P(i,j+1)-P(i,j)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j+1)-yv(i,j))) + &
               & Sx_north(i,j) / (dabs(xv(i,j+1)-xv(i,j))))


            if (choice_lat_BC .eq. 2) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            elseif (choice_lat_BC .eq. 1) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
                  & ( Sx_east(i,j) / dabs((xn(N,j)-xv(i,j)) + &
                  & (xv(1,j) - xn(1,j)) ) + &
                  & Sy_east(i,j) / dabs( yv(i,j)-yv(1,j)))
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            endif

            !=========== Mass flux due to gravity force =================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


            !====================== heat flux =======================================
            if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
               Fconds(i,j) = k_eff_south(i,j) * (Tground-T(i,j,1)) * &
                  & (Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
                  & Sx_south(i,j) / (dabs(xv(i,j)-(xn(i+1,j)+xn(i,j))*0.5_dp)))
            else
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif

            Fcondn(i,j) = k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) * &
               & (Sy_north(i,j) / (dabs(yv(i,j)-yv(i,j+1))) + &
               & Sx_north(i,j) / (dabs(xv(i,j)-xv(i,j+1))))
            Fconde(i,j) = 0.0_dp
            Fcondw(i,j) = k_eff_west(i,j)*(T(i-1,j,1)-T(i,j,1)) * &
               & (Sy_west(i,j) / (dabs(yv(i,j)-yv(i-1,j))) + &
               & Sx_west(i,j) / (dabs(xv(i,j)-xv(i-1,j))))


      else   ! flat snowpack (beta == 0)

            !=========== Mass flux due to water pressure gradients ================================
            Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) / (dabs(xv(i,j)-xv(i-1,j)))
            Fpress(i,j) = 0.0_dp
            Fpresn(i,j) = Sy_north(i,j) * K_north(i,j) * (P(i,j+1)-P(i,j)) / (dabs(yv(i,j+1)-yv(i,j)))

            if (choice_lat_BC .eq. 2) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            elseif (choice_lat_BC .eq. 1) then
               Fprese(i,j) = 0.0_dp
               Fgrave(i,j) = 0.0_dp
            elseif (choice_lat_BC .eq. 0) then
               Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
                  & ( Sx_east(i,j) / dabs((xn(N,j)-xv(i,j)) + &
                  & (xv(1,j) - xn(1,j)) ) )
               Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
            endif

            !=========== Mass flux due to gravity force ================================
            Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
            Fgravn(i,j) = K_north(i,j) * Sy_north(i,j)
            Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

            !===================== heat flux ============================
            if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
               Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (Tground-T(i,j,1)) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
            else
               Fconds(i,j) = Sy_south(i,j) * Qground
            endif
            Fcondn(i,j) = Sy_south(i,j) * k_eff_north(i,j) * (T(i,j+1,1)-T(i,j,1)) &
               & /(dabs(yv(i,j)-yv(i,j+1)))
            Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * &
                  & (T(1,j,1)-T(i,j,1)) / dabs(xv(i,j)-xv(i-1,j))
            Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) &
               & /(dabs(xv(i,j)-xv(i-1,j)))

      endif
   endif




   !================================= Upper right cell ====================================
   i=N-1
   if (M .gt. 2) then ! more than two horizontal snow layers
      j=M-1
      if (beta .gt. 0.0001) then ! sloping snowpack

         !=========== Mass flux due to water pressure gradients ================================
         Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
            & (Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)) + &
            & Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)))
         Fpress(i,j) = K_south(i,j) * (P(i,j-1)-P(i,j)) * &
            & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
            & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
         Fpresn(i,j) = 0.0_dp


         if (choice_lat_BC .eq. 2) then
            Fprese(i,j) = 0.0_dp
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
         elseif (choice_lat_BC .eq. 1) then
            Fprese(i,j) = 0.0_dp
            Fgrave(i,j) = 0.0_dp
         elseif (choice_lat_BC .eq. 0) then
            Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ) + &
               & Sy_east(i,j) / dabs( yv(i,j)-yv(1,j)) )
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
         endif
         !=========== Mass flux due to gravity force =================================
         Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
         Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
         Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


         !====================== heat flux ================================
         Fconds(i,j) = k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) * &
            & (Sy_south(i,j) / dabs(yv(i,j)-yv(i,j-1)) + &
            & Sx_south(i,j) / dabs(xv(i,j)-xv(i,j-1)))
         Fcondn(i,j) = k_eff_north(i,j) * (Tss(i)-T(i,j,1)) * &
            & (Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
            & Sx_north(i,j) / dabs(xv(i,j)-(xn(i,j+1)+xn(i+1,j+1))*0.5_dp))
         Fconde(i,j) = 0._dp
         Fcondw(i,j) = k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) * &
            & (Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)) + &
            & Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)))


      else   ! flat snowpack (beta == 0)

         !=========== Mass flux due to water pressure gradients ================================
         Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) / &
            & dabs(xv(i,j)-xv(i-1,j))
         Fpress(i,j) = Sy_south(i,j) * K_south(i,j) * (P(i,j-1)-P(i,j)) / &
            & dabs(yv(i,j)-yv(i,j-1))
         Fpresn(i,j) = 0.0_dp


         if (choice_lat_BC .eq. 2) then
            Fprese(i,j) = 0.0_dp
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
         elseif (choice_lat_BC .eq. 1) then
            Fprese(i,j) = 0.0_dp
            Fgrave(i,j) = 0.0_dp
         elseif (choice_lat_BC .eq. 0) then
            Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ) )
            Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
         endif

         !=========== Mass flux due to gravity force ================================
         Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
         Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
         Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

         !====================== heat flux =====================================
         Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (T(i,j-1,1)-T(i,j,1)) / &
            & dabs(yv(i,j)-yv(i,j-1))
         Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (Tss(i)-T(i,j,1)) / &
            &(dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
         Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * &
                  & (T(1,j,1)-T(i,j,1)) / dabs(xv(i,j)-xv(i-1,j))
         Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) / &
            & dabs(xv(i,j)-xv(i-1,j))

      endif
   endif


   !==================================================================================
   !======================= One snow layer ===========================================
   !==================================================================================


   !===================== Left cell ===================================
   if (M .eq. 2) then
      j=1
      i=1

      !=========== Mass flux due to water pressure gradients ================================
      if (beta .gt. 0.0001) then ! sloping snowpack
         Fprese(i,j) = K_east(i,j) * (P(i+1,j)-P(i,j)) * &
            & (Sx_east(i,j) / dabs(xv(i+1,j)-xv(i,j)) + &
            & Sy_east(i,j) / dabs(yv(i+1,j)-yv(i,j)))
         Fpress(i,j) = 0.0_dp
         Fpresn(i,j) = 0.0_dp
      else ! flat snowpack
         Fprese(i,j) = Sx_east(i,j) * K_east(i,j) * (P(i+1,j)-P(i,j)) / dabs(xv(i+1,j)-xv(i,j))
         Fpress(i,j) = 0.0_dp
         Fpresn(i,j) = 0.0_dp
      endif


      if (choice_lat_BC .eq. 2 .or. choice_lat_BC .eq. 1) then
         Fpresw(i,j) = 0.0_dp
         Fgravw(i,j) = 0.0_dp
      elseif (choice_lat_BC .eq. 0) then
         if (beta .gt. 0.0001) then ! sloping snowpack
            Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
               & ( Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ) + &
               & Sy_west(i,j) / dabs( yv(N-1,j)-yv(1,j)) )
         else
            Fpresw(i,j) = K_west(i,j) * (P(N-1,j)-P(i,j)) * &
               & ( Sx_west(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ) )
         endif
        Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
      endif

      !=========== Mass flux due to gravity force ================================
      Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
      Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
      Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)

      !====================== heat flux ========================================
      if (beta .gt. 0.0001) then ! sloping snowpack
         if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
            Fconds(i,j) = k_eff_south(i,j) * (Tground-T(i,j,1)) * &
               & (Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
               & Sx_south(i,j) / (dabs(xv(i,j)-(xn(i+1,j)+xn(i,j))*0.5_dp)))
         else
            Fconds(i,j) = Sy_south(i,j) * Qground
         endif
         Fcondn(i,j) = k_eff_north(i,j) * (Tss(i)-T(i,j,1)) * &
            & (Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
            & Sx_north(i,j) / (dabs(xv(i,j)-(xn(i,j+1)+xn(i+1,j+1))*0.5_dp)))
         Fconde(i,j) = k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) * &
            & (Sy_east(i,j) / dabs(yv(i,j)-yv(i+1,j)) + &
            & Sx_east(i,j) / dabs(xv(i,j)-xv(i+1,j)))
         Fcondw(i,j) = 0.0_dp

      else
         if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
            Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (Tground-T(i,j,1)) / &
               & (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
         else
            Fconds(i,j) = Sy_south(i,j) * Qground
         endif
         Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (Tss(i)-T(i,j,1)) / &
            & (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
         Fconde(i,j) = Sx_east(i,j) * k_eff_east(i,j) * (T(i+1,j,1)-T(i,j,1)) / &
            & (dabs(xv(i,j)-xv(i+1,j)))
         Fcondw(i,j)=0.0_dp
      endif



   !=========================== Middle cells ==================================

      j=1
      !=========== Mass flux due to water pressure gradients =================================
      if (beta .gt. 0.0001) then ! sloping snowpack
         Fpresw(2:N-2,j) = K_west(2:N-2,j) * (P(1:N-3,j)-P(2:N-2,j)) * &
            & (Sx_west(2:N-2,j) / dabs(xv(2:N-2,j)-xv(1:N-3,j)) + &
            & Sy_west(2:N-2,j) / dabs(yv(2:N-2,j)-yv(1:N-3,j)))
         Fprese(2:N-2,j) = K_east(2:N-2,j) * (P(3:N-1,j)-P(2:N-2,j)) * &
            & (Sx_east(2:N-2,j) / dabs(xv(3:N-1,j)-xv(2:N-2,j)) + &
            & Sy_east(2:N-2,j) / dabs(yv(3:N-1,j)-yv(2:N-2,j)))
         Fpress(2:N-2,j) = 0.0_dp
         Fpresn(2:N-2,j) = 0.0_dp
      else ! flat snowpack
         Fpresw(2:N-2,j) = Sx_west(2:N-2,j) * K_west(2:N-2,j) * (P(1:N-3,j)-P(2:N-2,j)) / &
            & dabs(xv(2:N-2,j)-xv(1:N-3,j))
         Fprese(2:N-2,j) = Sx_east(2:N-2,j) * K_east(2:N-2,j) * (P(3:N-1,j)-P(2:N-2,j)) / &
            & dabs(xv(3:N-1,j)-xv(2:N-2,j))
         Fpress(2:N-2,j) = 0.0_dp
         Fpresn(2:N-2,j) = 0.0_dp
      endif

      !=========== Mass flux due to gravity force =================================
      Fgravw(2:N-2,j) = K_west(2:N-2,j) * Sy_west(2:N-2,j)
      Fgrave(2:N-2,j) = -K_east(2:N-2,j) * Sy_east(2:N-2,j)
      Fgravn(2:N-2,j) = (Qin(2:N-2)+rain(2:N-2)) * Sy_north(2:N-2,j)
      Fgravs(2:N-2,j) = -K_south(2:N-2,j) * Sy_south(2:N-2,j)



      !====================== heat flux =======================================
      if (beta.gt. 0.00001) then ! sloping snowpack
         if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
            Fconds(2:N-2,j) = k_eff_south(2:N-2,j) * (Tground-T(2:N-2,j,1)) * &
               & (Sy_south(2:N-2,j) / (dabs(yn(2:N-2,j+1)-yn(2:N-2,j))*0.5_dp) + &
               & Sx_south(2:N-2,j) / (dabs(xv(2:N-2,j)-(xn(3:N-1,j)+xn(2:N-2,j))*0.5_dp)))
         else
            Fconds(2:N-2,j) = Sy_south(2:N-2,j) * Qground
         endif
         Fcondn(2:N-2,j) = k_eff_north(2:N-2,j) * (Tss(i)-T(2:N-2,j,1)) * &
            & (Sy_north(2:N-2,j) / (dabs(yn(2:N-2,j+1)-yn(2:N-2,j))*0.5_dp) + &
            & Sx_north(2:N-2,j) / (dabs(xv(2:N-2,j)-(xn(2:N-2,j+1)+xn(3:N-1,j+1)) *0.5_dp)))
         Fconde(2:N-2,j) = k_eff_east(2:N-2,j) * (T(3:N-1,j,1)-T(2:N-2,j,1)) * &
            & (Sy_east(2:N-2,j) / dabs(yv(2:N-2,j)-yv(3:N-1,j)) + &
            & Sx_east(2:N-2,j) / dabs(xv(2:N-2,j)-xv(3:N-1,j)))
         Fcondw(2:N-2,j) = k_eff_west(2:N-2,j) * (T(1:N-3,j,1)-T(2:N-2,j,1)) * &
            & (Sy_west(2:N-2,j) / dabs(yv(2:N-2,j)-yv(1:N-3,j)) + &
            & Sx_west(2:N-2,j) / dabs(xv(2:N-2,j)-xv(1:N-3,j)))
      else
         if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
            Fconds(2:N-2,j) = Sy_south(2:N-2,j) * k_eff_south(2:N-2,j) * (Tground-T(2:N-2,j,1)) / &
            & (dabs(yn(2:N-2,j+1)-yn(2:N-2,j))*0.5_dp)
         else
            Fconds(2:N-2,j) = Sy_south(2:N-2,j) * Qground
         endif
         Fcondn(2:N-2,j) = Sy_north(2:N-2,j) * k_eff_north(2:N-2,j) * (Tss(2:N-2)-T(2:N-2,j,1)) / &
            & (dabs(yn(2:N-2,j+1)-yn(2:N-2,j))*0.5_dp)
         Fconde(2:N-2,j) = Sx_east(2:N-2,j) * k_eff_east(2:N-2,j) * (T(3:N-1,j,1)-T(2:N-2,j,1)) / &
            & dabs(xv(2:N-2,j)-xv(3:N-1,j))
         Fcondw(2:N-2,j) = Sx_west(2:N-2,j) * k_eff_west(2:N-2,j) * (T(1:N-3,j,1)-T(2:N-2,j,1)) / &
            & dabs(xv(2:N-2,j)-xv(1:N-3,j))
      endif


   !===================================== Right cell ==========================
      i=N-1
      j=1
      !=========== Mass flux due to water pressure gradients ================================
      if (beta .gt. 0.0001)then ! sloping snowpack
         Fpresw(i,j) = K_west(i,j) * (P(i-1,j)-P(i,j)) * &
            & (Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)) + &
            & Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)))
         Fpress(i,j) = 0.0_dp
         Fpresn(i,j) = 0.0_dp
      else
         Fpresw(i,j) = Sx_west(i,j) * K_west(i,j) * (P(i-1,j)-P(i,j)) / &
            & dabs(xv(i,j)-xv(i-1,j))
         Fpress(i,j) = 0.0_dp
         Fpresn(i,j) = 0.0_dp
      endif


      if (choice_lat_BC .eq. 2) then
         Fprese(i,j) = 0.0_dp
         Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
      elseif (choice_lat_BC .eq. 1) then
         Fprese(i,j) = 0.0_dp
         Fgrave(i,j) = 0.0_dp
      elseif (choice_lat_BC .eq. 0) then
         if (beta .gt. 0.0001) then ! sloping snowpack
            Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ) + &
               & Sy_east(i,j) / dabs( yv(N+1,j)-yv(1,j)) )
         else
            Fprese(i,j) = K_east(i,j) * (P(1,j)-P(i,j)) * &
               & ( Sx_east(i,j) / dabs((xn(N,j)-xv(N-1,j)) + &
               & (xv(1,j) - xn(1,j)) ))

         endif
         Fgrave(i,j) = -K_east(i,j) * Sy_east(i,j)
      endif

      !=========== Mass flux due to gravity force ================================
      Fgravw(i,j) = K_west(i,j) * Sy_west(i,j)
      Fgravn(i,j) = (Qin(i)+rain(i)) * Sy_north(i,j)
      Fgravs(i,j) = -K_south(i,j) * Sy_south(i,j)


      !====================== heat flux ========================================
      if (beta .gt. 0.0001) then ! sloping snowpack
         if (index_ground .eq. 0) then ! use temperature at the interface snow-soil
            Fconds(i,j) = k_eff_south(i,j) * (Tground-T(i,j,1)) * &
               & (Sy_south(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
               & Sx_south(i,j) / (dabs(xv(i,j)-(xn(i+1,j)+xn(i,j))*0.5_dp)))
         else
            Fconds(i,j) = Sy_south(i,j) * Qground
         endif
         Fcondn(i,j) = k_eff_north(i,j) * (Tss(i)-T(i,j,1)) * &
            & (Sy_north(i,j) / (dabs(yn(i,j+1)-yn(i,j))*0.5_dp) + &
            & Sx_north(i,j) / (dabs(xv(i,j)-(xn(i,j+1)+xn(i+1,j+1))*0.5_dp)))
         Fconde(i,j) = 0.0_dp
         Fcondw(i,j) = k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) * &
            & (Sy_west(i,j) / dabs(yv(i,j)-yv(i-1,j)) + &
            & Sx_west(i,j) / dabs(xv(i,j)-xv(i-1,j)))
      else ! flat snowpack
         Fconds(i,j) = Sy_south(i,j) * k_eff_south(i,j) * (Tground-T(i,j,1)) / &
            & (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
         Fcondn(i,j) = Sy_north(i,j) * k_eff_north(i,j) * (Tss(i)-T(i,j,1)) / &
            & (dabs(yn(i,j+1)-yn(i,j))*0.5_dp)
         Fconde(i,j) = 0.0_dp
         Fcondw(i,j) = Sx_west(i,j) * k_eff_west(i,j) * (T(i-1,j,1)-T(i,j,1)) / &
            & dabs(xv(i,j)-xv(i-1,j))
      endif

   endif




return
end
