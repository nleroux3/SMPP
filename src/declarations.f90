! Declaration of all the variables used in the model
module declarations

integer, parameter :: dp = selected_real_kind(15,307)

real(dp), parameter :: rho_i = 917._dp, &        ! Density of ice [kg/m3]
                   &   Cpi = 2110._dp, &         ! Heat capacity of ice [J/kg/K]
                   &   Cpa = 1003.5_dp, &        ! Heat capacity of air [J/kg/K] 
                   &   rho_a = 1._dp , &         ! Air density [kg/m3]
                   &   Cpw = 4181.3_dp, &        ! Heat capacity of water [J/kg/K]
                   &   rho_w = 1000._dp, &       ! Density of water [kg/m3]
                   &   Lf = 338000._dp, &        ! latent heat of fusion [J/kg]
                   &   irreducible = 0.024_dp    ! Irreducible water content [m3/m3]
integer :: N, &                  ! Number of horizontal nodes
       &   M, &                  ! Number of vertical nodes
       &   i, &                  ! Dummy variable
       &   j, &                  ! Dummy variable
       &   iterations, &         ! Number of iterations for plotting and read met data
       &   index_ground, &       ! Choice between temperature or heat flux at the interface snow-soil
       &   index_energy, &       ! Choice between constant heat flux at surface of met data 
       &   choice_lat_BC, &      ! Choice for the lateral boundary conditions
       &   choice_turb_corr, &   ! Choice between Richardson number of Monin-Obukhov length for atmospheric stability
       &   choice_thermal_eq     ! Choice between thermal equilibrium  (1) or not (0)
real(dp) :: L1, &                ! Horitontal length of the snowpack [m]
        &   L2, &                ! Depth of the snowpack [m]
        &   tf, &                ! Final time specified in input file [s]
        &   Dt, &                ! Time step [s]
        &   beta, &              ! Snowpack angle [rad]
        &   Tground, &           ! Temperature at snow-soil interface [C]
        &   Dt_ini, &            ! Initial time step specified in input file [s]
        &   sum_fluxes, &        ! Sum of all the mass fluxes at the boudnaries [m3]
        &   Qground              ! Ground heat flux [W/m2]
real(dp),allocatable :: xn(:,:), &               ! X position of nodes [m]
                    &   yn(:,:), &               ! Y position of nodes [m]
                    &   xv(:,:), &               ! X position of middle of cells [m]
                    &   yv(:,:), &               ! Y position of middle of cells [m]
                    &   Vol(:,:), &              ! Volume of each cell [m3]
                    &   S(:,:), &                ! Effective water saturation [m3/m3]
                    &   Fgrav(:,:), &            ! Mass flux due to gravity forces [m3/s]
                    &   Fpres(:,:), &            ! Mass flux due to pressure forces [m3/s]
                    &   P(:,:), &                ! Head pressure [m]
                    &   Hn(:), &                 ! Heat flux at snow surface [W/m2]
                    &   Sx_west(:,:), &          ! Surface west boundary in x direction [m2]
                    &   Sy_west(:,:), &          ! Surface west boundary in y direction [m2]
                    &   Sx_east(:,:), &          ! Surface east boundary in x direction [m2]
                    &   Sy_east(:,:), &          ! Surface east boundary in y direction [m2]
                    &   Sx_north(:,:), &         ! Surface north boundary in x direction [m2]
                    &   Sy_north(:,:), &         ! Surface north boundary in y direction [m2]
                    &   Sx_south(:,:), &         ! Surface south boundary in x direction [m2]
                    &   Sy_south(:,:), &         ! Surface south boundary in y direction [m2]
                    &   density(:,:), &          ! Snow bulk density [kg/m3]
                    &   water_content(:,:,:), &  ! Volumetric water content [m3]
                    &   grain_classic(:,:), &    ! Classic grain diameter [m]
                    &   grain_opt(:,:), &        ! Optical grain diameter [m]
                    &   porosity(:,:), &         ! Snow porosity = air content + water content [-]
                    &   K(:,:), &                ! Hydraulic conductivity [m/s]
                    &   K_south(:,:), &          ! Hydraulic conductivity at south boundary [m/s]
                    &   K_north(:,:), &          ! Hydraulic conductivity at north boundary [m/s]
                    &   K_east(:,:), &           ! Hydraulic conductivity at east boundary [m/s]
                    &   K_west(:,:), &           ! Hydraulic conductivity at west boundary [m/s]
                    &   Ks(:,:), &               ! Saturated hydraulic conductivity [m/s]
                    &   mm(:,:), &               ! Parameter in van Genuchten model 
                    &   nn(:,:), &               ! Parameter in van Genuchten model 
                    &   alpha(:,:), &            ! Parameter in van Genuchten model for drainage
                    &   T(:,:,:), &              ! Snow temperature [C]
                    &   k_eff(:,:), &            ! Effective thermal conductivity [
                    &   k_eff_west(:,:), &       ! Effective thermal conductivity at west boundary [W/(mK)]
                    &   k_eff_east(:,:), &       ! Effective thermal conductivity at east boundary [W/(mK)]
                    &   k_eff_south(:,:), &      ! Effective thermal conductivity at south boundary [W/(mK)]
                    &   k_eff_north(:,:), &      ! Effective thermal conductivity at north boundary [W/(mK)]
                    &   Tss(:), &                ! Snow surface temperature [C]
                    &   Qin(:), &                ! Input water flux estimated from melt [m/s]
                    &   dry_density(:,:), &      ! Dry snow density [kg/m3]
                    &   outflow(:), &            ! Snowpack runoff [m3]
                    &   Fprese(:,:), &           ! Mass flux due to pressure forces  at east boundary [m3/s]
                    &   Fpresw(:,:), &           ! Mass flux due to pressure forces  at west boundary [m3/s]
                    &   Fpresn(:,:), &           ! Mass flux due to pressure forces  at north boundary [m3/s]
                    &   Fpress(:,:), &           ! Mass flux due to pressure forces  at south boundary [m3/s]
                    &   Fgrave(:,:), &           ! Mass flux due to gravity forces  at east boundary [m3/s]
                    &   Fgravw(:,:), &           ! Mass flux due to gravity forces  at west boundary [m3/s]
                    &   Fgravn(:,:), &           ! Mass flux due to gravity forces  at north boundary [m3/s]
                    &   Fgravs(:,:), &           ! Mass flux due to gravity forces  at south boundary [m3/s]
                    &   Fconv(:,:), &            ! Heat convection between air and ice in ventilation layer [W/m2]
                    &   Fcond(:,:), &            ! Heat flux due to conduction [W/m2]
                    &   Fconde(:,:), &           ! Heat flux due to conduction at east boundary [W/m2]
                    &   Fcondw(:,:), &           ! Heat flux due to conduction at west boundary [W/m2]
                    &   Fcondn(:,:), &           ! Heat flux due to conduction at north boundary [W/m2]
                    &   Fconds(:,:), &           ! Heat flux due to conduction at south boundary [W/m2]
                    &   sphericity(:,:), &       ! Sphericity of snow grains
                    &   dendricity(:,:), &       ! Dendricity of snow grains
                    &   rain(:), &               ! Rain flux [m/s]
                    &   alpha_wetting(:,:), &    ! Parameter in van Genuchten model for imbibition [1/m]
                    &   Pid(:,:), &              ! Heat pressure at the reversal point imbibition/drainage [m]
                    &   Pdi(:,:), &              ! Heat pressure at the reversal point draingage/imbibition [m]
                    &   theta_di(:,:,:), &       ! Water content at the reversal point imbibition/drainage [-]
                    &   theta_id(:,:,:), &       ! Water content at the reversal point drainage/imbibition [-]
                    &   theta_s(:,:,:), &        ! Variable used to compute scanning curves
                    &   theta_r(:,:,:)           ! Variable used to compute scanning curves
 
integer, allocatable :: wrc(:,:) , &             ! Variable to know which water retention curve is used
                    &   order_scanning(:,:)      ! Order of the scanning curve (max of 90)

end module
