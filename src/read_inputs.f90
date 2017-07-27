! Read inputs from file "inputs.txt" 
subroutine read_inputs
use Declarations 
implicit none
integer :: nb_layers, &                  ! Number of snow layers
       &   sss                           ! Variable
real(dp) :: thickness, &                 ! Thickness of each snow layer
        &   pert_grain, &                ! Perturbation in grain size [decimal value]
        &   pert_dens, &                 ! Perturbation in density [decimal value]
        &   rr, &                        ! Random variable
        &   rrr, &                       ! Random variable
        &   grain_classic_ini(1,M-1), &  ! Initial grain sizes from input file
        &   dry_density_ini(1,M-1)       ! Initial densities from input file
real(dp), allocatable :: cumu_layers(:)  ! Cumulative thicknesses


!========== Open file "snow_layers.txt" ==============
open(7,file='inputs.txt') 
read(7,*) 
read(7,*) 
read(7,*)   
read(7,*) 
read(7,*) pert_grain
read(7,*) pert_dens
read(7,*) beta  ! ground slope angle [rad]
read(7,*) nb_layers
read(7,*) tf ! final time [s]
read(7,*) Dt_ini ! time step [s]
read(7,*) iterations


!================================== Read energy fluxes =======================================
read(7,*) rain(1) ! [m/hr]

rain(2:N-1) = rain(1)

read(7,*) Tss(1)
Tss(2:N-1)=Tss(1)

read(7,*) Hn(1)
Hn(2:N-1) = Hn(1)

read(7,*) index_ground ! 0: temperature [C] at soil-snow ; 1: heat flux from soil (>0) [W/m2]
read(7,*) Tground
Qground = Tground

read(7,*) choice_lat_BC
read(7,*)


allocate (cumu_layers(nb_layers))

read(7,*) thickness
cumu_layers(1) = thickness
do sss=2,nb_layers
   read(7,*) thickness
   cumu_layers(sss) = cumu_layers(sss-1) + thickness
enddo

!================= Read layer properties ===================================
read(7,*)
read(7,*) dry_density_ini(1,1),water_content(1,1,1),grain_classic_ini(1,1),sphericity(1,1),dendricity(1,1)
sss=1
do j = 2,M-1
   if (yn(1,j)/dcos(beta) .lt. cumu_layers(sss)) then
      dry_density_ini(1,j) = dry_density_ini(1,j-1)
      grain_classic_ini(1,j) = grain_classic_ini(1,j-1)
      water_content(1,j,1) = water_content(1,j-1,1)
      sphericity(1,j) = sphericity(1,j-1)
      dendricity(1,j) = dendricity(1,j-1)
   else
      read(7,*) dry_density_ini(1,j),water_content(1,j,1),grain_classic_ini(1,j),sphericity(1,j),dendricity(1,j)
      sss=sss+1

   endif
enddo



do j = 1,M-1
   do i = 1,N-1

      call random_number(rr)
      call random_number(rrr)   
      grain_classic(i,j) = grain_classic_ini(1,j)*dabs(1.0_dp+pert_grain*dsqrt(-2._dp * dlog(rr)) * dcos(2._dp*3.14_dp*rrr))

      call random_number(rr)
      call random_number(rrr)   
      dry_density(i,j) = dry_density_ini(1,j)*dabs(1.0_dp+pert_dens*dsqrt(-2._dp * dlog(rr)) * dcos(2._dp*3.14_dp*rrr))

      water_content(i,j,1) = water_content(1,j,1)
      sphericity(i,j) = sphericity(1,j)
      dendricity(i,j) = dendricity(1,j)
      
    enddo
enddo


do j = 1,M-1
   do i = 1,N-1
      porosity(i,j) = 1._dp-dry_density(i,j)/rho_i
      S(i,j) = max((water_content(i,j,1)-irreducible)/(porosity(i,j)-irreducible),0._dp)  
      density(i,j) = dry_density(i,j)+rho_w*water_content(i,j,1) 

      alpha(i,j) = 4.4e6_dp * (grain_classic(i,j)/dry_density(i,j))**(0.98_dp)    ! (Yamaguchi et al., 2012)
      nn(i,j) = 1.0_dp + 2.7e-3_dp * (grain_classic(i,j)/dry_density(i,j))**(-0.61_dp)
      mm(i,j) = 1.0_dp - 1.0_dp/nn(i,j)
      alpha_wetting(i,j) = 2._dp * alpha(i,j)
      
      if (water_content(i,j,1) .gt. irreducible) then ! On main wetting curve
         wrc(i,j) = 1
         P(i,j) = -((S(i,j)**(-1.0_dp/mm(i,j)) - 1.0_dp) **(1.0_dp/nn(i,j))) / alpha_wetting(i,j)   
      else  ! on WEP
         wrc(i,j) = 0
         P(i,j) = -(0.0437_dp / grain_classic(i,j) * 1e-3_dp + 0.01074_dp)
      endif

  

   enddo
enddo

water_content(:,:,2) = water_content(:,:,1)


read(7,*)
sss = 1
read(7,*) T(1,1,1)

do j = 2,M-1
   if (yv(1,j).lt.cumu_layers(sss)) then
      T(1,j,1) = T(1,j-1,1)
   else
      read(7,*) T(1,j,1)
      sss=sss+1
   endif
enddo



do j=1,M-1
   do i=1,N-1
      T(i,j,1) = T(1,j,1)
    enddo
enddo
close (7)


T(:,:,2) = T(:,:,1)


deallocate(cumu_layers)


return
end
