! Writes outputs in .VTS (for Paraview) and .dat formats
subroutine write_outputs(time,mass_conservation)
use declarations
implicit none
real(dp), intent(in) :: time, &               ! Time since beginning of simulation [s]
                    &   mass_conservation     ! Mass conservation coefficient (=1 if mass conserved)


if (time .eq. 0._dp) then

   CALL VTSWriter(0,N,M,xn,yn,water_content(:,:,1),'theta_w','ini')
   CALL VTSWriter(0,N,M,xn,yn,T(:,:,1),'tempera','ini')
   CALL VTSWriter(0,N,M,xn,yn,dry_density(:,:),'density','ini')
   CALL VTSWriter(0,N,M,xn,yn,grain_opt(:,:),'grain_d','ini')
    

   open(18,file='./outputs/outflow.dat')
   open(19,file='./outputs/Tss.dat')
   open(20,file='./outputs/mass_conservation.dat')

elseif (mod(int(time/Dt_ini),iterations) .eq. 0 .and. time .lt. tf )  then
   if (int((time-Dt)/Dt_ini) .ne. int(time/Dt_ini)) then
   
      CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,water_content(:,:,1),'theta_w','int')
      CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,T(:,:,1),'tempera','int')
      CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,dry_density(:,:),'density','int')
      CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,grain_opt(:,:),'grain_d','int')

      write(18,*) (outflow(i),i=1,N-1)
      write(19,*) (Tss(i),i=1,N-1)
      write(20,*) mass_conservation

      write(*,*) 'time step = ', int(time/Dt_ini), '   snowdepth [m] = ', L2, '   Tss [C] = ', Tss(1)
   endif

elseif (time .ge. tf ) then

   CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,water_content(:,:,1),'theta_w','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,T(:,:,1),'tempera','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,dry_density(:,:),'density','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations,N,M,xn,yn,grain_opt(:,:),'grain_d','end')


   write(18,*) (outflow(i),i=1,N-1)
   close(18)
   write(19,*) (Tss(i),i=1,N-1)
   close(19)
   write(20,*) mass_conservation
   close(20)

   write(*,*) 'Final time [s] = ', time, 'outflow [m3] = ',sum(outflow(:)), 'snowdepth [m] = ', L2

   deallocate (Fgrav,P,alpha,nn,mm)
   deallocate (K,K_south,Fpres)
   deallocate (K_east,K_west,K_north,Ks,Fcond)
   deallocate (k_eff,k_eff_west,k_eff_east,k_eff_south,k_eff_north)
   deallocate (Qin,outflow)
   deallocate (Fgravw,Fgravn,Fgravs,grain_opt)
   deallocate (Fprese,Fpresw,Fpresn,Fpress,Fgrave)
   deallocate (Fconde,Fcondw,Fcondn,Fconds)
   deallocate (Tss,alpha_wetting)
   deallocate (water_content,grain_classic,porosity,T)
   deallocate (S,density,dry_density,Hn)
   deallocate (wrc,sphericity,dendricity)
   deallocate (xn,yn,xv,yv, Vol, order_scanning, theta_s, theta_r)
   deallocate (Sx_west,Sy_west,Sx_east,Sy_east)
   deallocate (Sx_north,Sy_north,Sx_south,Sy_south)
   deallocate (Pid, Pdi, theta_di, theta_id)

   stop

elseif (sum(Vol(:,:)) .eq. 0._dp) then ! Everything melted, stop model

   CALL VTSWriter(int(time/Dt_ini)/iterations+1,N,M,xn,yn,water_content(:,:,1),'theta_w','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations+1,N,M,xn,yn,T(:,:,1),'tempera','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations+1,N,M,xn,yn,dry_density(:,:),'density','end')
   CALL VTSWriter(int(time/Dt_ini)/iterations+1,N,M,xn,yn,grain_opt(:,:),'grain_d','end')

   write(18,*) (outflow(i),i=1,N-1)
   close(18)
   
   write(19,*) (Tss(i),i=1,N-1)
   close(19)
   write(20,*) mass_conservation
   close(20)

   write(*,*) 'Final time [s] = ', time, 'outflow [m3] = ',sum(outflow(:)), 'snowdepth [m] = ', L2

   deallocate (Fgrav,P,alpha,nn,mm)
   deallocate (K,K_south,Fpres)
   deallocate (K_east,K_west,K_north,Ks,Fcond)
   deallocate (k_eff,k_eff_west,k_eff_east,k_eff_south,k_eff_north)
   deallocate (Qin,outflow)
   deallocate (Fgravw,Fgravn,Fgravs,grain_opt)
   deallocate (Fprese,Fpresw,Fpresn,Fpress,Fgrave)
   deallocate (Fconde,Fcondw,Fcondn,Fconds)
   deallocate (Tss,alpha_wetting)
   deallocate (water_content,grain_classic,porosity,T)
   deallocate (S,density,dry_density,Hn)
   deallocate (wrc,sphericity,dendricity)
   deallocate (xn,yn,xv,yv, Vol, order_scanning, theta_s, theta_r)
   deallocate (Sx_west,Sy_west,Sx_east,Sy_east)
   deallocate (Sx_north,Sy_north,Sx_south,Sy_south)
   deallocate (Pid, Pdi, theta_di, theta_id)
   stop
endif


end subroutine
