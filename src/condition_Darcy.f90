! Makes sure that water does not flow from dry to wet cell due to negative pressure gradient
subroutine condition_Darcy
use Declarations
implicit none


!======== Make sure that water does not flow from dry to wet cells ===================
! e.g. Fgravs(i,j) + Fpress(i,j) < 0 means water leaves the cell through south face (i,j)

if (M .gt. 2) then ! more than one snow layer

         do j = 1,M-1
            do i = 1,N-1
         !============= Middle cells ========================
               if (i .ne. 1 .and. i .ne. N-1 .and. j .ne. 1 .and. j .ne. M-1) then
                  if (water_content(i,j,1) .le. irreducible) then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp)     then    
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     end if

                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp
                     end if
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     end if
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     end if
                  endif
               endif


               if (i .ne. 1 .and. i .ne. N-1 .and. j .eq. 1) then         
                  !============= South cells except corners (j = 1) ==============
                  if  (water_content(i,j,1) .le. irreducible) then
                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp
                     end if
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     end if
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     end if
                  end if
               end if    


               if (i .ne. 1 .and. i .ne. N-1 .and. j .eq. M-1) then         
                  !============== North cells (except corners) (j = M-1) =========
                  if  (water_content(i,j,1) .le. irreducible)  then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp) then
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     end if
                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp
                     end if
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     end if
                  end if
               end if


               if (j .ne. 1 .and. j .ne. M-1 .and. i .eq. N-1) then         
                  !============== East cells (except upper and lower cells) (i = N-1) =======================
                  if (water_content(i,j,1) .le. irreducible) then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp) then
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     end if
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     end if
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     end if
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                           Fprese(i,j) = 0.0_dp
                           Fgrave(i,j) = 0.0_dp

                           Fpresw(1,j) = 0.0_dp
                           Fgravw(1,j) = 0.0_dp
                        end if
                     end if
                  end if
               endif  
      
               if (j .ne. 1 .and. j .ne. M-1 .and. i .eq. 1) then         
                  !============ West cells (except upper and lower cells) (i=1) ==================
                  if  (water_content(i,j,1) .le. irreducible) then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp) then
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     end if
                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp
                     end if
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     end if
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                           Fprese(N-1,j) = 0.0_dp
                           Fgrave(N-1,j) = 0.0_dp

                           Fpresw(i,j) = 0.0_dp
                           Fgravw(i,j) = 0.0_dp
                        end if
                     end if
                  end if
               endif


               if (j .eq. M-1 .and. i .eq. N-1) then         
                  !============ upper right cell (i = N-1; j = M-1) ===========================
                  if  (water_content(i,j,1) .le. irreducible) then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp) then
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     endif
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     endif
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                           Fprese(i,j) = 0.0_dp
                           Fgrave(i,j) = 0.0_dp

                           Fpresw(1,j) = 0.0_dp
                           Fgravw(1,j) = 0.0_dp
                         end if
                     end if
                  endif
               endif

               if (j .eq. 1 .and. i .eq. N-1) then         
                  !==================== lower right cell (i=N-1; j = 1) ============================
                  if  (water_content(i,j,1) .le. irreducible)  then
                     if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                        Fpresw(i,j) = 0.0_dp
                        Fgravw(i,j) = 0.0_dp

                        Fprese(i-1,j) = 0.0_dp
                        Fgrave(i-1,j) = 0.0_dp
                     endif
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     endif
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                           Fprese(i,j) = 0.0_dp
                           Fgrave(i,j) = 0.0_dp

                           Fpresw(1,j) = 0.0_dp
                           Fgravw(1,j) = 0.0_dp
                         end if
                     end if
                  endif
               endif

               if (j .eq. 1 .and. i .eq. 1) then         
                  !============= Lower left cell (i=1; j = 1) ==================================
                  if  (water_content(i,j,1) .le. irreducible)  then
                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp 
                     endif
                     if ((Fgravn(i,j)+Fpresn(i,j)) .le. 0._dp) then
                        Fpresn(i,j) = 0.0_dp
                        Fgravn(i,j) = 0.0_dp

                        Fpress(i,j+1) = 0.0_dp
                        Fgravs(i,j+1) = 0.0_dp
                     endif
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                           Fprese(N-1,j) = 0.0_dp
                           Fgrave(N-1,j) = 0.0_dp

                           Fpresw(i,j) = 0.0_dp
                           Fgravw(i,j) = 0.0_dp
                         end if
                     end if
                  endif
               endif

               if (j .eq. M-1 .and. i .eq. 1) then         
                  !========== Upper left cell (i = 1; j = M-1) =============================
                  if  (water_content(i,j,1) .le. irreducible)  then
                     if ((Fgravs(i,j)+Fpress(i,j)) .le. 0._dp) then
                        Fpress(i,j) = 0.0_dp
                        Fgravs(i,j) = 0.0_dp

                        Fpresn(i,j-1) = 0.0_dp
                        Fgravn(i,j-1) = 0.0_dp
                     endif
                     if ((Fgrave(i,j)+Fprese(i,j)) .le. 0._dp) then
                        Fprese(i,j) = 0.0_dp
                        Fgrave(i,j) = 0.0_dp

                        Fpresw(i+1,j) = 0.0_dp
                        Fgravw(i+1,j) = 0.0_dp
                     endif
                     if (choice_lat_BC .eq. 0) then
                        if ((Fgravw(i,j)+Fpresw(i,j)) .le. 0._dp) then
                           Fprese(N-1,j) = 0.0_dp
                           Fgrave(N-1,j) = 0.0_dp

                           Fpresw(i,j) = 0.0_dp
                           Fgravw(i,j) = 0.0_dp
                         end if
                     end if
                  endif
               endif
            enddo
         enddo

else  ! M = 2, one snow layer, j = 1

         do i=2,N-2
            !============= Middle cells ==========================
            if  (water_content(i,1,1) .le. irreducible) then
               if ((Fgrave(i,1)+Fprese(i,1)) .le. 0._dp) then
                  Fprese(i,1) = 0.0_dp
                  Fgrave(i,1) = 0.0_dp

                  Fpresw(i+1,1) = 0.0_dp
                  Fgravw(i+1,1) = 0.0_dp
               end if
               if ((Fgravw(i,1)+Fpresw(i,1)) .le. 0._dp) then
                  Fpresw(i,1) = 0.0_dp
                  Fgravw(i,1) = 0.0_dp

                  Fprese(i-1,1) = 0.0_dp
                  Fgrave(i-1,1) = 0.0_dp
               end if
            endif

         enddo

         !============ Right cell (i = N-1) ========================
         if  (water_content(N-1,1,1) .le. irreducible) then
            if ((Fgravw(N-1,1)+Fpresw(N-1,1)) .le. 0._dp) then
               Fpresw(N-1,1) = 0.0_dp
               Fgravw(N-1,1) = 0.0_dp

               Fprese(N-2,1) = 0.0_dp
               Fgrave(N-2,1) = 0.0_dp
            endif
         endif
         if (choice_lat_BC .eq. 0) then
            if ((Fgrave(N-1,1)+Fprese(N-1,1)) .le. 0._dp) then
               Fprese(N-1,1) = 0.0_dp
               Fgrave(N-1,1) = 0.0_dp

               Fpresw(1,1) = 0.0_dp
               Fgravw(1,1) = 0.0_dp
             end if
         endif

         if (choice_lat_BC .eq. 0) then
            if ((Fgrave(N-1,1)+Fprese(N-1,1)) .ge. 0._dp) then
               Fprese(N-1,1) = 0.0_dp
               Fgrave(N-1,1) = 0.0_dp

               Fpresw(1,1) = 0.0_dp
               Fgravw(1,1) = 0.0_dp
             end if
         endif

         !============= Left cell (i = 1) =================
         if  (water_content(1,1,1) .le. irreducible)  then
            if ((Fgrave(1,1)+Fprese(1,1)) .le. 0._dp) then
               Fprese(1,1) = 0.0_dp
               Fgrave(1,1) = 0.0_dp

               Fpresw(2,1) = 0.0_dp
               Fgravw(2,1) = 0.0_dp
            endif
         endif
         if (choice_lat_BC .eq. 0) then
            if ((Fgravw(1,1)+Fpresw(1,1)) .le. 0._dp) then
               Fprese(N-1,1) = 0.0_dp
               Fgrave(N-1,1) = 0.0_dp

               Fpresw(1,1) = 0.0_dp
               Fgravw(1,1) = 0.0_dp
             end if
         endif

         if (choice_lat_BC .eq. 0) then
            if ((Fgravw(1,1)+Fpresw(1,1)) .ge. 0._dp) then
               Fprese(N-1,1) = 0.0_dp
               Fgrave(N-1,1) = 0.0_dp

               Fpresw(1,1) = 0.0_dp
               Fgravw(1,1) = 0.0_dp
             end if
         endif
endif




end subroutine
