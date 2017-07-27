! Converts from classical grain size to optical grain size (Vionnet et al., 2012)
subroutine conversion_grain
use declarations, only : dendricity, &    ! Dendricity of the grain
                     &   sphericity, &    ! Sphericity of the grain
                     &   M, &             ! Number of vertical nodes
                     &   grain_opt, &     ! Optical grain diameter [m] 
                     &   grain_classic, & ! Grain diameter [m]
                     &   dp               ! Precision of float numbers
implicit none

where (dendricity(:,1:M-1) .gt. 0.00001 ) ! dendritic case
   grain_opt(:,1:M-1) = 1e-4_dp*(dendricity(:,1:M-1)+(1.0_dp-dendricity(:,1:M-1)) &
      & *(4.0_dp-sphericity(:,1:M-1)))
elsewhere ! non-dendritic case
   grain_opt(:,1:M-1) = grain_classic(:,1:M-1)*sphericity(:,1:M-1) + &
      & (1.0_dp-sphericity(:,1:M-1))*max(4e-4_dp,grain_classic(:,1:M-1)*0.5_dp)
end where


return
end subroutine
