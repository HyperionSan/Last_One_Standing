module horizon_eth

contains

subroutine eth1 (nn, output, input, spin, e1)
!-----------------------------------------------------------------------
! wrapper for the d1 routine, calls it once for each patch
!-----------------------------------------------------------------------

   use null_eth
   use null_interp
   implicit none

   integer,                             intent (in)  :: nn
   double complex, dimension (nn,nn,2), intent (out) :: output
   double complex, dimension (nn,nn,2), intent (in)  :: input
   integer,                             intent (in)  :: spin, e1

   call null_d1 (output(:,:,1), input(:,:,1), spin, e1)
   call null_d1 (output(:,:,2), input(:,:,2), spin, e1)

   call null_cnsint (output(:,:,1), output(:,:,2), spin + e1)
   call null_cnsint (output(:,:,2), output(:,:,1), spin + e1)

end subroutine eth1

subroutine eth2 (nn, output, input, spin, e1, e2)
!-----------------------------------------------------------------------
! wrapper for the d2 routine, calls it once for each patch
!-----------------------------------------------------------------------

   use null_eth
   use null_interp
   implicit none

   integer,                             intent (in)  :: nn
   double complex, dimension (nn,nn,2), intent (out) :: output
   double complex, dimension (nn,nn,2), intent (in)  :: input
   integer,                             intent (in)  :: spin, e1, e2

   call null_d2 (output(:,:,1), input(:,:,1), spin, e1, e2)
   call null_d2 (output(:,:,2), input(:,:,2), spin, e1, e2)

   call null_cnsint (output(:,:,1), output(:,:,2), spin + e1 + e2)
   call null_cnsint (output(:,:,2), output(:,:,1), spin + e1 + e2)

end subroutine eth2

end module horizon_eth
