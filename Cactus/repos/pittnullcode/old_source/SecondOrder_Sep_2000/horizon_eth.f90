module horizon_eth

   use eth
   public eth1, eth2, edgeintc
   private

contains

subroutine eth1 (nn, output, input, spin, e1)
!-----------------------------------------------------------------------
! wrapper for the d1 routine, calls it once for each patch
!-----------------------------------------------------------------------

   implicit none

   integer,                             intent (in)  :: nn
   double complex, dimension (nn,nn,2), intent (out) :: output
   double complex, dimension (nn,nn,2), intent (in)  :: input
   integer,                             intent (in)  :: spin, e1

   call d1 (nn, output(1:nn,1:nn,1), input(1:nn,1:nn,1), spin, e1)
   call d1 (nn, output(1:nn,1:nn,2), input(1:nn,1:nn,2), spin, e1)

   call cnsint (nn, output(1:nn,1:nn,2), output(1:nn,1:nn,1), spin + e1)
   call cnsint (nn, output(1:nn,1:nn,1), output(1:nn,1:nn,2), spin + e1)

end subroutine eth1

subroutine eth2 (nn, output, input, spin, e1, e2)
!-----------------------------------------------------------------------
! wrapper for the d2 routine, calls it once for each patch
!-----------------------------------------------------------------------

   implicit none

   integer,                             intent (in)  :: nn
   double complex, dimension (nn,nn,2), intent (out) :: output
   double complex, dimension (nn,nn,2), intent (in)  :: input
   integer,                             intent (in)  :: spin, e1, e2

   call d2 (nn, output(1:nn,1:nn,1), input(1:nn,1:nn,1), spin, e1, e2)
   call d2 (nn, output(1:nn,1:nn,2), input(1:nn,1:nn,2), spin, e1, e2)

   call cnsint (nn, output(1:nn,1:nn,2), output(1:nn,1:nn,1), spin + e1 + e2)
   call cnsint (nn, output(1:nn,1:nn,1), output(1:nn,1:nn,2), spin + e1 + e2)

end subroutine eth2

subroutine edgeintc (nn, field, spin)
!-----------------------------------------------------------------------
! wrapper for the cnsint routine, calls it once for each patch
!-----------------------------------------------------------------------

   implicit none

   integer,                             intent (in)    :: nn, spin
   double complex, dimension (nn,nn,2), intent (inout) :: field

   call cnsint (nn, field(1:nn,1:nn,2), field(1:nn,1:nn,1), spin)
   call cnsint (nn, field(1:nn,1:nn,1), field(1:nn,1:nn,2), spin)

end subroutine edgeintc

end module horizon_eth
