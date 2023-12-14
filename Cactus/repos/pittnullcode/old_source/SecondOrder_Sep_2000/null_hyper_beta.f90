module null_hyper_beta

   implicit none

   double precision, dimension (:,:), allocatable, save, private :: dx_k
   double complex,   dimension (:,:), allocatable, save, private :: dx_j
   integer,          dimension (:,:), allocatable, save, private :: mask

contains

subroutine null_hyper_beta_allocate

   use null_grid, only : nn

   allocate (dx_j(nn,nn) ,dx_k(nn,nn), mask(nn,nn))

end subroutine null_hyper_beta_allocate

subroutine null_beta (i, jns, bns, masks)

   use null_grid
   use null_mask

   implicit none

   integer,                             intent (in)    :: i
   double complex,   dimension (:,:,:), intent (in)    :: jns
   double precision, dimension (:,:,:), intent (inout) :: bns 
   integer,          dimension (:,:),   intent (in)    :: masks

   double precision xhere

   call null_remask (i, masks+1, mask)

   xhere = xh(i-1)

   dx_j = (jns(:,:,i) - jns(:,:,i-1)) / dx

   dx_k = ( sqrt(1. + jns(:,:,i) * conjg(jns(:,:,i))) - sqrt(1. +&
        & jns(:,:,i-1) * conjg(jns(:,:,i-1))) ) / dx

   bns(:,:,i) = bns(:,:,i) * (1 - mask) + mask * (bns(:,:,i-1) + dx&
        & * xhere * (1. - xhere) / 8. * (dx_j * conjg(dx_j) - dx_k *&
        & dx_k) )

end subroutine null_beta

end module null_hyper_beta
